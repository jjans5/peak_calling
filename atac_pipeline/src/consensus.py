"""
Consensus Peak Calling
======================

Functions for creating consensus peaks across multiple samples.
Implements iterative peak filtering to resolve overlaps.
"""

from __future__ import annotations
import logging
import sys
import os
import pandas as pd
import pyranges as pr
import numpy as np
from typing import Dict, Optional, Union, List
import warnings

# Suppress FutureWarning messages from pyranges
warnings.simplefilter(action='ignore', category=FutureWarning)

# =============================================================================
# LOGGING SETUP
# =============================================================================

log = logging.getLogger("ConsensusPeaks")
if not log.handlers:
    handler = logging.StreamHandler(stream=sys.stdout)
    handler.setFormatter(logging.Formatter("%(asctime)s %(name)-12s %(levelname)-8s %(message)s"))
    log.addHandler(handler)
    log.setLevel(logging.INFO)

# =============================================================================
# CONSTANTS
# =============================================================================

NARROWPEAK_COLS = [
    "Chromosome", "Start", "End", "Name", "Score", "Strand",
    "SignalValue", "PValue", "QValue", "Summit"
]


# =============================================================================
# MAIN FUNCTIONS
# =============================================================================

def get_consensus_peaks(
    narrow_peaks_dict: Dict[str, pr.PyRanges],
    peak_half_width: int,
    chromsizes: Optional[Union[pr.PyRanges, pd.DataFrame]] = None,
    path_to_blacklist: Optional[str] = None,
) -> pr.PyRanges:
    """
    Generate consensus peaks from multiple samples.
    
    Parameters
    ----------
    narrow_peaks_dict : dict
        Dictionary mapping sample names to PyRanges of narrowPeak data
    peak_half_width : int
        Half-width for extending peaks from summit (total width = 2 * half_width)
    chromsizes : PyRanges or DataFrame, optional
        Chromosome sizes for filtering peaks that extend beyond chromosomes
    path_to_blacklist : str, optional
        Path to blacklist BED file for filtering problematic regions
    
    Returns
    -------
    PyRanges
        Consensus peaks after merging and filtering
    """
    if isinstance(chromsizes, pd.DataFrame):
        chromsizes = chromsizes.loc[:, ["Chromosome", "Start", "End"]]
        chromsizes = pr.PyRanges(chromsizes)

    log.info("Extending and merging peaks per class")
    
    center_extended_peaks = []
    for x in narrow_peaks_dict.keys():
        # 1. Extend and Intersect
        extended = calculate_peaks_and_extend(
            narrow_peaks_dict[x], peak_half_width, chromsizes, path_to_blacklist
        )
        
        # Safety check
        if len(extended) == 0:
            log.warning(f"⚠️  Sample '{x}' has 0 peaks after extending/intersecting! "
                       "Check chromosome naming (chr1 vs 1).")
            continue 

        # 2. Filter locally
        filtered = iterative_peak_filtering(extended)
        center_extended_peaks.append(filtered.df)

    if not center_extended_peaks:
        raise ValueError("❌ All samples were empty after intersection. Check your chromsizes format!")

    log.info("Normalizing peak scores (CPM)")
    center_extended_peaks_norm = [cpm(pr.PyRanges(x), "Score") for x in center_extended_peaks]
    
    # Combine all samples
    combined_pr = pr.PyRanges(pd.concat([x.df for x in center_extended_peaks_norm], axis=0, sort=False))

    log.info("Merging peaks globally")
    consensus_peaks = iterative_peak_filtering(combined_pr)
    
    log.info(f"Done! Found {len(consensus_peaks)} consensus peaks.")
    return consensus_peaks


def load_narrowpeaks(
    peak_dir: str,
    q_value_threshold: float = 0.05,
    min_peaks_per_sample: int = 5000,
    add_chr_prefix: bool = True,
) -> Dict[str, pr.PyRanges]:
    """
    Load and filter narrowPeak files from a directory.
    
    Parameters
    ----------
    peak_dir : str
        Directory containing narrowPeak files
    q_value_threshold : float
        Q-value threshold for filtering (keep peaks with q-value < threshold)
    min_peaks_per_sample : int
        Minimum number of peaks required to include a sample
    add_chr_prefix : bool
        Whether to add 'chr' prefix to chromosome names if missing
    
    Returns
    -------
    dict
        Dictionary mapping sample names to PyRanges
    """
    narrow_peaks_dict = {}
    files = [f for f in os.listdir(peak_dir) if f.endswith(".narrowPeak")]
    
    log.info(f"📂 Found {len(files)} peak files. Loading and filtering...")
    
    min_score = -np.log10(q_value_threshold)
    
    for f in files:
        sample_name = f.replace("_peaks.narrowPeak", "")
        file_path = os.path.join(peak_dir, f)
        
        df = pd.read_csv(
            file_path, sep="\t", header=None, names=NARROWPEAK_COLS,
            comment='#', dtype={0: str}
        )
        
        # Normalize chromosome names
        if add_chr_prefix and not df["Chromosome"].iloc[0].startswith("chr"):
            df["Chromosome"] = "chr" + df["Chromosome"]
        
        # Filter by q-value
        df_filtered = df[df["QValue"] >= min_score].copy()
        
        n_peaks = len(df_filtered)
        if n_peaks < min_peaks_per_sample:
            log.warning(f"   ⚠️  Skipping {sample_name}: Too few peaks ({n_peaks})")
            continue
        
        log.info(f"   ✅ {sample_name}: {n_peaks} peaks loaded")
        narrow_peaks_dict[sample_name] = pr.PyRanges(df_filtered)
    
    if len(narrow_peaks_dict) == 0:
        raise ValueError("No samples passed filtering!")
    
    return narrow_peaks_dict


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def cpm(x: pr.PyRanges, column: str) -> pr.PyRanges:
    """Normalize score column to Counts Per Million."""
    total_score = x.df[column].sum()
    if total_score > 0:
        x.df[column] = (x.df[column] / total_score) * 1e6
    return x


def calculate_peaks_and_extend(
    narrow_peaks: pr.PyRanges,
    peak_half_width: int,
    chromsizes: Optional[pr.PyRanges] = None,
    path_to_blacklist: Optional[str] = None
) -> pr.PyRanges:
    """
    Center peaks on summit and extend by half_width.
    
    Parameters
    ----------
    narrow_peaks : PyRanges
        Input peaks with Summit column
    peak_half_width : int
        Half-width for extension
    chromsizes : PyRanges, optional
        Chromosome sizes for filtering
    path_to_blacklist : str, optional
        Path to blacklist BED file
    
    Returns
    -------
    PyRanges
        Extended and filtered peaks
    """
    # Ensure Summit column exists
    if "Summit" not in narrow_peaks.columns:
        narrow_peaks.Summit = ((narrow_peaks.End - narrow_peaks.Start) / 2).astype(int)

    starts = narrow_peaks.Start + narrow_peaks.Summit - peak_half_width
    ends = narrow_peaks.Start + narrow_peaks.Summit + peak_half_width 
    
    df = narrow_peaks.df.copy()
    df["Start"] = starts
    df["End"] = ends
    
    # Filter out negative starts
    df = df[df["Start"] >= 0]
    
    center_extended_peaks = pr.PyRanges(df)
    
    if chromsizes is not None:
        center_extended_peaks = center_extended_peaks.intersect(chromsizes, how="containment")
    
    if path_to_blacklist is not None:
        if isinstance(path_to_blacklist, str):
            blacklist = pr.read_bed(path_to_blacklist)
        center_extended_peaks = center_extended_peaks.overlap(blacklist, invert=True)
        
    return center_extended_peaks


def iterative_peak_filtering(peaks: pr.PyRanges) -> pr.PyRanges:
    """
    Resolve overlaps by greedily keeping the highest scoring peak.
    
    Parameters
    ----------
    peaks : PyRanges
        Input peaks with Score column
    
    Returns
    -------
    PyRanges
        Non-overlapping peaks
    """
    # Sanitize input - remove old cluster columns
    df_clean = peaks.df
    cols_to_drop = [c for c in ["Cluster", "Count"] if c in df_clean.columns]
    if cols_to_drop:
        df_clean = df_clean.drop(columns=cols_to_drop)
    peaks = pr.PyRanges(df_clean)

    # Cluster overlapping peaks
    peaks = peaks.sort()
    peaks_with_cluster = peaks.cluster(count=True) 
    
    df = peaks_with_cluster.df
    
    # Separate non-overlapping peaks
    keepers = [df[df["Count"] == 1].copy()] 
    
    # Process clusters with overlaps
    complex_clusters = df[df["Count"] > 1]
    
    if not complex_clusters.empty:
        candidates = pr.PyRanges(complex_clusters)
        
        while len(candidates) > 0:
            # For every Cluster, find the peak with the Max Score
            best_indices = candidates.df.groupby("Cluster", observed=True)["Score"].idxmax()
            
            # Extract the best peaks
            best_peaks_df = candidates.df.loc[best_indices].copy()
            best_peaks = pr.PyRanges(best_peaks_df)
            
            keepers.append(best_peaks_df)
            
            # Remove candidates that overlap with the winners
            candidates = candidates.overlap(best_peaks, invert=True)
            
            if len(candidates) > 0:
                # Sanitize before re-clustering
                cand_df = candidates.df
                cols_to_drop = [c for c in ["Cluster", "Count"] if c in cand_df.columns]
                if cols_to_drop:
                    cand_df = cand_df.drop(columns=cols_to_drop)
                
                candidates = pr.PyRanges(cand_df)
                candidates = candidates.cluster(count=True)

    # Concatenate and clean up
    final_df = pd.concat(keepers).sort_values(["Chromosome", "Start"])
    
    cols_to_drop = [c for c in ["Cluster", "Count"] if c in final_df.columns]
    if cols_to_drop:
        final_df = final_df.drop(columns=cols_to_drop)

    return pr.PyRanges(final_df)


def harmonize_chromosomes(
    narrow_peaks_dict: Dict[str, pr.PyRanges],
    chromsizes: Union[pr.PyRanges, pd.DataFrame]
) -> tuple:
    """
    Ensure consistent chromosome naming (UCSC format with 'chr' prefix).
    
    Parameters
    ----------
    narrow_peaks_dict : dict
        Dictionary of sample name -> PyRanges
    chromsizes : PyRanges or DataFrame
        Chromosome sizes
    
    Returns
    -------
    tuple
        (harmonized_peaks_dict, harmonized_chromsizes)
    """
    log.info("🔄 Checking and harmonizing chromosome names...")

    # Fix chromsizes
    if isinstance(chromsizes, pr.PyRanges):
        cs_df = chromsizes.df
    else:
        cs_df = chromsizes.copy()
    
    example_chrom = str(cs_df["Chromosome"].iloc[0])
    if not example_chrom.startswith("chr"):
        log.info(f"   🛠️  Fixing chromsizes: '{example_chrom}' -> 'chr{example_chrom}'")
        cs_df["Chromosome"] = "chr" + cs_df["Chromosome"].astype(str)
    
    chromsizes_fixed = pr.PyRanges(cs_df)

    # Fix peak dictionaries
    narrow_peaks_fixed = {}
    
    for sample, peaks in narrow_peaks_dict.items():
        df = peaks.df
        example_chrom = str(df["Chromosome"].iloc[0])
        
        if not example_chrom.startswith("chr"):
            if len(narrow_peaks_fixed) == 0:
                log.info(f"   🛠️  Fixing peaks: '{example_chrom}' -> 'chr{example_chrom}'")
            
            df["Chromosome"] = "chr" + df["Chromosome"].astype(str)
            narrow_peaks_fixed[sample] = pr.PyRanges(df)
        else:
            narrow_peaks_fixed[sample] = peaks
            
    log.info("✅ Harmonization complete. All data uses 'chr' prefix.")
    return narrow_peaks_fixed, chromsizes_fixed
