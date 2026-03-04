"""
Cross-Species Peak Analysis Pipeline
=====================================

Pipeline for comparing ATAC-seq peaks across species using UCSC liftOver:

1. Liftover non-human species peaks to human (hg38)
2. Merge ALL species (including human) into unified consensus, tracking species origin
   (summit-based clustering with fixed-width output)
3. *(skipped — human-specific classification via orthology in Step 7)*
4. Liftover unified consensus peaks back to each species
   - 4b. (optional) Filter liftback peaks by size
   - 4c. (optional) Reciprocal liftover check (round-trip NHP→hg38 overlap)
   - 4d. (optional) Per-peak liftover conservation scoring (match_ratio 0–1)
5. Identify species-specific peaks for each non-human species
6. Annotate all peaks with closest gene from GTF
7. Build master annotation — reclassify peaks using orthology
   (unified vs human_specific based on liftback success)
8. Cross-map species-specific peaks via direct inter-species chains

Author: J. Janssens
"""

import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple, Union
from collections import defaultdict

import pandas as pd
import numpy as np

from .liftover import liftover_peaks, liftover_two_step, get_chain_file, DEFAULT_CHAIN_DIR, compute_liftover_similarity

# Reverse chain files (hg38 -> species)
REVERSE_CHAIN_FILES = {
    "Gorilla": "hg38ToGorGor4.over.chain",
    "Chimpanzee": "hg38ToPanTro5.over.chain",
    "Bonobo": "hg38ToPanPan2.over.chain",
    "Macaque": "hg38ToRheMac10.over.chain",
    "Marmoset_step1": "hg38ToCalJac4.over.chain",
    "Marmoset_step2": "calJac4ToCalJac1.over.chain",
}

# Default GTF files for closest-gene annotation
# NOTE: Cellranger reference directories are named using Ensembl conventions
# (panPan1 = panpan1.1, panTro3 = Pan_tro_3.0) but the underlying assemblies
# are the same as UCSC panPan2 and panTro5 respectively (confirmed via NCBI
# accessions GCA_000258655.2 and GCA_000001515.5). No coordinate mismatch.
DEFAULT_GTF_FILES = {
    "Human": "/cluster/home/jjanssens/jjans/analysis/cerebellum/genomes_new/homo_sapiens/gencode.v48.basic.annotation.gtf.gz",
    "Bonobo": "/cluster/work/treutlein/jjans/data/intestine/nhp_atlas/genomes/reference_/panPan1/genes/genes.gtf",       # Ensembl panpan1.1 = UCSC panPan2
    "Gorilla": "/cluster/work/treutlein/jjans/data/intestine/nhp_atlas/genomes/reference_/gorGor4/genes/genes.gtf",
    "Macaque": "/cluster/work/treutlein/jjans/data/intestine/nhp_atlas/genomes/reference_/Mmul10/genes/genes.gtf",
    "Marmoset": "/cluster/work/treutlein/jjans/data/intestine/nhp_atlas/genomes/reference_/calJac1_mito/genes/genes.gtf",
    "Chimpanzee": "/cluster/work/treutlein/jjans/data/intestine/nhp_atlas/genomes/reference_/panTro3/genes/genes.gtf.gz", # Ensembl Pan_tro_3.0 = UCSC panTro5
}


def get_reverse_chain_file(species: str, chain_dir: str = DEFAULT_CHAIN_DIR) -> str:
    """Get path to reverse chain file (hg38 -> species)."""
    if species not in REVERSE_CHAIN_FILES:
        raise ValueError(f"Unknown species: {species}. Available: {list(REVERSE_CHAIN_FILES.keys())}")
    return os.path.join(chain_dir, REVERSE_CHAIN_FILES[species])


def compute_species_similarity(
    input_bed: str,
    species_list: Optional[List[str]] = None,
    chain_dir: str = DEFAULT_CHAIN_DIR,
    liftover_path: str = "liftOver",
    output_tsv: Optional[str] = None,
    verbose: bool = True,
    ncpu: int = 1,
) -> pd.DataFrame:
    """
    Compute per-peak liftover similarity for multiple species.

    For each species, performs a round-trip liftover
    (hg38 → species → hg38) and measures how many bases survive.
    Returns a combined DataFrame with a ``species`` column.

    Parameters
    ----------
    input_bed : str
        BED file in hg38 coordinates (e.g. unified consensus).
    species_list : list of str, optional
        Species to compute similarity for.
        Default: ``["Bonobo", "Chimpanzee", "Gorilla", "Macaque", "Marmoset"]``
    chain_dir : str
        Directory containing chain files.
    liftover_path : str
        Path to UCSC ``liftOver`` executable.
    output_tsv : str, optional
        If provided, save combined results to this TSV file.
    verbose : bool
        Print progress.
    ncpu : int
        Workers for parallel liftover.

    Returns
    -------
    pd.DataFrame
        Combined DataFrame with columns from ``compute_liftover_similarity``
        plus a ``species`` column.
    """
    if species_list is None:
        species_list = ["Bonobo", "Chimpanzee", "Gorilla", "Macaque", "Marmoset"]

    all_dfs = []
    for species in species_list:
        if verbose:
            print(f"\n{'='*60}")
            print(f"  {species}")
            print(f"{'='*60}")

        # Resolve chain files
        if species == "Marmoset":
            chain_fwd = get_chain_file("Marmoset_step1", chain_dir)
            chain_fwd_2 = get_chain_file("Marmoset_step2", chain_dir)
            chain_rev = get_reverse_chain_file("Marmoset_step1", chain_dir)
            chain_rev_2 = get_reverse_chain_file("Marmoset_step2", chain_dir)
            # Forward: hg38 → calJac4 → calJac1  (reverse of species→hg38)
            # So forward = reverse chains, reverse = forward chains
            df = compute_liftover_similarity(
                input_bed=input_bed,
                chain_forward=chain_rev,      # hg38 → calJac4
                chain_reverse=chain_fwd,      # calJac4 → hg38
                chain_forward_2=chain_rev_2,  # calJac4 → calJac1
                chain_reverse_2=chain_fwd_2,  # calJac1 → calJac4
                liftover_path=liftover_path,
                auto_chr=True, verbose=verbose, ncpu=ncpu,
            )
        else:
            chain_rev = get_reverse_chain_file(species, chain_dir)
            chain_fwd = get_chain_file(species, chain_dir)
            df = compute_liftover_similarity(
                input_bed=input_bed,
                chain_forward=chain_rev,   # hg38 → species
                chain_reverse=chain_fwd,   # species → hg38
                liftover_path=liftover_path,
                auto_chr=True, verbose=verbose, ncpu=ncpu,
            )

        df['species'] = species
        all_dfs.append(df)

    combined = pd.concat(all_dfs, ignore_index=True)

    if output_tsv:
        os.makedirs(os.path.dirname(output_tsv) or '.', exist_ok=True)
        combined.to_csv(output_tsv, sep='\t', index=False)
        if verbose:
            print(f"\n💾 Saved combined similarity to: {output_tsv}")

    return combined


# =============================================================================
# MERGE WITH SPECIES TRACKING
# =============================================================================

def merge_with_species_tracking(
    species_beds: Dict[str, str],
    output_bed: str,
    merge_distance: int = 0,
    verbose: bool = True,
) -> Dict[str, Any]:
    """
    Merge BED files from multiple species, tracking which species contributed
    to each merged peak.

    Uses bedtools merge -c 4 -o distinct to collapse species tags.

    Output BED format: chr, start, end, species_list (comma-separated)

    Args:
        species_beds: Dict mapping species name to BED file path
        output_bed: Path to output merged BED file
        merge_distance: Maximum distance between features to merge
        verbose: Print progress

    Returns:
        Dictionary with merge statistics
    """
    if verbose:
        print(f"Merging {len(species_beds)} species BED files with species tracking...")

    os.makedirs(os.path.dirname(output_bed), exist_ok=True)

    # Count input peaks per species
    total_input = 0
    species_counts = {}
    for species, bed in species_beds.items():
        with open(bed) as f:
            count = sum(1 for line in f if line.strip() and not line.startswith('#'))
        species_counts[species] = count
        total_input += count
        if verbose:
            print(f"   {species}: {count:,} peaks")

    if verbose:
        print(f"   Total input peaks: {total_input:,}")

    with tempfile.TemporaryDirectory() as tmpdir:
        # Create tagged BED: chr, start, end, species_name
        tagged_file = os.path.join(tmpdir, "tagged.bed")
        with open(tagged_file, 'w') as fout:
            for species, bed in species_beds.items():
                with open(bed) as fin:
                    for line in fin:
                        if line.strip() and not line.startswith('#'):
                            parts = line.strip().split('\t')
                            if len(parts) >= 3:
                                fout.write(f"{parts[0]}\t{parts[1]}\t{parts[2]}\t{species}\n")

        # Sort
        sorted_file = os.path.join(tmpdir, "sorted.bed")
        sort_cmd = f"sort -k1,1 -k2,2n {tagged_file} > {sorted_file}"
        subprocess.run(sort_cmd, shell=True, check=True)

        # Merge with species tracking
        merge_cmd = (
            f"bedtools merge -i {sorted_file} -d {merge_distance} "
            f"-c 4 -o distinct > {output_bed}"
        )
        result = subprocess.run(merge_cmd, shell=True, capture_output=True, text=True)

        if result.returncode != 0:
            if verbose:
                print("   bedtools merge failed, using Python fallback...")
            _python_merge_with_tracking(sorted_file, output_bed, merge_distance)

    # Count output and compute species distribution
    output_count = 0
    species_dist = defaultdict(int)
    with open(output_bed) as f:
        for line in f:
            if line.strip():
                output_count += 1
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    for sp in parts[3].split(','):
                        species_dist[sp.strip()] += 1

    if verbose:
        print(f"   Merged peaks: {output_count:,}")
        print(f"   Reduction: {(1 - output_count / total_input) * 100:.1f}%")
        print(f"   Species distribution in merged peaks:")
        for sp in sorted(species_dist.keys()):
            print(f"      {sp}: present in {species_dist[sp]:,} peaks")

    return {
        "status": "success",
        "input_peaks": total_input,
        "merged_peaks": output_count,
        "species_counts": dict(species_counts),
        "species_distribution": dict(species_dist),
        "output_file": output_bed,
    }


def _python_merge_with_tracking(input_bed: str, output_bed: str, merge_distance: int = 0):
    """Python fallback for merge with species tracking."""
    regions = defaultdict(list)
    with open(input_bed) as f:
        for line in f:
            if line.strip():
                parts = line.strip().split('\t')
                chrom = parts[0]
                start, end = int(parts[1]), int(parts[2])
                species = parts[3] if len(parts) >= 4 else "unknown"
                regions[chrom].append((start, end, species))

    with open(output_bed, 'w') as fout:
        for chrom in sorted(regions.keys()):
            intervals = sorted(regions[chrom])
            if not intervals:
                continue
            merged_start, merged_end = intervals[0][0], intervals[0][1]
            merged_species = {intervals[0][2]}
            for start, end, species in intervals[1:]:
                if start <= merged_end + merge_distance:
                    merged_end = max(merged_end, end)
                    merged_species.add(species)
                else:
                    fout.write(f"{chrom}\t{merged_start}\t{merged_end}\t{','.join(sorted(merged_species))}\n")
                    merged_start, merged_end = start, end
                    merged_species = {species}
            fout.write(f"{chrom}\t{merged_start}\t{merged_end}\t{','.join(sorted(merged_species))}\n")


# =============================================================================
# SUMMIT-BASED CONSENSUS MERGING (v2)
# =============================================================================

def _normalize_scores_by_species(
    summits: List[Dict[str, Any]],
    verbose: bool = False,
) -> None:
    """
    Normalize peak scores within each species to [0, 1] using percentile rank.

    This ensures species with different sequencing depths or cell numbers
    contribute comparable weights to the score-weighted consensus.
    Modifies summits in-place by adding a 'norm_score' key.

    Args:
        summits: List of summit dicts, each with 'species' and 'score' keys
        verbose: Print per-species score ranges
    """
    by_species = defaultdict(list)
    for s in summits:
        by_species[s['species']].append(s)

    for species, peaks in by_species.items():
        scores = [p['score'] for p in peaks]
        n = len(scores)

        if n <= 1:
            for p in peaks:
                p['norm_score'] = 1.0
            continue

        # Rank-based normalization: ties get average rank
        sorted_indices = sorted(range(n), key=lambda i: scores[i])
        ranks = [0.0] * n
        i = 0
        while i < n:
            # Find ties
            j = i
            while j < n and scores[sorted_indices[j]] == scores[sorted_indices[i]]:
                j += 1
            avg_rank = (i + j - 1) / 2.0
            for k in range(i, j):
                ranks[sorted_indices[k]] = avg_rank
            i = j

        # Scale to [0, 1]
        max_rank = n - 1
        for idx, p in enumerate(peaks):
            p['norm_score'] = ranks[idx] / max_rank if max_rank > 0 else 0.5

        if verbose:
            raw = [p['score'] for p in peaks]
            norm = [p['norm_score'] for p in peaks]
            print(f"      {species}: raw [{min(raw):.0f}, {max(raw):.0f}] "
                  f"-> norm [{min(norm):.3f}, {max(norm):.3f}]")


def summit_based_merge(
    species_beds: Dict[str, str],
    output_bed: str,
    summit_window: int = 500,
    cluster_distance: int = 250,
    max_peak_size: int = 2000,
    min_peak_size: int = 100,
    score_col: int = 4,
    normalize_scores: bool = True,
    verbose: bool = True,
) -> Dict[str, Any]:
    """
    Summit-based cross-species consensus peak merging.

    Instead of merging overlapping intervals (which creates enormous peaks
    when lifted peaks from different species overlap due to structural
    rearrangements), this approach:

    1. Filters out peaks distorted by liftover (too large or too small)
    2. Extracts 1bp summit (center) from each peak
    3. Clusters nearby summits using ``bedtools cluster -d <distance>``
    4. Computes score-weighted average position within each cluster
    5. Expands to fixed-width windows (``summit_window`` bp)
    6. Tracks species origin per consensus peak

    Output BED format: chr, start, end, species_list (comma-separated)
    (same as merge_with_species_tracking for downstream compatibility)

    Args:
        species_beds: Dict mapping species name to BED file path (hg38 coords).
            Expects narrowPeak or at least BED5 format with a score column.
        output_bed: Path to output merged BED file
        summit_window: Width of output consensus peaks (default: 500bp)
        cluster_distance: Max distance between summits to cluster (default: 250bp).
            Two peaks whose centers are within this distance get merged.
        max_peak_size: Maximum peak size to include (filters liftover artifacts)
        min_peak_size: Minimum peak size to include
        score_col: 0-indexed column number for peak scores (default: 4, BED score).
            Use 6 for signalValue (fold enrichment) in narrowPeak format.
        normalize_scores: Whether to normalize scores across species using
            rank-based normalization (recommended when species have different
            sequencing depths or cell counts)
        verbose: Print progress

    Returns:
        Dictionary with merge statistics (compatible with existing pipeline)
    """
    if verbose:
        print(f"Summit-based merge of {len(species_beds)} species")
        print(f"   Parameters: window={summit_window}bp, cluster_d={cluster_distance}bp, "
              f"size_filter=[{min_peak_size}, {max_peak_size}]bp")

    os.makedirs(os.path.dirname(output_bed), exist_ok=True)

    # ── Step 1: Read peaks, filter by size, extract summits ──────────────
    summits = []
    species_input = {}
    species_filtered = {}

    for species, bed_path in species_beds.items():
        n_input = 0
        n_size_filtered = 0

        with open(bed_path) as f:
            for line in f:
                if not line.strip() or line.startswith('#'):
                    continue
                n_input += 1
                parts = line.strip().split('\t')
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                peak_size = end - start

                # Size filter: remove liftover-distorted peaks
                if peak_size < min_peak_size or peak_size > max_peak_size:
                    n_size_filtered += 1
                    continue

                # Extract score
                score = float(parts[score_col]) if len(parts) > score_col else 1.0
                if score <= 0:
                    score = 1.0

                # Summit = center of peak
                summit_pos = (start + end) // 2

                summits.append({
                    'chrom': chrom,
                    'summit': summit_pos,
                    'score': score,
                    'species': species,
                })

        species_input[species] = n_input
        species_filtered[species] = n_size_filtered

        if verbose:
            kept = n_input - n_size_filtered
            print(f"   {species}: {n_input:,} input, "
                  f"{n_size_filtered:,} size-filtered, {kept:,} kept")

    total_input = sum(species_input.values())
    total_kept = len(summits)

    if not summits:
        with open(output_bed, 'w') as f:
            pass
        return {
            "status": "error",
            "message": "No summits after filtering",
            "input_peaks": total_input,
            "merged_peaks": 0,
        }

    if verbose:
        print(f"   Total: {total_input:,} input -> {total_kept:,} summits after size filter")

    # ── Step 2: Normalize scores ─────────────────────────────────────────
    if normalize_scores:
        _normalize_scores_by_species(summits, verbose=verbose)
        if verbose:
            print(f"   Scores normalized (rank-based) across {len(species_beds)} species")
    else:
        for s in summits:
            s['norm_score'] = s['score']

    # ── Step 3: Cluster summits with bedtools ────────────────────────────
    with tempfile.TemporaryDirectory() as tmpdir:
        summits_bed = os.path.join(tmpdir, "summits.bed")
        sorted_bed = os.path.join(tmpdir, "summits_sorted.bed")
        clustered_bed = os.path.join(tmpdir, "summits_clustered.bed")

        # Write: chr, summit, summit+1, index, norm_score, species
        with open(summits_bed, 'w') as f:
            for i, s in enumerate(summits):
                f.write(f"{s['chrom']}\t{s['summit']}\t{s['summit'] + 1}\t"
                        f"{i}\t{s['norm_score']:.6f}\t{s['species']}\n")

        # Sort
        subprocess.run(
            f"sort -k1,1 -k2,2n {summits_bed} > {sorted_bed}",
            shell=True, check=True,
        )

        # Cluster: bedtools cluster adds a cluster ID as the last column
        subprocess.run(
            f"bedtools cluster -i {sorted_bed} -d {cluster_distance} > {clustered_bed}",
            shell=True, check=True,
        )

        # ── Step 4: Parse clusters, compute weighted centers ─────────────
        clusters = defaultdict(list)
        with open(clustered_bed) as f:
            for line in f:
                parts = line.strip().split('\t')
                chrom = parts[0]
                pos = int(parts[1])
                score = float(parts[4])
                species = parts[5]
                cluster_id = parts[6]  # bedtools cluster appends this

                clusters[(chrom, cluster_id)].append({
                    'pos': pos,
                    'score': score,
                    'species': species,
                })

    # ── Step 5: Score-weighted center → fixed-width peaks ────────────────
    half_window = summit_window // 2
    consensus_peaks = []

    for (chrom, _cid), members in clusters.items():
        total_weight = sum(m['score'] for m in members)

        if total_weight > 0:
            center = int(
                round(sum(m['pos'] * m['score'] for m in members) / total_weight)
            )
        else:
            # Fallback: unweighted mean
            center = int(round(sum(m['pos'] for m in members) / len(members)))

        species_set = sorted(set(m['species'] for m in members))

        start = max(0, center - half_window)
        end = center + half_window

        consensus_peaks.append((chrom, start, end, ','.join(species_set)))

    # Sort by genomic position
    consensus_peaks.sort(key=lambda x: (x[0], x[1]))

    # ── Step 6: Write output ─────────────────────────────────────────────
    with open(output_bed, 'w') as f:
        for chrom, start, end, species_list in consensus_peaks:
            f.write(f"{chrom}\t{start}\t{end}\t{species_list}\n")

    # ── Statistics ───────────────────────────────────────────────────────
    n_output = len(consensus_peaks)
    sizes = [p[2] - p[1] for p in consensus_peaks] if consensus_peaks else [0]
    species_dist = defaultdict(int)
    for _, _, _, sp_list in consensus_peaks:
        for sp in sp_list.split(','):
            species_dist[sp.strip()] += 1

    if verbose:
        print(f"\n   Clustered {total_kept:,} summits into {n_output:,} consensus peaks")
        print(f"   Output peak sizes: min={min(sizes)}, median={sorted(sizes)[len(sizes)//2]}, "
              f"max={max(sizes)}")
        reduction = (1 - n_output / total_input) * 100 if total_input > 0 else 0
        print(f"   Reduction: {reduction:.1f}%")
        print(f"   Species distribution in consensus peaks:")
        for sp in sorted(species_dist.keys()):
            print(f"      {sp}: present in {species_dist[sp]:,} peaks")

    return {
        "status": "success",
        "input_peaks": total_input,
        "merged_peaks": n_output,
        "species_counts": dict(species_input),
        "species_distribution": dict(species_dist),
        "output_file": output_bed,
        "summit_stats": {
            "total_summits": total_kept,
            "size_filtered": sum(species_filtered.values()),
            "clusters": n_output,
            "peak_size_range": (min(sizes), max(sizes)),
        },
    }


# =============================================================================
# SUMMIT-BASED LIFTOVER WITH CONSERVATION SCORING (v2)
# =============================================================================

def liftover_summits_with_conservation(
    input_bed: str,
    output_bed: str,
    chain_file: str,
    liftover_path: str = "liftOver",
    min_match: float = 0.1,
    summit_window: int = 500,
    min_conservation: float = 0.0,
    auto_chr: bool = True,
    verbose: bool = True,
) -> Dict[str, Any]:
    """
    Lift peaks by their 1bp summit, then re-expand to fixed-width windows.

    This avoids the classic liftover distortion problem where a 500bp peak
    in species A can map to a 50bp or 50,000bp region in species B due to
    indels, inversions, or assembly gaps. Instead:

    1. Extract 1bp summit (center) from each input peak
    2. LiftOver the 1bp coordinate with a *relaxed* ``min_match`` (since
       a single base either maps or doesn't — the stringent ``min_match``
       is not meaningful for 1bp intervals)
    3. Re-expand each lifted summit to a fixed ``summit_window`` bp window
    4. Optionally compute a conservation score: liftover the full peak
       interval and compute ``lifted_size / original_size``. Peaks below
       ``min_conservation`` are discarded.

    Use this as an alternative liftover strategy when you want all output
    peaks to be a uniform size, regardless of structural variation.

    Args:
        input_bed: Input BED file (≥BED3; column 4 = name if present)
        output_bed: Output BED file with fixed-width lifted peaks
        chain_file: Path to liftOver chain file (source → target)
        liftover_path: Path to liftOver binary
        min_match: Min ratio for the full-interval liftover used to compute
            conservation score. The summit liftover always uses 0.1.
        summit_window: Width of output peaks (default: 500bp)
        min_conservation: Minimum conservation score (lifted_size / original_size)
            to keep a peak. Set to 0.0 to keep all successfully lifted summits.
        auto_chr: Auto-detect and harmonize chr prefix
        verbose: Print progress

    Returns:
        Dict with statistics: n_input, n_lifted, n_conservation_filtered, etc.
    """
    if verbose:
        print(f"  Summit-based liftover (window={summit_window}bp, "
              f"min_conservation={min_conservation})")

    os.makedirs(os.path.dirname(output_bed), exist_ok=True)
    half = summit_window // 2

    with tempfile.TemporaryDirectory() as tmpdir:
        # --- Step 1: Read input, extract summits ---
        peaks = []  # list of (chrom, start, end, name)
        summit_bed = os.path.join(tmpdir, "summits.bed")

        with open(input_bed) as fin, open(summit_bed, 'w') as fout:
            for i, line in enumerate(fin):
                if not line.strip() or line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                name = parts[3] if len(parts) > 3 else f"peak_{i}"
                peaks.append((chrom, start, end, name))

                # 1bp summit
                center = (start + end) // 2
                fout.write(f"{chrom}\t{center}\t{center + 1}\t{name}\n")

        n_input = len(peaks)
        if n_input == 0:
            with open(output_bed, 'w') as f:
                pass
            return {"status": "success", "n_input": 0, "n_lifted": 0,
                    "n_conservation_filtered": 0, "output_file": output_bed}

        # --- Step 2: Liftover 1bp summits (relaxed minMatch) ---
        summit_lifted = os.path.join(tmpdir, "summits_lifted.bed")
        summit_unmapped = os.path.join(tmpdir, "summits_unmapped.bed")

        cmd_summit = (
            f"{liftover_path} {summit_bed} {chain_file} "
            f"{summit_lifted} {summit_unmapped} -minMatch=0.1"
        )
        subprocess.run(cmd_summit, shell=True, check=True,
                        capture_output=True, text=True)

        # Parse lifted summits: name -> (target_chrom, target_pos)
        lifted_summits = {}
        with open(summit_lifted) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    lifted_summits[parts[3]] = (parts[0], int(parts[1]))

        n_summit_lifted = len(lifted_summits)

        # --- Step 3: Optionally compute conservation score ---
        conservation = {}  # name -> score
        n_conservation_filtered = 0

        if min_conservation > 0:
            # Liftover full intervals to measure size distortion
            full_bed = os.path.join(tmpdir, "full_intervals.bed")
            full_lifted = os.path.join(tmpdir, "full_lifted.bed")
            full_unmapped = os.path.join(tmpdir, "full_unmapped.bed")

            with open(full_bed, 'w') as f:
                for chrom, start, end, name in peaks:
                    f.write(f"{chrom}\t{start}\t{end}\t{name}\n")

            cmd_full = (
                f"{liftover_path} {full_bed} {chain_file} "
                f"{full_lifted} {full_unmapped} -minMatch={min_match}"
            )
            subprocess.run(cmd_full, shell=True, check=True,
                            capture_output=True, text=True)

            # Build map: peak_name -> original_size
            original_sizes = {name: end - start for chrom, start, end, name in peaks}

            # Parse lifted full intervals
            with open(full_lifted) as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 4:
                        name = parts[3]
                        lifted_size = int(parts[2]) - int(parts[1])
                        orig_size = original_sizes.get(name, summit_window)
                        # Conservation = how much of the interval survived
                        conservation[name] = min(lifted_size / orig_size, 2.0) if orig_size > 0 else 0.0

        # --- Step 4: Build output with fixed-width windows ---
        results = []
        for name, (target_chrom, target_pos) in lifted_summits.items():
            # Conservation filter
            if min_conservation > 0:
                score = conservation.get(name, 0.0)
                if score < min_conservation:
                    n_conservation_filtered += 1
                    continue

            new_start = max(0, target_pos - half)
            new_end = target_pos + half
            results.append((target_chrom, new_start, new_end, name))

        # Sort by genomic position
        results.sort(key=lambda x: (x[0], x[1]))

        with open(output_bed, 'w') as f:
            for chrom, start, end, name in results:
                f.write(f"{chrom}\t{start}\t{end}\t{name}\n")

    n_output = len(results)
    n_unmapped = n_input - n_summit_lifted

    if verbose:
        print(f"    Input: {n_input:,}")
        print(f"    Summit lifted: {n_summit_lifted:,}, unmapped: {n_unmapped:,}")
        if min_conservation > 0:
            print(f"    Conservation filtered: {n_conservation_filtered:,}")
        print(f"    Output: {n_output:,} peaks ({summit_window}bp each)")

    return {
        "status": "success",
        "n_input": n_input,
        "n_summit_lifted": n_summit_lifted,
        "n_unmapped": n_unmapped,
        "n_conservation_filtered": n_conservation_filtered,
        "n_output": n_output,
        "output_file": output_bed,
    }


# =============================================================================
# HUMAN-SPECIFIC PEAK DETECTION
# =============================================================================

def find_human_specific_peaks(
    human_bed: str,
    nonhuman_lifted_beds: List[str],
    output_bed: str,
    merge_distance: int = 0,
    verbose: bool = True,
) -> Dict[str, Any]:
    """
    .. deprecated::
        This function uses **detection-based** overlap to find human peaks
        not overlapping any NHP called peak.  This conflates "no peak called"
        with "no orthologous sequence".  Use the orthology-based
        classification in ``build_master_annotation()`` instead, which
        identifies human-specific peaks as those where liftback to every
        NHP species failed (``n_species_orth == 1``).

    Find human peaks that do not overlap any non-human lifted peak.

    **Important**: These are peaks where no other species had a *called peak*
    in the orthologous region — i.e., the regulatory element appears to be
    active only in human. This does NOT mean the underlying sequence cannot
    be lifted to other genomes (many of these regions will successfully
    liftOver). For peaks whose *sequence* is unique to human (no orthologous
    region in any other species), check liftback failure status in the
    master annotation (``{Species}_orth == 0`` for all NHP species).

    Args:
        human_bed: Path to human consensus peaks (hg38)
        nonhuman_lifted_beds: List of non-human species peaks lifted to hg38
        output_bed: Output BED file for human-specific peaks
        merge_distance: Slop for overlap (0 = strict overlap)
        verbose: Print progress

    Returns:
        Dictionary with statistics
    """
    import warnings
    warnings.warn(
        "find_human_specific_peaks() is deprecated. Use orthology-based "
        "classification in build_master_annotation() instead.",
        DeprecationWarning, stacklevel=2,
    )

    if verbose:
        print("Finding human-specific peaks...")

    os.makedirs(os.path.dirname(output_bed), exist_ok=True)

    with tempfile.TemporaryDirectory() as tmpdir:
        # Merge all non-human lifted peaks into one BED
        nonhuman_merged = os.path.join(tmpdir, "nonhuman_merged.bed")

        # Concatenate non-human BEDs (BED3 only)
        cat_file = os.path.join(tmpdir, "nonhuman_cat.bed")
        with open(cat_file, 'w') as fout:
            for bed in nonhuman_lifted_beds:
                if os.path.exists(bed):
                    with open(bed) as fin:
                        for line in fin:
                            if line.strip() and not line.startswith('#'):
                                parts = line.strip().split('\t')
                                if len(parts) >= 3:
                                    fout.write(f"{parts[0]}\t{parts[1]}\t{parts[2]}\n")

        # Sort and merge non-human peaks
        sorted_file = os.path.join(tmpdir, "nonhuman_sorted.bed")
        subprocess.run(f"sort -k1,1 -k2,2n {cat_file} > {sorted_file}", shell=True, check=True)
        subprocess.run(f"bedtools merge -i {sorted_file} > {nonhuman_merged}", shell=True, check=True)

        # Extract BED3 from human peaks (may have extra columns)
        human_bed3 = os.path.join(tmpdir, "human_bed3.bed")
        with open(human_bed) as fin, open(human_bed3, 'w') as fout:
            for line in fin:
                if line.strip() and not line.startswith('#'):
                    parts = line.strip().split('\t')
                    fout.write(f"{parts[0]}\t{parts[1]}\t{parts[2]}\n")

        # Sort human peaks
        human_sorted = os.path.join(tmpdir, "human_sorted.bed")
        subprocess.run(f"sort -k1,1 -k2,2n {human_bed3} > {human_sorted}", shell=True, check=True)

        # Find human peaks NOT overlapping any non-human peak
        subtract_cmd = f"bedtools intersect -v -a {human_sorted} -b {nonhuman_merged} > {output_bed}"
        subprocess.run(subtract_cmd, shell=True, check=True)

    # Count results
    with open(human_bed) as f:
        total_human = sum(1 for l in f if l.strip() and not l.startswith('#'))
    with open(output_bed) as f:
        n_specific = sum(1 for l in f if l.strip())

    if verbose:
        print(f"   Total human peaks: {total_human:,}")
        print(f"   Human-specific peaks: {n_specific:,} ({n_specific / total_human * 100:.1f}%)")

    return {
        "status": "success",
        "total_human": total_human,
        "human_specific": n_specific,
        "pct_specific": n_specific / total_human * 100 if total_human > 0 else 0,
        "output_file": output_bed,
    }


# =============================================================================
# SPECIES-SPECIFIC PEAK DETECTION
# =============================================================================

def find_species_specific_peaks(
    original_bed: str,
    liftback_bed: str,
    output_bed: str,
    species: str,
    verbose: bool = True,
) -> Dict[str, Any]:
    """
    Find species-specific peaks: original species peaks that do NOT overlap
    any lifted-back unified peak.

    These represent peaks unique to that species that could not be
    lifted over to human and back.

    Args:
        original_bed: Original species consensus peaks (in species coordinates)
        liftback_bed: Unified peaks lifted back to species coordinates
        output_bed: Output BED file for species-specific peaks
        species: Species name (for peak ID prefix)
        verbose: Print progress

    Returns:
        Dictionary with statistics
    """
    if verbose:
        print(f"   Finding {species}-specific peaks...")

    os.makedirs(os.path.dirname(output_bed), exist_ok=True)

    with tempfile.TemporaryDirectory() as tmpdir:
        # Extract BED3 from original peaks
        orig_bed3 = os.path.join(tmpdir, "original_bed3.bed")
        with open(original_bed) as fin, open(orig_bed3, 'w') as fout:
            for line in fin:
                if line.strip() and not line.startswith('#'):
                    parts = line.strip().split('\t')
                    fout.write(f"{parts[0]}\t{parts[1]}\t{parts[2]}\n")

        orig_sorted = os.path.join(tmpdir, "original_sorted.bed")
        subprocess.run(f"sort -k1,1 -k2,2n {orig_bed3} > {orig_sorted}", shell=True, check=True)

        # Extract BED3 from liftback peaks and merge
        liftback_bed3 = os.path.join(tmpdir, "liftback_bed3.bed")
        with open(liftback_bed) as fin, open(liftback_bed3, 'w') as fout:
            for line in fin:
                if line.strip() and not line.startswith('#'):
                    parts = line.strip().split('\t')
                    fout.write(f"{parts[0]}\t{parts[1]}\t{parts[2]}\n")

        # Detect and fix chr prefix mismatch between original and liftback
        def _first_chrom(path: str) -> str:
            with open(path) as fh:
                for line in fh:
                    if line.strip() and not line.startswith('#'):
                        return line.split('\t')[0]
            return ""

        orig_chrom = _first_chrom(orig_bed3)
        liftback_chrom = _first_chrom(liftback_bed3)
        orig_has_chr = orig_chrom.startswith("chr")
        liftback_has_chr = liftback_chrom.startswith("chr")

        if orig_has_chr != liftback_has_chr:
            # Harmonize liftback to match original
            liftback_fixed = os.path.join(tmpdir, "liftback_fixed.bed")
            with open(liftback_bed3) as fin, open(liftback_fixed, 'w') as fout:
                for line in fin:
                    if line.strip():
                        if orig_has_chr and not liftback_has_chr:
                            fout.write("chr" + line)
                        else:
                            fout.write(line[3:] if line.startswith("chr") else line)
            liftback_bed3 = liftback_fixed
            if verbose:
                print(f"      (harmonized chr prefix: liftback {'added' if orig_has_chr else 'stripped'} chr)")

        liftback_sorted = os.path.join(tmpdir, "liftback_sorted.bed")
        liftback_merged = os.path.join(tmpdir, "liftback_merged.bed")
        subprocess.run(f"sort -k1,1 -k2,2n {liftback_bed3} > {liftback_sorted}", shell=True, check=True)
        subprocess.run(f"bedtools merge -i {liftback_sorted} > {liftback_merged}", shell=True, check=True)

        # Species-specific = original NOT overlapping liftback
        specific_raw = os.path.join(tmpdir, "specific_raw.bed")
        subprocess.run(
            f"bedtools intersect -v -a {orig_sorted} -b {liftback_merged} > {specific_raw}",
            shell=True, check=True,
        )

        # Add species-specific peak IDs
        sp_lower = species.lower()
        with open(specific_raw) as fin, open(output_bed, 'w') as fout:
            for i, line in enumerate(fin, 1):
                parts = line.strip().split('\t')
                peak_id = f"{sp_lower}_peak_{i:06d}"
                fout.write(f"{parts[0]}\t{parts[1]}\t{parts[2]}\t{peak_id}\n")

    # Count
    with open(original_bed) as f:
        total_original = sum(1 for l in f if l.strip() and not l.startswith('#'))
    with open(output_bed) as f:
        n_specific = sum(1 for l in f if l.strip())

    if verbose:
        print(f"      Original peaks: {total_original:,}")
        print(f"      {species}-specific: {n_specific:,} ({n_specific / total_original * 100:.1f}%)")

    return {
        "status": "success",
        "species": species,
        "total_original": total_original,
        "species_specific": n_specific,
        "pct_specific": n_specific / total_original * 100 if total_original > 0 else 0,
        "output_file": output_bed,
    }


# =============================================================================
# GENE ANNOTATION
# =============================================================================

def _parse_gene_name_from_gtf_attrs(attrs: str) -> str:
    """Extract gene_name from GTF attribute string."""
    match = re.search(r'gene_name\s+"([^"]+)"', attrs)
    if match:
        return match.group(1)
    match = re.search(r"gene_name\s+'([^']+)'", attrs)
    if match:
        return match.group(1)
    # Fall back to gene_id
    match = re.search(r'gene_id\s+"([^"]+)"', attrs)
    if match:
        return match.group(1)
    return "unknown"


def extract_gene_bed_from_gtf(
    gtf_path: str,
    output_bed: str,
    feature_type: str = "gene",
    verbose: bool = True,
) -> str:
    """
    Extract gene TSSs from a GTF file into a sorted BED6 file
    for use with bedtools closest.

    Output BED: chr, start, end, gene_name, 0, strand

    Args:
        gtf_path: Path to GTF file (plain or .gz)
        output_bed: Path to output BED file
        feature_type: Feature type to extract (default: "gene")
        verbose: Print progress

    Returns:
        Path to output BED file
    """
    os.makedirs(os.path.dirname(output_bed) or '.', exist_ok=True)

    if gtf_path.endswith('.gz'):
        import gzip
        opener = gzip.open
    else:
        opener = open

    genes = []
    with opener(gtf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            if parts[2] != feature_type:
                continue

            chrom = parts[0]
            start = int(parts[3]) - 1  # GTF 1-based -> BED 0-based
            end = int(parts[4])
            strand = parts[6]
            gene_name = _parse_gene_name_from_gtf_attrs(parts[8])

            # TSS position
            if strand == '+':
                tss = start
            else:
                tss = end - 1

            genes.append((chrom, tss, tss + 1, gene_name, 0, strand))

    genes.sort(key=lambda x: (x[0], x[1]))

    with open(output_bed, 'w') as fout:
        for chrom, start, end, name, score, strand in genes:
            fout.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n")

    if verbose:
        print(f"   Extracted {len(genes):,} gene TSSs from {os.path.basename(gtf_path)}")

    return output_bed


def annotate_with_closest_gene(
    peaks_bed: str,
    gene_bed: str,
    output_file: str,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Annotate peaks with the closest gene using bedtools closest.

    Args:
        peaks_bed: BED file with peaks (at least chr, start, end, peak_id)
        gene_bed: BED file with gene TSSs (from extract_gene_bed_from_gtf)
        output_file: Path for output TSV annotation file
        verbose: Print progress

    Returns:
        DataFrame with peak_id, closest_gene, distance_to_gene
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        # Detect chr prefix in peaks vs gene BED and harmonize
        def _first_chrom(path: str) -> str:
            with open(path) as fh:
                for line in fh:
                    if line.strip() and not line.startswith('#'):
                        return line.split('\t')[0]
            return ""

        peak_chrom = _first_chrom(peaks_bed)
        gene_chrom = _first_chrom(gene_bed)
        peaks_have_chr = peak_chrom.startswith("chr")
        genes_have_chr = gene_chrom.startswith("chr")

        # If mismatch, create a harmonized gene BED with matching chr convention
        if peaks_have_chr and not genes_have_chr:
            gene_bed_fixed = os.path.join(tmpdir, "genes_chr.bed")
            with open(gene_bed) as fin, open(gene_bed_fixed, 'w') as fout:
                for line in fin:
                    if line.strip():
                        fout.write("chr" + line)
            gene_bed = gene_bed_fixed
            if verbose:
                print(f"   (added chr prefix to gene BED for matching)")
        elif not peaks_have_chr and genes_have_chr:
            gene_bed_fixed = os.path.join(tmpdir, "genes_nochr.bed")
            with open(gene_bed) as fin, open(gene_bed_fixed, 'w') as fout:
                for line in fin:
                    if line.strip():
                        fout.write(line[3:] if line.startswith("chr") else line)
            gene_bed = gene_bed_fixed
            if verbose:
                print(f"   (stripped chr prefix from gene BED for matching)")

        # Ensure peaks are sorted
        peaks_sorted = os.path.join(tmpdir, "peaks_sorted.bed")
        subprocess.run(f"sort -k1,1 -k2,2n {peaks_bed} > {peaks_sorted}", shell=True, check=True)

        # Ensure genes are sorted
        genes_sorted = os.path.join(tmpdir, "genes_sorted.bed")
        subprocess.run(f"sort -k1,1 -k2,2n {gene_bed} > {genes_sorted}", shell=True, check=True)

        # Run bedtools closest (-d for distance, -t first to break ties)
        closest_out = os.path.join(tmpdir, "closest.tsv")
        cmd = f"bedtools closest -a {peaks_sorted} -b {genes_sorted} -d -t first > {closest_out}"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

        if result.returncode != 0:
            print(f"   bedtools closest failed: {result.stderr}")
            return pd.DataFrame()

        # Detect number of peak columns from input
        with open(peaks_sorted) as f:
            first_line = f.readline().strip()
            n_peak_cols = len(first_line.split('\t')) if first_line else 4

        # Parse output
        # Format: [peak_cols...] [gene_cols(6)] [distance]
        # gene_name is at index n_peak_cols + 3
        rows = []
        with open(closest_out) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < n_peak_cols + 7:
                    continue

                peak_id = parts[3] if n_peak_cols >= 4 else f"{parts[0]}:{parts[1]}-{parts[2]}"
                gene_name = parts[n_peak_cols + 3]
                distance = int(parts[-1]) if parts[-1] != '.' else -1

                rows.append({
                    "peak_id": peak_id,
                    "closest_gene": gene_name,
                    "distance_to_gene": distance,
                })

    df = pd.DataFrame(rows)

    if not df.empty:
        df.to_csv(output_file, sep='\t', index=False)
        if verbose:
            print(f"   Annotated {len(df):,} peaks with closest gene")
            median_dist = df["distance_to_gene"].median()
            at_tss = (df["distance_to_gene"] == 0).sum()
            print(f"   Median distance to TSS: {median_dist:,.0f} bp")
            print(f"   Peaks at TSS (distance=0): {at_tss:,}")

    return df


# =============================================================================
# PEAK ID FUNCTIONS
# =============================================================================

def add_peak_ids(
    input_bed: str,
    output_bed: str,
    prefix: str = "peak",
    verbose: bool = True,
) -> Dict[str, Any]:
    """
    Add unique peak IDs to a BED file.
    Creates BED format: chr, start, end, peak_id, [extra cols...]
    """
    peak_count = 0

    with open(input_bed) as fin, open(output_bed, 'w') as fout:
        for line in fin:
            if line.strip() and not line.startswith('#'):
                parts = line.strip().split('\t')
                chrom, start, end = parts[0], parts[1], parts[2]
                peak_count += 1
                peak_id = f"{prefix}_{peak_count:06d}"
                # Keep extra columns (e.g., species list) after the peak_id
                extra = parts[3:] if len(parts) > 3 else []
                if extra:
                    fout.write(f"{chrom}\t{start}\t{end}\t{peak_id}\t" + "\t".join(extra) + "\n")
                else:
                    fout.write(f"{chrom}\t{start}\t{end}\t{peak_id}\n")

    if verbose:
        print(f"Added IDs to {peak_count:,} peaks ({prefix}_000001 to {prefix}_{peak_count:06d})")

    return {
        "status": "success",
        "peak_count": peak_count,
        "output_file": output_bed,
    }


def add_peak_ids_to_human_specific(
    input_bed: str,
    output_bed: str,
    verbose: bool = True,
) -> Dict[str, Any]:
    """Add human-specific peak IDs (human_peak_NNNNNN).

    .. deprecated::
        No longer used in the pipeline.  Human-specific IDs are assigned
        in ``build_master_annotation()`` after orthology-based classification.
    """
    return add_peak_ids(input_bed, output_bed, prefix="human_peak", verbose=verbose)


# =============================================================================
# LIFTBACK FUNCTION
# =============================================================================

def liftback_peaks(
    input_bed: str,
    output_bed: str,
    species: str,
    chain_dir: str = DEFAULT_CHAIN_DIR,
    liftover_path: str = "liftOver",
    min_match: float = 0.95,
    min_blocks: Optional[float] = None,
    multiple: bool = False,
    auto_chr: bool = True,
    verbose: bool = True,
    ncpu: int = 1,
) -> Dict[str, Any]:
    """
    Liftover peaks from hg38 back to a species genome.
    Preserves peak_id (4th column) if present.
    For Marmoset, performs two-step liftover (hg38 -> calJac4 -> calJac1).
    """
    if verbose:
        print(f"\n   Lifting back to {species}...")

    if species == "Marmoset":
        chain1 = get_reverse_chain_file("Marmoset_step1", chain_dir)
        chain2 = get_reverse_chain_file("Marmoset_step2", chain_dir)
        result = liftover_two_step(
            input_bed=input_bed, output_bed=output_bed,
            chain_file_1=chain1, chain_file_2=chain2,
            liftover_path=liftover_path, min_match=min_match,
            min_blocks=min_blocks, multiple=multiple,
            auto_chr=auto_chr, verbose=verbose, ncpu=ncpu,
        )
    else:
        chain_file = get_reverse_chain_file(species, chain_dir)
        result = liftover_peaks(
            input_bed=input_bed, output_bed=output_bed,
            chain_file=chain_file, liftover_path=liftover_path,
            min_match=min_match, min_blocks=min_blocks,
            multiple=multiple, auto_chr=auto_chr,
            verbose=verbose, ncpu=ncpu,
        )

    return result


# =============================================================================
# POST-LIFTBACK SIZE FILTERING
# =============================================================================

def filter_liftback_by_size(
    liftback_dir: str,
    output_dir: str,
    species_list: List[str],
    min_liftback_size: int = 100,
    max_liftback_size: int = 5000,
    file_pattern: str = "unified_consensus_{species}.bed",
    verbose: bool = True,
) -> Dict[str, Any]:
    """
    Filter liftback BED files by peak size, discarding peaks that became
    too large or too small during the liftover round-trip.

    UCSC liftOver's ``-minMatch`` parameter controls what fraction of input
    bases must remap, but does NOT constrain the size of the output region.
    Structural rearrangements and lineage-specific insertions in the target
    genome can cause a 500 bp input peak to span hundreds of kilobases in
    the output.  This function removes such artefacts.

    For every species, writes:
      - ``{output_dir}/{stem}_filtered.bed``  — peaks that pass the size filter
      - ``{output_dir}/{stem}_size_rejected.bed`` — peaks that failed

    Parameters
    ----------
    liftback_dir : str
        Directory containing liftback BED files (from Step 4).
    output_dir : str
        Directory for filtered output files.
    species_list : list of str
        NHP species names (e.g. ["Bonobo", "Chimpanzee", ...]).
    min_liftback_size : int
        Minimum allowed peak size in target genome (default 100 bp).
    max_liftback_size : int
        Maximum allowed peak size in target genome (default 5000 bp).
    file_pattern : str
        Filename pattern with ``{species}`` placeholder.
    verbose : bool
        Print per-species filtering statistics.

    Returns
    -------
    dict
        ``{species: {"input": N, "kept": N, "rejected": N,
                     "output_file": path, "rejected_file": path}}``
    """
    os.makedirs(output_dir, exist_ok=True)
    results = {}

    if verbose:
        print(f"\nFiltering liftback peaks by size "
              f"[{min_liftback_size:,} – {max_liftback_size:,} bp]")
        print(f"  Input dir:  {liftback_dir}")
        print(f"  Output dir: {output_dir}")
        print(f"  {'Species':<15s} {'Input':>8s} {'Kept':>8s} "
              f"{'Rejected':>8s} {'%Kept':>8s}")
        print("  " + "-" * 55)

    for species in species_list:
        fname = file_pattern.format(species=species)
        input_bed = os.path.join(liftback_dir, fname)
        if not os.path.exists(input_bed):
            if verbose:
                print(f"  {species:<15s} (file not found, skipping)")
            continue

        stem = fname.replace(".bed", "")
        filtered_bed = os.path.join(output_dir, f"{stem}_filtered.bed")
        rejected_bed = os.path.join(output_dir, f"{stem}_size_rejected.bed")

        n_input = n_kept = n_rejected = 0
        with open(input_bed) as fin, \
             open(filtered_bed, 'w') as fout, \
             open(rejected_bed, 'w') as frej:
            for line in fin:
                if not line.strip() or line.startswith('#'):
                    continue
                n_input += 1
                parts = line.strip().split('\t')
                size = int(parts[2]) - int(parts[1])
                # Write only chr, start, end, peak_id (4 columns)
                bed4 = f"{parts[0]}\t{parts[1]}\t{parts[2]}\t{parts[3]}\n"
                if min_liftback_size <= size <= max_liftback_size:
                    fout.write(bed4)
                    n_kept += 1
                else:
                    frej.write(bed4)
                    n_rejected += 1

        pct = (n_kept / n_input * 100) if n_input > 0 else 0
        results[species] = {
            "input": n_input,
            "kept": n_kept,
            "rejected": n_rejected,
            "output_file": filtered_bed,
            "rejected_file": rejected_bed,
        }

        if verbose:
            print(f"  {species:<15s} {n_input:>8,} {n_kept:>8,} "
                  f"{n_rejected:>8,} {pct:>7.1f}%")

    return results


def update_master_after_filtering(
    master_file: str,
    filter_results: Dict[str, Dict],
    output_file: Optional[str] = None,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Update the master annotation table after liftback size filtering.

    For each species where peaks were rejected, sets ``{Species}_orth = 0``
    and clears the liftback coordinates.  Then recalculates
    ``n_species_orth`` and reclassifies peaks whose orthology count dropped
    to 1 as ``human_specific``.

    Parameters
    ----------
    master_file : str
        Path to the existing master_annotation.tsv (from Step 7).
    filter_results : dict
        Output of :func:`filter_liftback_by_size`.
    output_file : str, optional
        Where to save the updated table.  If *None*, overwrites
        ``master_file`` with ``_filtered`` suffix.
    verbose : bool
        Print summary of changes.

    Returns
    -------
    pd.DataFrame
        Updated master annotation.
    """
    df = pd.read_csv(master_file, sep="\t", index_col="peak_id")
    original_types = df["peak_type"].value_counts().to_dict()

    if output_file is None:
        output_file = master_file.replace(".tsv", "_filtered.tsv")

    # Collect rejected peak IDs per species from the rejected BED files
    rejected_ids = {}  # species -> set of peak_ids
    for species, info in filter_results.items():
        rej_file = info.get("rejected_file")
        if rej_file and os.path.exists(rej_file):
            ids = set()
            with open(rej_file) as f:
                for line in f:
                    if line.strip() and not line.startswith('#'):
                        parts = line.strip().split('\t')
                        if len(parts) >= 4:
                            ids.add(parts[3])
            rejected_ids[species] = ids

    # Clear orthology for rejected peaks
    total_cleared = 0
    for species, ids in rejected_ids.items():
        orth_col = f"{species}_orth"
        chr_col = f"{species}_chr"
        start_col = f"{species}_start"
        end_col = f"{species}_end"

        mask = df.index.isin(ids)
        n_cleared = mask.sum()
        total_cleared += n_cleared

        if orth_col in df.columns:
            df.loc[mask, orth_col] = 0
        for col in [chr_col, start_col, end_col]:
            if col in df.columns:
                df.loc[mask, col] = np.nan

        if verbose and n_cleared > 0:
            print(f"  {species}: cleared orthology for {n_cleared:,} "
                  f"size-rejected peaks")

    # Recalculate n_species_orth
    orth_cols = [c for c in df.columns if c.endswith("_orth")]
    df["n_species_orth"] = df[orth_cols].sum(axis=1)

    # Reclassify: unified peaks that lost all NHP orthologs -> human_specific
    mask_reclassify = (
        (df["peak_type"] == "unified") &
        (df["n_species_orth"] == 1)
    )
    n_reclassified = mask_reclassify.sum()
    df.loc[mask_reclassify, "peak_type"] = "human_specific"

    # Reassign IDs for newly reclassified peaks
    if n_reclassified > 0:
        # Find current max human_peak number
        existing_hs = [idx for idx in df.index if idx.startswith("human_peak_")]
        if existing_hs:
            max_hs_num = max(int(idx.split("_")[-1]) for idx in existing_hs)
        else:
            max_hs_num = 0

        new_mapping = {}
        counter = max_hs_num
        for idx in df.index[mask_reclassify]:
            if not idx.startswith("human_peak_"):
                counter += 1
                new_mapping[idx] = f"human_peak_{counter:06d}"

        if new_mapping:
            df = df.rename(index=new_mapping)

    # Save
    df.to_csv(output_file, sep="\t")

    new_types = df["peak_type"].value_counts().to_dict()

    if verbose:
        print(f"\n  Master annotation updated:")
        print(f"    Before: {original_types}")
        print(f"    After:  {new_types}")
        if n_reclassified > 0:
            print(f"    Reclassified {n_reclassified:,} unified -> human_specific "
                  f"(lost all NHP orthologs after size filter)")
        print(f"    Saved: {output_file}")

    return df


def filter_and_update_master(
    output_dir: str,
    species_list: List[str],
    min_liftback_size: int = 100,
    max_liftback_size: int = 5000,
    liftback_subdir: str = "04_lifted_back",
    master_subdir: str = "07_master_annotation",
    filter_subdir: str = "04_lifted_back_filtered",
    verbose: bool = True,
) -> Dict[str, Any]:
    """
    One-call convenience: filter all liftback files + update master annotation.

    This is the **recommended way** to apply post-liftback size filtering
    to an existing pipeline run without re-running the full pipeline.

    Creates files in ``{output_dir}/{filter_subdir}/`` and saves the updated
    master annotation as ``master_annotation_filtered.tsv``.

    Parameters
    ----------
    output_dir : str
        Top-level pipeline output directory (e.g. ``cross_species_consensus_v2``).
    species_list : list of str
        NHP species (e.g. ``["Bonobo", "Chimpanzee", ...]``).
    min_liftback_size : int
        Minimum peak size in target genome (default 100 bp).
    max_liftback_size : int
        Maximum peak size in target genome (default 5000 bp).
    liftback_subdir : str
        Subdirectory with liftback BED files (default ``04_lifted_back``).
    master_subdir : str
        Subdirectory with master annotation (default ``07_master_annotation``).
    filter_subdir : str
        Subdirectory for filtered output (default ``04_lifted_back_filtered``).
    verbose : bool
        Print filtering stats and master annotation changes.

    Returns
    -------
    dict
        ``{"filter_results": {species: stats}, "master_df": DataFrame,
          "master_file": path, "filter_dir": path}``

    Example
    -------
    >>> from src.cross_species import filter_and_update_master
    >>> result = filter_and_update_master(
    ...     output_dir="/path/to/cross_species_consensus_v2",
    ...     species_list=["Bonobo", "Chimpanzee", "Gorilla", "Macaque", "Marmoset"],
    ...     min_liftback_size=100,
    ...     max_liftback_size=5000,
    ... )
    """
    liftback_dir = os.path.join(output_dir, liftback_subdir)
    master_file = os.path.join(output_dir, master_subdir, "master_annotation.tsv")
    filter_dir = os.path.join(output_dir, filter_subdir)

    if verbose:
        print("=" * 70)
        print("POST-LIFTBACK SIZE FILTERING")
        print("=" * 70)

    # 1. Filter liftback BEDs
    filter_results = filter_liftback_by_size(
        liftback_dir=liftback_dir,
        output_dir=filter_dir,
        species_list=species_list,
        min_liftback_size=min_liftback_size,
        max_liftback_size=max_liftback_size,
        verbose=verbose,
    )

    # 2. Update master annotation
    filtered_master_file = os.path.join(
        output_dir, master_subdir, "master_annotation_filtered.tsv"
    )
    master_df = update_master_after_filtering(
        master_file=master_file,
        filter_results=filter_results,
        output_file=filtered_master_file,
        verbose=verbose,
    )

    # 3. Regenerate combined BED files with filtered liftback coords
    combined_dir = os.path.join(output_dir, "09_combined_filtered")
    os.makedirs(combined_dir, exist_ok=True)

    if verbose:
        print(f"\n  Regenerating combined BED files in {combined_dir}/")

    # Unified peaks hg38
    unified_hg38 = os.path.join(combined_dir, "unified_peaks_hg38.bed")
    n_uni = 0
    with open(unified_hg38, 'w') as fout:
        for pid, row in master_df[master_df["peak_type"] == "unified"].iterrows():
            fout.write(f"{row['Human_chr']}\t{int(row['Human_start'])}\t"
                       f"{int(row['Human_end'])}\t{pid}\n")
            n_uni += 1

    # Human-specific peaks
    hs_hg38 = os.path.join(combined_dir, "human_specific_peaks_hg38.bed")
    n_hs = 0
    with open(hs_hg38, 'w') as fout:
        for pid, row in master_df[master_df["peak_type"] == "human_specific"].iterrows():
            fout.write(f"{row['Human_chr']}\t{int(row['Human_start'])}\t"
                       f"{int(row['Human_end'])}\t{pid}\n")
            n_hs += 1

    # Per-species unified (from filtered liftback)
    for species in species_list:
        sp_info = filter_results.get(species, {})
        filtered_bed = sp_info.get("output_file")
        sp_out = os.path.join(combined_dir, f"unified_peaks_{species}.bed")
        n = 0
        if filtered_bed and os.path.exists(filtered_bed):
            with open(filtered_bed) as fin, open(sp_out, 'w') as fout:
                for line in fin:
                    if line.strip() and not line.startswith('#'):
                        parts = line.strip().split('\t')
                        fout.write(f"{parts[0]}\t{parts[1]}\t{parts[2]}\t{parts[3]}\n")
                        n += 1
        if verbose:
            print(f"    unified_peaks_{species}.bed: {n:,}")

    # Per-species all_peaks (filtered liftback + species-specific)
    sp_specific_dir = os.path.join(output_dir, "05_species_specific")
    for species in species_list:
        sp_all = os.path.join(combined_dir, f"all_peaks_{species}.bed")
        n = 0
        with open(sp_all, 'w') as fout:
            # Filtered liftback
            sp_info = filter_results.get(species, {})
            filtered_bed = sp_info.get("output_file")
            if filtered_bed and os.path.exists(filtered_bed):
                with open(filtered_bed) as fin:
                    for line in fin:
                        if line.strip() and not line.startswith('#'):
                            parts = line.strip().split('\t')
                            fout.write(f"{parts[0]}\t{parts[1]}\t{parts[2]}\t{parts[3]}\n")
                            n += 1
            # Species-specific
            sp_bed = os.path.join(sp_specific_dir, f"{species}_specific_peaks.bed")
            if os.path.exists(sp_bed):
                with open(sp_bed) as fin:
                    for line in fin:
                        if line.strip() and not line.startswith('#'):
                            parts = line.strip().split('\t')
                            fout.write(f"{parts[0]}\t{parts[1]}\t{parts[2]}\t{parts[3]}\n")
                            n += 1
        if verbose:
            print(f"    all_peaks_{species}.bed: {n:,}")

    # Human all_peaks
    human_all = os.path.join(combined_dir, "all_peaks_Human.bed")
    n_human = 0
    with open(human_all, 'w') as fout:
        for pid, row in master_df[
            master_df["peak_type"].isin(["unified", "human_specific"])
        ].iterrows():
            fout.write(f"{row['Human_chr']}\t{int(row['Human_start'])}\t"
                       f"{int(row['Human_end'])}\t{pid}\n")
            n_human += 1
    if verbose:
        print(f"    all_peaks_Human.bed: {n_human:,}")

    if verbose:
        print(f"\n  unified_peaks_hg38.bed:        {n_uni:,}")
        print(f"  human_specific_peaks_hg38.bed: {n_hs:,}")
        print(f"\n  Done. Use 09_combined_filtered/ for downstream analysis.")

    return {
        "filter_results": filter_results,
        "master_df": master_df,
        "master_file": filtered_master_file,
        "filter_dir": filter_dir,
        "combined_dir": combined_dir,
    }


# =============================================================================
# PEAK RESIZING
# =============================================================================

def resize_peaks(
    input_bed: str,
    output_bed: str,
    target_size: int = 500,
    narrowpeak_file: Optional[str] = None,
    chromsizes_file: Optional[str] = None,
    verbose: bool = True,
) -> Dict[str, Any]:
    """
    Resize peaks to a uniform fixed width.

    For every region in *input_bed*, determine a centre point and expand
    symmetrically to *target_size*.  Centre selection priority:

    1. **Summit from overlapping narrowPeak** — if *narrowpeak_file* is
       supplied and has a summit column (col 10), the summit position of
       the best-overlapping narrowPeak entry is used.
    2. **Centre of overlapping narrowPeak** — if the narrowPeak file has
       no summit column (or summit is -1), the midpoint of the
       overlapping narrowPeak region is used.
    3. **Own centre** — if no overlap is found, the midpoint of the
       input region itself is used.

    Parameters
    ----------
    input_bed : str
        BED file with regions to resize (≥4 columns: chr, start, end, id).
    output_bed : str
        Path for output BED file (chr, start, end, id).
    target_size : int
        Desired uniform peak width (default 500 bp).
    narrowpeak_file : str, optional
        NarrowPeak file (10 columns) or BED file with peak coordinates.
        Used to find summits or centres for anchoring.
    chromsizes_file : str, optional
        Tab-delimited ``chrom<TAB>size`` file.  If provided, peaks are
        clipped to chromosome boundaries.
    verbose : bool
        Print summary statistics.

    Returns
    -------
    dict
        ``{"input": N, "output": N, "summit_anchored": N,
          "narrowpeak_center_anchored": N, "self_center_anchored": N,
          "clipped": N, "output_file": path}``
    """
    os.makedirs(os.path.dirname(output_bed) or '.', exist_ok=True)
    half = target_size // 2

    # Load chromsizes if provided
    chrom_max = {}
    if chromsizes_file and os.path.exists(chromsizes_file):
        with open(chromsizes_file) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    chrom_max[parts[0]] = int(parts[1])

    # Load narrowPeak data if provided
    np_summits = {}   # (chrom, summit_pos) indexed by interval tree
    np_entries = []    # list of (chrom, start, end, summit_abs) tuples
    has_summit_col = False

    if narrowpeak_file and os.path.exists(narrowpeak_file):
        with open(narrowpeak_file) as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.strip().split('\t')
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                # NarrowPeak col 10 (0-indexed 9) is summit offset from start
                summit_abs = None
                if len(parts) >= 10:
                    try:
                        summit_offset = int(parts[9])
                        if summit_offset >= 0:
                            summit_abs = start + summit_offset
                            has_summit_col = True
                    except (ValueError, IndexError):
                        pass
                np_entries.append((chrom, start, end, summit_abs))

        if verbose:
            print(f"  NarrowPeak file: {len(np_entries):,} entries"
                  f" ({'with' if has_summit_col else 'without'} summits)")

    # Build simple interval lookup: chrom -> sorted list of (start, end, summit)
    np_by_chrom = defaultdict(list)
    for chrom, start, end, summit in np_entries:
        np_by_chrom[chrom].append((start, end, summit))
    for chrom in np_by_chrom:
        np_by_chrom[chrom].sort()

    def _find_best_overlap(chrom, qstart, qend):
        """Find the narrowPeak entry with the most overlap."""
        best_ovl = 0
        best_entry = None
        for (s, e, summit) in np_by_chrom.get(chrom, []):
            if s >= qend:
                break
            if e <= qstart:
                continue
            ovl = min(qend, e) - max(qstart, s)
            if ovl > best_ovl:
                best_ovl = ovl
                best_entry = (s, e, summit)
        return best_entry if best_ovl > 0 else None

    # Process peaks
    n_input = 0
    n_summit = 0
    n_np_center = 0
    n_self_center = 0
    n_clipped = 0

    with open(input_bed) as fin, open(output_bed, 'w') as fout:
        for line in fin:
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            peak_id = parts[3] if len(parts) > 3 else f"{chrom}:{start}-{end}"
            n_input += 1

            # Determine centre
            center = None
            overlap = _find_best_overlap(chrom, start, end) if np_entries else None

            if overlap is not None:
                np_s, np_e, np_summit = overlap
                if np_summit is not None:
                    center = np_summit
                    n_summit += 1
                else:
                    center = (np_s + np_e) // 2
                    n_np_center += 1
            else:
                center = (start + end) // 2
                n_self_center += 1

            # Expand around centre
            new_start = center - half
            new_end = center + half
            if target_size % 2 != 0:
                new_end += 1  # handle odd sizes

            # Clip to chromosome boundaries
            if new_start < 0:
                new_start = 0
                new_end = target_size
                n_clipped += 1
            if chrom_max and chrom in chrom_max:
                cmax = chrom_max[chrom]
                if new_end > cmax:
                    new_end = cmax
                    new_start = max(0, cmax - target_size)
                    n_clipped += 1

            fout.write(f"{chrom}\t{new_start}\t{new_end}\t{peak_id}\n")

    if verbose:
        print(f"  Resized {n_input:,} peaks to {target_size} bp:")
        print(f"    Summit-anchored:        {n_summit:,}")
        print(f"    NarrowPeak-center:      {n_np_center:,}")
        print(f"    Self-center:            {n_self_center:,}")
        if n_clipped:
            print(f"    Clipped to chrom edge:  {n_clipped:,}")

    return {
        "input": n_input,
        "output": n_input,
        "summit_anchored": n_summit,
        "narrowpeak_center_anchored": n_np_center,
        "self_center_anchored": n_self_center,
        "clipped": n_clipped,
        "output_file": output_bed,
    }


# =============================================================================
# RECIPROCAL LIFTOVER CHECK
# =============================================================================

def reciprocal_liftover_check(
    hg38_bed: str,
    liftback_dir: str,
    output_dir: str,
    species_list: List[str],
    chain_dir: str = DEFAULT_CHAIN_DIR,
    liftover_path: str = "liftOver",
    min_match: Union[float, Dict[str, float]] = 0.95,
    min_blocks: Optional[float] = None,
    multiple: bool = False,
    max_distance: int = 500,
    file_pattern: str = "unified_consensus_{species}.bed",
    verbose: bool = True,
    ncpu: int = 1,
) -> Dict[str, Any]:
    """
    Verify that unified peaks survive a reciprocal (round-trip) liftover.

    For each NHP species the liftback BED (hg38 → NHP from Step 4) is
    lifted *back* to hg38 using the **forward** chain (NHP → hg38).  Each
    peak passes only if the round-trip hg38 coordinates overlap the
    original hg38 coordinates (or are within *max_distance* bp).

    This catches peaks where the liftover lands in a non-orthologous
    location (e.g. a paralogous repeat on a different chromosome or
    a far-away region on the same chromosome).

    For every species, writes:
      - ``{species}_reciprocal_pass.bed``  — peaks that survive round-trip
      - ``{species}_reciprocal_fail.bed``  — peaks that fail

    Parameters
    ----------
    hg38_bed : str
        Original unified peaks in hg38 (chr, start, end, peak_id).
    liftback_dir : str
        Directory with liftback BED files (``unified_consensus_{species}.bed``).
    output_dir : str
        Directory for reciprocal-check output.
    species_list : list of str
        NHP species to check.
    chain_dir : str
        Directory containing chain files.
    liftover_path : str
        Path to UCSC ``liftOver`` executable.
    min_match : float or dict
        Minimum match ratio for the reciprocal liftover.
    min_blocks : float, optional
        Minimum block ratio.
    multiple : bool
        Allow multiple mappings (default False).
    max_distance : int
        Maximum distance (bp) between original and round-trip start
        positions to consider a pass (default 500).
    file_pattern : str
        Filename pattern with ``{species}`` placeholder for liftback BEDs
        (default ``unified_consensus_{species}.bed``).
    verbose : bool
        Print per-species statistics.
    ncpu : int
        Parallel workers for liftover.

    Returns
    -------
    dict
        ``{species: {"input": N, "passed": N, "failed": N, "unmapped": N,
                     "pass_file": path, "fail_file": path}}``
    """
    os.makedirs(output_dir, exist_ok=True)

    def _get_mm(species):
        if isinstance(min_match, dict):
            return min_match.get(species, 0.95)
        return min_match

    # Load original hg38 coordinates: peak_id -> (chr, start, end)
    hg38_coords = {}
    with open(hg38_bed) as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                parts = line.strip().split('\t')
                hg38_coords[parts[3]] = (parts[0], int(parts[1]), int(parts[2]))

    results = {}

    if verbose:
        print(f"\nReciprocal liftover check (max_distance={max_distance} bp)")
        print(f"  Original hg38 peaks: {len(hg38_coords):,}")
        print(f"  {'Species':<15s} {'Input':>8s} {'Lifted':>8s} "
              f"{'Pass':>8s} {'Fail':>8s} {'Unmap':>8s} {'%Pass':>8s}")
        print("  " + "-" * 65)

    for species in species_list:
        lb_file = os.path.join(liftback_dir,
                               file_pattern.format(species=species))
        if not os.path.exists(lb_file):
            if verbose:
                print(f"  {species:<15s} (liftback file not found, skipping)")
            continue

        # Clean the liftback to BED4 (in case it has extra columns)
        with tempfile.TemporaryDirectory() as tmpdir:
            lb_bed4 = os.path.join(tmpdir, "liftback_bed4.bed")
            n_input = 0
            with open(lb_file) as fin, open(lb_bed4, 'w') as fout:
                for line in fin:
                    if line.strip() and not line.startswith('#'):
                        parts = line.strip().split('\t')
                        fout.write(f"{parts[0]}\t{parts[1]}\t{parts[2]}\t{parts[3]}\n")
                        n_input += 1

            # Lift NHP coords back to hg38 using forward chain
            round_trip_bed = os.path.join(tmpdir, "round_trip_hg38.bed")

            if species == "Marmoset":
                # Reverse of the 2-step: calJac1 → calJac4 → hg38
                chain1 = get_chain_file("Marmoset_step1", chain_dir)  # calJac1→calJac4
                chain2 = get_chain_file("Marmoset_step2", chain_dir)  # calJac4→hg38
                liftover_two_step(
                    input_bed=lb_bed4, output_bed=round_trip_bed,
                    chain_file_1=chain1, chain_file_2=chain2,
                    liftover_path=liftover_path,
                    min_match=_get_mm(species),
                    min_blocks=min_blocks, multiple=multiple,
                    auto_chr=True, verbose=False, ncpu=ncpu,
                )
            else:
                chain_fwd = get_chain_file(species, chain_dir)  # NHP → hg38
                liftover_peaks(
                    input_bed=lb_bed4, output_bed=round_trip_bed,
                    chain_file=chain_fwd, liftover_path=liftover_path,
                    min_match=_get_mm(species),
                    min_blocks=min_blocks, multiple=multiple,
                    auto_chr=True, verbose=False, ncpu=ncpu,
                )

            # Read round-trip hg38 coordinates
            rt_coords = {}  # peak_id -> (chr, start, end)
            if os.path.exists(round_trip_bed):
                with open(round_trip_bed) as f:
                    for line in f:
                        if line.strip() and not line.startswith('#'):
                            parts = line.strip().split('\t')
                            rt_coords[parts[3]] = (parts[0], int(parts[1]), int(parts[2]))

        # Compare original vs round-trip
        pass_file = os.path.join(output_dir, f"{species}_reciprocal_pass.bed")
        fail_file = os.path.join(output_dir, f"{species}_reciprocal_fail.bed")

        n_pass = n_fail = n_unmapped = 0

        # Also read original liftback coords for writing pass/fail in NHP coords
        lb_peaks = {}  # peak_id -> original liftback line (BED4)
        with open(lb_file) as f:
            for line in f:
                if line.strip() and not line.startswith('#'):
                    parts = line.strip().split('\t')
                    lb_peaks[parts[3]] = f"{parts[0]}\t{parts[1]}\t{parts[2]}\t{parts[3]}\n"

        with open(pass_file, 'w') as fp, open(fail_file, 'w') as ff:
            for peak_id, lb_line in lb_peaks.items():
                if peak_id not in rt_coords:
                    ff.write(lb_line)
                    n_unmapped += 1
                    continue

                orig = hg38_coords.get(peak_id)
                if orig is None:
                    ff.write(lb_line)
                    n_fail += 1
                    continue

                rt_chr, rt_start, rt_end = rt_coords[peak_id]
                orig_chr, orig_start, orig_end = orig

                # Check: same chromosome and overlapping or within max_distance
                if (rt_chr == orig_chr and
                        rt_start <= orig_end + max_distance and
                        rt_end >= orig_start - max_distance):
                    fp.write(lb_line)
                    n_pass += 1
                else:
                    ff.write(lb_line)
                    n_fail += 1

        pct = (n_pass / n_input * 100) if n_input > 0 else 0
        results[species] = {
            "input": n_input,
            "passed": n_pass,
            "failed": n_fail,
            "unmapped": n_unmapped,
            "pass_file": pass_file,
            "fail_file": fail_file,
        }

        if verbose:
            print(f"  {species:<15s} {n_input:>8,} {n_pass + n_fail + n_unmapped:>8,} "
                  f"{n_pass:>8,} {n_fail:>8,} {n_unmapped:>8,} {pct:>7.1f}%")

    return results


# =============================================================================
# FINALIZE PEAK SETS (reciprocal + size filter + resize)
# =============================================================================

def finalize_peak_sets(
    output_dir: str,
    species_list: List[str],
    target_size: int = 500,
    narrowpeak_files: Optional[Dict[str, str]] = None,
    chromsizes_files: Optional[Dict[str, str]] = None,
    min_liftback_size: int = 100,
    max_liftback_size: int = 5000,
    reciprocal_check: bool = True,
    chain_dir: str = DEFAULT_CHAIN_DIR,
    liftover_path: str = "liftOver",
    min_match: Union[float, Dict[str, float]] = 0.95,
    max_distance: int = 500,
    ncpu: int = 1,
    species_beds: Optional[Dict[str, str]] = None,
    human_bed: Optional[str] = None,
    liftback_subdir: str = "04_lifted_back",
    master_subdir: str = "07_master_annotation",
    merged_subdir: str = "02_merged_consensus",
    species_specific_subdir: str = "05_species_specific",
    final_subdir: str = "10_final",
    verbose: bool = True,
) -> Dict[str, Any]:
    """
    Complete post-pipeline finalization: reciprocal check → size filter →
    resize → write final uniform-width BED files.

    This is designed to run on existing pipeline output **without**
    re-running the full pipeline.

    Processing order for unified peaks per species
    -----------------------------------------------
    1. **Size filter** — discard liftback peaks outside
       ``[min_liftback_size, max_liftback_size]``.
    2. **Reciprocal liftover check** — lift NHP coords back to hg38 and
       verify overlap with original hg38 coords.
    3. **Update master annotation** — recalculate ``n_species_orth`` and
       reclassify unified → human_specific where needed.
    4. **Resize** — expand/shrink every surviving peak to *target_size*
       using summit info from *narrowpeak_files*.

    Output in ``{output_dir}/{final_subdir}/``
    -------------------------------------------
    - ``all_peaks_{Species}.bed``  — unified + species-specific, all
      same *target_size* width
    - ``unified_peaks_hg38.bed``   — unified peaks in hg38
    - ``unified_peaks_{Species}.bed`` — unified peaks in species coords
    - ``human_specific_peaks_hg38.bed``
    - ``loss_summary.tsv``         — per-species, per-step peak counts
    - ``{species}_rescued_peaks.bed`` — original peaks lost during
      filtering, rescued back into the final set
    - Intermediate files in ``{final_subdir}/intermediates/``

    Parameters
    ----------
    output_dir : str
        Top-level pipeline output directory.
    species_list : list of str
        NHP species names.
    target_size : int
        Uniform peak width for final output (default 500 bp).
    narrowpeak_files : dict, optional
        ``{species: path}`` with species narrowPeak / peak files for
        summit-based centring.  Key ``"Human"`` for human peaks.
    chromsizes_files : dict, optional
        ``{species: path}`` to chrom.sizes files for boundary clipping.
    min_liftback_size : int
        Minimum liftback peak size (default 100).
    max_liftback_size : int
        Maximum liftback peak size (default 5000).
    chain_dir : str
        Chain file directory.
    liftover_path : str
        Path to ``liftOver`` executable.
    min_match : float or dict
        Min match ratio for reciprocal liftover.
    max_distance : int
        Max distance (bp) for reciprocal overlap test (default 500).
    ncpu : int
        Parallel workers.
    species_beds : dict, optional
        ``{species: path}`` mapping NHP species name to its original
        consensus peaks BED file.  When provided a rescue step (Step F)
        compares final output with original input and re-adds any peaks
        that were lost during filtering, labelled
        ``{species}_rescued_NNNNNN``.
    human_bed : str, optional
        Path to the original human consensus peaks BED file.
        Used together with *species_beds* for the rescue step.
    liftback_subdir, master_subdir, merged_subdir,
    species_specific_subdir, final_subdir : str
        Subdirectory names within *output_dir*.
    verbose : bool
        Print progress and loss summaries.

    Returns
    -------
    dict
        ``{"loss_summary": DataFrame, "master_df": DataFrame,
          "final_dir": path, ...}``
    """
    liftback_dir = os.path.join(output_dir, liftback_subdir)
    master_file = os.path.join(output_dir, master_subdir, "master_annotation.tsv")
    merged_bed = os.path.join(output_dir, merged_subdir,
                              "unified_consensus_hg38_bed4.bed")
    # Fallback if bed4 doesn't exist
    if not os.path.exists(merged_bed):
        merged_bed = os.path.join(output_dir, merged_subdir,
                                  "unified_consensus_hg38_with_ids.bed")
    sp_specific_dir = os.path.join(output_dir, species_specific_subdir)
    final_dir = os.path.join(output_dir, final_subdir)
    intermediates = os.path.join(final_dir, "intermediates")
    os.makedirs(intermediates, exist_ok=True)

    if narrowpeak_files is None:
        narrowpeak_files = {}
    if chromsizes_files is None:
        chromsizes_files = {}

    loss_rows = []  # for loss_summary.tsv

    if verbose:
        print("=" * 70)
        print("FINALIZE PEAK SETS")
        print(f"  target_size = {target_size} bp")
        print(f"  size filter = [{min_liftback_size}, {max_liftback_size}] bp")
        print(f"  reciprocal max_distance = {max_distance} bp")
        print("=" * 70)

    # ------------------------------------------------------------------
    # Load original hg38 peak counts (before any filtering)
    # ------------------------------------------------------------------
    hg38_bed = merged_bed
    # Prepare a clean BED4 from the merged file for reciprocal check
    hg38_bed4 = os.path.join(intermediates, "unified_hg38_bed4.bed")
    n_hg38_total = 0
    with open(hg38_bed) as fin, open(hg38_bed4, 'w') as fout:
        for line in fin:
            if line.strip() and not line.startswith('#'):
                parts = line.strip().split('\t')
                fout.write(f"{parts[0]}\t{parts[1]}\t{parts[2]}\t{parts[3]}\n")
                n_hg38_total += 1

    if verbose:
        print(f"\n  Total unified hg38 peaks (input): {n_hg38_total:,}")

    # ------------------------------------------------------------------
    # STEP A: Size-filter liftback BEDs
    # ------------------------------------------------------------------
    if verbose:
        print(f"\n--- Step A: Size filter [{min_liftback_size:,}–"
              f"{max_liftback_size:,} bp] ---")

    size_filter_dir = os.path.join(intermediates, "size_filtered")
    size_results = filter_liftback_by_size(
        liftback_dir=liftback_dir,
        output_dir=size_filter_dir,
        species_list=species_list,
        min_liftback_size=min_liftback_size,
        max_liftback_size=max_liftback_size,
        verbose=verbose,
    )

    for species, info in size_results.items():
        loss_rows.append({
            "species": species, "step": "A_size_filter",
            "input": info["input"], "kept": info["kept"],
            "removed": info["rejected"],
        })

    # ------------------------------------------------------------------
    # STEP B: Reciprocal liftover check (on size-filtered liftback)
    # ------------------------------------------------------------------
    reciprocal_results = {}
    if reciprocal_check:
        if verbose:
            print(f"\n--- Step B: Reciprocal liftover check ---")

        reciprocal_dir = os.path.join(intermediates, "reciprocal")
        reciprocal_results = reciprocal_liftover_check(
            hg38_bed=hg38_bed4,
            liftback_dir=size_filter_dir,
            output_dir=reciprocal_dir,
            species_list=species_list,
            chain_dir=chain_dir,
            liftover_path=liftover_path,
            min_match=min_match,
            max_distance=max_distance,
            file_pattern="unified_consensus_{species}_filtered.bed",
            verbose=verbose,
            ncpu=ncpu,
        )

        for species, info in reciprocal_results.items():
            loss_rows.append({
                "species": species, "step": "B_reciprocal_check",
                "input": info["input"],
                "kept": info["passed"],
                "removed": info["failed"] + info["unmapped"],
            })
    else:
        if verbose:
            print(f"\n--- Step B: Reciprocal liftover check (SKIPPED) ---")

    # ------------------------------------------------------------------
    # STEP C: Update master annotation with combined filtering
    # ------------------------------------------------------------------
    if verbose:
        print(f"\n--- Step C: Update master annotation ---")

    master_df = pd.read_csv(master_file, sep="\t", index_col="peak_id")
    original_types = master_df["peak_type"].value_counts().to_dict()

    # Collect peak IDs that PASSED all filters per species
    passed_ids = {}
    if reciprocal_check:
        # Use reciprocal pass files
        for species, info in reciprocal_results.items():
            pf = info.get("pass_file")
            if pf and os.path.exists(pf):
                ids = set()
                with open(pf) as f:
                    for line in f:
                        if line.strip() and not line.startswith('#'):
                            parts = line.strip().split('\t')
                            if len(parts) >= 4:
                                ids.add(parts[3])
                passed_ids[species] = ids
    else:
        # Use size-filtered files (no reciprocal check)
        for species, info in size_results.items():
            sf = info.get("output_file")
            if sf and os.path.exists(sf):
                ids = set()
                with open(sf) as f:
                    for line in f:
                        if line.strip() and not line.startswith('#'):
                            parts = line.strip().split('\t')
                            if len(parts) >= 4:
                                ids.add(parts[3])
                passed_ids[species] = ids

    # Clear orthology for peaks that did NOT pass
    for species in species_list:
        orth_col = f"{species}_orth"
        chr_col = f"{species}_chr"
        start_col = f"{species}_start"
        end_col = f"{species}_end"

        if orth_col not in master_df.columns:
            continue

        # Peaks that currently have orth=1 but didn't pass
        sp_passed = passed_ids.get(species, set())
        currently_orth = master_df.index[master_df[orth_col] == 1]
        failed = set(currently_orth) - sp_passed

        if failed:
            mask = master_df.index.isin(failed)
            master_df.loc[mask, orth_col] = 0
            for col in [chr_col, start_col, end_col]:
                if col in master_df.columns:
                    master_df.loc[mask, col] = np.nan

            if verbose:
                print(f"  {species}: cleared orthology for "
                      f"{len(failed):,} peaks (size + reciprocal filter)")

    # Recalculate n_species_orth
    orth_cols = [c for c in master_df.columns if c.endswith("_orth")]
    master_df["n_species_orth"] = master_df[orth_cols].sum(axis=1)

    # Reclassify
    mask_reclassify = (
        (master_df["peak_type"] == "unified") &
        (master_df["n_species_orth"] == 1)
    )
    n_reclassified = mask_reclassify.sum()
    master_df.loc[mask_reclassify, "peak_type"] = "human_specific"

    # Reassign IDs for newly reclassified peaks
    if n_reclassified > 0:
        existing_hs = [idx for idx in master_df.index
                       if idx.startswith("human_peak_")]
        max_hs_num = (max(int(idx.split("_")[-1]) for idx in existing_hs)
                      if existing_hs else 0)
        new_mapping = {}
        counter = max_hs_num
        for idx in master_df.index[mask_reclassify]:
            if not idx.startswith("human_peak_"):
                counter += 1
                new_mapping[idx] = f"human_peak_{counter:06d}"
        if new_mapping:
            master_df = master_df.rename(index=new_mapping)
            # Update passed_ids with new names so downstream is consistent
            reverse_map = {v: k for k, v in new_mapping.items()}
            for species in passed_ids:
                passed_ids[species] = {
                    new_mapping.get(pid, pid) for pid in passed_ids[species]
                }

    new_types = master_df["peak_type"].value_counts().to_dict()
    if verbose:
        print(f"\n  Peak types before: {original_types}")
        print(f"  Peak types after:  {new_types}")
        if n_reclassified:
            print(f"  Reclassified {n_reclassified:,} unified → human_specific")

    # Save updated master
    master_out = os.path.join(final_dir, "master_annotation_final.tsv")
    master_df.to_csv(master_out, sep="\t")

    loss_rows.append({
        "species": "ALL", "step": "C_reclassify",
        "input": sum(original_types.values()),
        "kept": sum(new_types.values()),
        "removed": n_reclassified,
    })

    # ------------------------------------------------------------------
    # STEP D: Write filtered BEDs, then resize to uniform width
    # ------------------------------------------------------------------
    if verbose:
        print(f"\n--- Step D: Resize all peaks to {target_size} bp ---")

    # D1. Unified peaks hg38 (before resize, for reference)
    uni_hg38_pre = os.path.join(intermediates, "unified_peaks_hg38_pre_resize.bed")
    n_uni = 0
    with open(uni_hg38_pre, 'w') as fout:
        for pid, row in master_df[master_df["peak_type"] == "unified"].iterrows():
            fout.write(f"{row['Human_chr']}\t{int(row['Human_start'])}\t"
                       f"{int(row['Human_end'])}\t{pid}\n")
            n_uni += 1

    # D2. Human-specific peaks (before resize)
    hs_hg38_pre = os.path.join(intermediates, "human_specific_hg38_pre_resize.bed")
    n_hs = 0
    with open(hs_hg38_pre, 'w') as fout:
        for pid, row in master_df[master_df["peak_type"] == "human_specific"].iterrows():
            fout.write(f"{row['Human_chr']}\t{int(row['Human_start'])}\t"
                       f"{int(row['Human_end'])}\t{pid}\n")
            n_hs += 1

    # D3. Per-species unified (from reciprocal pass or size-filtered, before resize)
    sp_uni_pre = {}
    for species in species_list:
        if reciprocal_check:
            info = reciprocal_results.get(species, {})
            pf = info.get("pass_file")
            if pf and os.path.exists(pf):
                sp_uni_pre[species] = pf
        else:
            info = size_results.get(species, {})
            sf = info.get("output_file")
            if sf and os.path.exists(sf):
                sp_uni_pre[species] = sf

    # D4. Per-species species-specific (before resize)
    sp_specific_pre = {}
    for species in species_list:
        sp_bed = os.path.join(sp_specific_dir, f"{species}_specific_peaks.bed")
        if os.path.exists(sp_bed):
            sp_specific_pre[species] = sp_bed

    # --- Resize everything ---
    # Human unified + human_specific (hg38)
    human_np = narrowpeak_files.get("Human")
    human_cs = chromsizes_files.get("Human")

    resize_results = {}

    # Unified hg38
    uni_hg38_final = os.path.join(final_dir, "unified_peaks_hg38.bed")
    if verbose:
        print(f"\n  Resizing unified hg38 ({n_uni:,} peaks):")
    resize_results["unified_hg38"] = resize_peaks(
        uni_hg38_pre, uni_hg38_final, target_size,
        narrowpeak_file=human_np, chromsizes_file=human_cs,
        verbose=verbose,
    )

    # Human-specific hg38
    hs_hg38_final = os.path.join(final_dir, "human_specific_peaks_hg38.bed")
    if verbose:
        print(f"\n  Resizing human-specific hg38 ({n_hs:,} peaks):")
    resize_results["hs_hg38"] = resize_peaks(
        hs_hg38_pre, hs_hg38_final, target_size,
        narrowpeak_file=human_np, chromsizes_file=human_cs,
        verbose=verbose,
    )

    # Per-species unified (NHP coords)
    for species, pf in sp_uni_pre.items():
        sp_final = os.path.join(final_dir, f"unified_peaks_{species}.bed")
        sp_np = narrowpeak_files.get(species)
        sp_cs = chromsizes_files.get(species)
        if verbose:
            print(f"\n  Resizing unified {species}:")
        resize_results[f"unified_{species}"] = resize_peaks(
            pf, sp_final, target_size,
            narrowpeak_file=sp_np, chromsizes_file=sp_cs,
            verbose=verbose,
        )

    # Per-species species-specific
    for species, sp_bed in sp_specific_pre.items():
        sp_final = os.path.join(final_dir, f"{species}_specific_peaks.bed")
        sp_np = narrowpeak_files.get(species)
        sp_cs = chromsizes_files.get(species)
        if verbose:
            print(f"\n  Resizing {species}-specific:")
        resize_results[f"specific_{species}"] = resize_peaks(
            sp_bed, sp_final, target_size,
            narrowpeak_file=sp_np, chromsizes_file=sp_cs,
            verbose=verbose,
        )

    # ------------------------------------------------------------------
    # STEP E: Combine into final all_peaks BEDs per species
    # ------------------------------------------------------------------
    if verbose:
        print(f"\n--- Step E: Write final all_peaks BEDs ---")

    # Human: unified + human_specific
    human_all = os.path.join(final_dir, "all_peaks_Human.bed")
    n_human = 0
    with open(human_all, 'w') as fout:
        for src in [uni_hg38_final, hs_hg38_final]:
            if os.path.exists(src):
                with open(src) as fin:
                    for line in fin:
                        if line.strip():
                            fout.write(line)
                            n_human += 1
    if verbose:
        print(f"  all_peaks_Human.bed: {n_human:,}")

    # NHP: unified + species-specific
    for species in species_list:
        sp_all = os.path.join(final_dir, f"all_peaks_{species}.bed")
        n = 0
        with open(sp_all, 'w') as fout:
            uni_f = os.path.join(final_dir, f"unified_peaks_{species}.bed")
            if os.path.exists(uni_f):
                with open(uni_f) as fin:
                    for line in fin:
                        if line.strip():
                            fout.write(line)
                            n += 1
            spc_f = os.path.join(final_dir, f"{species}_specific_peaks.bed")
            if os.path.exists(spc_f):
                with open(spc_f) as fin:
                    for line in fin:
                        if line.strip():
                            fout.write(line)
                            n += 1
        if verbose:
            print(f"  all_peaks_{species}.bed: {n:,}")

    # ------------------------------------------------------------------
    # STEP F: Rescue missing peaks from original input
    # ------------------------------------------------------------------
    rescued_counts = {}
    if species_beds or human_bed:
        if verbose:
            print(f"\n--- Step F: Rescue peaks lost during filtering ---")

        # Build dict of all original inputs and all final outputs
        original_beds = {}  # species -> original input BED
        final_beds = {}     # species -> final all_peaks BED

        if human_bed and os.path.exists(human_bed):
            original_beds["Human"] = human_bed
            final_beds["Human"] = os.path.join(final_dir, "all_peaks_Human.bed")

        if species_beds:
            for sp, bed in species_beds.items():
                if os.path.exists(bed):
                    original_beds[sp] = bed
                    final_beds[sp] = os.path.join(final_dir, f"all_peaks_{sp}.bed")

        with tempfile.TemporaryDirectory() as tmpdir:
            for sp in original_beds:
                orig = original_beds[sp]
                final = final_beds.get(sp)
                if not final or not os.path.exists(final):
                    continue

                # Extract BED3 from original and sort
                orig_bed3 = os.path.join(tmpdir, f"{sp}_orig_bed3.bed")
                with open(orig) as fin, open(orig_bed3, 'w') as fout:
                    for line in fin:
                        if line.strip() and not line.startswith('#'):
                            parts = line.strip().split('\t')
                            fout.write(f"{parts[0]}\t{parts[1]}\t{parts[2]}\n")

                orig_sorted = os.path.join(tmpdir, f"{sp}_orig_sorted.bed")
                subprocess.run(
                    f"sort -k1,1 -k2,2n {orig_bed3} > {orig_sorted}",
                    shell=True, check=True,
                )

                # Extract BED3 from final and sort/merge
                final_bed3 = os.path.join(tmpdir, f"{sp}_final_bed3.bed")
                with open(final) as fin, open(final_bed3, 'w') as fout:
                    for line in fin:
                        if line.strip() and not line.startswith('#'):
                            parts = line.strip().split('\t')
                            fout.write(f"{parts[0]}\t{parts[1]}\t{parts[2]}\n")

                final_sorted = os.path.join(tmpdir, f"{sp}_final_sorted.bed")
                final_merged = os.path.join(tmpdir, f"{sp}_final_merged.bed")
                subprocess.run(
                    f"sort -k1,1 -k2,2n {final_bed3} > {final_sorted}",
                    shell=True, check=True,
                )
                subprocess.run(
                    f"bedtools merge -i {final_sorted} > {final_merged}",
                    shell=True, check=True,
                )

                # Find original peaks NOT overlapping any final peak
                missing_raw = os.path.join(tmpdir, f"{sp}_missing.bed")
                subprocess.run(
                    f"bedtools intersect -v -a {orig_sorted} -b {final_merged} > {missing_raw}",
                    shell=True, check=True,
                )

                n_missing = sum(1 for line in open(missing_raw) if line.strip())
                if n_missing == 0:
                    if verbose:
                        print(f"  {sp}: no missing peaks — nothing to rescue")
                    rescued_counts[sp] = 0
                    continue

                # Resize rescued peaks to target_size and assign IDs
                sp_lower = sp.lower()
                rescued_bed = os.path.join(final_dir, f"{sp}_rescued_peaks.bed")
                sp_cs = chromsizes_files.get(sp)
                cs_lookup = {}
                if sp_cs and os.path.exists(sp_cs):
                    with open(sp_cs) as f:
                        for line in f:
                            parts = line.strip().split('\t')
                            if len(parts) >= 2:
                                cs_lookup[parts[0]] = int(parts[1])

                with open(missing_raw) as fin, open(rescued_bed, 'w') as fout:
                    for i, line in enumerate(fin, 1):
                        parts = line.strip().split('\t')
                        chrom = parts[0]
                        start = int(parts[1])
                        end = int(parts[2])
                        mid = (start + end) // 2
                        new_start = max(0, mid - target_size // 2)
                        new_end = new_start + target_size
                        # Clip to chromosome bounds
                        if chrom in cs_lookup:
                            if new_end > cs_lookup[chrom]:
                                new_end = cs_lookup[chrom]
                                new_start = max(0, new_end - target_size)
                        peak_id = f"{sp_lower}_rescued_{i:06d}"
                        fout.write(f"{chrom}\t{new_start}\t{new_end}\t{peak_id}\n")

                # Append rescued peaks to the all_peaks file
                with open(final, 'a') as fout:
                    with open(rescued_bed) as fin:
                        for line in fin:
                            if line.strip():
                                fout.write(line)

                rescued_counts[sp] = n_missing
                loss_rows.append({
                    "species": sp, "step": "F_rescued",
                    "input": n_missing, "kept": n_missing,
                    "removed": 0,
                })

                if verbose:
                    with open(final) as f:
                        n_final_total = sum(1 for l in f if l.strip())
                    print(f"  {sp}: rescued {n_missing:,} peaks "
                          f"(all_peaks now {n_final_total:,})")

    # ------------------------------------------------------------------
    # Loss summary
    # ------------------------------------------------------------------
    loss_df = pd.DataFrame(loss_rows)
    loss_file = os.path.join(final_dir, "loss_summary.tsv")
    loss_df.to_csv(loss_file, sep="\t", index=False)

    if verbose:
        print(f"\n{'=' * 70}")
        print("LOSS SUMMARY")
        print("=" * 70)
        print(loss_df.to_string(index=False))
        print(f"\n  Saved: {loss_file}")
        print(f"  Final BEDs in: {final_dir}")
        print(f"  Master annotation: {master_out}")

    return {
        "loss_summary": loss_df,
        "master_df": master_df,
        "master_file": master_out,
        "final_dir": final_dir,
        "size_filter_results": size_results,
        "reciprocal_results": reciprocal_results,
        "resize_results": resize_results,
        "rescued_counts": rescued_counts,
    }


# =============================================================================
# ANNOTATION FILE CREATION
# =============================================================================

def create_peak_annotation(
    unified_bed: str,
    human_specific_bed: str,
    species_specific_beds: Dict[str, str],
    gtf_files: Dict[str, str],
    output_dir: str,
    verbose: bool = True,
) -> Dict[str, Any]:
    """
    Create a unified peak annotation file for all peak types.

    Columns: peak_id, chr, start, end, peak_type, species_detected,
             closest_gene, distance_to_gene

    Args:
        unified_bed: Unified consensus BED (chr, start, end, peak_id, species_list)
        human_specific_bed: Human-specific peaks BED (chr, start, end, peak_id)
        species_specific_beds: Dict of species -> species-specific BED
        gtf_files: Dict of species -> GTF file path
        output_dir: Output directory for annotation files

    Returns:
        Dictionary with file paths and statistics
    """
    if verbose:
        print("\n" + "=" * 70)
        print("CREATING PEAK ANNOTATION")
        print("=" * 70)

    os.makedirs(output_dir, exist_ok=True)
    gene_bed_dir = os.path.join(output_dir, "gene_beds")
    os.makedirs(gene_bed_dir, exist_ok=True)

    results: Dict[str, Any] = {"annotations": {}}

    # --- 1. Annotate unified peaks (hg38 coordinates) ---
    if verbose:
        print("\n1. Annotating unified peaks (hg38)...")

    unified_rows = []
    with open(unified_bed) as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                parts = line.strip().split('\t')
                chrom, start, end = parts[0], parts[1], parts[2]
                peak_id = parts[3]
                species_list = parts[4] if len(parts) > 4 else "Human"
                unified_rows.append({
                    "peak_id": peak_id,
                    "chr": chrom,
                    "start": int(start),
                    "end": int(end),
                    "peak_type": "unified",
                    "species_detected": species_list,
                })

    # Get closest gene for unified peaks (hg38 = Human GTF)
    human_gene_bed = os.path.join(gene_bed_dir, "Human_genes.bed")
    if "Human" in gtf_files and os.path.exists(gtf_files["Human"]):
        extract_gene_bed_from_gtf(gtf_files["Human"], human_gene_bed, verbose=verbose)

        unified_gene_tsv = os.path.join(output_dir, "unified_gene_annotation.tsv")
        gene_df = annotate_with_closest_gene(unified_bed, human_gene_bed, unified_gene_tsv, verbose=verbose)

        if not gene_df.empty:
            gene_map = dict(zip(gene_df["peak_id"], zip(gene_df["closest_gene"], gene_df["distance_to_gene"])))
            for row in unified_rows:
                if row["peak_id"] in gene_map:
                    row["closest_gene"] = gene_map[row["peak_id"]][0]
                    row["distance_to_gene"] = gene_map[row["peak_id"]][1]
                else:
                    row["closest_gene"] = "."
                    row["distance_to_gene"] = -1
        else:
            for row in unified_rows:
                row["closest_gene"] = "."
                row["distance_to_gene"] = -1
    else:
        if verbose:
            print("   WARNING: Human GTF not found, skipping gene annotation for unified peaks")
        for row in unified_rows:
            row["closest_gene"] = "."
            row["distance_to_gene"] = -1

    # --- 2. Annotate human-specific peaks ---
    if verbose:
        print("\n2. Annotating human-specific peaks...")

    human_specific_rows = []
    if os.path.exists(human_specific_bed):
        with open(human_specific_bed) as f:
            for line in f:
                if line.strip() and not line.startswith('#'):
                    parts = line.strip().split('\t')
                    human_specific_rows.append({
                        "peak_id": parts[3] if len(parts) > 3 else f"{parts[0]}:{parts[1]}-{parts[2]}",
                        "chr": parts[0],
                        "start": int(parts[1]),
                        "end": int(parts[2]),
                        "peak_type": "human_specific",
                        "species_detected": "Human",
                    })

        # Closest gene (same Human GTF)
        if os.path.exists(human_gene_bed):
            hs_gene_tsv = os.path.join(output_dir, "human_specific_gene_annotation.tsv")
            gene_df = annotate_with_closest_gene(human_specific_bed, human_gene_bed, hs_gene_tsv, verbose=verbose)

            if not gene_df.empty:
                gene_map = dict(zip(gene_df["peak_id"], zip(gene_df["closest_gene"], gene_df["distance_to_gene"])))
                for row in human_specific_rows:
                    if row["peak_id"] in gene_map:
                        row["closest_gene"] = gene_map[row["peak_id"]][0]
                        row["distance_to_gene"] = gene_map[row["peak_id"]][1]
                    else:
                        row["closest_gene"] = "."
                        row["distance_to_gene"] = -1
            else:
                for row in human_specific_rows:
                    row["closest_gene"] = "."
                    row["distance_to_gene"] = -1
        else:
            for row in human_specific_rows:
                row["closest_gene"] = "."
                row["distance_to_gene"] = -1

    # --- 3. Annotate species-specific peaks ---
    if verbose:
        print("\n3. Annotating species-specific peaks...")

    species_specific_rows = {}
    for species, sp_bed in species_specific_beds.items():
        species_specific_rows[species] = []
        if not os.path.exists(sp_bed):
            continue

        rows = []
        with open(sp_bed) as f:
            for line in f:
                if line.strip() and not line.startswith('#'):
                    parts = line.strip().split('\t')
                    rows.append({
                        "peak_id": parts[3] if len(parts) > 3 else f"{parts[0]}:{parts[1]}-{parts[2]}",
                        "chr": parts[0],
                        "start": int(parts[1]),
                        "end": int(parts[2]),
                        "peak_type": f"{species.lower()}_specific",
                        "species_detected": species,
                    })

        # Closest gene for this species
        if species in gtf_files and os.path.exists(gtf_files[species]):
            sp_gene_bed = os.path.join(gene_bed_dir, f"{species}_genes.bed")
            extract_gene_bed_from_gtf(gtf_files[species], sp_gene_bed, verbose=verbose)

            sp_gene_tsv = os.path.join(output_dir, f"{species}_specific_gene_annotation.tsv")
            gene_df = annotate_with_closest_gene(sp_bed, sp_gene_bed, sp_gene_tsv, verbose=verbose)

            if not gene_df.empty:
                gene_map = dict(zip(gene_df["peak_id"], zip(gene_df["closest_gene"], gene_df["distance_to_gene"])))
                for row in rows:
                    if row["peak_id"] in gene_map:
                        row["closest_gene"] = gene_map[row["peak_id"]][0]
                        row["distance_to_gene"] = gene_map[row["peak_id"]][1]
                    else:
                        row["closest_gene"] = "."
                        row["distance_to_gene"] = -1
            else:
                for row in rows:
                    row["closest_gene"] = "."
                    row["distance_to_gene"] = -1
        else:
            if verbose:
                print(f"   WARNING: GTF not found for {species}, skipping gene annotation")
            for row in rows:
                row["closest_gene"] = "."
                row["distance_to_gene"] = -1

        species_specific_rows[species] = rows

    # --- 4. Combine into master annotation file ---
    if verbose:
        print("\n4. Writing master annotation file...")

    all_rows = unified_rows + human_specific_rows
    for species, rows in species_specific_rows.items():
        all_rows.extend(rows)

    annotation_df = pd.DataFrame(all_rows, columns=[
        "peak_id", "chr", "start", "end", "peak_type", "species_detected",
        "closest_gene", "distance_to_gene",
    ])

    annotation_file = os.path.join(output_dir, "peak_annotation.tsv")
    annotation_df.to_csv(annotation_file, sep='\t', index=False)

    if verbose:
        print(f"\n   Master annotation: {annotation_file}")
        print(f"   Total peaks annotated: {len(annotation_df):,}")
        print(f"   Peak types:")
        for pt, count in annotation_df["peak_type"].value_counts().items():
            print(f"      {pt}: {count:,}")

    results["annotation_file"] = annotation_file
    results["total_peaks"] = len(annotation_df)
    results["type_counts"] = annotation_df["peak_type"].value_counts().to_dict()

    return results


# =============================================================================
# MASTER ANNOTATION TABLE
# =============================================================================

def build_master_annotation(
    unified_bed: str,
    human_specific_bed: str,
    species_specific_beds: Dict[str, str],
    liftback_dir: str,
    gene_bed_dir: str,
    human_gene_annotation_tsv: str,
    species_list: List[str],
    output_file: str,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Build a comprehensive master annotation table with one row per peak.

    Peak classification (determined by **orthology**, not detection):

      - ``unified``: the orthologous genomic region exists in **≥2 species**
        (``n_species_orth >= 2``).  The peak may be *detected* (called) in
        only one species — that indicates regulatory divergence, not absence
        of the sequence.
      - ``human_specific``: the region exists **only in human**
        (``n_species_orth == 1``).  Liftback to every NHP species failed.
        This is symmetric with NHP ``{species}_specific`` peaks, which are
        original NHP peaks that could not be lifted to hg38.
      - ``{species}_specific``: original NHP peaks that could not be lifted
        to hg38 and therefore have no orthologous region in the human genome.

    Columns:
      - peak_id (index)
      - peak_type: ``'unified'``, ``'human_specific'``, ``'{species}_specific'``
      - Human_chr, Human_start, Human_end (hg38 coords for unified/hs peaks)
      - {Species}_chr, {Species}_start, {Species}_end (liftback coords)
      - Human_gene, Human_gene_dist
      - {Species}_gene, {Species}_gene_dist (from bedtools closest on liftback)
      - {Species}_det: binary 0/1 — a peak was *called* in that species
      - {Species}_orth: binary 0/1 — the orthologous genomic region *exists*
        in that species (liftback from hg38 succeeded)
      - n_species_det: count of species with ``_det=1``
      - n_species_orth: count of species with ``_orth=1``

    For species-specific peaks, only the native species coords are populated.

    Note:
        The ``human_specific_bed`` parameter is accepted for API
        compatibility but is **not used**.  Human-specific peaks are
        identified from the unified set based on liftback orthology
        (``n_species_orth == 1``), not from a separate BED file.

        Input peaks typically have temporary IDs (``temp_NNNNNN``).
        This function assigns **final** IDs (``unified_NNNNNN``,
        ``human_peak_NNNNNN``) based on the orthology classification.
        These final IDs are never changed downstream.
    """
    if verbose:
        print("Building master annotation table...")

    all_rows = []

    # ------------------------------------------------------------------
    # 1. Unified peaks -- hg38 coords + species detection
    # ------------------------------------------------------------------
    if verbose:
        print("  1. Loading unified peaks + human gene annotation...")

    # Load human gene annotation
    human_gene_map = {}
    if os.path.exists(human_gene_annotation_tsv):
        hg = pd.read_csv(human_gene_annotation_tsv, sep="\t")
        human_gene_map = dict(zip(
            hg["peak_id"],
            zip(hg["closest_gene"], hg["distance_to_gene"]),
        ))

    with open(unified_bed) as f:
        for line in f:
            if line.strip() and not line.startswith("#"):
                parts = line.strip().split("\t")
                pid = parts[3]
                species_detected = parts[4] if len(parts) > 4 else "Human"
                row = {
                    "peak_id": pid,
                    "peak_type": "unified",
                    "Human_chr": parts[0],
                    "Human_start": int(parts[1]),
                    "Human_end": int(parts[2]),
                    "species_detected": species_detected,
                }
                if pid in human_gene_map:
                    row["Human_gene"] = human_gene_map[pid][0]
                    row["Human_gene_dist"] = int(human_gene_map[pid][1])
                all_rows.append(row)

    # ------------------------------------------------------------------
    # 2. Human-specific classification (orthology-based)
    # ------------------------------------------------------------------
    # Human-specific peaks are NOT loaded from a separate BED file.
    # Instead, they are identified from the unified set in step 6
    # after liftback coordinates reveal which species have orthologous
    # regions.  Unified peaks with n_species_orth == 1 (liftback to
    # every NHP species failed) are reclassified as human_specific.
    # This is symmetric with NHP species-specific peaks, which are
    # original NHP peaks that could not be lifted to hg38.
    if verbose:
        print("  2. (human_specific peaks will be identified from "
              "unified set after liftback in step 6)")

    # ------------------------------------------------------------------
    # 3. Liftback coordinates for unified peaks
    # ------------------------------------------------------------------
    if verbose:
        print("  3. Loading liftback coordinates per species...")

    for species in species_list:
        lb_file = os.path.join(liftback_dir, f"unified_consensus_{species}.bed")
        if not os.path.exists(lb_file):
            continue

        # Build map: peak_id -> (chr, start, end) in species coords
        lb_map = {}
        with open(lb_file) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 4:
                    lb_map[parts[3]] = (parts[0], int(parts[1]), int(parts[2]))

        # Merge into rows
        for row in all_rows:
            if row["peak_type"] == "unified" and row["peak_id"] in lb_map:
                c, s, e = lb_map[row["peak_id"]]
                row[f"{species}_chr"] = c
                row[f"{species}_start"] = s
                row[f"{species}_end"] = e

        if verbose:
            print(f"    {species}: {len(lb_map):,} liftback coords")

    # ------------------------------------------------------------------
    # 4. Species gene annotation on liftback coords
    # ------------------------------------------------------------------
    if verbose:
        print("  4. Annotating liftback peaks with species genes...")

    for species in species_list:
        lb_file = os.path.join(liftback_dir, f"unified_consensus_{species}.bed")
        gene_bed = os.path.join(gene_bed_dir, f"{species}_genes.bed")
        if not os.path.exists(lb_file) or not os.path.exists(gene_bed):
            continue

        with tempfile.TemporaryDirectory() as tmpdir:
            # Harmonize chr prefix
            with open(lb_file) as f:
                lb_chr = f.readline().split("\t")[0]
            with open(gene_bed) as f:
                gene_chr = f.readline().split("\t")[0]

            gene_bed_use = gene_bed
            if lb_chr.startswith("chr") and not gene_chr.startswith("chr"):
                gene_bed_use = os.path.join(tmpdir, "genes_chr.bed")
                with open(gene_bed) as fin, open(gene_bed_use, "w") as fout:
                    for line in fin:
                        if line.strip():
                            fout.write("chr" + line)
            elif not lb_chr.startswith("chr") and gene_chr.startswith("chr"):
                gene_bed_use = os.path.join(tmpdir, "genes_nochr.bed")
                with open(gene_bed) as fin, open(gene_bed_use, "w") as fout:
                    for line in fin:
                        if line.strip():
                            fout.write(line[3:] if line.startswith("chr") else line)

            lb_sorted = os.path.join(tmpdir, "lb_sorted.bed")
            gene_sorted = os.path.join(tmpdir, "gene_sorted.bed")
            subprocess.run(f"sort -k1,1 -k2,2n {lb_file} > {lb_sorted}", shell=True, check=True)
            subprocess.run(f"sort -k1,1 -k2,2n {gene_bed_use} > {gene_sorted}", shell=True, check=True)

            closest_out = os.path.join(tmpdir, "closest.tsv")
            subprocess.run(
                f"bedtools closest -a {lb_sorted} -b {gene_sorted} -d -t first > {closest_out}",
                shell=True, check=True,
            )

            with open(lb_sorted) as f:
                n_lb_cols = len(f.readline().strip().split("\t"))

            sp_gene_map = {}
            with open(closest_out) as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) < n_lb_cols + 7:
                        continue
                    pid = parts[3]
                    gene_name = parts[n_lb_cols + 3]
                    dist = int(parts[-1]) if parts[-1] != "." else -1
                    sp_gene_map[pid] = (gene_name, dist)

        for row in all_rows:
            if row["peak_id"] in sp_gene_map:
                row[f"{species}_gene"] = sp_gene_map[row["peak_id"]][0]
                row[f"{species}_gene_dist"] = sp_gene_map[row["peak_id"]][1]

        if verbose:
            print(f"    {species}: {len(sp_gene_map):,} gene annotations")

    # ------------------------------------------------------------------
    # 5. Species-specific peaks
    # ------------------------------------------------------------------
    if verbose:
        print("  5. Loading species-specific peaks + gene annotation...")

    for species, sp_bed in species_specific_beds.items():
        if not os.path.exists(sp_bed):
            continue

        sp_gene_tsv = os.path.join(
            os.path.dirname(human_gene_annotation_tsv),
            f"{species}_specific_gene_annotation.tsv",
        )
        sp_gene_map = {}
        if os.path.exists(sp_gene_tsv):
            sg = pd.read_csv(sp_gene_tsv, sep="\t")
            sp_gene_map = dict(zip(
                sg["peak_id"],
                zip(sg["closest_gene"], sg["distance_to_gene"]),
            ))

        with open(sp_bed) as f:
            for line in f:
                if line.strip() and not line.startswith("#"):
                    parts = line.strip().split("\t")
                    pid = parts[3] if len(parts) > 3 else f"{parts[0]}:{parts[1]}-{parts[2]}"
                    row = {
                        "peak_id": pid,
                        "peak_type": f"{species.lower()}_specific",
                        f"{species}_chr": parts[0],
                        f"{species}_start": int(parts[1]),
                        f"{species}_end": int(parts[2]),
                        "species_detected": species,
                    }
                    if pid in sp_gene_map:
                        row[f"{species}_gene"] = sp_gene_map[pid][0]
                        row[f"{species}_gene_dist"] = int(sp_gene_map[pid][1])
                    all_rows.append(row)

    # ------------------------------------------------------------------
    # 6. Build DataFrame and add detection + orthology columns
    # ------------------------------------------------------------------
    if verbose:
        print("  6. Building DataFrame with detection + orthology columns...")

    all_species = ["Human"] + species_list
    df = pd.DataFrame(all_rows)

    # --- Binary detection columns from species_detected ---
    # These reflect which species had a CALLED PEAK overlapping this
    # consensus peak during the cross-species merge step.
    for sp in all_species:
        col = f"{sp}_det"
        df[col] = df["species_detected"].str.contains(sp, na=False).astype(int)

    df["n_species_det"] = df[[f"{sp}_det" for sp in all_species]].sum(axis=1)

    # --- Binary orthology columns from liftback coordinates ---
    # These reflect whether the orthologous region EXISTS in each species
    # (i.e., the liftback from hg38 succeeded). A peak with _orth=1 but
    # _det=0 means the genomic region is conserved but no peak was called
    # there in that species (regulatory divergence).
    # Human_orth=1 for all peaks that are defined in hg38 (currently all
    # have peak_type="unified"; human_specific is assigned below).
    # For NHP species, _orth=1 iff liftback coords exist.
    df["Human_orth"] = (df["peak_type"] == "unified").astype(int)

    for sp in species_list:
        chr_col = f"{sp}_chr"
        df[f"{sp}_orth"] = df[chr_col].notna().astype(int) if chr_col in df.columns else 0

    orth_cols = [f"{sp}_orth" for sp in all_species]
    df["n_species_orth"] = df[[c for c in orth_cols if c in df.columns]].sum(axis=1)

    # --- Reclassify unified peaks with no NHP ortholog as human_specific ---
    # A unified peak with n_species_orth == 1 means the liftback to EVERY
    # NHP species failed — the genomic region only exists in the human
    # genome.  This is the exact analogue of NHP species-specific peaks
    # (original NHP peaks that could not be lifted to hg38).
    mask_human_only = (df["peak_type"] == "unified") & (df["n_species_orth"] == 1)
    n_reclassified = mask_human_only.sum()
    df.loc[mask_human_only, "peak_type"] = "human_specific"

    if verbose:
        print(f"    Reclassified {n_reclassified:,} unified peaks as "
              f"human_specific (n_species_orth == 1)")
        print(f"    Remaining unified: "
              f"{(df['peak_type'] == 'unified').sum():,} "
              f"(n_species_orth >= 2)")

    # --- Assign final peak IDs based on classification ---
    # Intermediate files use temp_NNNNNN IDs (or whatever peak_prefix was).
    # Now that we know the final classification, assign clean sequential
    # IDs that will NEVER change downstream:
    #   unified        -> unified_000001, unified_000002, ...
    #   human_specific -> human_peak_000001, human_peak_000002, ...
    #   {sp}_specific  -> keep existing IDs (already correctly named)
    id_mapping = {}  # old_id -> new_id
    unified_counter = 0
    hs_counter = 0

    for idx in df.index:
        old_id = df.at[idx, "peak_id"]
        pt = df.at[idx, "peak_type"]
        if pt == "unified":
            unified_counter += 1
            id_mapping[old_id] = f"unified_{unified_counter:06d}"
        elif pt == "human_specific":
            hs_counter += 1
            id_mapping[old_id] = f"human_peak_{hs_counter:06d}"
        else:
            id_mapping[old_id] = old_id  # species-specific: keep as-is

    df["peak_id"] = df["peak_id"].map(id_mapping)

    # Save mapping so downstream notebook can update BED files
    mapping_file = output_file.replace(".tsv", "_id_mapping.tsv")
    pd.DataFrame(
        list(id_mapping.items()), columns=["old_id", "new_id"],
    ).to_csv(mapping_file, sep="\t", index=False)

    if verbose:
        n_changed = sum(1 for k, v in id_mapping.items() if k != v)
        print(f"    Renamed {n_changed:,} peak IDs "
              f"({unified_counter:,} unified, {hs_counter:,} human_specific)")
        print(f"    Saved ID mapping: {os.path.basename(mapping_file)}")

    if verbose:
        # Show the discrepancy between det and orth
        for sp in all_species:
            det_col = f"{sp}_det"
            orth_col = f"{sp}_orth"
            if det_col in df.columns and orth_col in df.columns:
                n_det = df[det_col].sum()
                n_orth = df[orth_col].sum()
                n_orth_not_det = ((df[orth_col] == 1) & (df[det_col] == 0)).sum()
                print(f"    {sp}: det={n_det:,}, orth={n_orth:,} "
                      f"(orth but not det: {n_orth_not_det:,})")

    # Reorder columns for clarity
    base_cols = ["peak_id", "peak_type", "species_detected",
                 "n_species_det", "n_species_orth"]
    det_cols = [f"{sp}_det" for sp in all_species]
    orth_cols_ordered = [f"{sp}_orth" for sp in all_species]

    coord_cols = []
    gene_cols = []
    for sp in all_species:
        coord_cols.extend([f"{sp}_chr", f"{sp}_start", f"{sp}_end"])
        gene_cols.extend([f"{sp}_gene", f"{sp}_gene_dist"])

    ordered = base_cols + det_cols + orth_cols_ordered + coord_cols + gene_cols
    # Only keep columns that actually exist
    ordered = [c for c in ordered if c in df.columns]
    remaining = [c for c in df.columns if c not in ordered]
    df = df[ordered + remaining]

    df = df.set_index("peak_id")

    # ------------------------------------------------------------------
    # 7. Save
    # ------------------------------------------------------------------
    df.to_csv(output_file, sep="\t")

    if verbose:
        print(f"\n  Master annotation saved: {output_file}")
        print(f"  Total peaks: {len(df):,}")
        print(f"  Columns: {len(df.columns)}")
        print(f"  Peak types: {df['peak_type'].value_counts().to_dict()}")
        print(f"\n  Detection summary (peak called = _det):")
        for sp in all_species:
            col = f"{sp}_det"
            if col in df.columns:
                print(f"    {sp}: {df[col].sum():,} peaks")
        print(f"\n  Orthology summary (liftback exists = _orth):")
        for sp in all_species:
            col = f"{sp}_orth"
            if col in df.columns:
                print(f"    {sp}: {df[col].sum():,} regions")

    return df


# =============================================================================
# CROSS-MAPPING SPECIES-SPECIFIC PEAKS VIA DIRECT CHAINS
# =============================================================================

# Direct inter-species liftover routes (avoiding hg38, since species-specific
# peaks already failed the hg38 round-trip).
# Each route is a list of chain file basenames to apply sequentially.
# Assembly mapping: Bonobo=panPan2, Chimpanzee=panTro5, Gorilla=gorGor4,
#                   Macaque=rheMac10, Marmoset=calJac1
CROSS_SPECIES_ROUTES = {
    # --- Bonobo (panPan2) -> others ---
    ("Bonobo", "Chimpanzee"):  ["panPan2ToPanTro4.over.chain", "panTro4ToPanTro5.over.chain"],
    ("Bonobo", "Gorilla"):     ["panPan2ToGorGor5.over.chain", "gorGor5ToGorGor4.over.chain"],
    ("Bonobo", "Macaque"):     ["panPan2ToPanTro4.over.chain", "panTro4ToPanTro5.over.chain",
                                "panTro5ToRheMac8.over.chain", "rheMac8ToRheMac10.over.chain"],
    ("Bonobo", "Marmoset"):    ["panPan2ToPanTro4.over.chain", "panTro4ToPanTro5.over.chain",
                                "panTro5ToRheMac8.over.chain", "rheMac8ToRheMac10.over.chain",
                                "rheMac10ToCalJac4.over.chain", "calJac4ToCalJac1.over.chain"],

    # --- Chimpanzee (panTro5) -> others ---
    ("Chimpanzee", "Bonobo"):  ["panTro5ToPanTro4.over.chain", "panTro4ToPanPan2.over.chain"],
    ("Chimpanzee", "Gorilla"): ["panTro5ToGorGor5.over.chain", "gorGor5ToGorGor4.over.chain"],
    ("Chimpanzee", "Macaque"): ["panTro5ToRheMac8.over.chain", "rheMac8ToRheMac10.over.chain"],
    ("Chimpanzee", "Marmoset"):["panTro5ToRheMac8.over.chain", "rheMac8ToRheMac10.over.chain",
                                "rheMac10ToCalJac4.over.chain", "calJac4ToCalJac1.over.chain"],

    # --- Gorilla (gorGor4) -> others ---
    ("Gorilla", "Bonobo"):     ["gorGor4ToGorGor5.over.chain", "gorGor5ToPanPan2.over.chain"],
    ("Gorilla", "Chimpanzee"): ["gorGor4ToGorGor5.over.chain", "gorGor5ToPanTro5.over.chain"],
    ("Gorilla", "Macaque"):    ["gorGor4ToGorGor5.over.chain", "gorGor5ToPanTro5.over.chain",
                                "panTro5ToRheMac8.over.chain", "rheMac8ToRheMac10.over.chain"],
    ("Gorilla", "Marmoset"):   ["gorGor4ToGorGor5.over.chain", "gorGor5ToPanTro5.over.chain",
                                "panTro5ToRheMac8.over.chain", "rheMac8ToRheMac10.over.chain",
                                "rheMac10ToCalJac4.over.chain", "calJac4ToCalJac1.over.chain"],

    # --- Macaque (rheMac10) -> others ---
    ("Macaque", "Chimpanzee"): ["rheMac10ToPanTro6.over.chain", "panTro6ToPanTro5.over.chain"],
    ("Macaque", "Bonobo"):     ["rheMac10ToPanTro6.over.chain", "panTro6ToPanTro5.over.chain",
                                "panTro5ToPanTro4.over.chain", "panTro4ToPanPan2.over.chain"],
    ("Macaque", "Gorilla"):    ["rheMac10ToPanTro6.over.chain", "panTro6ToPanTro5.over.chain",
                                "panTro5ToGorGor5.over.chain", "gorGor5ToGorGor4.over.chain"],
    ("Macaque", "Marmoset"):   ["rheMac10ToCalJac4.over.chain", "calJac4ToCalJac1.over.chain"],

    # --- Marmoset (calJac1) -> others ---
    ("Marmoset", "Macaque"):   ["calJac1ToCalJac4.over.chain", "calJac4ToRheMac10.over.chain"],
    ("Marmoset", "Chimpanzee"):["calJac1ToCalJac4.over.chain", "calJac4ToRheMac10.over.chain",
                                "rheMac10ToPanTro6.over.chain", "panTro6ToPanTro5.over.chain"],
    ("Marmoset", "Bonobo"):    ["calJac1ToCalJac4.over.chain", "calJac4ToRheMac10.over.chain",
                                "rheMac10ToPanTro6.over.chain", "panTro6ToPanTro5.over.chain",
                                "panTro5ToPanTro4.over.chain", "panTro4ToPanPan2.over.chain"],
    ("Marmoset", "Gorilla"):   ["calJac1ToCalJac4.over.chain", "calJac4ToRheMac10.over.chain",
                                "rheMac10ToPanTro6.over.chain", "panTro6ToPanTro5.over.chain",
                                "panTro5ToGorGor5.over.chain", "gorGor5ToGorGor4.over.chain"],
}


def _multi_step_liftover(
    input_bed: str,
    output_bed: str,
    chain_files: List[str],
    liftover_path: str,
    min_match: float = 0.5,
) -> Dict[str, Any]:
    """
    Chain multiple sequential liftOver calls.

    Each step's output becomes the next step's input.
    Returns the final lifted count and a per-step summary.
    """
    n_input = sum(1 for line in open(input_bed) if line.strip() and not line.startswith("#"))
    current_bed = input_bed
    step_summary = []

    with tempfile.TemporaryDirectory() as tmpdir:
        for i, chain in enumerate(chain_files):
            is_last = (i == len(chain_files) - 1)
            step_out = output_bed if is_last else os.path.join(tmpdir, f"step_{i}.bed")
            step_unmapped = os.path.join(tmpdir, f"step_{i}.unmapped.bed")

            result = liftover_peaks(
                input_bed=current_bed,
                output_bed=step_out,
                chain_file=chain,
                liftover_path=liftover_path,
                min_match=min_match,
                auto_chr=True,
                verbose=False,
                ncpu=1,
            )

            lifted = result.get("lifted", 0)
            step_summary.append({
                "step": i + 1,
                "chain": os.path.basename(chain),
                "input": sum(1 for line in open(current_bed) if line.strip() and not line.startswith("#")) if os.path.exists(current_bed) else 0,
                "lifted": lifted,
            })

            if lifted == 0:
                # Nothing left to lift; ensure output exists
                if not is_last:
                    Path(output_bed).touch()
                break

            current_bed = step_out

    n_output = sum(1 for line in open(output_bed) if line.strip() and not line.startswith("#")) if os.path.exists(output_bed) else 0
    return {
        "status": "success",
        "input": n_input,
        "lifted": n_output,
        "steps": step_summary,
    }


def _cross_map_one_pair(
    source: str,
    target: str,
    src_bed: str,
    n_source: int,
    chain_paths: List[str],
    chain_basenames: List[str],
    target_orig: str,
    liftover_path: str,
    min_match_val: float,
) -> Dict[str, Any]:
    """
    Cross-map one (source, target) pair: multi-step liftover + intersection.

    Designed to be called in parallel via ProcessPoolExecutor.
    All paths are resolved before submission so no closures are needed.

    Returns a dict with the cross-mapping result row.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        lifted_bed = os.path.join(tmpdir, f"lifted_to_{target}.bed")
        lift_result = _multi_step_liftover(
            input_bed=src_bed,
            output_bed=lifted_bed,
            chain_files=chain_paths,
            liftover_path=liftover_path,
            min_match=min_match_val,
        )

        n_lifted = lift_result.get("lifted", 0)
        step_summary = lift_result.get("steps", [])

        if n_lifted == 0:
            return {
                "source": source,
                "target": target,
                "source_specific": n_source,
                "n_hops": len(chain_basenames),
                "route": " -> ".join(chain_basenames),
                "lifted_to_target": 0,
                "overlap_target_peaks": 0,
                "pct_overlap": 0.0,
                "steps": step_summary,
            }

        # Harmonize chr prefix between lifted and target
        def _first_chr(path):
            with open(path) as fh:
                for ln in fh:
                    if ln.strip():
                        return ln.split("\t")[0]
            return ""

        lifted_chr = _first_chr(lifted_bed)
        target_chr = _first_chr(target_orig)

        target_use = target_orig
        if lifted_chr.startswith("chr") != target_chr.startswith("chr"):
            target_fixed = os.path.join(tmpdir, f"target_{target}_fixed.bed")
            with open(target_orig) as fin, open(target_fixed, "w") as fout:
                for line in fin:
                    if line.strip() and not line.startswith("#"):
                        if lifted_chr.startswith("chr") and not line.startswith("chr"):
                            fout.write("chr" + line)
                        elif not lifted_chr.startswith("chr") and line.startswith("chr"):
                            fout.write(line[3:])
                        else:
                            fout.write(line)
            target_use = target_fixed

        # Sort and intersect
        lifted_sorted = os.path.join(tmpdir, "lifted_sorted.bed")
        target_sorted = os.path.join(tmpdir, "target_sorted.bed")
        subprocess.run(f"sort -k1,1 -k2,2n {lifted_bed} > {lifted_sorted}",
                       shell=True, check=True, capture_output=True)
        subprocess.run(f"sort -k1,1 -k2,2n {target_use} > {target_sorted}",
                       shell=True, check=True, capture_output=True)

        overlap_out = os.path.join(tmpdir, "overlap.bed")
        subprocess.run(
            f"bedtools intersect -a {lifted_sorted} -b {target_sorted} -u > {overlap_out}",
            shell=True, check=True, capture_output=True,
        )
        n_overlap = sum(1 for _ in open(overlap_out))
        pct = n_overlap / n_source * 100 if n_source > 0 else 0

        return {
            "source": source,
            "target": target,
            "source_specific": n_source,
            "n_hops": len(chain_basenames),
            "route": " -> ".join(chain_basenames),
            "lifted_to_target": n_lifted,
            "overlap_target_peaks": n_overlap,
            "pct_overlap": pct,
            "steps": step_summary,
        }


def cross_map_species_specific_peaks(
    species_specific_beds: Dict[str, str],
    species_beds: Dict[str, str],
    chain_dir: str,
    liftover_path: str,
    output_dir: str,
    min_match: Union[float, Dict[str, float]] = 0.5,
    verbose: bool = True,
    ncpu: int = 1,
) -> pd.DataFrame:
    """
    Map each species' specific peaks to every other species via direct
    inter-species chain files (multi-step, avoiding hg38).

    Species-specific peaks are those that failed the hg38 round-trip,
    so routing them through hg38 again is unlikely to succeed. Instead,
    we use direct assembly-to-assembly chains (e.g. panPan2 -> gorGor5 ->
    gorGor4) downloaded from UCSC.

    The routes are defined in CROSS_SPECIES_ROUTES. Each route is a list
    of chain files applied sequentially. Great apes are closely related
    and typically need only 2 hops. Distant pairs (e.g. Gorilla->Marmoset)
    require more hops via intermediate assemblies.

    For each pair (source, target):
      1. Apply the multi-step liftover route
      2. Intersect the result with the target's original consensus peaks

    Args:
        species_specific_beds: Dict of species -> path to species-specific peaks
        species_beds: Dict of species -> path to original consensus peaks
        chain_dir: Chain file directory
        liftover_path: liftOver binary path
        output_dir: Where to write cross-map outputs
        min_match: Min match for liftover (float or per-species dict)
        verbose: Print progress
        ncpu: Number of parallel workers (default 1 = sequential).
              Each source-target pair is processed independently.

    Returns:
        DataFrame with cross-mapping summary.
    """
    from concurrent.futures import ProcessPoolExecutor, as_completed

    os.makedirs(output_dir, exist_ok=True)

    def _get_mm(species: str) -> float:
        if isinstance(min_match, dict):
            return min_match.get(species, 0.5)
        return min_match

    all_species = sorted(species_specific_beds.keys())

    if verbose:
        print("=" * 70)
        print("CROSS-MAPPING SPECIES-SPECIFIC PEAKS (DIRECT INTER-SPECIES CHAINS)")
        print("=" * 70)
        print("Strategy: direct multi-step liftover between assemblies")
        print("          (avoiding hg38, since these peaks failed the hg38 round-trip)")
        if ncpu > 1:
            print(f"⚙️  Using {ncpu} parallel workers")
        print()

    # ----- Build list of jobs -----
    jobs = []       # (source, target, src_bed, n_source, chain_paths, chain_basenames, target_orig, mm)
    skip_rows = []  # results for pairs we can skip immediately

    source_counts = {}
    for source in all_species:
        src_bed = species_specific_beds[source]
        if not os.path.exists(src_bed):
            continue
        n_source = sum(1 for line in open(src_bed) if line.strip() and not line.startswith("#"))
        if n_source == 0:
            continue
        source_counts[source] = n_source

        for target in all_species:
            if target == source:
                continue

            route_key = (source, target)
            if route_key not in CROSS_SPECIES_ROUTES:
                skip_rows.append({
                    "source": source, "target": target,
                    "source_specific": n_source, "n_hops": 0,
                    "route": "none", "lifted_to_target": 0,
                    "overlap_target_peaks": 0, "pct_overlap": 0.0,
                })
                continue

            chain_basenames = CROSS_SPECIES_ROUTES[route_key]
            chain_paths = [os.path.join(chain_dir, c) for c in chain_basenames]

            missing = [c for c in chain_paths if not os.path.exists(c)]
            if missing:
                if verbose:
                    print(f"  {source} -> {target}: MISSING chain files: "
                          f"{[os.path.basename(m) for m in missing]}")
                skip_rows.append({
                    "source": source, "target": target,
                    "source_specific": n_source,
                    "n_hops": len(chain_basenames),
                    "route": " -> ".join(chain_basenames),
                    "lifted_to_target": 0, "overlap_target_peaks": 0,
                    "pct_overlap": 0.0,
                })
                continue

            target_orig = species_beds.get(target, "")
            if not os.path.exists(target_orig):
                continue

            jobs.append((
                source, target, src_bed, n_source,
                chain_paths, chain_basenames, target_orig, _get_mm(source),
            ))

    if verbose:
        print(f"  {len(jobs)} liftover jobs to run, {len(skip_rows)} skipped")

    # ----- Execute jobs (parallel or sequential) -----
    results_rows = list(skip_rows)

    if ncpu > 1 and len(jobs) > 1:
        n_workers = min(ncpu, len(jobs))
        if verbose:
            print(f"  🚀 Submitting {len(jobs)} jobs to {n_workers} workers...")

        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            futures = {}
            for src, tgt, src_bed, n_src, c_paths, c_names, tgt_orig, mm in jobs:
                fut = executor.submit(
                    _cross_map_one_pair,
                    src, tgt, src_bed, n_src,
                    c_paths, c_names, tgt_orig,
                    liftover_path, mm,
                )
                futures[fut] = (src, tgt)

            for fut in as_completed(futures):
                source, target = futures[fut]
                try:
                    result = fut.result()
                    results_rows.append(result)
                    if verbose:
                        n_src = result["source_specific"]
                        n_hop = result["n_hops"]
                        n_lft = result["lifted_to_target"]
                        n_ovl = result["overlap_target_peaks"]
                        pct = result["pct_overlap"]
                        steps_str = ""
                        for s in result.get("steps", []):
                            steps_str += (f"\n      step {s['step']}: "
                                          f"{s['chain']}  {s['input']}->{s['lifted']}")
                        print(f"  {source} -> {target} ({n_hop} hops): "
                              f"{n_lft:,} lifted, {n_ovl:,} overlap ({pct:.1f}%)"
                              f"{steps_str}")
                except Exception as exc:
                    print(f"  {source} -> {target}: ERROR: {exc}")
                    results_rows.append({
                        "source": source, "target": target,
                        "source_specific": source_counts.get(source, 0),
                        "n_hops": 0, "route": "error",
                        "lifted_to_target": 0, "overlap_target_peaks": 0,
                        "pct_overlap": 0.0,
                    })
    else:
        # Sequential execution
        for src, tgt, src_bed_, n_src_, c_paths_, c_names_, tgt_orig_, mm_ in jobs:
            source, target = src, tgt
            if verbose:
                print(f"  {source} -> {target} ...", end="", flush=True)

            result = _cross_map_one_pair(
                src, tgt, src_bed_, n_src_,
                c_paths_, c_names_, tgt_orig_,
                liftover_path, mm_,
            )
            results_rows.append(result)

            if verbose:
                n_hop = result["n_hops"]
                n_lft = result["lifted_to_target"]
                n_ovl = result["overlap_target_peaks"]
                pct = result["pct_overlap"]
                steps_str = ""
                for s in result.get("steps", []):
                    steps_str += (f"\n      step {s['step']}: "
                                  f"{s['chain']}  {s['input']}->{s['lifted']}")
                print(f" ({n_hop} hops): "
                      f"{n_lft:,} lifted, {n_ovl:,} overlap ({pct:.1f}%)"
                      f"{steps_str}")

    # ----- Build summary -----
    # Remove internal 'steps' key before making DataFrame
    for row in results_rows:
        row.pop("steps", None)

    cross_df = pd.DataFrame(results_rows)

    # Save
    out_file = os.path.join(output_dir, "species_specific_cross_mapping.tsv")
    cross_df.to_csv(out_file, sep="\t", index=False)

    # Also create a summary matrix
    if len(cross_df) > 0:
        matrix = cross_df.pivot_table(
            index="source", columns="target",
            values="pct_overlap", fill_value=0,
        )
        matrix_file = os.path.join(output_dir, "cross_mapping_matrix_pct.tsv")
        matrix.to_csv(matrix_file, sep="\t")

        overlap_matrix = cross_df.pivot_table(
            index="source", columns="target",
            values="overlap_target_peaks", fill_value=0,
        )
        overlap_file = os.path.join(output_dir, "cross_mapping_matrix_counts.tsv")
        overlap_matrix.to_csv(overlap_file, sep="\t")

        if verbose:
            print(f"\n{'='*70}")
            print(f"CROSS-MAPPING SUMMARY (% of source-specific peaks overlapping target)")
            print(f"{'='*70}")
            print(matrix.round(1).to_string())
            print(f"\nSaved: {out_file}")
            print(f"Saved: {matrix_file}")
            print(f"Saved: {overlap_file}")

    return cross_df


# =============================================================================
# LEGACY: merge_bed_files (kept for backwards compatibility)
# =============================================================================

def merge_bed_files(
    input_beds: List[str],
    output_bed: str,
    merge_distance: int = 0,
    verbose: bool = True,
) -> Dict[str, Any]:
    """Legacy simple merge without species tracking.

    .. deprecated::
        Use ``merge_with_species_tracking()`` or ``summit_based_merge()``
        instead.  Kept for backwards compatibility only.
    """
    import warnings
    warnings.warn(
        "merge_bed_files() is deprecated. Use merge_with_species_tracking() "
        "or summit_based_merge() instead.",
        DeprecationWarning, stacklevel=2,
    )
    if verbose:
        print(f"Merging {len(input_beds)} BED files...")

    total_input = 0
    for bed in input_beds:
        with open(bed) as f:
            total_input += sum(1 for line in f if line.strip() and not line.startswith('#'))

    os.makedirs(os.path.dirname(output_bed), exist_ok=True)

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_concat = os.path.join(tmpdir, "concat.bed")
        with open(tmp_concat, 'w') as fout:
            for bed in input_beds:
                with open(bed) as f:
                    for line in f:
                        if line.strip() and not line.startswith('#'):
                            parts = line.strip().split('\t')
                            if len(parts) >= 3:
                                fout.write(f"{parts[0]}\t{parts[1]}\t{parts[2]}\n")

        tmp_sorted = os.path.join(tmpdir, "sorted.bed")
        subprocess.run(f"sort -k1,1 -k2,2n {tmp_concat} > {tmp_sorted}", shell=True, check=True)
        subprocess.run(f"bedtools merge -i {tmp_sorted} -d {merge_distance} > {output_bed}", shell=True, check=True)

    with open(output_bed) as f:
        output_count = sum(1 for line in f if line.strip())

    return {
        "status": "success",
        "input_peaks": total_input,
        "merged_peaks": output_count,
        "output_file": output_bed,
    }


# =============================================================================
# MAIN PIPELINE
# =============================================================================

def cross_species_consensus_pipeline(
    species_beds: Dict[str, str],
    human_bed: str,
    output_dir: str,
    chain_dir: str = DEFAULT_CHAIN_DIR,
    liftover_path: str = "liftOver",
    min_match: Union[float, Dict[str, float]] = 0.95,
    min_blocks: Optional[float] = None,
    multiple: bool = False,
    merge_distance: int = 0,
    merge_method: str = "summit",
    summit_window: int = 500,
    cluster_distance: int = 250,
    max_peak_size: int = 2000,
    min_peak_size: int = 100,
    score_col: int = 4,
    normalize_scores: bool = True,
    peak_prefix: str = "temp",
    max_liftback_size: Optional[int] = None,
    min_liftback_size: Optional[int] = None,
    reciprocal_check: bool = False,
    max_reciprocal_distance: int = 500,
    compute_conservation: bool = False,
    gtf_files: Optional[Dict[str, str]] = None,
    pre_lifted_beds: Optional[Dict[str, str]] = None,
    verbose: bool = True,
    ncpu: int = 1,
) -> Dict[str, Any]:
    """
    Complete cross-species peak comparison pipeline with species tracking.

    Steps:
    1. Liftover all non-human species peaks to hg38 (or load pre-lifted)
    2. Merge lifted peaks + human peaks (assigned temp IDs)
    3. (skipped — human-specific classification via orthology in Step 7)
    4. Liftover temp consensus peaks back to each species
       - 4b. (optional) Size-filter liftback peaks
       - 4c. (optional) Reciprocal liftover check
       - 4d. (optional) Per-peak liftover conservation scoring
    5. Identify species-specific peaks for each non-human species
    6. Annotate all peaks with closest gene
    7. Build master annotation — reclassify peaks using orthology:
       ``temp_NNNNNN`` → ``unified_NNNNNN`` (n_species_orth ≥ 2)
       or ``human_peak_NNNNNN`` (liftback to all NHP failed)
    8. Cross-map species-specific peaks via direct inter-species chains

    Peak IDs:
    - Intermediate files use ``temp_NNNNNN`` (or custom ``peak_prefix``).
    - Final IDs are assigned in Step 7 and are **never changed** after that.
    - Species-specific peaks use ``{species}_peak_NNNNNN`` throughout.

    Args:
        species_beds: Dict mapping non-human species name to input BED file path
        human_bed: Path to human consensus peaks BED file (hg38)
        output_dir: Output directory for all results
        chain_dir: Directory containing chain files
        liftover_path: Path to liftOver executable
        min_match: Minimum match ratio for liftover. Either a single float
            applied to all species, or a dict mapping species name to its
            match ratio (e.g. {'Bonobo': 0.9, 'Macaque': 0.8}).
        min_blocks: Minimum ratio of alignment blocks that must map
        multiple: Allow multiple output regions
        merge_distance: Distance for merging peaks (only for interval method)
        merge_method: ``"summit"`` (default) uses summit-based clustering
            with fixed-width output windows. ``"interval"`` uses the legacy
            bedtools merge approach on full intervals.
        summit_window: Width of output consensus peaks (summit method only)
        cluster_distance: Max distance between summits to cluster (summit only)
        max_peak_size: Filter out lifted peaks larger than this (summit only)
        min_peak_size: Filter out lifted peaks smaller than this (summit only)
        score_col: 0-indexed column for peak scores (summit only, default: 4)
        normalize_scores: Rank-normalize scores per species (summit only)
        peak_prefix: Prefix for temporary peak IDs used in intermediate
            files.  Final IDs (``unified_NNNNNN``, ``human_peak_NNNNNN``) are
            assigned in Step 7 after orthology-based classification.
        max_liftback_size: If set, discard liftback peaks larger than this
            (bp) after Step 4.  Recommended: 2000–5000.  *None* = no filter.
        min_liftback_size: If set, discard liftback peaks smaller than this
            (bp) after Step 4.  Recommended: 100.  *None* = no filter.
        reciprocal_check: If True, run a reciprocal liftover check
            (Step 4c) — lift NHP liftback coords back to hg38 and only
            keep peaks that overlap the original hg38 position.
        max_reciprocal_distance: Max distance (bp) between original and
            round-trip hg38 positions for a peak to pass the reciprocal
            check (default 500).
        compute_conservation: If True, run per-peak liftover conservation
            scoring (Step 4d) — round-trip liftover hg38 → species → hg38
            measuring match_ratio (fraction of bases recovered). Results
            are saved to ``10_similarity/`` and attached to ``results``.
        gtf_files: Dict mapping species name to GTF file path (optional)
        pre_lifted_beds: Optional dict mapping species name to path of
            pre-computed liftover BED files (already in hg38 coords).
            If provided, Step 1 is skipped for those species.
        verbose: Print detailed progress
        ncpu: Number of parallel workers

    Returns:
        Dictionary with all results and file paths
    """
    if gtf_files is None:
        gtf_files = DEFAULT_GTF_FILES.copy()

    results = {
        "lift_to_human": {},
        "merge": None,
        "lift_back": {},
        "species_specific": {},
        "output_files": {},
    }

    # Create output directories
    lift_human_dir = os.path.join(output_dir, "01_lifted_to_human")
    merged_dir = os.path.join(output_dir, "02_merged_consensus")
    lift_back_dir = os.path.join(output_dir, "04_lifted_back")
    species_specific_dir = os.path.join(output_dir, "05_species_specific")
    annotation_dir = os.path.join(output_dir, "06_annotation")

    for d in [lift_human_dir, merged_dir, lift_back_dir,
              species_specific_dir, annotation_dir]:
        os.makedirs(d, exist_ok=True)

    print("=" * 70)
    print(f"CROSS-SPECIES CONSENSUS PIPELINE (v2 - {merge_method} merge)")
    if ncpu > 1:
        print(f"Using {ncpu} parallel workers")
    print("=" * 70)

    # =========================================================================
    # STEP 1: Liftover all non-human species to hg38 (or load pre-lifted)
    # =========================================================================
    print("\nSTEP 1: Liftover non-human species -> hg38")
    print("-" * 50)

    lifted_beds = {}  # species -> lifted BED path

    # Helper: resolve per-species min_match
    def _get_min_match(species: str) -> float:
        if isinstance(min_match, dict):
            return min_match.get(species, 0.95)
        return min_match

    for species, input_bed in species_beds.items():
        # Check if pre-lifted file was provided
        if pre_lifted_beds and species in pre_lifted_beds:
            pre_lifted = pre_lifted_beds[species]
            if os.path.exists(pre_lifted):
                print(f"   {species}: using pre-lifted file ({pre_lifted})")
                lifted_beds[species] = pre_lifted
                n_lifted = sum(1 for _ in open(pre_lifted))
                n_original = sum(1 for l in open(input_bed) if l.strip() and not l.startswith('#'))
                results["lift_to_human"][species] = {
                    "status": "success",
                    "original": n_original,
                    "lifted": n_lifted,
                    "unmapped": n_original - n_lifted,
                    "source": "pre_lifted",
                }
                continue
            else:
                print(f"   {species}: pre-lifted file not found, running liftover")

        if not os.path.exists(input_bed):
            print(f"   Skipping {species} - file not found: {input_bed}")
            continue

        output_bed = os.path.join(lift_human_dir, f"{species}_hg38.bed")
        species_mm = _get_min_match(species)

        if verbose:
            print(f"\n   {species} (min_match={species_mm}):")

        if species == "Marmoset":
            chain1 = get_chain_file("Marmoset_step1", chain_dir)
            chain2 = get_chain_file("Marmoset_step2", chain_dir)
            result = liftover_two_step(
                input_bed=input_bed, output_bed=output_bed,
                chain_file_1=chain1, chain_file_2=chain2,
                liftover_path=liftover_path, min_match=species_mm,
                min_blocks=min_blocks, multiple=multiple,
                auto_chr=True, verbose=verbose, ncpu=ncpu,
            )
        else:
            chain_file = get_chain_file(species, chain_dir)
            result = liftover_peaks(
                input_bed=input_bed, output_bed=output_bed,
                chain_file=chain_file, liftover_path=liftover_path,
                min_match=species_mm, min_blocks=min_blocks,
                multiple=multiple, auto_chr=True, verbose=verbose, ncpu=ncpu,
            )

        results["lift_to_human"][species] = result

        if result["status"] == "success" and result["lifted"] > 0:
            lifted_beds[species] = output_bed

    if not lifted_beds:
        return {
            "status": "error",
            "message": "No peaks were successfully lifted to hg38",
            **results,
        }

    # =========================================================================
    # STEP 2: Merge all species (including human) with species tracking
    # =========================================================================
    print("\n" + "=" * 70)
    print(f"STEP 2: Merge all species with species tracking (method={merge_method})")
    print("-" * 50)

    # Build dict of all species BEDs in hg38 coordinates
    all_species_hg38 = {"Human": human_bed}
    all_species_hg38.update(lifted_beds)

    merged_bed = os.path.join(merged_dir, "unified_consensus_hg38_merged.bed")

    if merge_method == "summit":
        merge_result = summit_based_merge(
            species_beds=all_species_hg38,
            output_bed=merged_bed,
            summit_window=summit_window,
            cluster_distance=cluster_distance,
            max_peak_size=max_peak_size,
            min_peak_size=min_peak_size,
            score_col=score_col,
            normalize_scores=normalize_scores,
            verbose=verbose,
        )
    elif merge_method == "interval":
        merge_result = merge_with_species_tracking(
            species_beds=all_species_hg38,
            output_bed=merged_bed,
            merge_distance=merge_distance,
            verbose=verbose,
        )
    else:
        raise ValueError(f"Unknown merge_method: {merge_method!r}. "
                         f"Use 'summit' or 'interval'.")
    results["merge"] = merge_result

    # Add peak IDs (preserve species column)
    merged_with_ids = os.path.join(merged_dir, "unified_consensus_hg38_with_ids.bed")
    id_result = add_peak_ids(
        input_bed=merged_bed,
        output_bed=merged_with_ids,
        prefix=peak_prefix,
        verbose=verbose,
    )
    results["output_files"]["unified_consensus"] = merged_with_ids

    # =========================================================================
    # STEP 3: (skipped) Human-specific peaks identified in Step 7
    # =========================================================================
    # Human-specific peak classification is done in Step 7
    # (build_master_annotation) using liftback orthology: peaks where the
    # liftback to EVERY NHP species failed (n_species_orth == 1) are
    # classified as human_specific.  This is the correct definition —
    # regions that cannot be transferred to other genomes.
    #
    # The old Step 3 used detection-based overlap (bedtools intersect -v)
    # which conflated "no peak called" with "no orthologous sequence".
    print("\n" + "=" * 70)
    print("STEP 3: (skipped — human-specific identified via orthology in Step 7)")
    print("-" * 50)
    print("   Human-specific peaks will be identified after liftback (Step 7).")

    # =========================================================================
    # STEP 4: Liftover unified peaks back to each species
    # =========================================================================
    print("\n" + "=" * 70)
    print("STEP 4: Liftover unified peaks back to species genomes")
    print("-" * 50)

    # Create a BED4 (no species column) for liftback to avoid extra column issues
    unified_bed4 = os.path.join(merged_dir, "unified_consensus_hg38_bed4.bed")
    with open(merged_with_ids) as fin, open(unified_bed4, 'w') as fout:
        for line in fin:
            parts = line.strip().split('\t')
            # Write only chr, start, end, peak_id
            fout.write(f"{parts[0]}\t{parts[1]}\t{parts[2]}\t{parts[3]}\n")

    for species in species_beds.keys():
        if species not in results["lift_to_human"]:
            continue
        if results["lift_to_human"][species]["status"] != "success":
            continue

        output_bed = os.path.join(lift_back_dir, f"unified_consensus_{species}.bed")

        result = liftback_peaks(
            input_bed=unified_bed4,
            output_bed=output_bed,
            species=species,
            chain_dir=chain_dir,
            liftover_path=liftover_path,
            min_match=_get_min_match(species),
            min_blocks=min_blocks,
            multiple=multiple,
            auto_chr=True,
            verbose=verbose,
            ncpu=ncpu,
        )

        results["lift_back"][species] = result
        results["output_files"][f"liftback_{species}"] = output_bed

    # =========================================================================
    # STEP 4b: (optional) Filter liftback peaks by size
    # =========================================================================
    if max_liftback_size is not None or min_liftback_size is not None:
        _min = min_liftback_size if min_liftback_size is not None else 0
        _max = max_liftback_size if max_liftback_size is not None else 999_999_999
        print("\n" + "=" * 70)
        print(f"STEP 4b: Filter liftback peaks by size [{_min:,} – {_max:,} bp]")
        print("-" * 50)

        filter_dir = os.path.join(output_dir, "04_lifted_back_filtered")
        filter_results = filter_liftback_by_size(
            liftback_dir=lift_back_dir,
            output_dir=filter_dir,
            species_list=list(species_beds.keys()),
            min_liftback_size=_min,
            max_liftback_size=_max,
            verbose=verbose,
        )
        results["liftback_filter"] = filter_results

        # Update liftback output paths to use filtered files
        for species, info in filter_results.items():
            filtered_bed = info.get("output_file")
            if filtered_bed and os.path.exists(filtered_bed):
                results["output_files"][f"liftback_{species}"] = filtered_bed
                # Overwrite original liftback dir file so downstream steps
                # (species_specific, annotation, master) use filtered coords
                import shutil
                dst = os.path.join(lift_back_dir, f"unified_consensus_{species}.bed")
                shutil.copy2(filtered_bed, dst)

    # =========================================================================
    # STEP 4c: (optional) Reciprocal liftover check
    # =========================================================================
    if reciprocal_check:
        print("\n" + "=" * 70)
        print(f"STEP 4c: Reciprocal liftover check (max_distance={max_reciprocal_distance} bp)")
        print("-" * 50)

        reciprocal_dir = os.path.join(output_dir, "04_reciprocal_check")
        reciprocal_results = reciprocal_liftover_check(
            hg38_bed=unified_bed4,
            liftback_dir=lift_back_dir,
            output_dir=reciprocal_dir,
            species_list=list(species_beds.keys()),
            chain_dir=chain_dir,
            liftover_path=liftover_path,
            min_match=min_match,
            max_distance=max_reciprocal_distance,
            verbose=verbose,
            ncpu=ncpu,
        )
        results["reciprocal_check"] = reciprocal_results

        # Overwrite liftback files with reciprocal-pass peaks only
        for species, info in reciprocal_results.items():
            pass_file = info.get("pass_file")
            if pass_file and os.path.exists(pass_file):
                dst = os.path.join(lift_back_dir,
                                   f"unified_consensus_{species}.bed")
                shutil.copy2(pass_file, dst)
                results["output_files"][f"liftback_{species}"] = dst

    # =========================================================================
    # STEP 4d: (optional) Per-peak liftover conservation scoring
    # =========================================================================
    if compute_conservation:
        print("\n" + "=" * 70)
        print("STEP 4d: Per-peak liftover conservation scoring (round-trip match ratio)")
        print("-" * 50)

        similarity_dir = os.path.join(output_dir, "10_similarity")
        os.makedirs(similarity_dir, exist_ok=True)
        similarity_tsv = os.path.join(similarity_dir,
                                      "liftover_similarity_all_species.tsv")

        similarity_df = compute_species_similarity(
            input_bed=unified_bed4,
            species_list=list(species_beds.keys()),
            chain_dir=chain_dir,
            liftover_path=liftover_path,
            output_tsv=similarity_tsv,
            verbose=verbose,
            ncpu=ncpu,
        )

        # Pivot to a peak × species match-ratio matrix
        match_matrix = similarity_df.pivot(
            index="peak_id", columns="species", values="match_ratio"
        )
        match_matrix["Human"] = 1.0
        col_order = ["Human"] + [s for s in ["Chimpanzee", "Bonobo",
                     "Gorilla", "Macaque", "Marmoset"]
                     if s in match_matrix.columns]
        match_matrix = match_matrix[
            [c for c in col_order if c in match_matrix.columns]
        ]

        # Attach coordinates
        first_sp = similarity_df["species"].iloc[0]
        coords = similarity_df[similarity_df["species"] == first_sp][
            ["peak_id", "chrom", "start", "end", "original_size"]
        ].set_index("peak_id")
        match_matrix = coords.join(match_matrix)

        matrix_tsv = os.path.join(similarity_dir, "match_ratio_matrix.tsv")
        match_matrix.to_csv(matrix_tsv, sep="\t")

        results["conservation"] = {
            "similarity_tsv": similarity_tsv,
            "matrix_tsv": matrix_tsv,
            "n_peaks": len(match_matrix),
            "per_species_median": {
                sp: float(match_matrix[sp].median())
                for sp in list(species_beds.keys())
                if sp in match_matrix.columns
            },
        }
        results["output_files"]["similarity_matrix"] = matrix_tsv
        results["output_files"]["similarity_detail"] = similarity_tsv

        if verbose:
            print(f"\n   Conservation scoring complete:")
            print(f"   Peaks scored: {len(match_matrix):,}")
            for sp, med in results["conservation"]["per_species_median"].items():
                print(f"      {sp}: median match_ratio = {med:.3f}")
            print(f"   Saved: {matrix_tsv}")

    # =========================================================================
    # STEP 5: Identify species-specific peaks
    # =========================================================================
    print("\n" + "=" * 70)
    print("STEP 5: Identify species-specific peaks")
    print("-" * 50)

    species_specific_beds = {}
    for species, original_bed in species_beds.items():
        liftback_file = os.path.join(lift_back_dir, f"unified_consensus_{species}.bed")

        if not os.path.exists(liftback_file):
            print(f"   Skipping {species} - no liftback file")
            continue

        sp_specific_bed = os.path.join(species_specific_dir, f"{species}_specific_peaks.bed")

        sp_result = find_species_specific_peaks(
            original_bed=original_bed,
            liftback_bed=liftback_file,
            output_bed=sp_specific_bed,
            species=species,
            verbose=verbose,
        )

        results["species_specific"][species] = sp_result
        results["output_files"][f"species_specific_{species}"] = sp_specific_bed
        species_specific_beds[species] = sp_specific_bed

    # =========================================================================
    # STEP 6: Create peak annotation
    # =========================================================================
    print("\n" + "=" * 70)
    print("STEP 6: Annotate peaks with closest gene")
    print("-" * 50)

    annotation_result = create_peak_annotation(
        unified_bed=merged_with_ids,
        human_specific_bed="",  # not available yet; classified in Step 7
        species_specific_beds=species_specific_beds,
        gtf_files=gtf_files,
        output_dir=annotation_dir,
        verbose=verbose,
    )
    results["annotation"] = annotation_result
    results["output_files"]["annotation"] = os.path.join(annotation_dir, "peak_annotation.tsv")

    # =========================================================================
    # STEP 7: Master annotation table
    # =========================================================================
    print("\n" + "=" * 70)
    print("STEP 7: Building master annotation table")
    print("-" * 50)

    master_dir = os.path.join(output_dir, "07_master_annotation")
    os.makedirs(master_dir, exist_ok=True)
    master_file = os.path.join(master_dir, "master_annotation.tsv")

    master_df = build_master_annotation(
        unified_bed=merged_with_ids,
        human_specific_bed="",  # not used; orthology-based classification
        species_specific_beds=species_specific_beds,
        liftback_dir=lift_back_dir,
        gene_bed_dir=os.path.join(annotation_dir, "gene_beds"),
        human_gene_annotation_tsv=os.path.join(annotation_dir, "unified_gene_annotation.tsv"),
        species_list=list(species_beds.keys()),
        output_file=master_file,
        verbose=verbose,
    )
    results["master_annotation"] = master_df
    results["output_files"]["master_annotation"] = master_file

    # =========================================================================
    # STEP 8: Cross-map species-specific peaks
    # =========================================================================
    print("\n" + "=" * 70)
    print("STEP 8: Cross-mapping species-specific peaks (direct inter-species chains)")
    print("-" * 50)

    cross_map_dir = os.path.join(output_dir, "08_cross_mapping")
    cross_map_df = cross_map_species_specific_peaks(
        species_specific_beds=species_specific_beds,
        species_beds=species_beds,
        chain_dir=chain_dir,
        liftover_path=liftover_path,
        output_dir=cross_map_dir,
        min_match=min_match,
        verbose=verbose,
        ncpu=ncpu,
    )
    results["cross_mapping"] = cross_map_df
    results["output_files"]["cross_mapping"] = os.path.join(cross_map_dir, "species_specific_cross_mapping.tsv")

    # =========================================================================
    # SUMMARY (orthology-based classification from master annotation)
    # =========================================================================
    print("\n" + "=" * 70)
    print("PIPELINE SUMMARY")
    print("=" * 70)

    # Peak counts from master annotation (orthology-based)
    type_counts = master_df["peak_type"].value_counts()
    n_unified = type_counts.get("unified", 0)
    n_human_specific = type_counts.get("human_specific", 0)
    n_total_master = len(master_df)

    print(f"\n   Initial merged peaks (temp IDs): {merge_result['merged_peaks']:,}")
    print(f"   After orthology-based classification:")
    print(f"      Unified (n_species_orth >= 2):  {n_unified:,}")
    print(f"      Human-specific (no NHP orth):   {n_human_specific:,}")
    for pt, count in type_counts.items():
        if pt not in ("unified", "human_specific"):
            print(f"      {pt}: {count:,}")
    print(f"      Total: {n_total_master:,}")

    print(f"\n   Liftback results:")
    print(f"   {'Species':<15} {'Lifted':>10} {'Unmapped':>10} {'Success':>10}")
    print("   " + "-" * 50)
    for species, result in results["lift_back"].items():
        if "lifted" in result:
            lifted = result["lifted"]
            unmapped = result.get("unmapped", result.get("total_unmapped", 0))
            total = lifted + unmapped
            pct = (lifted / total * 100) if total > 0 else 0
            print(f"   {species:<15} {lifted:>10,} {unmapped:>10,} {pct:>9.1f}%")

    print(f"\n   Species-specific peaks (NHP, could not lift to hg38):")
    for species, sp_result in results["species_specific"].items():
        n = sp_result.get("species_specific", 0)
        pct = sp_result.get("pct_specific", 0)
        print(f"      {species}: {n:,} ({pct:.1f}%)")

    print(f"\n   Output files:")
    for key, filepath in results["output_files"].items():
        if filepath and os.path.exists(filepath):
            print(f"      {key}: {filepath}")

    results["status"] = "success"
    results["message"] = (
        f"Pipeline complete. "
        f"Unified: {n_unified:,}, "
        f"Human-specific (orthology): {n_human_specific:,}, "
        f"Total: {n_total_master:,}"
    )

    return results


# =============================================================================
# LEGACY: create_peak_matrix (kept for backwards compatibility)
# =============================================================================

def create_peak_matrix(
    unified_human_bed: str,
    species_beds: Dict[str, str],
    output_file: str,
    verbose: bool = True,
) -> Dict[str, Any]:
    """
    Create a peak presence/absence matrix across species.

    .. deprecated::
        Use the master annotation (``build_master_annotation()``) which provides
        per-species orthology columns.  Kept for backwards compatibility only.
    """
    import warnings
    warnings.warn(
        "create_peak_matrix() is deprecated. Use the master annotation "
        "(build_master_annotation()) which provides per-species orthology columns.",
        DeprecationWarning, stacklevel=2,
    )

    peaks = {}
    with open(unified_human_bed) as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    peak_id = parts[3]
                    peaks[peak_id] = {
                        "chr": parts[0],
                        "start": parts[1],
                        "end": parts[2],
                    }
                    if len(parts) >= 5:
                        peaks[peak_id]["species_detected"] = parts[4]

    species_list = sorted(species_beds.keys())
    for peak_id in peaks:
        for species in species_list:
            peaks[peak_id][species] = 0

    for species, bed_file in species_beds.items():
        if not os.path.exists(bed_file):
            continue
        with open(bed_file) as f:
            for line in f:
                if line.strip() and not line.startswith('#'):
                    parts = line.strip().split('\t')
                    if len(parts) >= 4:
                        peak_id = parts[3]
                        if peak_id in peaks:
                            peaks[peak_id][species] = 1

    with open(output_file, 'w') as fout:
        header = ["peak_id", "chr", "start", "end"] + species_list
        fout.write("\t".join(header) + "\n")
        for peak_id in sorted(peaks.keys()):
            row = [peak_id, peaks[peak_id]["chr"], peaks[peak_id]["start"], peaks[peak_id]["end"]]
            row += [str(peaks[peak_id][sp]) for sp in species_list]
            fout.write("\t".join(row) + "\n")

    n_peaks = len(peaks)
    conserved_all = sum(1 for p in peaks.values() if all(p[sp] == 1 for sp in species_list))
    conserved_any = sum(1 for p in peaks.values() if any(p[sp] == 1 for sp in species_list))

    if verbose:
        print(f"\nPeak Matrix Summary:")
        print(f"   Total unified peaks: {n_peaks:,}")
        print(f"   Present in all species: {conserved_all:,} ({conserved_all / n_peaks * 100:.1f}%)")
        print(f"   Present in >=1 species: {conserved_any:,} ({conserved_any / n_peaks * 100:.1f}%)")

    return {
        "status": "success",
        "total_peaks": n_peaks,
        "conserved_all": conserved_all,
        "conserved_any": conserved_any,
        "species": species_list,
        "output_file": output_file,
    }
