"""
Cross-Species Peak Analysis Pipeline
=====================================

Pipeline for comparing ATAC-seq peaks across species:
1. Liftover species peaks to human (hg38)
2. Merge into unified human consensus
3. Liftover back to original species with peak IDs for tracking

Author: J. Janssens
"""

import os
import subprocess
from pathlib import Path
from typing import Dict, List, Any, Optional
from collections import defaultdict

from .liftover import liftover_peaks, liftover_two_step, get_chain_file, DEFAULT_CHAIN_DIR

# Reverse chain files (hg38 → species)
# NOTE: These must match the assemblies used for the original consensus peaks
REVERSE_CHAIN_FILES = {
    "Gorilla": "hg38ToGorGor4.over.chain",      # hg38 → gorGor4
    "Chimpanzee": "hg38ToPanTro5.over.chain",   # hg38 → panTro5
    "Bonobo": "hg38ToPanPan2.over.chain",       # hg38 → panPan2 (matches panpan1.1)
    "Macaque": "hg38ToRheMac10.over.chain",     # hg38 → rheMac10
    "Marmoset_step1": "hg38ToCalJac4.over.chain",  # hg38 → calJac4
    "Marmoset_step2": "calJac4ToCalJac1.over.chain",  # calJac4 → calJac1
}


def get_reverse_chain_file(species: str, chain_dir: str = DEFAULT_CHAIN_DIR) -> str:
    """Get path to reverse chain file (hg38 → species)."""
    if species not in REVERSE_CHAIN_FILES:
        raise ValueError(f"Unknown species: {species}. Available: {list(REVERSE_CHAIN_FILES.keys())}")
    return os.path.join(chain_dir, REVERSE_CHAIN_FILES[species])


def merge_bed_files(
    input_beds: List[str],
    output_bed: str,
    merge_distance: int = 0,
    verbose: bool = True,
) -> Dict[str, Any]:
    """
    Merge multiple BED files into a single consensus BED file.
    
    Uses bedtools to concatenate and merge overlapping regions.
    
    Args:
        input_beds: List of input BED file paths
        output_bed: Path to output merged BED file
        merge_distance: Maximum distance between features to merge (default: 0 = overlapping only)
        verbose: Print progress
    
    Returns:
        Dictionary with merge statistics
    """
    import tempfile
    
    if verbose:
        print(f"🔗 Merging {len(input_beds)} BED files...")
    
    # Count input peaks
    total_input = 0
    for bed in input_beds:
        with open(bed) as f:
            total_input += sum(1 for line in f if line.strip() and not line.startswith('#'))
    
    if verbose:
        print(f"   Total input peaks: {total_input:,}")
    
    # Create output directory
    os.makedirs(os.path.dirname(output_bed), exist_ok=True)
    
    # Concatenate all files, sort, and merge
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as tmp:
        tmp_concat = tmp.name
        for bed in input_beds:
            with open(bed) as f:
                for line in f:
                    if line.strip() and not line.startswith('#'):
                        # Only take first 3 columns (chr, start, end)
                        parts = line.strip().split('\t')
                        if len(parts) >= 3:
                            tmp.write(f"{parts[0]}\t{parts[1]}\t{parts[2]}\n")
    
    try:
        # Sort
        with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as tmp_sorted_file:
            tmp_sorted = tmp_sorted_file.name
        
        sort_cmd = f"sort -k1,1 -k2,2n {tmp_concat} > {tmp_sorted}"
        subprocess.run(sort_cmd, shell=True, check=True)
        
        # Merge with bedtools
        merge_cmd = f"bedtools merge -i {tmp_sorted} -d {merge_distance} > {output_bed}"
        result = subprocess.run(merge_cmd, shell=True, capture_output=True, text=True)
        
        if result.returncode != 0:
            # Try without bedtools (simple Python merge)
            if verbose:
                print("   ⚠️  bedtools not found, using Python merge...")
            _python_merge(tmp_sorted, output_bed, merge_distance)
        
        # Count output peaks
        with open(output_bed) as f:
            output_count = sum(1 for line in f if line.strip())
        
        if verbose:
            print(f"   Merged peaks: {output_count:,}")
            print(f"   Reduction: {(1 - output_count/total_input)*100:.1f}%")
        
        return {
            "status": "success",
            "input_peaks": total_input,
            "merged_peaks": output_count,
            "output_file": output_bed,
            "message": f"✅ Merged {total_input:,} → {output_count:,} peaks"
        }
        
    finally:
        # Cleanup
        if os.path.exists(tmp_concat):
            os.unlink(tmp_concat)
        if 'tmp_sorted' in locals() and os.path.exists(tmp_sorted):
            os.unlink(tmp_sorted)


def _python_merge(input_bed: str, output_bed: str, merge_distance: int = 0):
    """Simple Python-based BED merge (fallback if bedtools unavailable)."""
    # Load all regions
    regions = defaultdict(list)
    with open(input_bed) as f:
        for line in f:
            if line.strip():
                parts = line.strip().split('\t')
                chrom, start, end = parts[0], int(parts[1]), int(parts[2])
                regions[chrom].append((start, end))
    
    # Merge overlapping regions per chromosome
    with open(output_bed, 'w') as fout:
        for chrom in sorted(regions.keys()):
            intervals = sorted(regions[chrom])
            merged = []
            
            for start, end in intervals:
                if merged and start <= merged[-1][1] + merge_distance:
                    # Overlaps with previous - extend
                    merged[-1] = (merged[-1][0], max(merged[-1][1], end))
                else:
                    # New interval
                    merged.append((start, end))
            
            for start, end in merged:
                fout.write(f"{chrom}\t{start}\t{end}\n")


def add_peak_ids(
    input_bed: str,
    output_bed: str,
    prefix: str = "peak",
    verbose: bool = True,
) -> Dict[str, Any]:
    """
    Add unique peak IDs to a BED file.
    
    Creates BED4 format: chr, start, end, peak_id
    
    Args:
        input_bed: Input BED file
        output_bed: Output BED file with peak IDs
        prefix: Prefix for peak IDs (default: "peak")
        verbose: Print progress
    
    Returns:
        Dictionary with statistics
    """
    peak_count = 0
    
    with open(input_bed) as fin, open(output_bed, 'w') as fout:
        for line in fin:
            if line.strip() and not line.startswith('#'):
                parts = line.strip().split('\t')
                chrom, start, end = parts[0], parts[1], parts[2]
                peak_count += 1
                peak_id = f"{prefix}_{peak_count:06d}"
                fout.write(f"{chrom}\t{start}\t{end}\t{peak_id}\n")
    
    if verbose:
        print(f"✅ Added IDs to {peak_count:,} peaks ({prefix}_000001 to {prefix}_{peak_count:06d})")
    
    return {
        "status": "success",
        "peak_count": peak_count,
        "output_file": output_bed,
        "message": f"✅ Added {peak_count:,} peak IDs"
    }


def liftback_peaks(
    input_bed: str,
    output_bed: str,
    species: str,
    chain_dir: str = DEFAULT_CHAIN_DIR,
    liftover_path: str = "liftOver",
    min_match: float = 0.95,
    min_blocks: float = None,
    multiple: bool = False,
    auto_chr: bool = True,
    verbose: bool = True,
    ncpu: int = 1,
) -> Dict[str, Any]:
    """
    Liftover peaks from hg38 back to a species genome.
    
    Preserves peak_id (4th column) if present.
    For Marmoset, performs two-step liftover (hg38 → calJac4 → calJac1).
    
    Args:
        input_bed: Input BED file (hg38 coordinates with peak IDs)
        output_bed: Output BED file (species coordinates)
        species: Target species name
        chain_dir: Directory containing chain files
        liftover_path: Path to liftOver executable
        min_match: Minimum match ratio
        min_blocks: Minimum ratio of alignment blocks that must map (optional)
        multiple: Allow multiple output regions for a single input (default False)
        auto_chr: Auto-fix chromosome prefix
        verbose: Print progress
        ncpu: Number of parallel workers (default 1)
    
    Returns:
        Dictionary with liftover results
    """
    if verbose:
        print(f"\n🔙 Lifting back to {species}...")
    
    if species == "Marmoset":
        # Two-step reverse: hg38 → calJac4 → calJac1
        chain1 = get_reverse_chain_file("Marmoset_step1", chain_dir)
        chain2 = get_reverse_chain_file("Marmoset_step2", chain_dir)
        
        result = liftover_two_step(
            input_bed=input_bed,
            output_bed=output_bed,
            chain_file_1=chain1,
            chain_file_2=chain2,
            liftover_path=liftover_path,
            min_match=min_match,
            min_blocks=min_blocks,
            multiple=multiple,
            auto_chr=auto_chr,
            verbose=verbose,
            ncpu=ncpu,
        )
    else:
        chain_file = get_reverse_chain_file(species, chain_dir)
        
        result = liftover_peaks(
            input_bed=input_bed,
            output_bed=output_bed,
            chain_file=chain_file,
            liftover_path=liftover_path,
            min_match=min_match,
            min_blocks=min_blocks,
            multiple=multiple,
            auto_chr=auto_chr,
            verbose=verbose,
            ncpu=ncpu,
        )
    
    return result


def cross_species_consensus_pipeline(
    species_beds: Dict[str, str],
    output_dir: str,
    chain_dir: str = DEFAULT_CHAIN_DIR,
    liftover_path: str = "liftOver",
    min_match: float = 0.95,
    min_blocks: float = None,
    multiple: bool = False,
    merge_distance: int = 0,
    peak_prefix: str = "unified",
    verbose: bool = True,
    ncpu: int = 1,
) -> Dict[str, Any]:
    """
    Complete cross-species peak comparison pipeline.
    
    Steps:
    1. Liftover all species peaks to hg38
    2. Merge into unified human consensus
    3. Add peak IDs for tracking
    4. Liftover back to each species
    
    Args:
        species_beds: Dict mapping species name to input BED file path
        output_dir: Output directory for all results
        chain_dir: Directory containing chain files
        liftover_path: Path to liftOver executable
        min_match: Minimum match ratio for liftover
        min_blocks: Minimum ratio of alignment blocks that must map (optional)
        multiple: Allow multiple output regions for a single input (default False)
        merge_distance: Distance for merging peaks (0 = overlapping only)
        peak_prefix: Prefix for peak IDs
        verbose: Print detailed progress
        ncpu: Number of parallel workers (default 1)
    
    Returns:
        Dictionary with all results and file paths
    """
    results = {
        "lift_to_human": {},
        "merge": None,
        "lift_back": {},
        "output_files": {},
    }
    
    # Create output directories
    lift_human_dir = os.path.join(output_dir, "01_lifted_to_human")
    merged_dir = os.path.join(output_dir, "02_merged_consensus")
    lift_back_dir = os.path.join(output_dir, "03_lifted_back")
    
    for d in [lift_human_dir, merged_dir, lift_back_dir]:
        os.makedirs(d, exist_ok=True)
    
    print("=" * 70)
    print("CROSS-SPECIES CONSENSUS PIPELINE")
    if ncpu > 1:
        print(f"⚙️  Using {ncpu} parallel workers")
    print("=" * 70)
    
    # =========================================================================
    # STEP 1: Liftover all species to human (hg38)
    # =========================================================================
    print("\n📍 STEP 1: Liftover species → hg38")
    print("-" * 50)
    
    lifted_beds = []
    
    for species, input_bed in species_beds.items():
        if not os.path.exists(input_bed):
            print(f"⏭️  Skipping {species} - file not found: {input_bed}")
            continue
        
        output_bed = os.path.join(lift_human_dir, f"{species}_hg38.bed")
        
        if verbose:
            print(f"\n🧬 {species}:")
        
        if species == "Marmoset":
            chain1 = get_chain_file("Marmoset_step1", chain_dir)
            chain2 = get_chain_file("Marmoset_step2", chain_dir)
            result = liftover_two_step(
                input_bed=input_bed,
                output_bed=output_bed,
                chain_file_1=chain1,
                chain_file_2=chain2,
                liftover_path=liftover_path,
                min_match=min_match,
                min_blocks=min_blocks,
                multiple=multiple,
                auto_chr=True,
                verbose=verbose,
                ncpu=ncpu,
            )
        else:
            chain_file = get_chain_file(species, chain_dir)
            result = liftover_peaks(
                input_bed=input_bed,
                output_bed=output_bed,
                chain_file=chain_file,
                liftover_path=liftover_path,
                min_match=min_match,
                min_blocks=min_blocks,
                multiple=multiple,
                auto_chr=True,
                verbose=verbose,
                ncpu=ncpu,
            )
        
        results["lift_to_human"][species] = result
        
        if result["status"] == "success" and result["lifted"] > 0:
            lifted_beds.append(output_bed)
    
    if not lifted_beds:
        return {
            "status": "error",
            "message": "❌ No peaks were successfully lifted to hg38",
            **results
        }
    
    # =========================================================================
    # STEP 2: Merge into unified human consensus
    # =========================================================================
    print("\n" + "=" * 70)
    print("📍 STEP 2: Merge into unified human consensus")
    print("-" * 50)
    
    merged_bed = os.path.join(merged_dir, "unified_consensus_hg38.bed")
    merge_result = merge_bed_files(
        input_beds=lifted_beds,
        output_bed=merged_bed,
        merge_distance=merge_distance,
        verbose=verbose,
    )
    results["merge"] = merge_result
    
    # =========================================================================
    # STEP 3: Add peak IDs
    # =========================================================================
    print("\n" + "=" * 70)
    print("📍 STEP 3: Add peak IDs for tracking")
    print("-" * 50)
    
    merged_with_ids = os.path.join(merged_dir, "unified_consensus_hg38_with_ids.bed")
    id_result = add_peak_ids(
        input_bed=merged_bed,
        output_bed=merged_with_ids,
        prefix=peak_prefix,
        verbose=verbose,
    )
    results["output_files"]["human_consensus"] = merged_with_ids
    
    # =========================================================================
    # STEP 4: Liftover back to each species
    # =========================================================================
    print("\n" + "=" * 70)
    print("📍 STEP 4: Liftover back to species genomes")
    print("-" * 50)
    
    for species in species_beds.keys():
        if species not in results["lift_to_human"]:
            continue
        if results["lift_to_human"][species]["status"] != "success":
            continue
        
        output_bed = os.path.join(lift_back_dir, f"unified_consensus_{species}.bed")
        
        result = liftback_peaks(
            input_bed=merged_with_ids,
            output_bed=output_bed,
            species=species,
            chain_dir=chain_dir,
            liftover_path=liftover_path,
            min_match=min_match,
            min_blocks=min_blocks,
            multiple=multiple,
            auto_chr=True,
            verbose=verbose,
            ncpu=ncpu,
        )
        
        results["lift_back"][species] = result
        results["output_files"][species] = output_bed
    
    # =========================================================================
    # SUMMARY
    # =========================================================================
    print("\n" + "=" * 70)
    print("PIPELINE SUMMARY")
    print("=" * 70)
    
    print(f"\n📊 Unified human consensus: {merge_result['merged_peaks']:,} peaks")
    print(f"   File: {merged_with_ids}")
    
    print(f"\n📊 Liftback results:")
    print(f"{'Species':<15} {'Lifted':>10} {'Unmapped':>10} {'Success':>10}")
    print("-" * 50)
    
    for species, result in results["lift_back"].items():
        if "lifted" in result:
            lifted = result["lifted"]
            unmapped = result.get("unmapped", result.get("total_unmapped", 0))
            total = lifted + unmapped
            pct = (lifted / total * 100) if total > 0 else 0
            print(f"{species:<15} {lifted:>10,} {unmapped:>10,} {pct:>9.1f}%")
    
    print("\n📁 Output files:")
    for key, filepath in results["output_files"].items():
        if os.path.exists(filepath):
            print(f"   {key}: {filepath}")
    
    results["status"] = "success"
    results["message"] = f"✅ Pipeline complete. Unified consensus: {merge_result['merged_peaks']:,} peaks"
    
    return results


def create_peak_matrix(
    unified_human_bed: str,
    species_beds: Dict[str, str],
    output_file: str,
    verbose: bool = True,
) -> Dict[str, Any]:
    """
    Create a peak presence/absence matrix across species.
    
    For each unified peak, check if it was successfully lifted to each species.
    
    Args:
        unified_human_bed: Human consensus BED with peak IDs (col 4)
        species_beds: Dict mapping species to lifted-back BED files
        output_file: Output TSV file for the matrix
    
    Returns:
        Dictionary with matrix statistics
    """
    # Load unified peaks
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
    
    # Initialize presence matrix
    species_list = sorted(species_beds.keys())
    for peak_id in peaks:
        for species in species_list:
            peaks[peak_id][species] = 0
    
    # Check presence in each species
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
    
    # Write matrix
    with open(output_file, 'w') as fout:
        # Header
        header = ["peak_id", "chr", "start", "end"] + species_list
        fout.write("\t".join(header) + "\n")
        
        # Data rows
        for peak_id in sorted(peaks.keys()):
            row = [
                peak_id,
                peaks[peak_id]["chr"],
                peaks[peak_id]["start"],
                peaks[peak_id]["end"],
            ]
            row += [str(peaks[peak_id][sp]) for sp in species_list]
            fout.write("\t".join(row) + "\n")
    
    # Calculate statistics
    n_peaks = len(peaks)
    conserved_all = sum(1 for p in peaks.values() if all(p[sp] == 1 for sp in species_list))
    conserved_any = sum(1 for p in peaks.values() if any(p[sp] == 1 for sp in species_list))
    
    if verbose:
        print(f"\n📊 Peak Matrix Summary:")
        print(f"   Total unified peaks: {n_peaks:,}")
        print(f"   Present in all species: {conserved_all:,} ({conserved_all/n_peaks*100:.1f}%)")
        print(f"   Present in ≥1 species: {conserved_any:,} ({conserved_any/n_peaks*100:.1f}%)")
        print(f"   Output: {output_file}")
    
    return {
        "status": "success",
        "total_peaks": n_peaks,
        "conserved_all": conserved_all,
        "conserved_any": conserved_any,
        "species": species_list,
        "output_file": output_file,
    }
