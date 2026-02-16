"""
Cross-Species Peak Analysis Pipeline
=====================================

Pipeline for comparing ATAC-seq peaks across species:
1. Liftover non-human species peaks to human (hg38)
2. Merge ALL species (including human) into unified consensus, tracking species origin
3. Identify human-specific peaks (not overlapping any non-human lifted peak)
4. Add peak IDs for tracking
5. Liftover unified peaks back to each species
6. Identify species-specific peaks (original species peaks not in unified set)
7. Annotate all peaks with closest gene from GTF

Author: J. Janssens
"""

import os
import re
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple
from collections import defaultdict

import pandas as pd

from .liftover import liftover_peaks, liftover_two_step, get_chain_file, DEFAULT_CHAIN_DIR

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
DEFAULT_GTF_FILES = {
    "Human": "/cluster/home/jjanssens/jjans/analysis/cerebellum/genomes_new/homo_sapiens/gencode.v48.basic.annotation.gtf.gz",
    "Bonobo": "/cluster/work/treutlein/jjans/data/intestine/nhp_atlas/genomes/reference_/panPan1/genes/genes.gtf",
    "Gorilla": "/cluster/work/treutlein/jjans/data/intestine/nhp_atlas/genomes/reference_/gorGor4/genes/genes.gtf",
    "Macaque": "/cluster/work/treutlein/jjans/data/intestine/nhp_atlas/genomes/reference_/Mmul10/genes/genes.gtf",
    "Marmoset": "/cluster/work/treutlein/jjans/data/intestine/nhp_atlas/genomes/reference_/calJac1_mito/genes/genes.gtf",
    # Chimpanzee: panTro5 assembly peaks, using panTro3 refGene as closest available
    "Chimpanzee": "/cluster/home/jjanssens/jjans/analysis/cerebellum/genomes_new/pan_troglodytes/panTro3.refGene.gtf.gz",
}


def get_reverse_chain_file(species: str, chain_dir: str = DEFAULT_CHAIN_DIR) -> str:
    """Get path to reverse chain file (hg38 -> species)."""
    if species not in REVERSE_CHAIN_FILES:
        raise ValueError(f"Unknown species: {species}. Available: {list(REVERSE_CHAIN_FILES.keys())}")
    return os.path.join(chain_dir, REVERSE_CHAIN_FILES[species])


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
    Find human peaks that do not overlap any non-human lifted peak.

    These are peaks detected in the human genome that could not be
    mapped from any other species (i.e., human-specific regulatory elements).

    Args:
        human_bed: Path to human consensus peaks (hg38)
        nonhuman_lifted_beds: List of non-human species peaks lifted to hg38
        output_bed: Output BED file for human-specific peaks
        merge_distance: Slop for overlap (0 = strict overlap)
        verbose: Print progress

    Returns:
        Dictionary with statistics
    """
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
    """Add human-specific peak IDs (human_peak_NNNNNN)."""
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
# LEGACY: merge_bed_files (kept for backwards compatibility)
# =============================================================================

def merge_bed_files(
    input_beds: List[str],
    output_bed: str,
    merge_distance: int = 0,
    verbose: bool = True,
) -> Dict[str, Any]:
    """Legacy simple merge without species tracking."""
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
    min_match: float = 0.95,
    min_blocks: Optional[float] = None,
    multiple: bool = False,
    merge_distance: int = 0,
    peak_prefix: str = "unified",
    gtf_files: Optional[Dict[str, str]] = None,
    verbose: bool = True,
    ncpu: int = 1,
) -> Dict[str, Any]:
    """
    Complete cross-species peak comparison pipeline with species tracking.

    Steps:
    1. Liftover all non-human species peaks to hg38
    2. Merge lifted peaks + human peaks, tracking species of origin
    3. Identify human-specific peaks (human peaks not overlapping any non-human)
    4. Add peak IDs for tracking
    5. Liftover unified peaks back to each species
    6. Identify species-specific peaks for each non-human species
    7. Annotate all peaks with closest gene

    Args:
        species_beds: Dict mapping non-human species name to input BED file path
        human_bed: Path to human consensus peaks BED file (hg38)
        output_dir: Output directory for all results
        chain_dir: Directory containing chain files
        liftover_path: Path to liftOver executable
        min_match: Minimum match ratio for liftover
        min_blocks: Minimum ratio of alignment blocks that must map
        multiple: Allow multiple output regions
        merge_distance: Distance for merging peaks
        peak_prefix: Prefix for unified peak IDs
        gtf_files: Dict mapping species name to GTF file path (optional)
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
        "human_specific": None,
        "lift_back": {},
        "species_specific": {},
        "output_files": {},
    }

    # Create output directories
    lift_human_dir = os.path.join(output_dir, "01_lifted_to_human")
    merged_dir = os.path.join(output_dir, "02_merged_consensus")
    human_specific_dir = os.path.join(output_dir, "03_human_specific")
    lift_back_dir = os.path.join(output_dir, "04_lifted_back")
    species_specific_dir = os.path.join(output_dir, "05_species_specific")
    annotation_dir = os.path.join(output_dir, "06_annotation")

    for d in [lift_human_dir, merged_dir, human_specific_dir, lift_back_dir,
              species_specific_dir, annotation_dir]:
        os.makedirs(d, exist_ok=True)

    print("=" * 70)
    print("CROSS-SPECIES CONSENSUS PIPELINE (v2 - with species tracking)")
    if ncpu > 1:
        print(f"Using {ncpu} parallel workers")
    print("=" * 70)

    # =========================================================================
    # STEP 1: Liftover all non-human species to hg38
    # =========================================================================
    print("\nSTEP 1: Liftover non-human species -> hg38")
    print("-" * 50)

    lifted_beds = {}  # species -> lifted BED path

    for species, input_bed in species_beds.items():
        if not os.path.exists(input_bed):
            print(f"   Skipping {species} - file not found: {input_bed}")
            continue

        output_bed = os.path.join(lift_human_dir, f"{species}_hg38.bed")

        if verbose:
            print(f"\n   {species}:")

        if species == "Marmoset":
            chain1 = get_chain_file("Marmoset_step1", chain_dir)
            chain2 = get_chain_file("Marmoset_step2", chain_dir)
            result = liftover_two_step(
                input_bed=input_bed, output_bed=output_bed,
                chain_file_1=chain1, chain_file_2=chain2,
                liftover_path=liftover_path, min_match=min_match,
                min_blocks=min_blocks, multiple=multiple,
                auto_chr=True, verbose=verbose, ncpu=ncpu,
            )
        else:
            chain_file = get_chain_file(species, chain_dir)
            result = liftover_peaks(
                input_bed=input_bed, output_bed=output_bed,
                chain_file=chain_file, liftover_path=liftover_path,
                min_match=min_match, min_blocks=min_blocks,
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
    print("STEP 2: Merge all species with species tracking")
    print("-" * 50)

    # Build dict of all species BEDs in hg38 coordinates
    all_species_hg38 = {"Human": human_bed}
    all_species_hg38.update(lifted_beds)

    merged_bed = os.path.join(merged_dir, "unified_consensus_hg38_merged.bed")
    merge_result = merge_with_species_tracking(
        species_beds=all_species_hg38,
        output_bed=merged_bed,
        merge_distance=merge_distance,
        verbose=verbose,
    )
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
    # STEP 3: Identify human-specific peaks
    # =========================================================================
    print("\n" + "=" * 70)
    print("STEP 3: Identify human-specific peaks")
    print("-" * 50)

    human_specific_raw = os.path.join(human_specific_dir, "human_specific_peaks_raw.bed")
    hs_result = find_human_specific_peaks(
        human_bed=human_bed,
        nonhuman_lifted_beds=list(lifted_beds.values()),
        output_bed=human_specific_raw,
        verbose=verbose,
    )

    # Add peak IDs to human-specific
    human_specific_bed = os.path.join(human_specific_dir, "human_specific_peaks.bed")
    add_peak_ids_to_human_specific(human_specific_raw, human_specific_bed, verbose=verbose)
    results["human_specific"] = hs_result
    results["output_files"]["human_specific"] = human_specific_bed

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
            min_match=min_match,
            min_blocks=min_blocks,
            multiple=multiple,
            auto_chr=True,
            verbose=verbose,
            ncpu=ncpu,
        )

        results["lift_back"][species] = result
        results["output_files"][f"liftback_{species}"] = output_bed

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
        human_specific_bed=human_specific_bed,
        species_specific_beds=species_specific_beds,
        gtf_files=gtf_files,
        output_dir=annotation_dir,
        verbose=verbose,
    )
    results["annotation"] = annotation_result
    results["output_files"]["annotation"] = os.path.join(annotation_dir, "peak_annotation.tsv")

    # =========================================================================
    # SUMMARY
    # =========================================================================
    print("\n" + "=" * 70)
    print("PIPELINE SUMMARY")
    print("=" * 70)

    print(f"\n   Unified consensus: {merge_result['merged_peaks']:,} peaks")
    print(f"   Human-specific: {hs_result['human_specific']:,} peaks")

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

    print(f"\n   Species-specific peaks:")
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
        f"Unified: {merge_result['merged_peaks']:,}, "
        f"Human-specific: {hs_result['human_specific']:,}"
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
    Uses the species_detected column from the unified BED if available,
    otherwise checks lifted-back BED files.
    """
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
