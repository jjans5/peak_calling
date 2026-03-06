"""
pipeline_steps.py
=================
Thin, self-contained wrapper functions for the cross-species consensus
pipeline.  Each function corresponds to one step in the step-by-step
notebook (``04_cross_species_pipeline_v3.ipynb``).

All heavy lifting is done by the functions in ``cross_species.py`` and
``liftover.py``; this module only adds:

  * A clear 1-function-per-step API
  * Consistent ``verbose`` logging
  * Automatic ID mapping when writing final BEDs
  * Rescue logic (original peaks missing from final set)
  * Distance classification (promoter / proximal / distal)

Usage
-----
>>> from src.pipeline_steps import (
...     lift_to_human,
...     merge_consensus,
...     lift_back_to_species,
...     filter_liftback,
...     extract_species_specific,
...     generate_master_annotation,
...     rescue_unmapped_peaks,
...     classify_distances,
...     export_final_beds,
... )
"""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from typing import Any, Dict, List, Optional, Union

import numpy as np
import pandas as pd

from .cross_species import (
    add_peak_ids,
    annotate_with_closest_gene,
    build_master_annotation,
    classify_peak_distance,
    create_peak_annotation,
    extract_gene_bed_from_gtf,
    filter_liftback_by_size,
    find_species_specific_peaks,
    liftback_peaks,
    reciprocal_liftover_check,
    resize_peaks,
    summit_based_merge,
    DEFAULT_GTF_FILES,
)
from .liftover import (
    get_chain_file,
    liftover_peaks,
    liftover_two_step,
    DEFAULT_CHAIN_DIR,
)


# =====================================================================
# Step 1 — Lift NHP peaks → hg38
# =====================================================================

def lift_to_human(
    species_beds: Dict[str, str],
    output_dir: str,
    chain_dir: str = DEFAULT_CHAIN_DIR,
    liftover_path: str = "liftOver",
    min_match: Union[float, Dict[str, float]] = 0.95,
    ncpu: int = 1,
    auto_chr: bool = True,
    verbose: bool = True,
) -> Dict[str, str]:
    """Lift each NHP species' peaks to hg38.

    Returns
    -------
    dict
        ``{species: path_to_lifted_bed}``
    """
    os.makedirs(output_dir, exist_ok=True)
    lifted = {}

    for species, input_bed in species_beds.items():
        output_bed = os.path.join(output_dir, f"{species}_hg38.bed")
        mm = min_match[species] if isinstance(min_match, dict) else min_match

        if verbose:
            print(f"  Lifting {species} → hg38  (minMatch={mm})")

        if species == "Marmoset":
            chain1 = get_chain_file("Marmoset_step1", chain_dir)
            chain2 = get_chain_file("Marmoset_step2", chain_dir)
            result = liftover_two_step(
                input_bed=input_bed, output_bed=output_bed,
                chain_file_1=chain1, chain_file_2=chain2,
                liftover_path=liftover_path, min_match=mm,
                auto_chr=auto_chr, verbose=verbose, ncpu=ncpu,
            )
        else:
            chain = get_chain_file(species, chain_dir)
            result = liftover_peaks(
                input_bed=input_bed, output_bed=output_bed,
                chain_file=chain, liftover_path=liftover_path,
                min_match=mm, auto_chr=auto_chr, verbose=verbose, ncpu=ncpu,
            )

        lifted[species] = output_bed
        if verbose:
            n = result.get("lifted", 0)
            u = result.get("unmapped", 0)
            print(f"    → {n:,} lifted, {u:,} unmapped")

    return lifted


# =====================================================================
# Step 2 — Summit-based merge → unified consensus
# =====================================================================

def merge_consensus(
    all_beds: Dict[str, str],
    output_dir: str,
    summit_window: int = 500,
    cluster_distance: int = 250,
    max_peak_size: int = 2000,
    min_peak_size: int = 100,
    verbose: bool = True,
) -> Dict[str, Any]:
    """Merge all species BEDs into a unified hg38 consensus.

    Returns
    -------
    dict
        Keys: ``merged_bed``, ``merged_with_ids``, ``bed4``,
        ``n_consensus``, ``merge_result``.
    """
    os.makedirs(output_dir, exist_ok=True)

    merged_bed = os.path.join(output_dir, "unified_consensus_hg38_merged.bed")
    merge_result = summit_based_merge(
        species_beds=all_beds,
        output_bed=merged_bed,
        summit_window=summit_window,
        cluster_distance=cluster_distance,
        max_peak_size=max_peak_size,
        min_peak_size=min_peak_size,
        normalize_scores=True,
        verbose=verbose,
    )

    # Temporary IDs (final IDs assigned in step 6)
    merged_with_ids = os.path.join(output_dir, "unified_consensus_hg38_with_ids.bed")
    add_peak_ids(merged_bed, merged_with_ids, prefix="temp")

    # BED4 for liftback
    bed4 = os.path.join(output_dir, "unified_consensus_hg38_bed4.bed")
    with open(merged_with_ids) as fin, open(bed4, "w") as fout:
        for line in fin:
            parts = line.strip().split("\t")
            fout.write(f"{parts[0]}\t{parts[1]}\t{parts[2]}\t{parts[3]}\n")

    n = sum(1 for _ in open(bed4))
    if verbose:
        print(f"  → {n:,} unified consensus peaks (temp IDs)")

    return {
        "merged_bed": merged_bed,
        "merged_with_ids": merged_with_ids,
        "bed4": bed4,
        "n_consensus": n,
        "merge_result": merge_result,
    }


# =====================================================================
# Step 3 — Liftback unified consensus → each species
# =====================================================================

def lift_back_to_species(
    unified_bed4: str,
    nhp_species: List[str],
    output_dir: str,
    chain_dir: str = DEFAULT_CHAIN_DIR,
    liftover_path: str = "liftOver",
    min_match: Union[float, Dict[str, float]] = 0.95,
    ncpu: int = 1,
    verbose: bool = True,
) -> Dict[str, Dict[str, Any]]:
    """Lift unified hg38 consensus back to each NHP species."""
    os.makedirs(output_dir, exist_ok=True)
    results = {}

    for species in nhp_species:
        output_bed = os.path.join(output_dir, f"unified_consensus_{species}.bed")
        mm = min_match[species] if isinstance(min_match, dict) else min_match

        result = liftback_peaks(
            input_bed=unified_bed4,
            output_bed=output_bed,
            species=species,
            chain_dir=chain_dir,
            liftover_path=liftover_path,
            min_match=mm,
            auto_chr=True,
            verbose=verbose,
            ncpu=ncpu,
        )
        results[species] = result
        if verbose:
            n = sum(1 for _ in open(output_bed)) if os.path.exists(output_bed) else 0
            print(f"  {species}: {n:,} peaks lifted back")

    return results


# =====================================================================
# Step 4 — Quality-filter liftback (size + reciprocal)
# =====================================================================

def filter_liftback(
    liftback_dir: str,
    unified_bed4: str,
    nhp_species: List[str],
    output_dir_size: str,
    output_dir_recip: str,
    chain_dir: str = DEFAULT_CHAIN_DIR,
    liftover_path: str = "liftOver",
    min_match: Union[float, Dict[str, float]] = 0.95,
    min_liftback_size: int = 100,
    max_liftback_size: int = 2000,
    max_recip_distance: int = 500,
    ncpu: int = 1,
    verbose: bool = True,
) -> Dict[str, Any]:
    """Size-filter + reciprocal liftover check.

    After filtering, the best-quality file per species is copied back to
    ``liftback_dir`` (overwriting the raw file) so that downstream steps
    read filtered coordinates automatically.

    Returns
    -------
    dict
        Keys: ``size_result``, ``recip_result``.
    """
    os.makedirs(output_dir_size, exist_ok=True)
    os.makedirs(output_dir_recip, exist_ok=True)

    # 4a — size filter
    size_result = filter_liftback_by_size(
        liftback_dir=liftback_dir,
        output_dir=output_dir_size,
        species_list=nhp_species,
        min_liftback_size=min_liftback_size,
        max_liftback_size=max_liftback_size,
        verbose=verbose,
    )

    # 4b — reciprocal check
    recip_result = reciprocal_liftover_check(
        liftback_dir=output_dir_size,
        hg38_bed=unified_bed4,
        output_dir=output_dir_recip,
        species_list=nhp_species,
        chain_dir=chain_dir,
        liftover_path=liftover_path,
        min_match=min_match,
        max_distance=max_recip_distance,
        file_pattern="unified_consensus_{species}_filtered.bed",
        verbose=verbose,
        ncpu=ncpu,
    )

    # Copy best file back to liftback_dir
    for species in nhp_species:
        pass_file = os.path.join(
            output_dir_recip,
            f"unified_consensus_{species}_filtered_reciprocal_pass.bed",
        )
        filt_file = os.path.join(
            output_dir_size,
            f"unified_consensus_{species}_filtered.bed",
        )
        dst = os.path.join(liftback_dir, f"unified_consensus_{species}.bed")

        if os.path.exists(pass_file) and os.path.getsize(pass_file) > 0:
            shutil.copy2(pass_file, dst)
        elif os.path.exists(filt_file):
            shutil.copy2(filt_file, dst)

        n = sum(1 for _ in open(dst)) if os.path.exists(dst) else 0
        if verbose:
            print(f"  {species} final liftback: {n:,}")

    return {"size_result": size_result, "recip_result": recip_result}


# =====================================================================
# Step 5 — Species-specific peaks
# =====================================================================

def extract_species_specific(
    species_beds: Dict[str, str],
    liftback_dir: str,
    output_dir: str,
    nhp_species: List[str],
    verbose: bool = True,
) -> Dict[str, str]:
    """Find peaks that exist only in one NHP species.

    Returns
    -------
    dict
        ``{species: path_to_species_specific_bed}``
    """
    os.makedirs(output_dir, exist_ok=True)
    specific = {}

    for species in nhp_species:
        original_bed = species_beds[species]
        liftback_file = os.path.join(liftback_dir, f"unified_consensus_{species}.bed")
        if not os.path.exists(liftback_file):
            if verbose:
                print(f"  {species}: no liftback file, skipping")
            continue

        sp_specific_bed = os.path.join(output_dir, f"{species}_specific_peaks.bed")
        find_species_specific_peaks(
            original_bed=original_bed,
            liftback_bed=liftback_file,
            output_bed=sp_specific_bed,
            species=species,
            verbose=verbose,
        )
        specific[species] = sp_specific_bed
        n = sum(1 for _ in open(sp_specific_bed)) if os.path.exists(sp_specific_bed) else 0
        if verbose:
            print(f"  {species}: {n:,} species-specific peaks")

    return specific


# =====================================================================
# Step 6 — Gene annotation + master annotation
# =====================================================================

def generate_master_annotation(
    merged_with_ids: str,
    species_specific_beds: Dict[str, str],
    liftback_dir: str,
    annotation_dir: str,
    master_dir: str,
    nhp_species: List[str],
    gtf_files: Optional[Dict[str, str]] = None,
    verbose: bool = True,
) -> pd.DataFrame:
    """Run gene annotation and build master annotation.

    Returns
    -------
    pd.DataFrame
        The master annotation table (index = ``peak_id``).
    """
    if gtf_files is None:
        gtf_files = DEFAULT_GTF_FILES.copy()

    os.makedirs(annotation_dir, exist_ok=True)
    os.makedirs(master_dir, exist_ok=True)

    # 6a — gene annotation (creates gene BEDs + per-peak annotation TSVs)
    create_peak_annotation(
        unified_bed=merged_with_ids,
        human_specific_bed="",
        species_specific_beds=species_specific_beds,
        gtf_files=gtf_files,
        output_dir=annotation_dir,
        verbose=verbose,
    )

    # 6b — master annotation (ID renaming happens here)
    annotation_tsv = os.path.join(annotation_dir, "peak_annotation.tsv")
    master_output = os.path.join(master_dir, "master_annotation.tsv")

    master = build_master_annotation(
        unified_bed=merged_with_ids,
        human_gene_annotation_tsv=annotation_tsv,
        liftback_dir=liftback_dir,
        species_specific_beds=species_specific_beds,
        gene_bed_dir=annotation_dir,
        output_file=master_output,
        species_list=nhp_species,
        gtf_files=gtf_files,
        verbose=verbose,
    )

    if verbose:
        print(f"\n  Peak types: {master['peak_type'].value_counts().to_dict()}")
        print(f"  Total peaks: {len(master):,}")

    return master


# =====================================================================
# Step 7 — Rescue missing peaks
# =====================================================================

def rescue_unmapped_peaks(
    species_beds: Dict[str, str],
    liftback_dir: str,
    species_specific_beds: Dict[str, str],
    nhp_species: List[str],
    output_dir: str,
    target_size: int = 500,
    verbose: bool = True,
) -> Dict[str, Any]:
    """Rescue original NHP peaks lost during filtering.

    For each species, finds original peaks that do NOT overlap any
    final liftback or species-specific peak.  These are resized to
    ``target_size`` and given ``{species}_rescued_NNNNNN`` IDs.

    Returns
    -------
    dict
        ``rescued_beds``: ``{species: path}``,
        ``rescued_rows``: list of dicts for appending to master.
    """
    os.makedirs(output_dir, exist_ok=True)
    rescued_beds: Dict[str, str] = {}
    rescued_rows: list = []

    for species in nhp_species:
        original_bed = species_beds[species]
        liftback_bed = os.path.join(liftback_dir, f"unified_consensus_{species}.bed")
        sp_specific_bed = species_specific_beds.get(species, "")

        if not os.path.exists(liftback_bed):
            continue

        with tempfile.TemporaryDirectory() as tmpdir:
            # Combine liftback + species-specific into "covered" BED
            covered_bed = os.path.join(tmpdir, "covered.bed")
            with open(covered_bed, "w") as fout:
                for src in [liftback_bed, sp_specific_bed]:
                    if src and os.path.exists(src):
                        with open(src) as fin:
                            for line in fin:
                                parts = line.strip().split("\t")
                                if len(parts) >= 3:
                                    fout.write(f"{parts[0]}\t{parts[1]}\t{parts[2]}\n")

            # Extract BED3 from original
            orig_bed3 = os.path.join(tmpdir, "orig.bed")
            with open(original_bed) as fin, open(orig_bed3, "w") as fout:
                for line in fin:
                    parts = line.strip().split("\t")
                    if len(parts) >= 3:
                        fout.write(f"{parts[0]}\t{parts[1]}\t{parts[2]}\n")

            # Harmonize chr prefix
            def _first_chr(path):
                with open(path) as f:
                    for line in f:
                        if line.strip():
                            return line.split("\t")[0].startswith("chr")
                return True

            orig_has_chr = _first_chr(orig_bed3)
            cov_has_chr = _first_chr(covered_bed)
            if orig_has_chr != cov_has_chr:
                fixed = os.path.join(tmpdir, "covered_fixed.bed")
                with open(covered_bed) as fin, open(fixed, "w") as fout:
                    for line in fin:
                        if orig_has_chr and not line.startswith("chr"):
                            fout.write("chr" + line)
                        elif not orig_has_chr and line.startswith("chr"):
                            fout.write(line[3:])
                        else:
                            fout.write(line)
                covered_bed = fixed

            # Sort & find missing
            orig_sorted = os.path.join(tmpdir, "orig_sorted.bed")
            cov_sorted = os.path.join(tmpdir, "cov_sorted.bed")
            subprocess.run(f"sort -k1,1 -k2,2n {orig_bed3} > {orig_sorted}",
                           shell=True, check=True)
            subprocess.run(f"sort -k1,1 -k2,2n {covered_bed} | bedtools merge -i - > {cov_sorted}",
                           shell=True, check=True)

            missing_bed = os.path.join(tmpdir, "missing.bed")
            subprocess.run(
                f"bedtools intersect -a {orig_sorted} -b {cov_sorted} -v > {missing_bed}",
                shell=True, check=True,
            )

            n_missing = sum(1 for _ in open(missing_bed))
            if n_missing == 0:
                if verbose:
                    print(f"  {species}: 0 peaks to rescue")
                continue

            # Resize & assign IDs
            rescued_out = os.path.join(output_dir, f"{species}_rescued_peaks.bed")
            sp_lower = species.lower()
            counter = 0
            with open(missing_bed) as fin, open(rescued_out, "w") as fout:
                for line in fin:
                    parts = line.strip().split("\t")
                    chrom, start, end = parts[0], int(parts[1]), int(parts[2])
                    center = (start + end) // 2
                    new_start = max(0, center - target_size // 2)
                    new_end = new_start + target_size
                    counter += 1
                    pid = f"{sp_lower}_rescued_{counter:06d}"
                    fout.write(f"{chrom}\t{new_start}\t{new_end}\t{pid}\n")

                    rescued_rows.append({
                        "peak_id": pid,
                        "peak_type": f"{sp_lower}_rescued",
                        f"{species}_chr": chrom,
                        f"{species}_start": new_start,
                        f"{species}_end": new_end,
                        "species_detected": species,
                        "n_species_det": 1,
                        "n_species_orth": 0,
                        f"{species}_det": 1,
                        f"{species}_orth": 0,
                    })

            rescued_beds[species] = rescued_out
            if verbose:
                print(f"  {species}: {counter:,} rescued peaks")

    if verbose:
        print(f"  Total rescued: {len(rescued_rows):,}")

    return {"rescued_beds": rescued_beds, "rescued_rows": rescued_rows}


# =====================================================================
# Step 8 — Distance classification
# =====================================================================

def classify_distances(
    master: pd.DataFrame,
    rescued_rows: list,
    master_dir: str,
    promoter_threshold: int = 200,
    proximal_threshold: int = 2000,
    verbose: bool = True,
) -> pd.DataFrame:
    """Append rescued peaks, classify all peaks by distance to TSS.

    Returns the final master DataFrame.
    """
    # Append rescued rows
    if rescued_rows:
        rescued_df = pd.DataFrame(rescued_rows).set_index("peak_id")
        master = pd.concat([master, rescued_df])
        if verbose:
            print(f"  Appended {len(rescued_df):,} rescued peaks to master")

    # Classify distance
    master = classify_peak_distance(
        master,
        promoter_threshold=promoter_threshold,
        proximal_threshold=proximal_threshold,
    )

    # Save
    out = os.path.join(master_dir, "master_annotation_final.tsv")
    master.to_csv(out, sep="\t")
    if verbose:
        # Show region counts for species that have them
        region_cols = [c for c in master.columns if c.endswith("_region")]
        for rc in region_cols:
            print(f"  {rc}: {master[rc].value_counts().to_dict()}")
        print(f"  Saved: {out}")

    return master


# =====================================================================
# Step 9 — Write final BED files (one per species)
# =====================================================================

def export_final_beds(
    master: pd.DataFrame,
    liftback_dir: str,
    species_specific_beds: Dict[str, str],
    rescued_beds: Dict[str, str],
    merged_with_ids: str,
    output_dir: str,
    nhp_species: List[str],
    master_dir: str,
    target_size: int = 500,
    verbose: bool = True,
) -> Dict[str, str]:
    """Write one BED per species with unified + specific + rescued peaks,
    all resized to ``target_size``.

    Peak IDs in the BEDs match the master annotation (final IDs, not
    temp_ IDs).

    Returns
    -------
    dict
        ``{species: path_to_final_bed}``
    """
    os.makedirs(output_dir, exist_ok=True)

    # Load ID mapping (temp_NNNNNN → unified_NNNNNN / human_peak_NNNNNN)
    mapping_file = os.path.join(master_dir, "master_annotation_id_mapping.tsv")
    id_map: Dict[str, str] = {}
    if os.path.exists(mapping_file):
        mdf = pd.read_csv(mapping_file, sep="\t")
        id_map = dict(zip(mdf["old_id"], mdf["new_id"]))
        if verbose:
            print(f"  Loaded ID mapping: {len(id_map):,} entries")

    all_species = ["Human"] + nhp_species
    final_beds: Dict[str, str] = {}

    for species in all_species:
        out_bed = os.path.join(output_dir, f"all_peaks_{species}.bed")
        entries: list = []  # (chrom, start, end, peak_id)

        if species == "Human":
            # Unified peaks in hg38 (from master)
            mask = master["peak_type"].isin(["unified", "human_specific"])
            for pid, row in master.loc[mask].iterrows():
                if pd.notna(row.get("Human_chr")):
                    entries.append((
                        str(row["Human_chr"]),
                        int(row["Human_start"]),
                        int(row["Human_end"]),
                        str(pid),
                    ))
        else:
            # a) Unified peaks from liftback (filtered, in species coords)
            lb_file = os.path.join(liftback_dir, f"unified_consensus_{species}.bed")
            if os.path.exists(lb_file):
                with open(lb_file) as f:
                    for line in f:
                        parts = line.strip().split("\t")
                        if len(parts) >= 4:
                            raw_id = parts[3]
                            final_id = id_map.get(raw_id, raw_id)
                            entries.append((parts[0], int(parts[1]), int(parts[2]), final_id))

            # b) Species-specific peaks
            sp_bed = species_specific_beds.get(species, "")
            if sp_bed and os.path.exists(sp_bed):
                with open(sp_bed) as f:
                    for line in f:
                        parts = line.strip().split("\t")
                        if len(parts) >= 4:
                            entries.append((parts[0], int(parts[1]), int(parts[2]), parts[3]))

            # c) Rescued peaks
            resc_bed = rescued_beds.get(species, "")
            if resc_bed and os.path.exists(resc_bed):
                with open(resc_bed) as f:
                    for line in f:
                        parts = line.strip().split("\t")
                        if len(parts) >= 4:
                            entries.append((parts[0], int(parts[1]), int(parts[2]), parts[3]))

        # Resize to uniform width
        half = target_size // 2
        with open(out_bed, "w") as fout:
            for chrom, start, end, pid in entries:
                center = (start + end) // 2
                ns = max(0, center - half)
                ne = ns + target_size
                fout.write(f"{chrom}\t{ns}\t{ne}\t{pid}\n")

        final_beds[species] = out_bed
        if verbose:
            print(f"  {species}: {len(entries):,} peaks → {out_bed}")

    return final_beds
