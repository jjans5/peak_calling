"""
Fragment-based count matrix construction & pseudobulk aggregation.

Functions for building single-cell peaks × barcodes count matrices from
fragment files, and for aggregating them into pseudobulk profiles.
"""

from __future__ import annotations

import os
import re
import glob
from collections import Counter

import numpy as np
import pandas as pd
import scipy.sparse as sp_sparse

import polars as pl

from pycisTopic.fragments import (
    read_fragments_to_polars_df,
    read_bed_to_polars_df,
    filter_fragments_by_cb,
)
from pycisTopic.genomic_ranges import intersection as gr_intersection


# ---------------------------------------------------------------------------
# Barcode re-indexing helpers
# ---------------------------------------------------------------------------

def reindex_nhp(cell_data: pd.DataFrame) -> pd.DataFrame:
    """Re-index NHP metadata so barcodes match ``{barcode}-{sample}`` format.

    Original index:  ``Sample_086_AAACAGCCAAGGAATC-1``
    Re-indexed:      ``AAACAGCCAAGGAATC-1-Sample_086``
    """
    cell_data = cell_data.copy()
    cell_data["old_barcode"] = cell_data.index
    idx = cell_data.index.to_series()
    idx = idx.str.replace("_plus_resequenced", "", regex=False)
    idx = idx.str.replace("Sample_", "Sample", regex=False)
    idx = idx.apply(lambda x: f"{x.split('_')[1]}-{x.split('_')[0]}")
    idx = idx.str.replace("Sample", "Sample_", regex=False)
    idx = idx.str.replace("MSN", "Sample_", regex=False)
    cell_data.index = idx
    return cell_data


def reindex_human(cell_data: pd.DataFrame) -> pd.DataFrame:
    """Re-index Human metadata so barcodes match ``{barcode}-{sample}`` format.

    Original index:  ``B006-A-002#AAACAGCCAACTAGCC-1``
    Re-indexed:      ``AAACAGCCAACTAGCC-1-B006_A_002``
    """
    cell_data = cell_data.copy()
    idx = cell_data.index.to_series()
    idx = idx.str.replace("_plus_resequenced", "", regex=False)
    idx = idx.str.replace("Sample_", "Sample", regex=False)
    idx = idx.apply(lambda x: f"{x.split('#')[1]}-{x.split('#')[0].replace('-', '_')}")
    cell_data.index = idx
    for col in ["orig.ident", "Sample.ID", "SampleNameOnly"]:
        if col in cell_data.columns:
            cell_data[col] = cell_data[col].str.replace("-", "_", regex=False)
    return cell_data


# ---------------------------------------------------------------------------
# Metadata & fragment-dict loading
# ---------------------------------------------------------------------------

def load_species_data(
    species: str,
    meta_dir: str,
    qc_dir: str,
    human_frag_dir: str,
    nhp_frag_dir: str,
    nhp_samples: dict[str, list[str]],
    verbose: bool = True,
) -> tuple[pd.DataFrame, dict[str, str]]:
    """Load metadata + fragment paths for one species, filtered to Adult + QC.

    Parameters
    ----------
    species : str
        Species name (e.g. "Human", "Bonobo").
    meta_dir : str
        Directory containing ``{Species}_RNA.tsv`` metadata files.
    qc_dir : str
        Directory containing ``qc_{Species}/{sample}.cbs_for_otsu_thresholds.tsv``.
    human_frag_dir : str
        Directory with Human fragment files (``*_atac_fragments.tsv.gz``).
    nhp_frag_dir : str
        Directory with NHP fragment files (``{Sample}_atac.fragments.tsv.gz``).
    nhp_samples : dict
        Mapping species → list of sample IDs for NHP species.
    verbose : bool
        Print diagnostic messages.

    Returns
    -------
    cell_data : pd.DataFrame
        Barcode-indexed metadata, filtered to Adult + QC-passing.
    fragments_dict : dict[str, str]
        sample_id → path to fragment file.
    """
    meta_path = os.path.join(meta_dir, f"{species}_RNA.tsv")
    cell_data = pd.read_csv(meta_path, sep="\t", index_col=0)

    # Re-index
    if species == "Human":
        cell_data = reindex_human(cell_data)
    else:
        cell_data = reindex_nhp(cell_data)

    # Filter to Adult
    cell_data = cell_data[cell_data["Age"] == "Adult"].copy()

    # Build fragment dict
    if species == "Human":
        frag_files = glob.glob(os.path.join(human_frag_dir, "*_atac_fragments.tsv.gz"))
        fragments_dict = {
            os.path.basename(f).split("_atac_fragments")[0].replace("-", "_"): f
            for f in frag_files
        }
        valid_samples = set(cell_data["orig.ident"].unique())
        fragments_dict = {k: v for k, v in fragments_dict.items() if k in valid_samples}
    else:
        fragments_dict = {
            s: os.path.join(nhp_frag_dir, f"{s}_atac.fragments.tsv.gz")
            for s in nhp_samples.get(species, [])
        }

    # QC barcode filtering
    qc_species_dir = os.path.join(qc_dir, f"qc_{species}")
    samples = sorted(set(cell_data["orig.ident"]))
    all_qc_barcodes = []

    for sample in samples:
        qc_file = os.path.join(qc_species_dir, f"{sample}.cbs_for_otsu_thresholds.tsv")
        if not os.path.exists(qc_file):
            if verbose:
                print(f"    ⚠️  QC file not found: {qc_file}")
            continue
        qc_df = pd.read_csv(qc_file, sep="\t", header=None)
        barcodes = [f"{bc}-{sample}" for bc in qc_df[0].tolist()]
        all_qc_barcodes.extend(barcodes)

    keep = sorted(set(all_qc_barcodes) & set(cell_data.index))
    cell_data = cell_data.loc[keep]

    # Deduplicate
    if cell_data.index.nunique() != len(cell_data):
        n_dupes = len(cell_data) - cell_data.index.nunique()
        if verbose:
            print(f"    ⚠️  {n_dupes} duplicate barcodes — keeping first occurrence")
        cell_data = cell_data[~cell_data.index.duplicated(keep="first")]

    # Keep only samples that still have cells
    remaining_samples = sorted(cell_data["orig.ident"].unique())
    fragments_dict = {s: fragments_dict[s] for s in remaining_samples
                      if s in fragments_dict}

    if verbose and len(cell_data) > 0:
        print(f"  Barcode format: '{cell_data.index[0]}' "
              f"(sample='{cell_data.iloc[0]['orig.ident']}')")

    return cell_data, fragments_dict


# ---------------------------------------------------------------------------
# Polars helpers
# ---------------------------------------------------------------------------

def load_regions_as_polars(peak_file: str) -> pl.DataFrame:
    """Load a BED file with ``RegionID`` (coordinates) and ``PeakID`` (name).

    ``RegionID`` is always ``chr:start-end`` (used internally for intersection
    matching).  ``PeakID`` comes from column 4 of the BED file when present
    (e.g. ``unified_000001``, ``human_peak_000003``); otherwise it falls back
    to the coordinate string.
    """
    regions = read_bed_to_polars_df(peak_file, engine="pyarrow", min_column_count=3)
    # Coordinate-based ID for intersection matching
    regions = regions.with_columns(
        (
            pl.col("Chromosome").cast(pl.Utf8)
            + ":" + pl.col("Start").cast(pl.Utf8)
            + "-" + pl.col("End").cast(pl.Utf8)
        ).alias("RegionID")
    )
    # Use Name column (BED col 4) as immutable PeakID when available
    if "Name" in regions.columns:
        regions = regions.with_columns(pl.col("Name").cast(pl.Utf8).alias("PeakID"))
    else:
        regions = regions.with_columns(pl.col("RegionID").alias("PeakID"))
    return regions.select(["Chromosome", "Start", "End", "RegionID", "PeakID"])


def harmonize_chroms(
    regions_pl: pl.DataFrame,
    fragments_pl: pl.DataFrame,
) -> tuple[pl.DataFrame, pl.DataFrame, bool]:
    """Ensure Chromosome columns use the same chr-prefix convention.

    Returns ``(regions, fragments, regions_changed)``.
    """
    reg_chroms = set(regions_pl.get_column("Chromosome").unique().to_list())
    frag_chroms = set(fragments_pl.get_column("Chromosome").unique().to_list())

    if reg_chroms & frag_chroms:
        return regions_pl, fragments_pl, False

    reg_has_chr = any(str(c).startswith("chr") for c in reg_chroms)
    frag_has_chr = any(str(c).startswith("chr") for c in frag_chroms)

    if reg_has_chr and not frag_has_chr:
        regions_pl = regions_pl.with_columns(
            pl.col("Chromosome").cast(pl.Utf8).str.replace("^chr", "")
            .cast(pl.Categorical).alias("Chromosome")
        )
        regions_pl = regions_pl.with_columns(
            (
                pl.col("Chromosome").cast(pl.Utf8)
                + ":" + pl.col("Start").cast(pl.Utf8)
                + "-" + pl.col("End").cast(pl.Utf8)
            ).alias("RegionID")
        )
        print("    [chr harmonize] Stripped 'chr' from peaks to match fragments")
        return regions_pl, fragments_pl, True

    if frag_has_chr and not reg_has_chr:
        fragments_pl = fragments_pl.with_columns(
            pl.col("Chromosome").cast(pl.Utf8).str.replace("^chr", "")
            .cast(pl.Categorical).alias("Chromosome")
        )
        print("    [chr harmonize] Stripped 'chr' from fragments to match peaks")
        return regions_pl, fragments_pl, False

    return regions_pl, fragments_pl, False


# ---------------------------------------------------------------------------
# Core matrix builder
# ---------------------------------------------------------------------------

def build_fragment_matrix(
    species: str,
    cell_data: pd.DataFrame,
    fragments_dict: dict[str, str],
    peak_file: str,
    verbose: bool = True,
) -> tuple[sp_sparse.csr_matrix | None, list[str], list[str], pd.DataFrame]:
    """Build a sparse (peaks × barcodes) count matrix for one species.

    Parameters
    ----------
    species : str
        Species name (for logging only).
    cell_data : pd.DataFrame
        Barcode-indexed metadata with ``orig.ident`` column.
    fragments_dict : dict
        sample_id → fragment file path.
    peak_file : str
        Path to the BED file of peaks (e.g. ``all_peaks_{Species}.bed``).
    verbose : bool
        Print progress per sample.

    Returns
    -------
    mat : csr_matrix or None
        Peaks × barcodes sparse count matrix.
    barcodes : list[str]
        Barcode names (columns).
    peak_ids : list[str]
        Peak IDs (rows).  Uses the Name column from the BED file when
        available (e.g. ``unified_000001``); falls back to coordinates.
    meta_df : pd.DataFrame
        Barcode metadata (cell_type, donor, region, age, run where available).
    """
    regions_pl = load_regions_as_polars(peak_file)
    region_ids = regions_pl.get_column("RegionID").to_list()  # coords for matching
    peak_ids = regions_pl.get_column("PeakID").to_list()       # immutable names
    n_regions = len(peak_ids)

    # Build coord→row_index mapping.
    # Coordinates may collide after liftover, so we map to ALL rows.
    region_to_idxs: dict[str, list[int]] = {}
    for i, rid in enumerate(region_ids):
        region_to_idxs.setdefault(rid, []).append(i)

    all_barcodes: list[str] = []
    all_rows: list[int] = []
    all_cols: list[int] = []
    all_vals: list[int] = []
    meta_records: list[dict] = []

    has_individual = "Individual" in cell_data.columns
    has_run = "Sequencing.Run" in cell_data.columns
    has_celltype = "cell_type_initial" in cell_data.columns
    has_region = "Region" in cell_data.columns
    has_age = "Age" in cell_data.columns
    _first_sample = True

    for sample_idx, (sample, frag_path) in enumerate(sorted(fragments_dict.items())):
        if not os.path.exists(frag_path):
            if verbose:
                print(f"  ⚠️  Fragment file not found: {frag_path}")
            continue

        sample_suffix = f"-{sample}"
        sample_barcodes = cell_data[cell_data["orig.ident"] == sample].index.tolist()
        if not sample_barcodes:
            continue

        if verbose:
            print(f"  [{sample_idx+1}/{len(fragments_dict)}] {sample}: "
                  f"{len(sample_barcodes):,} barcodes", end="", flush=True)

        cbs_raw = []
        for bc in sample_barcodes:
            if bc.endswith(sample_suffix):
                cbs_raw.append(bc[: -len(sample_suffix)])
            else:
                cbs_raw.append(bc)
        cbs_series = pl.Series("CB", cbs_raw, dtype=pl.Categorical)

        fragments_pl = read_fragments_to_polars_df(
            fragments_bed_filename=frag_path, engine="pyarrow",
        )
        fragments_filt = filter_fragments_by_cb(fragments_pl, cbs_series)
        del fragments_pl

        if fragments_filt.height == 0:
            if verbose:
                print(" → 0 fragments after CB filter")
            continue

        regions_pl, fragments_filt, regions_changed = harmonize_chroms(
            regions_pl, fragments_filt,
        )
        if regions_changed:
            region_ids = regions_pl.get_column("RegionID").to_list()
            region_to_idxs = {}
            for i, rid in enumerate(region_ids):
                region_to_idxs.setdefault(rid, []).append(i)

        overlap_df = gr_intersection(
            regions1_df_pl=regions_pl,
            regions2_df_pl=fragments_filt,
            how="all",
            regions1_info=True, regions2_info=True,
            regions1_coord=False, regions2_coord=False,
            regions1_suffix="@1", regions2_suffix="@2",
        )

        if _first_sample and verbose:
            print(f"\n    [DEBUG] intersection columns: {overlap_df.columns}")
            print(f"    [DEBUG] intersection rows:    {overlap_df.height:,}")
            _first_sample = False

        if overlap_df.height == 0:
            if verbose:
                print(" → 0 intersections")
            continue

        region_col = next(
            (c for c in ["RegionID@1", "RegionID"] if c in overlap_df.columns), None,
        )
        cb_col = next(
            (c for c in ["Name@2", "CB@2", "Name", "CB"] if c in overlap_df.columns), None,
        )
        if region_col is None or cb_col is None:
            if verbose:
                print(f" → ⚠️ missing columns (RegionID={region_col}, CB={cb_col})")
                print(f"       available: {overlap_df.columns}")
            continue

        region_cb_counts = (
            overlap_df.lazy()
            .group_by([region_col, cb_col])
            .agg(pl.len().cast(pl.UInt32).alias("count"))
            .collect()
        )
        del overlap_df, fragments_filt

        unique_cbs_raw = sorted(region_cb_counts.get_column(cb_col).unique().to_list())
        local_bc_map: dict[str, int] = {}
        sample_barcodes_set = set(sample_barcodes)

        for raw_bc in unique_cbs_raw:
            full_bc = f"{raw_bc}{sample_suffix}"
            if full_bc in sample_barcodes_set:
                local_bc_map[raw_bc] = len(all_barcodes)
                all_barcodes.append(full_bc)

                rec = {"barcode": full_bc, "sample": sample}
                if full_bc in cell_data.index:
                    if has_celltype:
                        rec["cell_type"] = cell_data.at[full_bc, "cell_type_initial"]
                    if has_individual:
                        rec["donor"] = cell_data.at[full_bc, "Individual"]
                    if has_run:
                        rec["run"] = cell_data.at[full_bc, "Sequencing.Run"]
                    if has_region:
                        rec["region"] = cell_data.at[full_bc, "Region"]
                    if has_age:
                        rec["age"] = cell_data.at[full_bc, "Age"]
                meta_records.append(rec)

        region_id_vals = region_cb_counts.get_column(region_col).to_list()
        cbs_vals = region_cb_counts.get_column(cb_col).to_list()
        count_vals = region_cb_counts.get_column("count").to_numpy()

        for i in range(len(region_id_vals)):
            bc = cbs_vals[i]
            rid = region_id_vals[i]
            if bc in local_bc_map and rid in region_to_idxs:
                col_idx = local_bc_map[bc]
                cnt = int(count_vals[i])
                for row_idx in region_to_idxs[rid]:
                    all_rows.append(row_idx)
                    all_cols.append(col_idx)
                    all_vals.append(cnt)

        if verbose:
            print(f" → {len(local_bc_map):,} barcodes with signal")

    # Assemble sparse matrix
    n_barcodes = len(all_barcodes)
    if n_barcodes == 0:
        if verbose:
            print("  ⚠️  No barcodes found!")
        return None, [], peak_ids, pd.DataFrame()

    if len(set(all_barcodes)) != n_barcodes:
        dupes = [bc for bc, cnt in Counter(all_barcodes).items() if cnt > 1]
        if verbose:
            print(f"  ⚠️  {len(dupes)} duplicate barcodes! Examples: {dupes[:5]}")

    mat = sp_sparse.csr_matrix(
        (np.array(all_vals, dtype=np.uint32),
         (np.array(all_rows, dtype=np.int64),
          np.array(all_cols, dtype=np.int64))),
        shape=(n_regions, n_barcodes),
    )

    meta_df = pd.DataFrame(meta_records)
    if len(meta_df) > 0:
        meta_df = meta_df.set_index("barcode")

    if verbose:
        print(f"\n  ✅ Final matrix: {n_regions:,} peaks × {n_barcodes:,} barcodes, "
              f"{mat.nnz:,} nonzero entries")
        if peak_ids and not peak_ids[0].startswith(("chr", "Chr")):
            print(f"      Peak ID format: '{peak_ids[0]}' … '{peak_ids[-1]}'")

    return mat, all_barcodes, peak_ids, meta_df


# ---------------------------------------------------------------------------
# Pseudobulk aggregation
# ---------------------------------------------------------------------------

def create_pseudobulk(
    mat: sp_sparse.csr_matrix,
    barcodes: list[str],
    region_ids: list[str],
    meta_df: pd.DataFrame,
    group_by: list[str] | None = None,
    min_cells: int = 5,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Sum single-cell counts into pseudobulk profiles.

    Parameters
    ----------
    mat : csr_matrix
        Peaks × barcodes sparse matrix.
    barcodes : list
        Barcode names (columns of *mat*).
    region_ids : list
        Region / peak IDs (rows of *mat*).  Typically the immutable peak
        names returned by :func:`build_fragment_matrix`.
    meta_df : pd.DataFrame
        Barcode metadata with columns matching *group_by*.
    group_by : list of str, optional
        Columns to group barcodes by (default: ``["cell_type", "donor"]``).
    min_cells : int
        Minimum barcodes per group.

    Returns
    -------
    counts_df : pd.DataFrame
        Peaks × pseudobulk samples.
    info_df : pd.DataFrame
        Group metadata (label, n_cells, etc.).
    """
    if group_by is None:
        group_by = ["cell_type", "donor"]

    for col in group_by:
        if col not in meta_df.columns:
            raise ValueError(f"Column '{col}' not in metadata: {list(meta_df.columns)}")

    bc_to_idx = {bc: i for i, bc in enumerate(barcodes)}
    groups = meta_df.groupby(group_by, observed=True)

    pseudobulk: dict[str, np.ndarray] = {}
    group_info: list[dict] = []

    for group_key, group_df in groups:
        if len(group_df) < min_cells:
            continue

        col_indices = [bc_to_idx[bc] for bc in group_df.index if bc in bc_to_idx]
        if len(col_indices) < min_cells:
            continue

        sub = mat[:, col_indices]
        counts = np.asarray(sub.sum(axis=1)).ravel()

        if isinstance(group_key, tuple):
            label = "__".join(str(x) for x in group_key)
        else:
            label = str(group_key)

        label_clean = re.sub(r"[/\\:*?\"<>|]", "_", label)
        label_clean = re.sub(r"\s+", "_", label_clean)

        pseudobulk[label_clean] = counts
        rec = {col: val for col, val in zip(
            group_by,
            group_key if isinstance(group_key, tuple) else [group_key],
        )}
        rec["label"] = label_clean
        rec["n_cells"] = len(col_indices)
        group_info.append(rec)

    counts_df = pd.DataFrame(pseudobulk, index=region_ids)
    info_df = pd.DataFrame(group_info)
    return counts_df, info_df
