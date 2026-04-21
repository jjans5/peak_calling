"""
Evolutionary Heatmap Utilities
===============================

Functions for building evolutionary accessibility modules across cell types,
plotting blocked heatmaps, exporting BED files, and running SNC enrichment.

All functions are parameterized (no module-level globals) so they can be
called from any notebook or script.
"""

import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Patch
from pathlib import Path
from typing import Optional, Union

try:
    import pybedtools
    HAS_PYBEDTOOLS = True
except ImportError:
    HAS_PYBEDTOOLS = False

try:
    from scipy.stats import fisher_exact
    from statsmodels.stats.multitest import multipletests
    HAS_STATS = True
except ImportError:
    HAS_STATS = False


# =============================================================================
# LOW-LEVEL HELPERS
# =============================================================================

def is_true(series):
    """Coerce a mixed bool/str/int annotation column to a boolean mask."""
    if series.dtype == bool:
        return series.fillna(False)
    if series.dtype == object:
        return series.astype(str).str.strip().str.upper().isin(["TRUE", "1"])
    return series.fillna(0).astype(bool)


def row_zscore(df):
    """Row-wise z-score, preserving NaNs for non-orthologous cells."""
    valid = df.notna()
    mu = df.mean(axis=1, skipna=True)
    sd = df.std(axis=1, skipna=True, ddof=0)
    z = df.sub(mu, axis=0).div(sd.replace(0, np.nan), axis=0)
    bad_rows = sd.isna() | (sd == 0)
    if bad_rows.any():
        z.loc[bad_rows] = z.loc[bad_rows].where(~valid.loc[bad_rows], 0.0)
    return z.where(valid)


def sanitize_name(x):
    """Return a filesystem-safe version of a string (spaces/special chars → _)."""
    x = str(x).strip()
    x = re.sub(r"\s+", "_", x)
    x = re.sub(r"[^A-Za-z0-9_.+-]", "_", x)
    return x


def celltype_to_dirname(cell_type):
    """
    Convert a cell type display name to the R-style directory name used in
    the DESeq2 results (mirrors R's make.names(): non-alphanumeric → dot).
    e.g. 'Crypt Fibroblasts WNT2B+' → 'Crypt.Fibroblasts.WNT2B.'
    """
    return re.sub(r"[^A-Za-z0-9]", ".", cell_type)


# =============================================================================
# I/O
# =============================================================================

def load_robust(contrast, robust_dir, cell_type):
    """Load ultra-robust peak IDs for a given contrast and cell type."""
    ct_dir = celltype_to_dirname(cell_type)
    path = Path(robust_dir) / ct_dir / f"{contrast}_UltraRobust.csv"
    if not path.exists():
        print(f"  [MISSING robust] {cell_type} / {contrast}")
        return set()
    return set(pd.read_csv(path).iloc[:, 0].astype(str))


def load_branch_df(contrast, branch_dir, cell_type):
    """Load the full DESeq2 branch result DataFrame for a contrast + cell type."""
    ct_dir = celltype_to_dirname(cell_type)
    path = Path(branch_dir) / ct_dir / f"{contrast}.csv"
    if not path.exists():
        return None
    df = pd.read_csv(path)
    id_col = "peak_id" if "peak_id" in df.columns else df.columns[0]
    return df.set_index(id_col)


def load_branch_sig(contrast, branch_dir, cell_type, padj=0.05, lfc=1.0):
    """Load significant UP peak IDs from a branch DESeq2 result."""
    df = load_branch_df(contrast, branch_dir, cell_type)
    if df is None:
        print(f"  [MISSING branch] {cell_type} / {contrast}")
        return set()
    df = df.dropna(subset=["padj", "log2FoldChange"])
    sig = df[(df["padj"] < padj) & (df["log2FoldChange"] > lfc)]
    return set(sig.index.astype(str))


def load_pseudobulk_matrices(cell_type, species_list, frag_dir):
    """
    Load pseudobulk count matrices for a given cell type across all species.

    Returns
    -------
    acc_species_df : DataFrame  (peaks × species)
    acc_donor_df   : DataFrame  (peaks × Species|Donor)
    acc_donor_region_df : DataFrame  (peaks × Species|Donor|Region)
    """
    acc_species = {}
    acc_donor = {}
    acc_donor_region = {}

    for sp in species_list:
        counts_path = Path(frag_dir) / sp / "pseudobulk_counts.parquet"
        groups_path = Path(frag_dir) / sp / "pseudobulk_groups.parquet"
        if not counts_path.exists():
            counts_path = counts_path.with_suffix(".tsv")
            groups_path = groups_path.with_suffix(".tsv")

        if not groups_path.exists():
            print(f"  {sp}: groups file missing, skipping")
            continue

        groups = (pd.read_parquet(groups_path) if str(groups_path).endswith(".parquet")
                  else pd.read_csv(groups_path, sep="\t"))

        ct_groups = groups[groups["cell_type"] == cell_type].copy()
        if len(ct_groups) == 0:
            continue

        labels = ct_groups["label"].tolist()
        print(f"  {sp}: {len(labels)} {cell_type} pseudobulks ({ct_groups['n_cells'].sum()} cells)")

        if not counts_path.exists():
            print(f"  {sp}: counts file missing, skipping")
            continue

        df = (pd.read_parquet(counts_path) if str(counts_path).endswith(".parquet")
              else pd.read_csv(counts_path, sep="\t", index_col=0))
        if "region_id" in df.columns:
            df = df.set_index("region_id")

        valid_labels = [l for l in labels if l in df.columns]
        if not valid_labels:
            continue
        sub = df[valid_labels]

        # Species-level: mean across all pseudobulks
        mc = sub.mean(axis=1)
        tot = mc.sum()
        acc_species[sp] = np.log2((mc / tot * 1e6 if tot > 0 else mc) + 1)

        meta = ct_groups.set_index("label")[["donor", "region"]].astype(str)

        # Donor-level
        donor_map = {}
        for lbl in valid_labels:
            d = meta.at[lbl, "donor"] if lbl in meta.index else "Unknown"
            donor_map.setdefault(d, []).append(lbl)
        for donor, dlbls in donor_map.items():
            dc = sub[dlbls].mean(axis=1)
            dt = dc.sum()
            acc_donor[f"{sp}|{donor}"] = np.log2((dc / dt * 1e6 if dt > 0 else dc) + 1)

        # Donor-region-level
        dr_map = {}
        for lbl in valid_labels:
            if lbl in meta.index:
                d, r = meta.at[lbl, "donor"], meta.at[lbl, "region"]
            else:
                d, r = "Unknown", "Unknown"
            dr_map.setdefault(f"{d}|{r}", []).append(lbl)
        for key, klbls in dr_map.items():
            kc = sub[klbls].mean(axis=1)
            kt = kc.sum()
            acc_donor_region[f"{sp}|{key}"] = np.log2((kc / kt * 1e6 if kt > 0 else kc) + 1)

    return (pd.DataFrame(acc_species),
            pd.DataFrame(acc_donor),
            pd.DataFrame(acc_donor_region))


# =============================================================================
# FILTERING & RANKING
# =============================================================================

def ortho_filter(peaks, required_species, anno, orth_cols):
    """Keep only peaks with orthologous sequence in ALL required species."""
    mask = anno[[orth_cols[s] for s in required_species]].apply(
        lambda c: is_true(c)).all(axis=1)
    return set(peaks) & set(anno.index[mask])


def rank_peaks(peak_set, contrast, branch_dir, cell_type, n=None):
    """
    Rank peaks by |LFC| × −log10(padj) from the DESeq2 branch result.
    n=None returns the full ranked list.
    """
    df = load_branch_df(contrast, branch_dir, cell_type)
    if df is None:
        out = list(peak_set)
        return out if n is None else out[:n]
    df = df.loc[df.index.intersection(peak_set)].dropna(
        subset=["padj", "log2FoldChange"])
    df = df.copy()
    df["score"] = df["log2FoldChange"].abs() * (
        -np.log10(df["padj"].clip(lower=1e-300)))
    ranked = list(df.sort_values("score", ascending=False).index)
    return ranked if n is None else ranked[:n]


def rank_by_basemean(peak_set, branch_dir, cell_type, filter_acc_df, n=None):
    """
    Rank peaks by baseMean from an available branch contrast.
    Falls back to species-level accessibility, then arbitrary order.
    n=None returns the full ranked list.
    """
    for contrast in ["Div_Human_vs_Apes", "Node1_Human_vs_Pan"]:
        df = load_branch_df(contrast, branch_dir, cell_type)
        if df is not None:
            overlap = df.index.intersection(peak_set)
            if len(overlap) > 0:
                ranked = list(
                    df.loc[overlap].sort_values("baseMean", ascending=False).index)
                return ranked if n is None else ranked[:n]
    if "Human" in filter_acc_df.columns:
        overlap = filter_acc_df.index.intersection(peak_set)
        ranked = list(
            filter_acc_df.loc[overlap, "Human"].sort_values(ascending=False).index)
        return ranked if n is None else ranked[:n]
    out = list(peak_set)
    return out if n is None else out[:n]


def acc_waterfall(peak_set, up_species, down_species, filter_acc_df, anno, orth_cols):
    """
    Strict waterfall: every down-species that has an ortholog must have
    accessibility strictly below the minimum of up-species.
    """
    peaks = sorted(set(peak_set) & set(filter_acc_df.index))
    if not peaks:
        return set()
    sub = filter_acc_df.loc[peaks]
    up_min = sub[up_species].min(axis=1)
    keep = pd.Series(True, index=peaks)
    for sp in down_species:
        has_orth = is_true(anno.reindex(peaks)[orth_cols[sp]])
        sp_val = sub[sp].fillna(-999)
        keep &= ~(has_orth & (sp_val >= up_min))
    return set(pd.Index(peaks)[keep])


def posthoc_filter(peak_list, up_species, down_species, filter_acc_df, anno, orth_cols,
                   threshold=0.8):
    """
    Extra-strict clean-up: remove peaks where any down-species (with ortholog)
    reaches > threshold × up-species mean accessibility.

    threshold : float  fraction of up-species mean (default 0.8 = 80%)
    """
    if not peak_list:
        return peak_list
    peaks = [p for p in peak_list if p in filter_acc_df.index]
    sub = filter_acc_df.loc[peaks]
    up_mean = sub[up_species].mean(axis=1)
    keep = pd.Series(True, index=peaks)
    for sp in down_species:
        has_orth = is_true(anno.reindex(peaks)[orth_cols[sp]])
        sp_val = sub[sp].fillna(0)
        keep &= ~(has_orth & (sp_val > up_mean * threshold))
    return [p for p in peaks if keep[p]]


# =============================================================================
# HIGH-LEVEL PIPELINE
# =============================================================================

def build_regions(cell_type, anno, filter_acc_df, robust_dir, branch_dir,
                  orth_cols, det_cols, species, max_per_block=50,
                  padj_thresh_sig=0.1, lfc_thresh_sig=0.5,
                  posthoc_threshold=0.8):
    """
    Build all 7 evolutionary peak modules for a given cell type.

    Parameters
    ----------
    padj_thresh_sig   : float  padj threshold used when excluding peaks
                               significant at deeper phylogenetic nodes
                               (default 0.1)
    lfc_thresh_sig    : float  |LFC| threshold for the same exclusion step
                               (default 0.5)
    posthoc_threshold : float  max fraction of up-species mean accessibility
                               allowed in a down-species to pass posthoc_filter
                               (default 0.8 = 80%)

    Notes on ultra-robust pre-filtering (applied in R, not re-applied here):
      Stage 1 (DESeq2): padj < 0.05, |LFC| > 1.0 in the pseudobulk model.
      Stage 2 (complete separation, single-pass min/max across ALL pseudobulks):
        - min_pos_logCPM > max_neg_logCPM  (every positive pseudobulk exceeds
                                            every negative pseudobulk)
        - min_pos_logCPM > 1               (peak genuinely accessible: > 1
                                            log2CPM in all positive pseudobulks)
      This is NOT cross-validation — it is a single-pass min/max check.

    Returns
    -------
    regions          : dict  label → top-N peak list (for plotting)
    regions_all      : dict  label → full filtered peak list (for BED/enrichment)
    total_candidates : dict  label → final candidate count
    filter_stats     : dict  label → {n_candidates, n_after_ortho,
                                      n_after_waterfall, n_final}
    """
    other_sp = [s for s in species if s != "Human"]
    regions          = {}
    regions_all      = {}
    total_candidates = {}
    filter_stats     = {}

    # Pre-load ultra-robust sets
    robust = {}
    for c in [
        "Node4_OldWorld_vs_Marmoset",
        "Node3_GreatApes_vs_Macaque",
        "ILS_HumanGorilla_vs_Pan",
        "ILS_HumanChimp_vs_GorillaBonobo",
        "ILS_HumanBonobo_vs_ChimpGorilla",
        "Div_Human_vs_Apes",
        "Node1_Human_vs_Pan",
    ]:
        robust[c] = load_robust(c, robust_dir, cell_type)

    # Pre-load significant peaks for hierarchical subtraction
    sig_node4 = load_branch_sig("Node4_OldWorld_vs_Marmoset", branch_dir, cell_type,
                                padj=padj_thresh_sig, lfc=lfc_thresh_sig)
    sig_node3 = load_branch_sig("Node3_GreatApes_vs_Macaque", branch_dir, cell_type,
                                padj=padj_thresh_sig, lfc=lfc_thresh_sig)
    deeper_than_ils = (robust["Node3_GreatApes_vs_Macaque"]
                       | robust["Node4_OldWorld_vs_Marmoset"]
                       | sig_node3 | sig_node4)
    deeper_than_n3 = robust["Node4_OldWorld_vs_Marmoset"] | sig_node4

    # Tracks every peak assigned to a module so far — enforces mutual exclusivity.
    # Each peak is assigned to exactly one module (the first it qualifies for
    # in the hierarchy below).  Hierarchy: Human-specific DNA → Human UP vs All
    # → ILS (Gorilla, Chimp, Bonobo) → Great Apes vs Macaque → Old World vs Marmoset
    assigned_peaks: set = set()

    # ── 1. Human-Specific DNA ────────────────────────────────────────────────
    human_spec_mask = (
        is_true(anno[orth_cols["Human"]])
        & anno[[orth_cols[s] for s in other_sp]].apply(
            lambda c: ~is_true(c)).all(axis=1)
        & is_true(anno[det_cols["Human"]])
    )
    human_spec_all_ids = set(anno.index[human_spec_mask])
    ranked_all = rank_by_basemean(human_spec_all_ids, branch_dir, cell_type,
                                  filter_acc_df, n=None)
    regions_all["Human-specific DNA"] = ranked_all
    regions["Human-specific DNA"] = ranked_all[:max_per_block]
    total_candidates["Human-specific DNA"] = len(ranked_all)
    filter_stats["Human-specific DNA"] = {
        "n_candidates": len(human_spec_all_ids),
        "n_after_ortho": len(human_spec_all_ids),
        "n_after_waterfall": len(human_spec_all_ids),
        "n_after_dedup": len(ranked_all),
        "n_final": len(ranked_all),
    }
    assigned_peaks.update(ranked_all)
    print(f"  Human-specific DNA:           {len(ranked_all):>5} -> "
          f"{len(regions['Human-specific DNA'])} shown")

    # ── 2. Human UP vs All Others (Divergence) ───────────────────────────────
    n_robust_div = len(robust["Div_Human_vs_Apes"])
    human_up_all = ortho_filter(
        robust["Div_Human_vs_Apes"],
        ["Human", "Chimpanzee", "Bonobo", "Gorilla"],
        anno, orth_cols,
    )
    human_up_all -= assigned_peaks          # removes human-specific DNA
    n_after_ortho = len(human_up_all)
    human_up_all = acc_waterfall(
        human_up_all, ["Human"], other_sp, filter_acc_df, anno, orth_cols)
    n_after_waterfall = len(human_up_all)
    ranked_all = rank_peaks(human_up_all, "Div_Human_vs_Apes",
                            branch_dir, cell_type, n=None)
    ranked_all = posthoc_filter(
        ranked_all, ["Human"], other_sp, filter_acc_df, anno, orth_cols,
        threshold=posthoc_threshold)
    # Dedup: remove any peak already assigned (should be 0 here, belt-and-suspenders)
    n_before_dedup = len(ranked_all)
    ranked_all = [p for p in ranked_all if p not in assigned_peaks]
    regions_all["Human UP vs All"] = ranked_all
    regions["Human UP vs All"] = ranked_all[:max_per_block]
    total_candidates["Human UP vs All"] = len(ranked_all)
    filter_stats["Human UP vs All"] = {
        "n_candidates": n_robust_div,
        "n_after_ortho": n_after_ortho,
        "n_after_waterfall": n_after_waterfall,
        "n_after_dedup": len(ranked_all),
        "n_final": len(ranked_all),
    }
    assigned_peaks.update(ranked_all)
    print(f"  Human UP vs All:              {len(ranked_all):>5} -> "
          f"{len(regions['Human UP vs All'])} shown")

    # ── 3–5. ILS contrasts ───────────────────────────────────────────────────
    ils_defs = [
        ("ILS: Human + Gorilla UP", "ILS_HumanGorilla_vs_Pan",
         ["Human", "Gorilla"], ["Chimpanzee", "Bonobo"]),
        ("ILS: Human + Chimp UP", "ILS_HumanChimp_vs_GorillaBonobo",
         ["Human", "Chimpanzee"], ["Gorilla", "Bonobo"]),
        ("ILS: Human + Bonobo UP", "ILS_HumanBonobo_vs_ChimpGorilla",
         ["Human", "Bonobo"], ["Chimpanzee", "Gorilla"]),
    ]
    for label, contrast, up_sp, contrast_down_sp in ils_defs:
        n_robust_c = len(robust[contrast])
        filtered = ortho_filter(
            robust[contrast], up_sp + contrast_down_sp, anno, orth_cols)
        filtered -= deeper_than_ils
        filtered -= assigned_peaks          # removes modules 1, 2, and prior ILS
        n_after_ortho_c = len(filtered)
        all_down = contrast_down_sp + ["Macaque", "Marmoset"]
        filtered = acc_waterfall(
            filtered, up_sp, all_down, filter_acc_df, anno, orth_cols)
        n_after_waterfall_c = len(filtered)
        ranked_all = rank_peaks(filtered, contrast, branch_dir, cell_type, n=None)
        ranked_all = posthoc_filter(
            ranked_all, up_sp, all_down, filter_acc_df, anno, orth_cols,
            threshold=posthoc_threshold)
        n_before_dedup_c = len(ranked_all)
        ranked_all = [p for p in ranked_all if p not in assigned_peaks]
        regions_all[label] = ranked_all
        regions[label] = ranked_all[:max_per_block]
        total_candidates[label] = len(ranked_all)
        filter_stats[label] = {
            "n_candidates": n_robust_c,
            "n_after_ortho": n_after_ortho_c,
            "n_after_waterfall": n_after_waterfall_c,
            "n_after_dedup": len(ranked_all),
            "n_final": len(ranked_all),
        }
        assigned_peaks.update(ranked_all)
        print(f"  {label}:   {len(ranked_all):>5} -> {len(regions[label])} shown")

    # ── 6. Great Apes UP vs Macaque ──────────────────────────────────────────
    n_robust_n3 = len(robust["Node3_GreatApes_vs_Macaque"])
    n3_filtered = ortho_filter(
        robust["Node3_GreatApes_vs_Macaque"],
        ["Human", "Chimpanzee", "Bonobo", "Gorilla", "Macaque"],
        anno, orth_cols,
    )
    n3_filtered -= deeper_than_n3
    n3_filtered -= assigned_peaks           # removes modules 1–5
    n_after_ortho_n3 = len(n3_filtered)
    n3_filtered = acc_waterfall(
        n3_filtered,
        ["Human", "Chimpanzee", "Bonobo", "Gorilla"], ["Macaque", "Marmoset"],
        filter_acc_df, anno, orth_cols,
    )
    n_after_waterfall_n3 = len(n3_filtered)
    ranked_all = rank_peaks(n3_filtered, "Node3_GreatApes_vs_Macaque",
                            branch_dir, cell_type, n=None)
    ranked_all = posthoc_filter(
        ranked_all,
        ["Human", "Chimpanzee", "Bonobo", "Gorilla"], ["Macaque", "Marmoset"],
        filter_acc_df, anno, orth_cols,
        threshold=posthoc_threshold,
    )
    ranked_all = [p for p in ranked_all if p not in assigned_peaks]
    regions_all["Great Apes UP vs Macaque"] = ranked_all
    regions["Great Apes UP vs Macaque"] = ranked_all[:max_per_block]
    total_candidates["Great Apes UP vs Macaque"] = len(ranked_all)
    filter_stats["Great Apes UP vs Macaque"] = {
        "n_candidates": n_robust_n3,
        "n_after_ortho": n_after_ortho_n3,
        "n_after_waterfall": n_after_waterfall_n3,
        "n_after_dedup": len(ranked_all),
        "n_final": len(ranked_all),
    }
    assigned_peaks.update(ranked_all)
    print(f"  Great Apes UP vs Macaque:     {len(ranked_all):>5} -> "
          f"{len(regions['Great Apes UP vs Macaque'])} shown")

    # ── 7. Old World UP vs Marmoset ──────────────────────────────────────────
    n_robust_n4 = len(robust["Node4_OldWorld_vs_Marmoset"])
    n4_filtered = ortho_filter(
        robust["Node4_OldWorld_vs_Marmoset"], species, anno, orth_cols)
    n4_filtered -= assigned_peaks           # removes modules 1–6
    n_after_ortho_n4 = len(n4_filtered)
    n4_filtered = acc_waterfall(
        n4_filtered,
        ["Human", "Chimpanzee", "Bonobo", "Gorilla", "Macaque"], ["Marmoset"],
        filter_acc_df, anno, orth_cols,
    )
    n_after_waterfall_n4 = len(n4_filtered)
    ranked_all = rank_peaks(n4_filtered, "Node4_OldWorld_vs_Marmoset",
                            branch_dir, cell_type, n=None)
    ranked_all = posthoc_filter(
        ranked_all,
        ["Human", "Chimpanzee", "Bonobo", "Gorilla", "Macaque"], ["Marmoset"],
        filter_acc_df, anno, orth_cols,
        threshold=posthoc_threshold,
    )
    ranked_all = [p for p in ranked_all if p not in assigned_peaks]
    regions_all["Old World UP vs Marmoset"] = ranked_all
    regions["Old World UP vs Marmoset"] = ranked_all[:max_per_block]
    total_candidates["Old World UP vs Marmoset"] = len(ranked_all)
    filter_stats["Old World UP vs Marmoset"] = {
        "n_candidates": n_robust_n4,
        "n_after_ortho": n_after_ortho_n4,
        "n_after_waterfall": n_after_waterfall_n4,
        "n_after_dedup": len(ranked_all),
        "n_final": len(ranked_all),
    }
    assigned_peaks.update(ranked_all)
    print(f"  Old World UP vs Marmoset:     {len(ranked_all):>5} -> "
          f"{len(regions['Old World UP vs Marmoset'])} shown")

    # Summary
    print(f"\n{'='*62}")
    print(f"  {'Category':<35} {'Candidates':>10} {'Shown':>7}")
    print(f"  {'-'*58}")
    all_shown = set()
    for label, peaks in regions.items():
        all_shown |= set(peaks)
        print(f"  {label:<35} {total_candidates.get(label, len(peaks)):>10} "
              f"{len(peaks):>7}")
    print(f"  {'-'*58}")
    print(f"  {'Total unique peaks shown':<35} {'':>10} {len(all_shown):>7}")

    return regions, regions_all, total_candidates, filter_stats


def build_heatmap_matrix(regions, acc_df, anno, orth_cols, apply_row_zscore=False):
    """
    Build the heatmap data matrix from a regions dict and an accessibility df.

    Only peaks present in acc_df are included (so block_sizes exactly match
    the number of plotted rows). NaN is enforced where ortholog is absent.

    Returns
    -------
    heatmap_plot_df : DataFrame  (ordered peaks × columns)
    block_labels    : list of str
    block_sizes     : list of int  (actual rows per block)
    """
    acc_index = set(acc_df.index)
    ordered_peaks = []
    block_labels = []
    block_sizes = []

    for label, peaks in regions.items():
        valid = [p for p in peaks if p in acc_index]
        ordered_peaks.extend(valid)
        block_labels.append(label)
        block_sizes.append(len(valid))

    heatmap_df = acc_df.reindex(ordered_peaks).copy()

    for col in heatmap_df.columns:
        sp = col.split("|", 1)[0]
        if sp not in orth_cols:
            continue
        orth_mask = is_true(anno.reindex(ordered_peaks)[orth_cols[sp]])
        heatmap_df.loc[~orth_mask.values, col] = np.nan

    if apply_row_zscore:
        heatmap_plot_df = row_zscore(heatmap_df)
    else:
        heatmap_plot_df = heatmap_df.copy()

    return heatmap_plot_df, block_labels, block_sizes


def plot_evolutionary_heatmap(heatmap_plot_df, block_labels, block_sizes,
                               cell_type, title_mode, out_dir,
                               apply_row_zscore=False):
    """
    Plot and save a blocked evolutionary accessibility heatmap.

    Saves {cell_type}_evolutionary_heatmap.pdf and .png into out_dir.
    Returns the matplotlib Figure.
    """
    n_peaks = heatmap_plot_df.shape[0]
    fig_h = max(6, n_peaks * 0.04 + 2)
    fig, (ax, cax) = plt.subplots(
        1, 2, figsize=(9, fig_h),
        gridspec_kw={"width_ratios": [30, 1], "wspace": 0.35},
    )

    cmap = (plt.cm.RdBu_r if apply_row_zscore else plt.cm.YlOrRd).copy()
    cmap.set_bad(color="#e0e0e0")

    data = heatmap_plot_df.values.astype(float)

    if apply_row_zscore:
        abs_max = np.nanpercentile(np.abs(data[~np.isnan(data)]), 98)
        abs_max = max(abs_max, 1e-6)
        norm = mcolors.TwoSlopeNorm(vmin=-abs_max, vcenter=0.0, vmax=abs_max)
        im = ax.imshow(data, aspect="auto", cmap=cmap, norm=norm,
                       interpolation="nearest")
    else:
        vmin, vmax = np.nanpercentile(data[~np.isnan(data)], [2, 98])
        im = ax.imshow(data, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax,
                       interpolation="nearest")

    ax.set_xticks(range(len(heatmap_plot_df.columns)))
    ax.set_xticklabels(list(heatmap_plot_df.columns),
                       rotation=60, ha="right", fontsize=9)
    ax.xaxis.set_ticks_position("bottom")

    cum = 0
    label_positions = []
    for i, (label, size) in enumerate(zip(block_labels, block_sizes)):
        if i > 0:
            ax.axhline(cum - 0.5, color="white", linewidth=2.5)
            ax.axhline(cum - 0.5, color="black", linewidth=0.8)
        label_positions.append(cum + size / 2)
        cum += size

    ax.set_yticks(label_positions)
    ax.set_yticklabels(
        [f"{l} (n={s})" for l, s in zip(block_labels, block_sizes)],
        fontsize=9, fontweight="bold",
    )
    ax.tick_params(axis="y", length=0, pad=8)

    cb = fig.colorbar(im, cax=cax)
    cb.set_label("Row z-score" if apply_row_zscore else "log2(CPM + 1)",
                 fontsize=10)

    legend_elements = [
        Patch(facecolor="#e0e0e0", edgecolor="black", linewidth=0.5,
              label="No ortholog")
    ]
    ax.legend(handles=legend_elements, loc="upper right", fontsize=8,
              framealpha=0.9, edgecolor="gray")

    ax.set_title(
        f"{cell_type} — Evolutionary Accessibility Modules ({title_mode})",
        fontsize=12, fontweight="bold", pad=12,
    )

    out_dir = Path(out_dir)
    ct_clean = sanitize_name(cell_type)
    for ext in ["pdf", "png"]:
        fig.savefig(out_dir / f"{ct_clean}_evolutionary_heatmap.{ext}",
                    bbox_inches="tight", dpi=150)
    print(f"  Heatmap saved → {out_dir}/{ct_clean}_evolutionary_heatmap.pdf")
    plt.close(fig)
    return fig


def export_beds(regions_all, anno, out_dir, cell_type, coord_system="Human",
                export_tag="strict_waterfall", export_mode="full_candidates"):
    """
    Export each module in regions_all as a sorted BED file plus a manifest TSV.

    Coordinates are taken from {coord_system}_chr/start/end columns in anno.
    Returns list of (module_name, path, n_peaks) tuples.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    species_triplet = (f"{coord_system}_chr",
                       f"{coord_system}_start",
                       f"{coord_system}_end")
    generic_candidates = [
        ("chr", "start", "end"),
        ("chrom", "start", "end"),
        ("seqnames", "start", "end"),
    ]

    if all(c in anno.columns for c in species_triplet):
        chrom_col, start_col, end_col = species_triplet
        coord_label = coord_system
    else:
        found = next((cols for cols in generic_candidates
                      if all(c in anno.columns for c in cols)), None)
        if found is None:
            available = sorted(c[:-4] for c in anno.columns if c.endswith("_chr"))
            raise ValueError(
                f"No usable coordinate columns for '{coord_system}'. "
                f"Available species: {available}"
            )
        chrom_col, start_col, end_col = found
        coord_label = "generic"

    cell_clean = sanitize_name(cell_type)
    coord_clean = sanitize_name(coord_label)
    mode_clean = sanitize_name(export_mode)
    tag_clean = sanitize_name(export_tag)

    written = []
    for module_name, peak_ids in regions_all.items():
        if not peak_ids:
            continue
        peak_ids = [p for p in peak_ids if p in anno.index]
        if not peak_ids:
            print(f"  [SKIP] {module_name}: no peaks in annotation index")
            continue

        bed = anno.loc[peak_ids, [chrom_col, start_col, end_col]].copy()
        bed = bed.rename(columns={chrom_col: "chrom",
                                  start_col: "start",
                                  end_col: "end"})
        bed["name"] = peak_ids
        bed = bed[["chrom", "start", "end", "name"]].dropna()
        if bed.empty:
            print(f"  [SKIP] {module_name}: no non-null coordinates")
            continue
        bed["start"] = bed["start"].astype(int)
        bed["end"] = bed["end"].astype(int)
        bed = bed[bed["end"] > bed["start"]]
        bed = bed.sort_values(["chrom", "start", "end"])

        module_clean = sanitize_name(module_name)
        fname = (f"{cell_clean}__{tag_clean}__{mode_clean}"
                 f"__{coord_clean}__{module_clean}.bed")
        out_path = out_dir / fname
        bed.to_csv(out_path, sep="\t", header=False, index=False)
        written.append((module_name, out_path, len(bed)))
        print(f"  [WROTE] {module_name:<35} {len(bed):>5} peaks → {out_path}")

    if written:
        manifest = pd.DataFrame([
            {"module": m, "n_peaks": n, "coord_system": coord_label,
             "export_mode": export_mode, "bed_path": str(p)}
            for m, p, n in written
        ])
        mpath = (out_dir
                 / f"{cell_clean}__{tag_clean}__{mode_clean}"
                   f"__{coord_clean}__manifest.tsv")
        manifest.to_csv(mpath, sep="\t", index=False)
        print(f"  Manifest → {mpath}")

    return written


# =============================================================================
# SNC ENRICHMENT
# =============================================================================

def load_snc_bedtool(snc_path, filter_fixed=True):
    """
    Load the Prufer SNC BED file as a sorted pybedtools.BedTool.

    Parameters
    ----------
    snc_path     : str  path to Prufer_SNC_hg38_full.bed
    filter_fixed : bool  if True (default), keep only rows where
                   modern_human_specific_and_fixed_0.9 == 1
    """
    if not HAS_PYBEDTOOLS:
        raise ImportError("pybedtools is required for SNC enrichment analysis")
    snc_df = pd.read_csv(snc_path, sep="\t", header=0, low_memory=False,
                         dtype={"seqnames": str})
    print(f"  SNC file loaded: {len(snc_df):,} rows")
    if filter_fixed:
        snc_df = snc_df[snc_df["modern_human_specific_and_fixed_0.9"] == 1]
        print(f"  SNCs after fixed filter (>=0.9): {len(snc_df):,}")
    snc_df = snc_df[["seqnames", "start", "end"]].copy()
    snc_df.columns = ["chrom", "start", "end"]
    snc_df["start"] = snc_df["start"].astype(int)
    snc_df["end"]   = snc_df["end"].astype(int)
    bt = pybedtools.BedTool.from_dataframe(snc_df)
    return bt.sort()


def build_snc_backgrounds(anno, det_cols, species):
    """
    Return contrast-specific background peak sets based on which peaks were
    detectable in the species relevant to each DESeq2 contrast.

    Backgrounds:
      "Human-specific DNA"       → detected in Human only
                                   (expected negative control: no SNCs here)
      "Human UP vs All"          → detected in all 6 species
      ILS contrasts              → detected in all 4 great apes
                                   (Human, Bonobo, Chimpanzee, Gorilla)
      "Great Apes UP vs Macaque" → detected in all 5 (great apes + Macaque)
      "Old World UP vs Marmoset" → detected in all 6 species
    """
    def det_mask(sp_list):
        mask = pd.Series(True, index=anno.index)
        for sp in sp_list:
            mask &= (anno[det_cols[sp]] == 1)
        return list(anno.index[mask])

    great_apes = ["Human", "Bonobo", "Chimpanzee", "Gorilla"]

    return {
        "Human-specific DNA":       det_mask(["Human"]),
        "Human UP vs All":          det_mask(species),
        "ILS: Human + Gorilla UP":  det_mask(great_apes),
        "ILS: Human + Chimp UP":    det_mask(great_apes),
        "ILS: Human + Bonobo UP":   det_mask(great_apes),
        "Great Apes UP vs Macaque": det_mask(great_apes + ["Macaque"]),
        "Old World UP vs Marmoset": det_mask(species),
    }


def precompute_snc_counts(anno, snc_bt, out_path, coord_system="Human",
                          flank_bp=0):
    """
    Intersect ALL peaks in the annotation with the SNC BedTool, recording
    the count of SNCs per peak.  Result is cached as a parquet file; if the
    cache already exists it is loaded directly (pybedtools not required).

    The peak_id is stored as the BED 'name' column so pybedtools order is
    irrelevant — the result is always correctly aligned.

    Parameters
    ----------
    anno        : DataFrame  master annotation (peak_id index)
    snc_bt      : pybedtools.BedTool  sorted, filtered SNC intervals
    out_path    : str | Path  where to save / load the parquet cache
    coord_system: str  species coordinate system (default "Human")
    flank_bp    : int  bp to expand each peak on both sides before
                  intersecting (default 0 = exact peak coordinates).
                  Use different out_path values for different flank_bp
                  values so caches do not collide.

    Returns
    -------
    DataFrame indexed by peak_id, column "n_snc" (int ≥ 0).
    """
    out_path = Path(out_path)
    if out_path.exists():
        print(f"  Loading cached SNC counts from {out_path}")
        return pd.read_parquet(out_path)

    if not HAS_PYBEDTOOLS:
        raise ImportError("pybedtools is required for precompute_snc_counts")

    chrom_col = f"{coord_system}_chr"
    start_col = f"{coord_system}_start"
    end_col   = f"{coord_system}_end"

    # 4-column BED with peak_id as the name column
    bed_df = anno[[chrom_col, start_col, end_col]].copy().dropna()
    bed_df.columns = ["chrom", "start", "end"]
    bed_df.insert(3, "name", bed_df.index.astype(str))
    bed_df["start"] = bed_df["start"].astype(int)
    bed_df["end"]   = bed_df["end"].astype(int)
    bed_df = bed_df[bed_df["end"] > bed_df["start"]]

    if flank_bp > 0:
        bed_df["start"] = (bed_df["start"] - flank_bp).clip(lower=0)
        bed_df["end"]   =  bed_df["end"] + flank_bp

    bed_df = bed_df.sort_values(["chrom", "start"])

    flank_str = f" (flank ±{flank_bp} bp)" if flank_bp > 0 else ""
    print(f"  Building BedTool from {len(bed_df):,} peaks{flank_str}...")
    bt = pybedtools.BedTool.from_dataframe(bed_df)
    print(f"  Intersecting with {len(snc_bt):,} SNCs (c=True)...", flush=True)
    counted = bt.intersect(snc_bt, c=True)

    result = counted.to_dataframe(
        names=["chrom", "start", "end", "peak_id", "n_snc"])
    result = result[["peak_id", "n_snc"]].copy()
    result["n_snc"] = result["n_snc"].astype(int)
    result = result.set_index("peak_id")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_parquet(out_path)
    print(f"  SNC counts ({len(result):,} peaks) saved → {out_path}")
    return result


def _peaks_to_bedtool(peak_ids, anno, coord_system="Human"):
    """
    Convert a list of peak IDs to a sorted pybedtools.BedTool using
    the {coord_system}_chr/start/end columns from anno.
    Returns (BedTool, valid_peak_ids) — valid_peak_ids preserves order and
    excludes peaks with null coordinates.
    """
    chrom_col = f"{coord_system}_chr"
    start_col = f"{coord_system}_start"
    end_col   = f"{coord_system}_end"

    sub = anno.reindex(peak_ids)[[chrom_col, start_col, end_col]].dropna()
    sub = sub.rename(columns={chrom_col: "chrom",
                               start_col: "start",
                               end_col: "end"})
    sub["start"] = sub["start"].astype(int)
    sub["end"]   = sub["end"].astype(int)
    sub = sub[sub["end"] > sub["start"]]
    if sub.empty:
        return None, []
    bt = pybedtools.BedTool.from_dataframe(sub.reset_index(drop=True))
    return bt.sort(), list(sub.index)


def _get_snc_flags(peak_ids, anno, snc_bt, coord_system="Human"):
    """
    For a list of peak IDs, return a numpy bool array indicating which peaks
    overlap ≥1 SNC.  The array is aligned to the order of peak_ids; peaks with
    no valid coordinates are treated as non-overlapping (False).
    """
    bt, valid_ids = _peaks_to_bedtool(peak_ids, anno, coord_system)
    if bt is None:
        return np.zeros(len(peak_ids), dtype=bool)

    # Intersect and get the set of valid_ids that have a hit
    hit_bt = bt.intersect(snc_bt, u=True)
    if len(hit_bt) == 0:
        hit_set = set()
    else:
        hit_df = hit_bt.to_dataframe(names=["chrom", "start", "end"])
        # Reconstruct which peak_ids had hits by matching coordinates
        valid_coords = anno.reindex(valid_ids)[
            [f"{coord_system}_chr", f"{coord_system}_start", f"{coord_system}_end"]
        ].copy()
        valid_coords.columns = ["chrom", "start", "end"]
        valid_coords["start"] = valid_coords["start"].astype(int)
        valid_coords["end"]   = valid_coords["end"].astype(int)
        # merge on coordinates to find which peak_ids are hits
        hit_coords = set(zip(hit_df["chrom"], hit_df["start"], hit_df["end"]))
        hit_set = {
            pid for pid, row in valid_coords.iterrows()
            if (row["chrom"], row["start"], row["end"]) in hit_coords
        }

    flags = np.array([pid in hit_set for pid in peak_ids], dtype=bool)
    return flags


def _count_sncs_in_peaks(peak_ids, anno, snc_bt, coord_system="Human"):
    """
    Count total SNCs (not just presence/absence) across a set of peaks using
    pybedtools intersect(c=True).  Returns total SNC count.
    """
    bt, _ = _peaks_to_bedtool(peak_ids, anno, coord_system)
    if bt is None:
        return 0
    counted = bt.intersect(snc_bt, c=True)
    counts_df = counted.to_dataframe()
    return int(counts_df.iloc[:, -1].sum())


def compute_snc_enrichment(cell_type, regions_all, module_backgrounds,
                            snc_counts_df=None,
                            anno=None, snc_bt=None,
                            n_permutations=10_000,
                            coord_system="Human",
                            peak_width_bp=500):
    """
    Test SNC enrichment for each module using contrast-specific backgrounds.

    For each module computes:
      - Presence/absence overlap + Fisher's exact test (one-tailed)
      - SNC density (SNCs per kb, all peaks are peak_width_bp wide)
      - Empirical permutation p-value (random resampling from background)

    Parameters
    ----------
    cell_type          : str
    regions_all        : dict  label → full peak list
    module_backgrounds : dict  label → list of background peak IDs
                         (from build_snc_backgrounds)
    snc_counts_df      : DataFrame  peak_id-indexed, "n_snc" column
                         (from precompute_snc_counts) — preferred fast path
    anno               : DataFrame  master annotation; required only when
                         snc_counts_df is None (legacy pybedtools path)
    snc_bt             : pybedtools.BedTool; required only when
                         snc_counts_df is None (legacy pybedtools path)
    n_permutations     : int  number of permutation draws (default 10,000)
    coord_system       : str  used only for legacy pybedtools path
    peak_width_bp      : int  peak width in bp (default 500)

    Returns
    -------
    List of dicts with enrichment stats per module.
    """
    if snc_counts_df is None:
        if anno is None or snc_bt is None:
            raise ValueError(
                "Provide snc_counts_df (preferred) OR both anno + snc_bt")
        if not HAS_PYBEDTOOLS:
            raise ImportError("pybedtools required when snc_counts_df is None")
    if not HAS_STATS:
        raise ImportError("scipy and statsmodels required")

    rng = np.random.default_rng(seed=42)
    rows = []

    for module in module_backgrounds.keys():
        bg_peak_ids = module_backgrounds.get(module, [])
        n_bg = len(bg_peak_ids)

        fg_peak_ids = regions_all.get(module, [])
        n_fg = len(fg_peak_ids)

        empty_row = {
            "cell_type": cell_type, "module": module,
            "n_fg": n_fg, "n_fg_snc": 0, "pct_fg_snc": np.nan,
            "n_bg": n_bg, "n_bg_snc": np.nan, "pct_bg_snc": np.nan,
            "fold_enrichment": np.nan,
            "density_fg_snc_per_kb": np.nan,
            "density_bg_snc_per_kb": np.nan,
            "density_fold_enrichment": np.nan,
            "odds_ratio": np.nan,
            "fisher_pvalue": np.nan,
            "perm_pvalue": np.nan,
            "mean_snc_fg": np.nan,
            "mean_snc_bg": np.nan,
        }

        if n_fg == 0 or n_bg == 0:
            rows.append(empty_row)
            continue

        # ── Get SNC flags for the full background ────────────────────────────
        if snc_counts_df is not None:
            # Fast path: use precomputed per-peak counts (no pybedtools)
            bg_snc_counts = (snc_counts_df
                             .reindex(bg_peak_ids, fill_value=0)["n_snc"]
                             .to_numpy(dtype=int))
            bg_snc_flags = bg_snc_counts >= 1
        else:
            # Legacy path: pybedtools intersect for this background
            print(f"    {module}: intersecting background ({n_bg:,} peaks)...",
                  end=" ", flush=True)
            bg_snc_flags = _get_snc_flags(bg_peak_ids, anno, snc_bt, coord_system)
            print(f"{bg_snc_flags.sum()} with SNC")

        n_bg_snc = int(bg_snc_flags.sum())
        # Array used for the permutation: actual integer counts when available
        # (more sensitive than binary flags and genuinely different from Fisher's);
        # falls back to binary integers in the legacy pybedtools path.
        _perm_arr   = (bg_snc_counts
                       if snc_counts_df is not None
                       else bg_snc_flags.astype(int))
        mean_snc_bg = float(_perm_arr.mean())

        # ── Foreground subset ─────────────────────────────────────────────────
        bg_id_to_idx = {pid: i for i, pid in enumerate(bg_peak_ids)}
        fg_in_bg = [p for p in fg_peak_ids if p in bg_id_to_idx]

        if not fg_in_bg:
            rows.append(empty_row)
            continue

        fg_indices = np.array([bg_id_to_idx[p] for p in fg_in_bg])
        n_fg_valid = len(fg_indices)
        n_fg_snc = int(bg_snc_flags[fg_indices].sum())

        # ── Fisher's exact test ───────────────────────────────────────────────
        a = n_fg_snc
        b = n_fg_valid - n_fg_snc
        c = max(n_bg_snc - n_fg_snc, 0)
        d = max(n_bg - n_bg_snc - b, 0)
        odds_ratio, fisher_pvalue = fisher_exact([[a, b], [c, d]],
                                                 alternative="greater")
        # Use BG-only rate (excluding FG peaks) to match the Fisher's exact test,
        # which compares FG vs BG-only in its 2x2 table.
        n_bg_only = n_bg - n_fg_valid
        n_bg_only_snc = n_bg_snc - n_fg_snc
        fe = ((n_fg_snc / n_fg_valid) / (n_bg_only_snc / n_bg_only)
              if (n_bg_only > 0 and n_bg_only_snc > 0) else np.nan)

        # ── SNC density (SNCs per kb) ─────────────────────────────────────────
        if snc_counts_df is not None:
            total_fg_sncs = int(
                snc_counts_df.reindex(fg_in_bg, fill_value=0)["n_snc"].sum())
            total_bg_sncs = int(
                snc_counts_df.reindex(bg_peak_ids, fill_value=0)["n_snc"].sum())
        else:
            total_fg_sncs = _count_sncs_in_peaks(fg_in_bg, anno, snc_bt, coord_system)
            total_bg_sncs = _count_sncs_in_peaks(bg_peak_ids, anno, snc_bt, coord_system)

        kb = peak_width_bp / 1000
        density_fg = total_fg_sncs / (n_fg_valid * kb) if n_fg_valid > 0 else np.nan
        density_bg = total_bg_sncs / (n_bg * kb) if n_bg > 0 else np.nan
        # BG-only density (excluding FG) for fold enrichment, consistent with fe above.
        total_bg_only_sncs = total_bg_sncs - total_fg_sncs
        density_bg_only = (total_bg_only_sncs / (n_bg_only * kb)
                           if n_bg_only > 0 else np.nan)
        density_fe = (density_fg / density_bg_only
                      if (density_bg_only and density_bg_only > 0) else np.nan)

        # ── Permutation test (mean SNC count) ────────────────────────────────
        # Tests whether the mean n_snc per foreground peak exceeds what is
        # expected from randomly drawing the same number of peaks from the
        # background.  Using actual counts (not binary flags) is more sensitive
        # than presence/absence and is genuinely distinct from Fisher's exact,
        # which already tests the hypergeometric null for binary overlap.
        obs_mean_snc = float(_perm_arr[fg_indices].mean())
        perm_means = np.array([
            _perm_arr[rng.choice(n_bg, n_fg_valid, replace=False)].mean()
            for _ in range(n_permutations)
        ])
        perm_pvalue = max((perm_means >= obs_mean_snc).mean(), 1.0 / n_permutations)

        rows.append({
            "cell_type": cell_type,
            "module": module,
            "n_fg": n_fg_valid,
            "n_fg_snc": n_fg_snc,
            "pct_fg_snc": round(100 * n_fg_snc / n_fg_valid, 2),
            "n_bg": n_bg,
            "n_bg_snc": n_bg_snc,
            "pct_bg_snc": round(100 * n_bg_snc / n_bg, 2),
            "fold_enrichment": round(fe, 3) if not np.isnan(fe) else np.nan,
            "density_fg_snc_per_kb": round(density_fg, 4) if not np.isnan(density_fg) else np.nan,
            "density_bg_snc_per_kb": round(density_bg, 4) if not np.isnan(density_bg) else np.nan,
            "density_fold_enrichment": round(density_fe, 3) if not np.isnan(density_fe) else np.nan,
            "odds_ratio": round(odds_ratio, 3),
            "fisher_pvalue": fisher_pvalue,
            "perm_pvalue": perm_pvalue,
            "mean_snc_fg": round(obs_mean_snc, 4),
            "mean_snc_bg": round(mean_snc_bg, 4),
        })

    return rows


def plot_snc_enrichment_dotplot(enrich_df, out_path, padj_col="fisher_padj",
                                 fe_col="fold_enrichment", pval_threshold=0.05):
    """
    Plot a dot plot of SNC enrichment: rows = modules, columns = cell types.
    Dot color = log2(fold enrichment), dot size = −log10(padj).
    Gray dots = not significant (padj ≥ pval_threshold).
    Saves to out_path as PDF.
    """
    modules = enrich_df["module"].unique().tolist()
    cell_types = enrich_df["cell_type"].unique().tolist()

    fig, ax = plt.subplots(figsize=(max(6, len(cell_types) * 0.8 + 1),
                                    max(3, len(modules) * 0.6 + 1)))

    cmap = plt.cm.RdBu_r
    fe_vals = enrich_df[fe_col].replace([np.inf, -np.inf], np.nan).dropna()
    lfe_abs = np.log2(fe_vals.clip(lower=0.01)).abs()
    vmax = np.nanpercentile(lfe_abs, 95) if len(lfe_abs) > 0 else 2
    vmax = max(vmax, 0.5)
    norm = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)

    for mi, module in enumerate(modules):
        for ci, ct in enumerate(cell_types):
            sub = enrich_df[(enrich_df["module"] == module) &
                            (enrich_df["cell_type"] == ct)]
            if sub.empty:
                continue
            row = sub.iloc[0]
            fe = row[fe_col]
            p = row[padj_col] if padj_col in row.index else row["pvalue"]
            if np.isnan(fe) or np.isnan(p):
                continue
            lfe = np.log2(max(fe, 0.01))
            size = min(200, -np.log10(max(p, 1e-10)) * 15)
            color = cmap(norm(lfe)) if p < pval_threshold else "#cccccc"
            ax.scatter(ci, mi, s=size, c=[color], zorder=3, linewidths=0.3,
                       edgecolors="black" if p < pval_threshold else "none")

    ax.set_xticks(range(len(cell_types)))
    ax.set_xticklabels(cell_types, rotation=45, ha="right", fontsize=9)
    ax.set_yticks(range(len(modules)))
    ax.set_yticklabels(modules, fontsize=9)
    ax.set_xlim(-0.5, len(cell_types) - 0.5)
    ax.set_ylim(-0.5, len(modules) - 0.5)
    ax.grid(True, linewidth=0.3, alpha=0.4)
    ax.set_axisbelow(True)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, shrink=0.5, pad=0.02)
    cb.set_label("log2(fold enrichment)", fontsize=9)

    # Size legend
    for s_val, label in [(15, "p=0.1"), (50, "p=0.01"), (150, "p=0.001")]:
        ax.scatter([], [], s=s_val, c="gray", label=label, linewidths=0)
    ax.legend(title="−log10(padj) scale", fontsize=7, title_fontsize=7,
              loc="upper right", framealpha=0.8)

    ax.set_title("SNC Enrichment in Human-Upregulated Modules", fontsize=11,
                 fontweight="bold", pad=10)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"  Dotplot saved → {out_path}")
    return fig


_MODULE_TO_CONTRAST = {
    "Human UP vs All":          "Div_Human_vs_Apes",
    "ILS: Human + Gorilla UP":  "ILS_HumanGorilla_vs_Pan",
    "ILS: Human + Chimp UP":    "ILS_HumanChimp_vs_GorillaBonobo",
    "ILS: Human + Bonobo UP":   "ILS_HumanBonobo_vs_ChimpGorilla",
    "Great Apes UP vs Macaque": "Node3_GreatApes_vs_Macaque",
    "Old World UP vs Marmoset": "Node4_OldWorld_vs_Marmoset",
}


def plot_snc_vs_stats(regions_all_df, snc_counts_df, module_backgrounds,
                      branch_dir, out_dir, module_to_contrast=None):
    """
    Generate SNC-vs-DESeq2-stats plots in two sets.

    **Per-cell-type** (one PDF per cell type):
      rows = modules that have a DESeq2 contrast
      left panel  : boxplot of log2FoldChange grouped by n_snc category
      right panel : boxplot of sign(LFC)·−log10(padj) grouped by n_snc category

    **Per-module** (one PDF per module):
      facets = cell types
      each facet: boxplot of log2FoldChange grouped by n_snc category

    x-axis in each panel: n_snc per peak (0 / 1 / 2 / 3+ grouping).
    Shows whether SNC burden within foreground peaks predicts differential
    accessibility strength — i.e. do peaks with more SNCs show stronger signal?

    Parameters
    ----------
    regions_all_df     : DataFrame  long format (cell_type, module, peak_id, rank)
    snc_counts_df      : DataFrame  peak_id-indexed, "n_snc" column
    module_backgrounds : dict  (unused here, kept for API consistency)
    branch_dir         : str | Path  DESeq2 branch results directory
    out_dir            : str | Path  output root
    module_to_contrast : dict  override module-to-contrast-name mapping
    """
    MAX_SNV_CAT = 3   # groups: 0, 1, 2, 3+  (everything ≥ MAX_SNV_CAT is pooled)
    MIN_PER_CAT = 3   # skip a category if it has fewer peaks than this

    def _boxplot_groups(df, stat_col):
        """
        Split df rows by n_snc category, return (labels, data_arrays) for
        matplotlib boxplot.  Categories with < MIN_PER_CAT points are skipped.
        """
        labels, data = [], []
        for k in range(MAX_SNV_CAT + 1):
            if k < MAX_SNV_CAT:
                mask = df["n_snc"] == k
                lab  = str(k)
            else:
                mask = df["n_snc"] >= MAX_SNV_CAT
                lab  = f"{MAX_SNV_CAT}+"
            vals = df.loc[mask, stat_col].dropna().values
            if len(vals) >= MIN_PER_CAT:
                labels.append(f"{lab}\n(n={len(vals)})")
                data.append(vals)
        return labels, data

    if module_to_contrast is None:
        module_to_contrast = _MODULE_TO_CONTRAST

    out_dir     = Path(out_dir)
    per_ct_dir  = out_dir / "snc_vs_stats_by_celltype"
    per_mod_dir = out_dir / "snc_vs_stats_by_module"
    per_ct_dir.mkdir(parents=True, exist_ok=True)
    per_mod_dir.mkdir(parents=True, exist_ok=True)

    cell_types = regions_all_df["cell_type"].unique().tolist()
    modules    = [m for m in module_to_contrast
                  if m in regions_all_df["module"].unique()]

    # ── Precompute joint (DESeq2 stats + n_snc) for each (ct, module) pair ───
    joint_data = {}
    for ct in cell_types:
        ct_sub = regions_all_df[regions_all_df["cell_type"] == ct]
        for mod, contrast in module_to_contrast.items():
            fg_ids = ct_sub[ct_sub["module"] == mod]["peak_id"].tolist()
            if not fg_ids:
                continue
            branch_df = load_branch_df(contrast, branch_dir, ct)
            if branch_df is None:
                continue
            df = pd.DataFrame({"peak_id": fg_ids}).set_index("peak_id")
            df["n_snc"] = snc_counts_df.reindex(df.index, fill_value=0)["n_snc"]
            df["lfc"]    = branch_df.reindex(df.index)["log2FoldChange"]
            df["padj"]   = branch_df.reindex(df.index)["padj"]
            df["signed_stat"] = (
                np.sign(df["lfc"])
                * -np.log10(df["padj"].clip(lower=1e-300)))
            df = df.dropna(subset=["lfc"])
            if df.empty:
                continue
            joint_data[(ct, mod)] = df

    # ── Per-cell-type view ────────────────────────────────────────────────────
    for ct in cell_types:
        ct_modules = [m for m in modules if (ct, m) in joint_data]
        if not ct_modules:
            continue
        n_mods = len(ct_modules)
        fig, axes = plt.subplots(
            n_mods, 2,
            figsize=(10, max(3, n_mods * 2.5 + 0.5)),
            squeeze=False,
        )
        for row_i, mod in enumerate(ct_modules):
            df     = joint_data[(ct, mod)]
            ax_lfc = axes[row_i, 0]
            ax_sig = axes[row_i, 1]

            # Left: log2FC grouped by n_snc category
            labels_lfc, data_lfc = _boxplot_groups(df, "lfc")
            if data_lfc:
                ax_lfc.boxplot(data_lfc, labels=labels_lfc, showfliers=False,
                               medianprops={"color": "firebrick", "linewidth": 2})
                ax_lfc.axhline(0, color="gray", linewidth=0.8, linestyle="--")
            else:
                ax_lfc.text(0.5, 0.5, "insufficient data",
                            transform=ax_lfc.transAxes,
                            ha="center", va="center", fontsize=8, color="gray")
            ax_lfc.set_xlabel("SNCs per peak", fontsize=8)
            ax_lfc.set_ylabel("log2 Fold Change", fontsize=8)
            ax_lfc.set_title(mod, fontsize=8, fontweight="bold")
            ax_lfc.tick_params(labelsize=7)

            # Right: sign(LFC)·−log10(padj) grouped by n_snc category
            labels_sig, data_sig = _boxplot_groups(df, "signed_stat")
            if data_sig:
                ax_sig.boxplot(data_sig, labels=labels_sig, showfliers=False,
                               medianprops={"color": "steelblue", "linewidth": 2})
                ax_sig.axhline(0, color="gray", linewidth=0.8, linestyle="--")
            else:
                ax_sig.text(0.5, 0.5, "insufficient data",
                            transform=ax_sig.transAxes,
                            ha="center", va="center", fontsize=8, color="gray")
            ax_sig.set_xlabel("SNCs per peak", fontsize=8)
            ax_sig.set_ylabel("sign(LFC)·−log10(padj)", fontsize=8)
            ax_sig.set_title(mod, fontsize=8, fontweight="bold")
            ax_sig.tick_params(labelsize=7)

        fig.suptitle(f"{ct} — SNC count vs DESeq2 stats (foreground peaks)",
                     fontsize=10, fontweight="bold")
        fig.tight_layout(rect=[0, 0, 1, 0.97])
        out_path = per_ct_dir / f"{sanitize_name(ct)}_snc_vs_stats.pdf"
        fig.savefig(out_path, bbox_inches="tight", dpi=150)
        plt.close(fig)
        print(f"  Per-cell-type → {out_path}")

    # ── Per-module view ───────────────────────────────────────────────────────
    for mod in modules:
        ct_list = [ct for ct in cell_types if (ct, mod) in joint_data]
        if not ct_list:
            continue
        ncols = min(4, len(ct_list))
        nrows = (len(ct_list) + ncols - 1) // ncols
        fig, axes = plt.subplots(
            nrows, ncols,
            figsize=(ncols * 3.5, nrows * 3 + 0.5),
            squeeze=False,
        )
        for pi, ct in enumerate(ct_list):
            ax = axes[pi // ncols, pi % ncols]
            df = joint_data[(ct, mod)]
            labels_lfc, data_lfc = _boxplot_groups(df, "lfc")
            if data_lfc:
                ax.boxplot(data_lfc, labels=labels_lfc, showfliers=False,
                           medianprops={"color": "firebrick", "linewidth": 2})
                ax.axhline(0, color="gray", linewidth=0.8, linestyle="--")
            else:
                ax.text(0.5, 0.5, "insufficient data",
                        transform=ax.transAxes,
                        ha="center", va="center", fontsize=8, color="gray")
            ax.set_xlabel("SNCs per peak", fontsize=7)
            ax.set_ylabel("log2 Fold Change", fontsize=7)
            ax.set_title(ct, fontsize=8, fontweight="bold")
            ax.tick_params(labelsize=7)

        for pi in range(len(ct_list), nrows * ncols):
            axes[pi // ncols, pi % ncols].set_visible(False)

        fig.suptitle(f"{mod} — SNC count vs LFC (all cell types)",
                     fontsize=10, fontweight="bold")
        fig.tight_layout(rect=[0, 0, 1, 0.97])
        out_path = per_mod_dir / f"{sanitize_name(mod)}_snc_vs_stats.pdf"
        fig.savefig(out_path, bbox_inches="tight", dpi=150)
        plt.close(fig)
        print(f"  Per-module     → {out_path}")


def discover_cell_types(robust_dir, frag_dir, species_list):
    """
    Return a list of cell type display names that have both:
      - a robust results directory under robust_dir
      - at least one species with pseudobulk data for that cell type

    The mapping from dot-notation dir names back to display names is
    recovered by cross-referencing against the pseudobulk group files.
    """
    robust_dirs = {p.name for p in Path(robust_dir).iterdir() if p.is_dir()}

    # Collect all cell type display names from pseudobulk group files
    display_names = set()
    for sp in species_list:
        for ext in [".parquet", ".tsv"]:
            gpath = Path(frag_dir) / sp / f"pseudobulk_groups{ext}"
            if gpath.exists():
                g = (pd.read_parquet(gpath) if ext == ".parquet"
                     else pd.read_csv(gpath, sep="\t"))
                display_names.update(g["cell_type"].unique())
                break

    # Match display names to robust dir names
    available = []
    for name in sorted(display_names):
        if celltype_to_dirname(name) in robust_dirs:
            available.append(name)

    return available
