"""Utilities for building publication-oriented cross-species figure panels."""

from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DEFAULT_CELL_TYPES: Tuple[str, ...] = ("Macrophages", "T cells", "Enterocytes")
DEFAULT_CONTRASTS: Tuple[str, ...] = (
    "Pair_Human_vs_Chimp",
    "Pair_Human_vs_Gorilla",
    "Pair_Human_vs_Bonobo",
    "Div_Human_vs_Apes",
)

VOLCANO_COLORS = {
    "positive": "#c23b22",
    "negative": "#2563eb",
    "neutral": "#c7c7c7",
}


def celltype_to_dirname(cell_type: str) -> str:
    """Match the directory naming convention used by the DESeq2 outputs."""
    return re.sub(r"[^A-Za-z0-9]", ".", cell_type)


def sanitize_name(value: str) -> str:
    """Create a filesystem-safe name."""
    value = str(value).strip()
    value = re.sub(r"\s+", "_", value)
    return re.sub(r"[^A-Za-z0-9_.+-]", "_", value)


def parse_contrast_name(contrast: str) -> Dict[str, str]:
    """Return display metadata for a branch contrast name."""
    base = contrast
    for prefix in ("Pair_", "Div_", "Node", "ILS_"):
        if base.startswith(prefix):
            if prefix in {"Pair_", "Div_", "ILS_"}:
                base = base[len(prefix):]
            break

    if "_vs_" in base:
        lhs, rhs = base.split("_vs_", 1)
    else:
        lhs, rhs = base, "Others"

    lhs_display = lhs.replace("_", " ")
    rhs_display = rhs.replace("_", " ")
    return {
        "contrast": contrast,
        "lhs": lhs_display,
        "rhs": rhs_display,
        "display": f"{lhs_display} vs {rhs_display}",
        "positive_label": f"{lhs_display}-high",
        "negative_label": f"{rhs_display}-high",
    }


def load_peak_annotation(annotation_path: str) -> pd.DataFrame:
    """Load a peak annotation table and expose a single nearest-gene column."""
    annotation_df = pd.read_csv(annotation_path, sep="\t", low_memory=False)
    gene_candidates = [
        "Human_gene",
        "closest_gene",
        "gene",
    ]
    distance_candidates = [
        "Human_gene_dist",
        "distance_to_gene",
        "gene_distance",
    ]

    gene_col = next((column for column in gene_candidates if column in annotation_df.columns), None)
    distance_col = next(
        (column for column in distance_candidates if column in annotation_df.columns),
        None,
    )

    keep_columns = ["peak_id"]
    if gene_col is not None:
        keep_columns.append(gene_col)
    if distance_col is not None:
        keep_columns.append(distance_col)

    annotation_df = annotation_df[keep_columns].copy()
    rename_map = {}
    if gene_col is not None:
        rename_map[gene_col] = "nearest_gene"
    if distance_col is not None:
        rename_map[distance_col] = "gene_distance"
    annotation_df = annotation_df.rename(columns=rename_map)

    if "nearest_gene" not in annotation_df.columns:
        annotation_df["nearest_gene"] = pd.NA
    if "gene_distance" not in annotation_df.columns:
        annotation_df["gene_distance"] = pd.NA

    annotation_df["nearest_gene"] = (
        annotation_df["nearest_gene"]
        .replace({".": pd.NA, "": pd.NA})
        .astype("string")
    )
    return annotation_df.drop_duplicates(subset=["peak_id"])


def load_deseq_result(result_path: str) -> pd.DataFrame:
    """Load a DESeq2 result CSV with a normalized peak identifier column."""
    df = pd.read_csv(result_path)
    id_col = "peak_id"
    if "peak_id" not in df.columns:
        if "region_id" in df.columns:
            id_col = "region_id"
        else:
            id_col = df.columns[0]
    df = df.rename(columns={id_col: "peak_id"})
    return df


def prepare_volcano_dataframe(
    result_df: pd.DataFrame,
    annotation_df: pd.DataFrame,
    contrast_info: Dict[str, str],
    padj_thresh: float = 0.05,
    lfc_thresh: float = 1.0,
) -> pd.DataFrame:
    """Merge DESeq2 output with annotation and derive volcano plotting fields."""
    plot_df = result_df.merge(annotation_df, on="peak_id", how="left")
    plot_df = plot_df.dropna(subset=["padj", "log2FoldChange"]).copy()
    plot_df["padj_safe"] = plot_df["padj"].clip(lower=1e-300)
    plot_df["neg_log10_padj"] = -np.log10(plot_df["padj_safe"])
    plot_df["rank_score"] = plot_df["log2FoldChange"].abs() * plot_df["neg_log10_padj"]

    plot_df["direction"] = "Not significant"
    positive_mask = (plot_df["padj"] < padj_thresh) & (plot_df["log2FoldChange"] >= lfc_thresh)
    negative_mask = (plot_df["padj"] < padj_thresh) & (plot_df["log2FoldChange"] <= -lfc_thresh)
    plot_df.loc[positive_mask, "direction"] = contrast_info["positive_label"]
    plot_df.loc[negative_mask, "direction"] = contrast_info["negative_label"]

    plot_df["label_gene"] = plot_df["nearest_gene"].fillna(plot_df["peak_id"]).astype(str)
    plot_df.loc[plot_df["label_gene"].isin(["<NA>", "nan", "."]), "label_gene"] = plot_df["peak_id"]
    return plot_df


def select_top_labels(plot_df: pd.DataFrame, contrast_info: Dict[str, str], n_per_direction: int) -> pd.DataFrame:
    """Choose the highest-ranked annotated peaks on each side of the volcano."""
    label_frames: List[pd.DataFrame] = []
    for direction_key in (contrast_info["positive_label"], contrast_info["negative_label"]):
        direction_df = plot_df[plot_df["direction"] == direction_key].copy()
        if direction_df.empty:
            continue
        direction_df = direction_df.sort_values(
            ["rank_score", "baseMean"],
            ascending=[False, False],
        )
        direction_df = direction_df.drop_duplicates(subset=["label_gene"])
        label_frames.append(direction_df.head(n_per_direction))

    if not label_frames:
        return plot_df.iloc[0:0].copy()
    label_df = pd.concat(label_frames, ignore_index=False)
    return label_df.sort_values("neg_log10_padj", ascending=False)


def infer_axis_limits(plot_frames: Sequence[pd.DataFrame], lfc_thresh: float, padj_thresh: float) -> Tuple[float, float]:
    """Derive shared axis limits across a panel."""
    if not plot_frames:
        return max(2.0, lfc_thresh * 1.5), max(5.0, -math.log10(padj_thresh) * 1.5)

    merged = pd.concat(plot_frames, ignore_index=True)
    x_max = float(np.nanquantile(np.abs(merged["log2FoldChange"]), 0.995))
    y_max = float(np.nanquantile(merged["neg_log10_padj"], 0.995))
    x_max = max(x_max, lfc_thresh * 1.5, 2.0)
    y_max = max(y_max, -math.log10(padj_thresh) * 1.5, 5.0)
    return x_max * 1.05, y_max * 1.05


def _annotate_points(ax: plt.Axes, label_df: pd.DataFrame) -> None:
    """Annotate selected peaks with simple offset placement."""
    for idx, (_, row) in enumerate(label_df.iterrows()):
        x_offset = 10 if row["log2FoldChange"] >= 0 else -10
        y_offset = 6 + (idx % 5) * 5
        horizontal_alignment = "left" if row["log2FoldChange"] >= 0 else "right"
        ax.annotate(
            row["label_gene"],
            xy=(row["log2FoldChange"], row["neg_log10_padj"]),
            xytext=(x_offset, y_offset),
            textcoords="offset points",
            fontsize=7,
            fontstyle="italic",
            ha=horizontal_alignment,
            va="bottom",
            arrowprops={
                "arrowstyle": "-",
                "lw": 0.6,
                "color": "#666666",
                "alpha": 0.7,
            },
        )


def build_result_manifest(
    results_root: str,
    annotation_path: str,
    output_path: str,
    cell_types: Sequence[str],
    contrasts: Sequence[str],
    padj_thresh: float = 0.05,
    lfc_thresh: float = 1.0,
) -> pd.DataFrame:
    """Create a manifest describing available DA result sources and basic stats."""
    annotation_df = load_peak_annotation(annotation_path)
    records: List[Dict[str, object]] = []

    for cell_type in cell_types:
        cell_dir = Path(results_root) / celltype_to_dirname(cell_type)
        for contrast in contrasts:
            result_path = cell_dir / f"{contrast}.csv"
            contrast_info = parse_contrast_name(contrast)
            record: Dict[str, object] = {
                "panel": "differential_accessibility_volcano",
                "cell_type": cell_type,
                "contrast": contrast,
                "contrast_display": contrast_info["display"],
                "result_path": str(result_path),
                "exists": result_path.exists(),
            }
            if result_path.exists():
                result_df = load_deseq_result(str(result_path))
                plot_df = prepare_volcano_dataframe(
                    result_df=result_df,
                    annotation_df=annotation_df,
                    contrast_info=contrast_info,
                    padj_thresh=padj_thresh,
                    lfc_thresh=lfc_thresh,
                )
                record.update(
                    {
                        "n_rows": len(plot_df),
                        "n_positive": int((plot_df["direction"] == contrast_info["positive_label"]).sum()),
                        "n_negative": int((plot_df["direction"] == contrast_info["negative_label"]).sum()),
                        "min_padj": float(plot_df["padj"].min()) if not plot_df.empty else np.nan,
                        "max_abs_log2fc": float(plot_df["log2FoldChange"].abs().max()) if not plot_df.empty else np.nan,
                    }
                )
            records.append(record)

    manifest_df = pd.DataFrame.from_records(records)
    output_file = Path(output_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    manifest_df.to_csv(output_file, sep="\t", index=False)
    return manifest_df


def plot_volcano_panel(
    results_root: str,
    annotation_path: str,
    output_dir: str,
    contrast: str,
    cell_types: Sequence[str],
    padj_thresh: float = 0.05,
    lfc_thresh: float = 1.0,
    n_per_direction: int = 8,
) -> Dict[str, str]:
    """Build a multi-panel volcano figure for one contrast across selected cell types."""
    annotation_df = load_peak_annotation(annotation_path)
    contrast_info = parse_contrast_name(contrast)

    plot_frames: List[pd.DataFrame] = []
    label_tables: List[pd.DataFrame] = []
    panel_entries: List[Tuple[str, Optional[pd.DataFrame], Optional[pd.DataFrame]]] = []

    for cell_type in cell_types:
        result_path = Path(results_root) / celltype_to_dirname(cell_type) / f"{contrast}.csv"
        if not result_path.exists():
            panel_entries.append((cell_type, None, None))
            continue

        result_df = load_deseq_result(str(result_path))
        plot_df = prepare_volcano_dataframe(
            result_df=result_df,
            annotation_df=annotation_df,
            contrast_info=contrast_info,
            padj_thresh=padj_thresh,
            lfc_thresh=lfc_thresh,
        )
        label_df = select_top_labels(plot_df, contrast_info, n_per_direction)
        label_df = label_df.assign(cell_type=cell_type, contrast=contrast)
        plot_frames.append(plot_df)
        label_tables.append(label_df)
        panel_entries.append((cell_type, plot_df, label_df))

    x_lim, y_lim = infer_axis_limits(plot_frames, lfc_thresh, padj_thresh)
    figure_dir = Path(output_dir)
    figure_dir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, len(cell_types), figsize=(5.1 * len(cell_types), 5.4), sharex=True, sharey=True)
    if len(cell_types) == 1:
        axes = [axes]

    for ax, (cell_type, plot_df, label_df) in zip(axes, panel_entries):
        if plot_df is None:
            ax.axis("off")
            ax.text(0.5, 0.5, f"Missing result\n{cell_type}", ha="center", va="center", fontsize=11)
            continue

        neutral_df = plot_df[plot_df["direction"] == "Not significant"]
        positive_df = plot_df[plot_df["direction"] == contrast_info["positive_label"]]
        negative_df = plot_df[plot_df["direction"] == contrast_info["negative_label"]]

        ax.scatter(
            neutral_df["log2FoldChange"],
            neutral_df["neg_log10_padj"],
            s=8,
            c=VOLCANO_COLORS["neutral"],
            alpha=0.55,
            linewidths=0,
        )
        ax.scatter(
            positive_df["log2FoldChange"],
            positive_df["neg_log10_padj"],
            s=12,
            c=VOLCANO_COLORS["positive"],
            alpha=0.8,
            linewidths=0,
        )
        ax.scatter(
            negative_df["log2FoldChange"],
            negative_df["neg_log10_padj"],
            s=12,
            c=VOLCANO_COLORS["negative"],
            alpha=0.8,
            linewidths=0,
        )

        _annotate_points(ax, label_df)

        ax.axvline(lfc_thresh, linestyle="--", color="#555555", linewidth=0.9)
        ax.axvline(-lfc_thresh, linestyle="--", color="#555555", linewidth=0.9)
        ax.axhline(-math.log10(padj_thresh), linestyle="--", color="#555555", linewidth=0.9)
        ax.set_xlim(-x_lim, x_lim)
        ax.set_ylim(0, y_lim)
        ax.set_title(
            f"{cell_type}\n{len(positive_df)} {contrast_info['positive_label']} | {len(negative_df)} {contrast_info['negative_label']}",
            fontsize=11,
            fontweight="bold",
        )
        ax.set_xlabel("log2 fold change")
        ax.grid(alpha=0.15, linewidth=0.4)

    axes[0].set_ylabel("-log10 adjusted p-value")
    fig.suptitle(
        f"Differential accessibility: {contrast_info['display']}",
        fontsize=14,
        fontweight="bold",
        y=1.02,
    )
    fig.tight_layout()

    stem = sanitize_name(contrast)
    pdf_path = figure_dir / f"{stem}_volcano_panel.pdf"
    png_path = figure_dir / f"{stem}_volcano_panel.png"
    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    labels_path = figure_dir / f"{stem}_volcano_labels.tsv"
    if label_tables:
        labels_df = pd.concat(label_tables, ignore_index=True)
    else:
        labels_df = pd.DataFrame(
            columns=[
                "peak_id",
                "baseMean",
                "log2FoldChange",
                "padj",
                "nearest_gene",
                "gene_distance",
                "direction",
                "cell_type",
                "contrast",
            ]
        )
    labels_df.to_csv(labels_path, sep="\t", index=False)

    return {
        "pdf": str(pdf_path),
        "png": str(png_path),
        "labels": str(labels_path),
    }


def build_volcano_suite(
    results_root: str,
    annotation_path: str,
    output_dir: str,
    cell_types: Sequence[str] = DEFAULT_CELL_TYPES,
    contrasts: Sequence[str] = DEFAULT_CONTRASTS,
    padj_thresh: float = 0.05,
    lfc_thresh: float = 1.0,
    n_per_direction: int = 8,
) -> Dict[str, object]:
    """Create the manifest and all requested DA volcano panels."""
    output_root = Path(output_dir)
    output_root.mkdir(parents=True, exist_ok=True)

    manifest_path = output_root / "da_volcano_manifest.tsv"
    manifest_df = build_result_manifest(
        results_root=results_root,
        annotation_path=annotation_path,
        output_path=str(manifest_path),
        cell_types=cell_types,
        contrasts=contrasts,
        padj_thresh=padj_thresh,
        lfc_thresh=lfc_thresh,
    )

    outputs = []
    for contrast in contrasts:
        outputs.append(
            {
                "contrast": contrast,
                **plot_volcano_panel(
                    results_root=results_root,
                    annotation_path=annotation_path,
                    output_dir=str(output_root),
                    contrast=contrast,
                    cell_types=cell_types,
                    padj_thresh=padj_thresh,
                    lfc_thresh=lfc_thresh,
                    n_per_direction=n_per_direction,
                ),
            }
        )

    return {
        "manifest": str(manifest_path),
        "outputs": outputs,
        "manifest_rows": len(manifest_df),
    }