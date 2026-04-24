#!/usr/bin/env python3
"""
Create RNA pseudobulk count matrices from merged h5ad for cross-species DESeq2.

Pseudobulk unit: Individual × ct (cell type) × Region
Aggregation:     raw count sum
Adult cells only (Age == 'Adult')

Saves per-species parquet files in a format compatible with load_pseudobulk_data()-
equivalent R loading code, plus combined RDS-ready CSVs.
"""
import sys
import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path
from scipy.sparse import issparse

H5AD = Path("/links/groups/treutlein/USERS/jjans/analysis/adult_intestine/rna/integration"
            "/scvi_integration/nhp_atlas_merged_scvi_input_20250520_112026.h5ad")
OUT_DIR = Path("/links/groups/treutlein/USERS/jjans/analysis/adult_intestine/rna"
               "/pseudobulk_deseq2")
OUT_DIR.mkdir(parents=True, exist_ok=True)

SPECIES   = ["Human", "Bonobo", "Chimpanzee", "Gorilla", "Macaque", "Marmoset"]
GROUP_COLS = ["Individual", "ct", "Region"]   # pseudobulk grouping

# QC thresholds (same philosophy as ATAC pipeline)
MIN_COUNTS = 500    # min total counts per pseudobulk sample
MIN_CELLS  = 10     # min cells contributing to a pseudobulk sample

print("Loading h5ad …", flush=True)
adata = sc.read_h5ad(H5AD)
print(f"  Full atlas: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")
print(f"  X dtype: {adata.X.dtype}, sparse: {issparse(adata.X)}")

# Verify raw counts: spot-check a few cells
sample_sums = np.asarray(adata.X[:5].sum(axis=1)).flatten()
print(f"  Row sums (first 5 cells): {sample_sums}  ← should look like raw counts (100s–10,000s)")

# Filter adults
adata = adata[adata.obs["Age"] == "Adult"].copy()
print(f"  Adult cells: {adata.shape[0]:,}")
print(f"  Species: {adata.obs['Species'].value_counts().to_dict()}")
print(f"  Cell types (ct): {sorted(adata.obs['ct'].unique())}")


def make_pseudobulk(adata_sub, group_cols, min_counts, min_cells, species_label):
    """Sum raw counts per group, return (counts_df, meta_df)."""
    obs = adata_sub.obs[group_cols + ["Species"]].copy()
    obs["_gkey"] = obs[group_cols].astype(str).agg("__".join, axis=1)

    unique_keys = obs["_gkey"].unique()
    idx_map = {k: np.where(obs["_gkey"].values == k)[0] for k in unique_keys}

    X = adata_sub.X
    rows, meta_rows = [], []

    for gkey, idx in idx_map.items():
        if issparse(X):
            vec = np.asarray(X[idx].sum(axis=0)).flatten()
        else:
            vec = X[idx].sum(axis=0)
        n_counts = vec.sum()
        n_cells  = len(idx)
        if n_counts < min_counts or n_cells < min_cells:
            continue
        rows.append(vec)
        vals = obs.iloc[idx[0]][group_cols].to_dict()
        vals.update({
            "sample_id":  gkey,
            "n_cells":    n_cells,
            "n_counts":   n_counts,
            "species":    species_label,
        })
        meta_rows.append(vals)

    if not rows:
        return None, None

    counts_df = pd.DataFrame(
        np.vstack(rows),
        index=[r["sample_id"] for r in meta_rows],
        columns=adata_sub.var_names,
    ).T  # genes × samples (matches ATAC convention)

    meta_df = pd.DataFrame(meta_rows).set_index("sample_id")
    # Rename to match ATAC pipeline column conventions
    meta_df = meta_df.rename(columns={"Individual": "donor", "ct": "cell_type",
                                       "Region": "region"})
    return counts_df, meta_df


all_meta = []

for species in SPECIES:
    print(f"\n{'='*50}\n{species}", flush=True)
    sp_adata = adata[adata.obs["Species"] == species].copy()
    print(f"  {sp_adata.shape[0]:,} adult cells")

    counts_df, meta_df = make_pseudobulk(
        sp_adata, GROUP_COLS, MIN_COUNTS, MIN_CELLS, species
    )
    if counts_df is None:
        print(f"  WARNING: no pseudobulks passed QC for {species}")
        continue

    sp_dir = OUT_DIR / species
    sp_dir.mkdir(exist_ok=True)

    # Save counts (genes × samples) as parquet
    counts_df.to_parquet(sp_dir / "pseudobulk_counts.parquet")
    # Save metadata
    meta_df.to_parquet(sp_dir / "pseudobulk_meta.parquet")

    all_meta.append(meta_df.assign(species=species))

    ct_counts = meta_df["cell_type"].value_counts()
    print(f"  {counts_df.shape[0]:,} genes × {counts_df.shape[1]:,} samples (QC passed)")
    print(f"  Cell types: {dict(ct_counts)}")
    print(f"  Saved to: {sp_dir}")

# Combined metadata summary
combined_meta = pd.concat(all_meta)
combined_meta.to_csv(OUT_DIR / "pseudobulk_meta_all_species.csv")
print(f"\nTotal pseudobulk samples across all species: {len(combined_meta)}")
print(combined_meta.groupby(["species", "cell_type"]).size().unstack(fill_value=0).to_string())
print(f"\nAll files saved to: {OUT_DIR}")
