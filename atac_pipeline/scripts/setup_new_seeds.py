#!/usr/bin/env python3
"""
Setup SCENIC+ pipeline directories and downsampled h5ad files for new seeds.

Reads config from atac_pipeline/config/scenicplus_prep.yaml.
Skips seeds/species where the output h5ad and Snakemake directory already exist.
Appends new runs to the manifest and writes a run script for the new seeds only.

Usage:
    python setup_new_seeds.py [--config PATH]
"""
import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

# ── resolve src/ relative to this script ──────────────────────────────────────
SCRIPT_DIR = Path(__file__).resolve().parent
SRC_DIR = SCRIPT_DIR.parent / "src"
sys.path.insert(0, str(SRC_DIR))

from scenicplus_prep import (
    load_prep_config,
    compute_downsample_targets,
    build_run_name,
    build_sub_rna_name,
    generate_config_yaml,
    init_snakemake_dir,
    write_config,
    write_run_manifest,
)

# ── CLI ────────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument(
    "--config",
    default=str(SCRIPT_DIR.parent / "config" / "scenicplus_prep.yaml"),
    help="Path to scenicplus_prep.yaml",
)
args = parser.parse_args()

# ── Load config ────────────────────────────────────────────────────────────────
prep_cfg = load_prep_config(args.config)
PATHS = prep_cfg["paths"]
RUN_CFG = prep_cfg["run"]
SPECIES_CONFIG = prep_cfg["species_config"]

SPECIES_LIST = RUN_CFG["species_list"]
SEEDS = RUN_CFG["seeds"]
CELL_TYPE_COL = RUN_CFG["cell_type_col"]
DOWNSAMPLE_MODE = RUN_CFG["downsample_mode"]
CELL_TYPES_SUBSET = RUN_CFG.get("cell_type_subset")
N_CPU = RUN_CFG["n_cpu"]

INPUT_FILES_DIR = PATHS["input_files_dir"]
SCENICPLUS_DIR = PATHS["scenicplus_dir"]
METADATA_DIR = PATHS["metadata_dir"]
TEMP_DIR = PATHS["temp_dir"]
MANIFESTS_DIR = PATHS["manifests_dir"]
LOGS_DIR = PATHS["logs_dir"]
SCRIPTS_DIR = PATHS["scripts_dir"]

MANIFEST_PATH = os.path.join(MANIFESTS_DIR, "scenicplus_runs_dsMin.tsv")

def downsample_cells(meta_df, cell_type_col, target_per_ct, seed, available_barcodes=None, min_cells=50):
    rng = np.random.RandomState(seed)

    if available_barcodes is not None:
        meta_df = meta_df.loc[meta_df.index.intersection(list(available_barcodes))]

    if target_per_ct is None:
        if min_cells > 0:
            counts = meta_df[cell_type_col].value_counts()
            valid_cts = counts[counts >= min_cells].index
            return meta_df[meta_df[cell_type_col].isin(valid_cts)].index.tolist()
        return meta_df.index.tolist()

    selected = []
    for ct in target_per_ct.index:
        ct_cells = meta_df[meta_df[cell_type_col] == ct].index.tolist()
        if len(ct_cells) < min_cells:
            continue
        n_target = min(target_per_ct[ct], len(ct_cells))
        if n_target > 0:
            chosen = rng.choice(ct_cells, size=n_target, replace=False)
            selected.extend(chosen)
    return selected


# ── Load composition directly from existing h5ad obs ──────────────────────────
# The {species}_adata_rna_sub.h5ad files already have cell_type_initial in obs,
# so no need to load/reformat the separate metadata TSVs.
print("Loading cell type composition from existing h5ad files...")
composition = {}

for species in SPECIES_LIST:
    rna_path = os.path.join(INPUT_FILES_DIR, f"{species}_adata_rna_sub.h5ad")
    adata = sc.read_h5ad(rna_path, backed="r")
    if CELL_TYPE_COL not in adata.obs.columns:
        raise ValueError(f"{species}: column '{CELL_TYPE_COL}' not found in obs. Available: {list(adata.obs.columns)}")
    ct_counts = adata.obs[CELL_TYPE_COL].value_counts()
    composition[species] = ct_counts
    adata.file.close()
    print(f"  {species}: {sum(ct_counts)} cells, {len(ct_counts)} cell types")

comp_df = pd.DataFrame(composition).fillna(0)

if CELL_TYPES_SUBSET is not None:
    comp_df = comp_df.loc[comp_df.index.intersection(CELL_TYPES_SUBSET)]

target_per_ct, ds_label = compute_downsample_targets(
    comp_df=comp_df,
    mode=DOWNSAMPLE_MODE,
    downsample_n=RUN_CFG.get("downsample_n"),
)

# Global min-cells filter (mirror notebook cell 22)
MIN_CELLS_ACROSS_ALL = 50
valid_cell_types = comp_df[(comp_df >= MIN_CELLS_ACROSS_ALL).all(axis=1)].index
dropped = [ct for ct in comp_df.index if ct not in valid_cell_types]
if dropped:
    print(f"Dropping {len(dropped)} cell types (< {MIN_CELLS_ACROSS_ALL} cells in at least one species): {dropped}")
if target_per_ct is not None:
    target_per_ct = target_per_ct.loc[target_per_ct.index.isin(valid_cell_types)]

print(f"Proceeding with {len(valid_cell_types)} cell types, downsample mode: {ds_label}")

# ── Load existing manifest to detect already-done runs ────────────────────────
if os.path.exists(MANIFEST_PATH):
    existing_manifest = pd.read_csv(MANIFEST_PATH, sep="\t")
    existing_keys = set(zip(existing_manifest["species"].astype(str), existing_manifest["seed"].astype(int)))
else:
    existing_manifest = pd.DataFrame()
    existing_keys = set()

# ── Main loop: only new seeds ──────────────────────────────────────────────────
new_run_summary = []

for species in SPECIES_LIST:
    new_seeds = [s for s in SEEDS if (species, int(s)) not in existing_keys]
    if not new_seeds:
        print(f"\n{species}: all seeds already set up, skipping.")
        continue

    print(f"\n{'='*60}")
    print(f"  {species}  —  {len(new_seeds)} new seeds: {new_seeds}")
    print(f"{'='*60}")

    species_cfg = SPECIES_CONFIG[species]
    rna_path = os.path.join(INPUT_FILES_DIR, f"{species}_adata_rna_sub.h5ad")
    adata_rna = sc.read_h5ad(rna_path)
    print(f"  Loaded RNA: {adata_rna.shape[0]} cells × {adata_rna.shape[1]} genes")

    # cell_type_initial is already in obs
    adata_rna_with_ct = adata_rna

    cistopic_path = os.path.join(INPUT_FILES_DIR, f"{species}_cistopic_object_subset.pkl")

    for seed in new_seeds:
        seed = int(seed)
        print(f"\n  --- Seed {seed} ---")

        sub_rna_fname = build_sub_rna_name(species, seed, ds_label, CELL_TYPES_SUBSET)
        sub_rna_path = os.path.join(INPUT_FILES_DIR, sub_rna_fname)

        if os.path.exists(sub_rna_path):
            print(f"  h5ad already exists, loading: {sub_rna_path}")
            adata_sub = sc.read_h5ad(sub_rna_path)
        else:
            selected_cells = downsample_cells(
                meta_df=adata_rna_with_ct.obs,
                cell_type_col=CELL_TYPE_COL,
                target_per_ct=target_per_ct,
                seed=seed,
            )
            adata_sub = adata_rna_with_ct[selected_cells].copy()
            adata_sub.write_h5ad(sub_rna_path)
            print(f"  Saved RNA ({adata_sub.shape[0]} cells): {sub_rna_path}")

        n_cells = adata_sub.shape[0]
        ct_counts = adata_sub.obs[CELL_TYPE_COL].value_counts()

        run_name = build_run_name(species, seed, ds_label, CELL_TYPES_SUBSET)
        run_dir = os.path.join(SCENICPLUS_DIR, run_name)

        initialized = init_snakemake_dir(run_dir)
        if initialized:
            print(f"  Initialized Snakemake: {run_dir}")
        else:
            print(f"  Snakemake dir already exists: {run_dir}")

        snakemake_dir = os.path.join(run_dir, "Snakemake")
        config = generate_config_yaml(
            species=species,
            seed=seed,
            cistopic_path=cistopic_path,
            rna_path=sub_rna_path,
            species_cfg=species_cfg,
            temp_dir=TEMP_DIR,
            n_cpu=N_CPU,
        )
        config_path = os.path.join(snakemake_dir, "config", "config.yaml")
        write_config(config, config_path)
        print(f"  Wrote config: {config_path}")

        new_run_summary.append({
            "species": species,
            "seed": seed,
            "downsample": ds_label,
            "cell_types": len(ct_counts),
            "n_cells": n_cells,
            "run_dir": run_name,
            "run_dir_abs": run_dir,
            "snakemake_dir": snakemake_dir,
            "config_path": config_path,
            "rna_file": sub_rna_fname,
            "rna_path": sub_rna_path,
        })

    del adata_rna

# ── Update manifest ────────────────────────────────────────────────────────────
if new_run_summary:
    new_run_df = pd.DataFrame(new_run_summary)
    if not existing_manifest.empty:
        updated_manifest = pd.concat([existing_manifest, new_run_df], ignore_index=True)
    else:
        updated_manifest = new_run_df
    write_run_manifest(updated_manifest, MANIFEST_PATH)
    print(f"\nManifest updated: {MANIFEST_PATH} ({len(updated_manifest)} total runs)")

    # ── Generate run script for new seeds only ─────────────────────────────────
    run_script_path = os.path.join(SCRIPTS_DIR, "run_scenicplus_dsMin_new_seeds.sh")
    lines = [
        "#!/bin/bash",
        "set -euo pipefail",
        "",
        'BASE="/home/jjanssens/jjans/analysis/adult_intestine/scenicplus"',
        'LOG_DIR="${BASE}/logs"',
        "JOBS=40",
        "",
        "# Initialize Mamba for this script",
        'eval "$(mamba shell hook --shell bash)"',
        "mamba activate scenicplus_scenicplus10a2_py311",
        'export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"',
        "",
        "run_pipeline() {",
        '    local SPECIES="$1"',
        '    local SEED="$2"',
        '    local TAG="${SPECIES}_seed${SEED}_dsMin"',
        '    local WDIR="${BASE}/scplus_pipeline_${TAG}/Snakemake"',
        '    local SNAKEFILE="${WDIR}/workflow/Snakefile"',
        '    local CONFIG="${WDIR}/config/config.yaml"',
        '    local LOG="${LOG_DIR}/${SPECIES}_seed${SEED}.log"',
        "",
        '    echo "[RUN] ${SPECIES} seed=${SEED}"',
        "",
        '    if [ ! -d "$WDIR" ]; then',
        '        echo "  [SKIP] Directory not found: $WDIR"',
        "        return 0",
        "    fi",
        "",
        '    cd "$WDIR"',
        "",
        '    echo "  [UNLOCK] ${TAG}"',
        "    snakemake \\",
        '        --snakefile "$SNAKEFILE" \\',
        '        --configfile "$CONFIG" \\',
        "        --unlock \\",
        "        --use-conda \\",
        '        2>&1 | tee -a "$LOG" || true',
        "",
        '    echo "  [SNAKEMAKE] ${TAG}"',
        "    snakemake \\",
        '        --snakefile "$SNAKEFILE" \\',
        '        --configfile "$CONFIG" \\',
        '        -j "$JOBS" \\',
        "        --use-conda \\",
        "        --rerun-incomplete \\",
        '        2>&1 | tee -a "$LOG"',
        "",
        '    echo "  [DONE] ${TAG}"',
        "}",
        "",
        'mkdir -p "$LOG_DIR"',
        "",
    ]

    for _, row in new_run_df.iterrows():
        lines.append(f'run_pipeline {row["species"]:12s} {int(row["seed"])}')

    lines += ["", 'echo "[ALL DONE]"', ""]

    with open(run_script_path, "w") as f:
        f.write("\n".join(lines))
    os.chmod(run_script_path, 0o755)
    print(f"Run script written: {run_script_path}")
else:
    print("\nNo new runs to set up.")
