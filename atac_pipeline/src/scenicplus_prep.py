import os
import glob
import stat
import yaml
import subprocess
from copy import deepcopy
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

import pandas as pd


def _format_string(value: str, ctx: Dict[str, Any]) -> str:
    formatted = value
    for _ in range(3):
        new_value = formatted.format(**ctx)
        if new_value == formatted:
            return new_value
        formatted = new_value
    return formatted


def _resolve_templates(obj: Any, ctx: Dict[str, Any]) -> Any:
    if isinstance(obj, str):
        return _format_string(obj, ctx)
    if isinstance(obj, dict):
        return {k: _resolve_templates(v, ctx) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_resolve_templates(v, ctx) for v in obj]
    return obj


def load_prep_config(config_path: str) -> Dict[str, Any]:
    with open(config_path, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    if "paths" not in cfg or "run" not in cfg or "species_config" not in cfg:
        raise ValueError("Config must contain paths, run, and species_config sections")

    cfg = deepcopy(cfg)

    path_ctx = dict(cfg["paths"])
    resolved_paths = {}
    for key, val in path_ctx.items():
        if isinstance(val, str):
            resolved_paths[key] = _format_string(val, {**path_ctx, **resolved_paths})
        else:
            resolved_paths[key] = val

    cfg["paths"] = resolved_paths
    cfg["species_config"] = _resolve_templates(
        cfg["species_config"],
        {**resolved_paths, **{"base_dir": resolved_paths["base_dir"], "input_files_dir": resolved_paths["input_files_dir"]}},
    )

    _validate_prep_config(cfg)
    return cfg


def _validate_prep_config(cfg: Dict[str, Any]) -> None:
    required_run = ["species_list", "seeds", "n_cpu", "cell_type_col", "downsample_mode"]
    missing = [k for k in required_run if k not in cfg["run"]]
    if missing:
        raise ValueError(f"Missing run config keys: {missing}")

    species_required = [
        "biomart_species",
        "ctx_db",
        "dem_db",
        "motif_annotations",
        "region_set_folder",
        "genome_annotation",
        "chromsizes",
    ]
    for species in cfg["run"]["species_list"]:
        if species not in cfg["species_config"]:
            raise ValueError(f"Species {species} missing from species_config")
        sp_cfg = cfg["species_config"][species]
        missing_sp = [k for k in species_required if k not in sp_cfg]
        if missing_sp:
            raise ValueError(f"Species {species} missing keys: {missing_sp}")


def compute_downsample_targets(comp_df: pd.DataFrame, mode: Any, downsample_n: Optional[int] = None):
    if mode == "min":
        target_per_ct = comp_df.apply(lambda row: row[row > 0].min(), axis=1).astype(int)
        return target_per_ct, "dsMin"
    if isinstance(mode, int) or (isinstance(mode, str) and mode.isdigit()):
        n = int(mode)
        return pd.Series(n, index=comp_df.index), f"ds{n}"
    if downsample_n is not None:
        n = int(downsample_n)
        return pd.Series(n, index=comp_df.index), f"ds{n}"
    return None, "full"


def build_run_name(species: str, seed: int, ds_label: str, cell_types_subset: Optional[Iterable[str]]) -> str:
    ct_suffix = ""
    if cell_types_subset is not None:
        ct_suffix = f"_ct{len(list(cell_types_subset))}"
    return f"scplus_pipeline_{species}_seed{seed}_{ds_label}{ct_suffix}"


def build_sub_rna_name(species: str, seed: int, ds_label: str, cell_types_subset: Optional[Iterable[str]]) -> str:
    ct_suffix = ""
    if cell_types_subset is not None:
        ct_suffix = f"_ct{len(list(cell_types_subset))}"
    return f"{species}_adata_rna_sub_seed{seed}_{ds_label}{ct_suffix}.h5ad"


def generate_config_yaml(
    species: str,
    seed: int,
    cistopic_path: str,
    rna_path: str,
    species_cfg: Dict[str, Any],
    temp_dir: str,
    n_cpu: int,
) -> Dict[str, Any]:
    return {
        "input_data": {
            "cisTopic_obj_fname": cistopic_path,
            "GEX_anndata_fname": rna_path,
            "region_set_folder": species_cfg["region_set_folder"],
            "ctx_db_fname": species_cfg["ctx_db"],
            "dem_db_fname": species_cfg["dem_db"],
            "path_to_motif_annotations": species_cfg["motif_annotations"],
        },
        "output_data": {
            "combined_GEX_ACC_mudata": "ACC_GEX.h5mu",
            "dem_result_fname": "dem_results.hdf5",
            "ctx_result_fname": "ctx_results.hdf5",
            "output_fname_dem_html": "dem_results.html",
            "output_fname_ctx_html": "ctx_results.html",
            "cistromes_direct": "cistromes_direct.h5ad",
            "cistromes_extended": "cistromes_extended.h5ad",
            "tf_names": "tf_names.txt",
            "genome_annotation": species_cfg["genome_annotation"],
            "chromsizes": species_cfg["chromsizes"],
            "search_space": "search_space.tsv",
            "tf_to_gene_adjacencies": "tf_to_gene_adj.tsv",
            "region_to_gene_adjacencies": "region_to_gene_adj.tsv",
            "eRegulons_direct": "eRegulon_direct.tsv",
            "eRegulons_extended": "eRegulons_extended.tsv",
            "AUCell_direct": "AUCell_direct.h5mu",
            "AUCell_extended": "AUCell_extended.h5mu",
            "scplus_mdata": "scplusmdata.h5mu",
        },
        "params_general": {
            "temp_dir": temp_dir,
            "n_cpu": int(n_cpu),
            "seed": int(seed),
        },
        "params_data_preparation": {
            "bc_transform_func": '"lambda x: f\'{x}\'"',
            "is_multiome": True,
            "key_to_group_by": "",
            "nr_cells_per_metacells": 10,
            "direct_annotation": "Direct_annot",
            "extended_annotation": "Orthology_annot",
            "species": species_cfg["biomart_species"],
            "biomart_host": "http://www.ensembl.org",
            "search_space_upstream": "1000 150000",
            "search_space_downstream": "1000 150000",
            "search_space_extend_tss": "10 10",
        },
        "params_motif_enrichment": {
            "species": "homo_sapiens",
            "annotation_version": "v10nr_clust",
            "motif_similarity_fdr": 0.001,
            "orthologous_identity_threshold": 0.0,
            "annotations_to_use": "Direct_annot Orthology_annot",
            "fraction_overlap_w_dem_database": 0.4,
            "dem_max_bg_regions": 500,
            "dem_balance_number_of_promoters": True,
            "dem_promoter_space": 1000,
            "dem_adj_pval_thr": 0.05,
            "dem_log2fc_thr": 1.0,
            "dem_mean_fg_thr": 0.0,
            "dem_motif_hit_thr": 3.0,
            "fraction_overlap_w_ctx_database": 0.4,
            "ctx_auc_threshold": 0.005,
            "ctx_nes_threshold": 3.0,
            "ctx_rank_threshold": 0.05,
        },
        "params_inference": {
            "tf_to_gene_importance_method": "GBM",
            "region_to_gene_importance_method": "GBM",
            "region_to_gene_correlation_method": "SR",
            "order_regions_to_genes_by": "importance",
            "order_TFs_to_genes_by": "importance",
            "gsea_n_perm": 1000,
            "quantile_thresholds_region_to_gene": "0.85 0.90 0.95",
            "top_n_regionTogenes_per_gene": "5 10 15",
            "top_n_regionTogenes_per_region": "",
            "min_regions_per_gene": 0,
            "rho_threshold": 0.05,
            "min_target_genes": 10,
        },
    }


def init_snakemake_dir(run_dir: str) -> bool:
    snakemake_dir = os.path.join(run_dir, "Snakemake")
    if os.path.exists(snakemake_dir):
        return False
    subprocess.run(["scenicplus", "init_snakemake", "--out_dir", run_dir], check=True, capture_output=True, text=True)
    return True


def write_config(config: Dict[str, Any], config_path: str) -> None:
    os.makedirs(os.path.dirname(config_path), exist_ok=True)
    with open(config_path, "w", encoding="utf-8") as f:
        yaml.safe_dump(config, f, default_flow_style=False, sort_keys=False)


def write_run_manifest(run_df: pd.DataFrame, manifest_path: str) -> None:
    os.makedirs(os.path.dirname(manifest_path), exist_ok=True)
    run_df.to_csv(manifest_path, sep="\t", index=False)


def generate_snakemake_command(run_snakemake_dir: str, n_cpu: int, extra_args: str = "--use-conda") -> str:
    snakefile = os.path.join(run_snakemake_dir, "workflow", "Snakefile")
    configfile = os.path.join(run_snakemake_dir, "config", "config.yaml")
    parts = [
        f"cd {run_snakemake_dir}",
        f"snakemake --snakefile {snakefile} --configfile {configfile} -j {int(n_cpu)}",
    ]
    if extra_args:
        parts[-1] = f"{parts[-1]} {extra_args.strip()}"
    return " && ".join(parts)


def generate_batch_script(
    run_df: pd.DataFrame,
    script_path: str,
    n_cpu: int,
    extra_snakemake_args: str = "--use-conda",
    log_dir: Optional[str] = None,
    shared_post_commands: Optional[Dict[str, List[str]]] = None,
) -> str:
    os.makedirs(os.path.dirname(script_path), exist_ok=True)
    if log_dir:
        os.makedirs(log_dir, exist_ok=True)

    lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "",
        "echo \"Starting SCENIC+ batch run\"",
        "",
    ]

    for _, row in run_df.iterrows():
        species = row["species"]
        seed = int(row["seed"])
        run_snakemake_dir = os.path.join(row["run_dir_abs"], "Snakemake")
        cmd = generate_snakemake_command(run_snakemake_dir, n_cpu=n_cpu, extra_args=extra_snakemake_args)
        lines.append(f"echo \"[RUN] {species} seed={seed}\"")
        if log_dir:
            log_path = os.path.join(log_dir, f"{species}_seed{seed}.log")
            lines.append(f"{cmd} |& tee {log_path}")
        else:
            lines.append(cmd)
        lines.append("")

    if shared_post_commands:
        lines.append("echo \"Starting shared post-processing steps\"")
        lines.append("")
        for species, commands in shared_post_commands.items():
            lines.append(f"echo \"[SHARED] {species}\"")
            for cmd in commands:
                lines.append(cmd)
            lines.append("")

    with open(script_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")

    current_mode = os.stat(script_path).st_mode
    os.chmod(script_path, current_mode | stat.S_IXUSR | stat.S_IXGRP)
    return script_path


def _infer_key_columns(df: pd.DataFrame, kind: str) -> List[str]:
    lower_to_col = {c.lower(): c for c in df.columns}

    if kind == "tf_to_gene":
        tf_candidates = ["tf", "tf_name", "transcription_factor"]
        gene_candidates = ["target", "gene", "target_gene", "target_genes", "gene_name"]
        tf_col = next((lower_to_col[c] for c in tf_candidates if c in lower_to_col), None)
        gene_col = next((lower_to_col[c] for c in gene_candidates if c in lower_to_col), None)
        if tf_col and gene_col:
            return [tf_col, gene_col]

    if kind == "region_to_gene":
        region_candidates = ["region", "region_name", "region_id"]
        gene_candidates = ["target", "gene", "target_gene", "target_genes", "gene_name"]
        region_col = next((lower_to_col[c] for c in region_candidates if c in lower_to_col), None)
        gene_col = next((lower_to_col[c] for c in gene_candidates if c in lower_to_col), None)
        if region_col and gene_col:
            return [region_col, gene_col]

    fallback = [c for c in df.columns if df[c].dtype == object][:2]
    if len(fallback) < 2:
        raise ValueError(f"Unable to infer key columns for {kind}. Columns: {list(df.columns)}")
    return fallback


def _aggregate_edge_tables(paths: List[str], kind: str, min_seed_support: int) -> pd.DataFrame:
    tables = []
    for path in paths:
        seed_token = Path(path).parts[-3]
        seed_str = seed_token.split("_seed")[-1].split("_")[0]
        df = pd.read_csv(path, sep="\t")
        df["seed"] = int(seed_str)
        tables.append(df)

    if not tables:
        return pd.DataFrame()

    full = pd.concat(tables, ignore_index=True)
    key_cols = _infer_key_columns(full, kind=kind)

    agg = full.groupby(key_cols, dropna=False).agg(seed_support=("seed", "nunique"))
    agg = agg.reset_index()
    agg["seed_fraction"] = agg["seed_support"] / full["seed"].nunique()

    numeric_cols = [c for c in full.columns if c not in key_cols + ["seed"] and pd.api.types.is_numeric_dtype(full[c])]
    if numeric_cols:
        stats = full.groupby(key_cols, dropna=False)[numeric_cols].mean().reset_index()
        stats = stats.rename(columns={c: f"mean_{c}" for c in numeric_cols})
        agg = agg.merge(stats, on=key_cols, how="left")

    agg = agg[agg["seed_support"] >= int(min_seed_support)].copy()
    return agg.sort_values(["seed_support"], ascending=False)


def aggregate_links_for_species(
    scenicplus_dir: str,
    species: str,
    output_dir: str,
    min_seed_support: int = 2,
) -> Dict[str, str]:
    tf_glob = os.path.join(scenicplus_dir, f"scplus_pipeline_{species}_seed*", "Snakemake", "tf_to_gene_adj.tsv")
    rg_glob = os.path.join(scenicplus_dir, f"scplus_pipeline_{species}_seed*", "Snakemake", "region_to_gene_adj.tsv")

    tf_paths = sorted(glob.glob(tf_glob))
    rg_paths = sorted(glob.glob(rg_glob))

    os.makedirs(output_dir, exist_ok=True)
    out = {}

    tf_agg = _aggregate_edge_tables(tf_paths, kind="tf_to_gene", min_seed_support=min_seed_support)
    tf_out = os.path.join(output_dir, f"{species}_tf_to_gene_consensus.tsv")
    tf_agg.to_csv(tf_out, sep="\t", index=False)
    out["tf_to_gene"] = tf_out

    rg_agg = _aggregate_edge_tables(rg_paths, kind="region_to_gene", min_seed_support=min_seed_support)
    rg_out = os.path.join(output_dir, f"{species}_region_to_gene_consensus.tsv")
    rg_agg.to_csv(rg_out, sep="\t", index=False)
    out["region_to_gene"] = rg_out

    return out
