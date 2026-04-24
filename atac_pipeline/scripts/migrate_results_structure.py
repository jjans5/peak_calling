#!/usr/bin/env python3
"""
Migrate DESeq2 results from old analysis-type-first layout to new cell-type-first layout.

Old ATAC layout (under ATAC_BASE):
  differential_results/shared_peaks/{CT}/{Species}_vs_All_Others.csv
  differential_results/evolutionary_branches/{CT}/{Contrast}.csv
  differential_results/ultra_robust_branches/{CT}/{Contrast}_UltraRobust.csv
  differential_results/region_specific/{Region}/{CT}/{Contrast}.csv
  differential_results/enterocytes/...  (deep-dive, nested structure)
  differential_results/cross_celltype_*/
  plots/{CT}_Visualizations/

New layout (under ATAC_BASE):
  {CT}/atac_shared/{Species}_vs_rest_DA.csv
  {CT}/atac_evolutionary/{Contrast}_DA.csv
  {CT}/atac_ultra_robust/{Contrast}_DA.csv
  {CT}/atac_region/{Region}_{Contrast}_DA.csv
  {CT}/enrichment/
  {CT}/motif/
  {CT}/plots/
  _summary/  (RDS objects, cross-CT plots, outlier logs)

Old RNA layout (under RNA_BASE):
  rna_differential_results/{CT}/{Contrast}_RNA_DE.csv
  rna_differential_results/species_specific/{CT}/{Species}_vs_rest_RNA_DE.csv

New RNA layout (under RNA_BASE):
  {CT}/rna_evolutionary/{Contrast}_DE.csv
  {CT}/rna_species/{Species}_vs_rest_DE.csv
  _summary/

Everything that doesn't map to the new structure is moved to BACKUP/.
Run with --dry-run first to preview, then without to execute.
"""

import argparse
import os
import shutil
from pathlib import Path
from datetime import datetime


ATAC_BASE = Path("/links/groups/treutlein/USERS/jjans/analysis/adult_intestine"
                 "/peaks/cross_species_consensus_v3/13_deseq2_R_pseudobulk")
RNA_BASE  = Path("/links/groups/treutlein/USERS/jjans/analysis/adult_intestine"
                 "/rna/pseudobulk_deseq2")
BACKUP    = ATAC_BASE / f"BACKUP_OLD_STRUCTURE_{datetime.now().strftime('%Y%m%d')}"


def move(src: Path, dst: Path, dry_run: bool):
    if not src.exists():
        return
    if dry_run:
        print(f"  MOVE  {src.relative_to(ATAC_BASE.parent.parent.parent.parent.parent)}\n"
              f"     →  {dst.relative_to(ATAC_BASE.parent.parent.parent.parent.parent)}")
        return
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists():
        # Don't overwrite — send to backup instead
        backup_dst = BACKUP / src.relative_to(ATAC_BASE)
        backup_dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(src), str(backup_dst))
        print(f"  [CONFLICT] backed up {src.name}")
        return
    shutil.move(str(src), str(dst))


def backup(src: Path, dry_run: bool):
    if not src.exists():
        return
    rel = src.relative_to(ATAC_BASE)
    dst = BACKUP / rel
    if dry_run:
        print(f"  BACKUP {rel}")
        return
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(str(src), str(dst))


def migrate_atac(dry_run: bool):
    print("\n=== ATAC: shared_peaks ===")
    sp_dir = ATAC_BASE / "differential_results" / "shared_peaks"
    if sp_dir.exists():
        for ct_dir in sorted(sp_dir.iterdir()):
            if not ct_dir.is_dir():
                continue
            ct = ct_dir.name
            for f in ct_dir.glob("*.csv"):
                # {Species}_vs_All_Others.csv → {CT}/atac_shared/{Species}_vs_rest_DA.csv
                new_name = f.name.replace("_vs_All_Others.csv", "_vs_rest_DA.csv")
                move(f, ATAC_BASE / ct / "atac_shared" / new_name, dry_run)
        for rds in sp_dir.glob("*.rds"):
            move(rds, ATAC_BASE / "_summary" / rds.name, dry_run)
        for csv in sp_dir.glob("outlier*.csv"):
            move(csv, ATAC_BASE / "_summary" / csv.name, dry_run)

    print("\n=== ATAC: evolutionary_branches ===")
    evo_dir = ATAC_BASE / "differential_results" / "evolutionary_branches"
    if evo_dir.exists():
        for ct_dir in sorted(evo_dir.iterdir()):
            if not ct_dir.is_dir():
                continue
            ct = ct_dir.name
            for f in ct_dir.glob("*.csv"):
                # {Contrast}.csv → {CT}/atac_evolutionary/{Contrast}_DA.csv
                new_name = f.stem + "_DA.csv"
                move(f, ATAC_BASE / ct / "atac_evolutionary" / new_name, dry_run)
        for rds in evo_dir.glob("*.rds"):
            move(rds, ATAC_BASE / "_summary" / rds.name, dry_run)

    print("\n=== ATAC: ultra_robust_branches ===")
    ur_dir = ATAC_BASE / "differential_results" / "ultra_robust_branches"
    if ur_dir.exists():
        for ct_dir in sorted(ur_dir.iterdir()):
            if not ct_dir.is_dir():
                continue
            ct = ct_dir.name
            for f in ct_dir.glob("*.csv"):
                # {Contrast}_UltraRobust.csv → {CT}/atac_ultra_robust/{Contrast}_DA.csv
                new_name = f.name.replace("_UltraRobust.csv", "_DA.csv")
                move(f, ATAC_BASE / ct / "atac_ultra_robust" / new_name, dry_run)
        for rds in ur_dir.glob("*.rds"):
            move(rds, ATAC_BASE / "_summary" / rds.name, dry_run)

    print("\n=== ATAC: region_specific ===")
    reg_dir = ATAC_BASE / "differential_results" / "region_specific"
    if reg_dir.exists():
        for region_dir in sorted(reg_dir.iterdir()):
            if not region_dir.is_dir():
                continue
            region = region_dir.name
            for ct_dir in sorted(region_dir.iterdir()):
                if not ct_dir.is_dir():
                    continue
                ct = ct_dir.name
                for f in ct_dir.glob("*.csv"):
                    # {Contrast}.csv → {CT}/atac_region/{Region}_{Contrast}_DA.csv
                    new_name = f"{region}_{f.stem}_DA.csv"
                    move(f, ATAC_BASE / ct / "atac_region" / new_name, dry_run)
        for rds in reg_dir.glob("*.rds"):
            move(rds, ATAC_BASE / "_summary" / rds.name, dry_run)
        for csv in reg_dir.glob("outlier*.csv"):
            move(csv, ATAC_BASE / "_summary" / csv.name, dry_run)

    print("\n=== ATAC: enterocytes deep-dive ===")
    ent_dir = ATAC_BASE / "differential_results" / "enterocytes"
    if ent_dir.exists():
        ct = "Enterocytes"
        # Nested evolutionary_branches
        for f in (ent_dir / "differential_results" / "evolutionary_branches" / ct).glob("*.csv"):
            new_name = f.stem + "_DA.csv"
            move(f, ATAC_BASE / ct / "atac_evolutionary" / new_name, dry_run)
        # Nested ultra_robust_branches
        for f in (ent_dir / "differential_results" / "ultra_robust_branches" / ct).glob("*.csv"):
            new_name = f.name.replace("_UltraRobust.csv", "_DA.csv")
            move(f, ATAC_BASE / ct / "atac_ultra_robust" / new_name, dry_run)
        # Nested shared_peaks
        for f in (ent_dir / "differential_results" / "shared_peaks" / ct).glob("*.csv"):
            new_name = f.name.replace("_vs_All_Others.csv", "_vs_rest_DA.csv")
            move(f, ATAC_BASE / ct / "atac_shared" / new_name, dry_run)
        # Nested region_specific
        for region_dir in (ent_dir / "differential_results" / "region_specific").glob("*"):
            if not region_dir.is_dir():
                continue
            region = region_dir.name
            for f in (region_dir / ct).glob("*.csv"):
                new_name = f"{region}_{f.stem}_DA.csv"
                move(f, ATAC_BASE / ct / "atac_region" / new_name, dry_run)
        # Nested enrichment
        enr_src = ent_dir / "differential_results" / "enrichment"
        if enr_src.exists():
            for f in enr_src.iterdir():
                move(f, ATAC_BASE / ct / "enrichment" / f.name, dry_run)
        # motif_enrichment → motif
        motif_src = ent_dir / "motif_enrichment"
        if motif_src.exists():
            for f in motif_src.iterdir():
                move(f, ATAC_BASE / ct / "motif" / f.name, dry_run)
        # gene_annotation
        ga_src = ent_dir / "gene_annotation"
        if ga_src.exists():
            for f in ga_src.iterdir():
                move(f, ATAC_BASE / ct / "gene_annotation" / f.name, dry_run)
        # RDS files in nested differential_results
        for rds in (ent_dir / "differential_results").glob("*.rds"):
            move(rds, ATAC_BASE / "_summary" / f"Enterocytes_{rds.name}", dry_run)
        # backup the now-empty enterocytes dir skeleton
        backup(ent_dir, dry_run)

    print("\n=== ATAC: summary-level files ===")
    diff_dir = ATAC_BASE / "differential_results"
    if diff_dir.exists():
        for rds in diff_dir.glob("*.rds"):
            move(rds, ATAC_BASE / "_summary" / rds.name, dry_run)
        for csv in diff_dir.glob("*.csv"):
            move(csv, ATAC_BASE / "_summary" / csv.name, dry_run)
        # cross_celltype_modules + cross_celltype_clusters
        for d in ["cross_celltype_modules", "cross_celltype_clusters"]:
            src_d = diff_dir / d
            if src_d.exists():
                dst_d = ATAC_BASE / "_summary" / d
                if dry_run:
                    print(f"  MOVEDIR {d} → _summary/{d}")
                else:
                    dst_d.parent.mkdir(parents=True, exist_ok=True)
                    shutil.move(str(src_d), str(dst_d))

    print("\n=== ATAC: backup old cell-type dirs in differential_results ===")
    # Any remaining dirs under differential_results go to backup
    if diff_dir.exists():
        for d in sorted(diff_dir.iterdir()):
            if d.is_dir():
                backup(d, dry_run)
        # now backup the differential_results dir itself if empty
        if not dry_run:
            try:
                diff_dir.rmdir()
                print("  Removed now-empty differential_results/")
            except OSError:
                backup(diff_dir, dry_run)

    print("\n=== ATAC: plots ===")
    plots_dir = ATAC_BASE / "plots"
    if plots_dir.exists():
        for d in plots_dir.iterdir():
            # {CT}_Visualizations/ → {CT}/plots/
            if d.is_dir() and "_Visualizations" in d.name:
                ct = d.name.replace("_Visualizations", "")
                if dry_run:
                    print(f"  MOVEDIR plots/{d.name} → {ct}/plots/")
                else:
                    dst = ATAC_BASE / ct / "plots"
                    dst.mkdir(parents=True, exist_ok=True)
                    for f in d.iterdir():
                        move(f, dst / f.name, dry_run)
                    try:
                        d.rmdir()
                    except OSError:
                        pass
            elif d.is_dir():
                # Unrecognized subdir → _summary/plots/
                move(d, ATAC_BASE / "_summary" / "plots" / d.name, dry_run)
            else:
                move(d, ATAC_BASE / "_summary" / "plots" / d.name, dry_run)

    print("\n=== ATAC: backup unrelated dirs ===")
    for name in ["bigwig_regions", "region_beds", "ILS_UltraRobust"]:
        backup(ATAC_BASE / name, dry_run)
    # Also backup top-level loose PDFs/PNGs
    for f in ATAC_BASE.glob("*.pdf"):
        backup(f, dry_run)
    for f in ATAC_BASE.glob("*.png"):
        backup(f, dry_run)


def migrate_rna(dry_run: bool):
    print("\n=== RNA: evolutionary ===")
    evo_dir = RNA_BASE / "rna_differential_results"
    if evo_dir.exists():
        for ct_dir in sorted(evo_dir.iterdir()):
            if not ct_dir.is_dir() or ct_dir.name == "species_specific":
                continue
            ct = ct_dir.name
            for f in ct_dir.glob("*.csv"):
                # {Contrast}_RNA_DE.csv → {CT}/rna_evolutionary/{Contrast}_DE.csv
                new_name = f.name.replace("_RNA_DE.csv", "_DE.csv")
                move(f, RNA_BASE / ct / "rna_evolutionary" / new_name, dry_run)
            for rds in ct_dir.glob("*.rds"):
                move(rds, RNA_BASE / ct / "rna_evolutionary" / rds.name, dry_run)

        print("\n=== RNA: species_specific ===")
        sp_dir = evo_dir / "species_specific"
        if sp_dir.exists():
            for ct_dir in sorted(sp_dir.iterdir()):
                if not ct_dir.is_dir():
                    continue
                ct = ct_dir.name
                for f in ct_dir.glob("*.csv"):
                    # {Species}_vs_rest_RNA_DE.csv → {CT}/rna_species/{Species}_vs_rest_DE.csv
                    new_name = f.name.replace("_RNA_DE.csv", "_DE.csv")
                    move(f, RNA_BASE / ct / "rna_species" / new_name, dry_run)
            for rds in sp_dir.glob("*.rds"):
                move(rds, RNA_BASE / "_summary" / rds.name, dry_run)

        # top-level RDS in rna_differential_results
        for rds in evo_dir.glob("*.rds"):
            move(rds, RNA_BASE / "_summary" / rds.name, dry_run)
        for csv in evo_dir.glob("*.csv"):
            move(csv, RNA_BASE / "_summary" / csv.name, dry_run)

        # backup now-empty dirs
        if not dry_run:
            try:
                sp_dir.rmdir()
            except OSError:
                pass
            try:
                evo_dir.rmdir()
                print("  Removed now-empty rna_differential_results/")
            except OSError:
                backup(evo_dir, dry_run)

    print("\n=== RNA: summary-level loose files ===")
    for f in RNA_BASE.glob("RNA_DE_*.rds"):
        move(f, RNA_BASE / "_summary" / f.name, dry_run)
    for f in RNA_BASE.glob("RNA_DE_*.csv"):
        move(f, RNA_BASE / "_summary" / f.name, dry_run)
    for f in RNA_BASE.glob("*.pdf"):
        move(f, RNA_BASE / "_summary" / f.name, dry_run)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true",
                        help="Print planned moves without executing")
    args = parser.parse_args()

    if args.dry_run:
        print("=== DRY RUN — no files will be moved ===\n")
    else:
        print(f"=== LIVE RUN — backup dir: {BACKUP} ===\n")
        BACKUP.mkdir(parents=True, exist_ok=True)

    migrate_atac(args.dry_run)
    migrate_rna(args.dry_run)

    if args.dry_run:
        print("\n=== DRY RUN complete — run without --dry-run to execute ===")
    else:
        print(f"\n=== Migration complete ===")
        print(f"Backed-up old files: {BACKUP}")
        print("Delete BACKUP/ after verifying the new structure looks correct.")


if __name__ == "__main__":
    main()
