#!/usr/bin/env Rscript
# Generate per-region pseudobulk RDS files for downstream notebooks.
# Requires pb_data_shared.rds to already exist (run 10a first, or use existing).

suppressPackageStartupMessages({
  library(DESeq2)
  library(arrow)
  library(dplyr)
  library(tibble)
  library(readr)
})

SCRIPT_DIR <- "/mp/nfs4krb5/bs-nas02/bsse_group_treutlein/USERS/jjans/analysis/adult_intestine/peaks/peak_calling/atac_pipeline"
source(file.path(SCRIPT_DIR, "src", "deseq2_utils.R"))

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
BASE      <- "/links/groups/treutlein/USERS/jjans/analysis/adult_intestine/peaks"
QUANT_DIR <- file.path(BASE, "cross_species_consensus_v3/12_fragment_matrices")
SAVE_DIR  <- file.path(BASE, "cross_species_consensus_v3/13_deseq2_R_pseudobulk/pseudobulk")
SPECIES   <- c("Human", "Bonobo", "Chimpanzee", "Gorilla", "Macaque", "Marmoset")

MIN_COUNTS_REGION <- 30000
MIN_CELLS_REGION  <- 50

message("Loading existing shared-peak data to get peak universe...")
shared_data   <- readRDS(file.path(SAVE_DIR, "pb_data_shared.rds"))
shared_peak_ids <- rownames(shared_data$counts)
meta_shared     <- shared_data$meta
message("  Shared peaks: ", length(shared_peak_ids))
message("  Shared-peak samples: ", ncol(shared_data$counts))

# ---------------------------------------------------------------------------
# Load raw (non-aggregated) data
# ---------------------------------------------------------------------------
message("Loading raw pseudobulk data (non-aggregated)...")
raw <- load_pseudobulk_data(QUANT_DIR, SPECIES)
message("  Species loaded: ", paste(names(raw$all_counts), collapse = ", "))

# Inner join for shared-peak per-region data
message("Building inner-join matrix (shared peaks, non-aggregated)...")
shared_raw <- merge_pseudobulk(raw$all_counts, raw$all_meta, join_type = "inner")

# Filter Adults
meta_preg   <- shared_raw$meta[shared_raw$meta$age == "Adult", ]
counts_preg <- shared_raw$counts[, rownames(meta_preg), drop = FALSE]
message("  Adult samples (inner): ", ncol(counts_preg))

# Region-level QC
meta_preg$total_counts <- colSums(counts_preg)
keep_preg <- meta_preg$total_counts >= MIN_COUNTS_REGION
if ("n_cells" %in% colnames(meta_preg)) {
  keep_preg <- keep_preg & (meta_preg$n_cells >= MIN_CELLS_REGION)
}
meta_preg   <- meta_preg[keep_preg, ]
counts_preg <- counts_preg[, rownames(meta_preg), drop = FALSE]
message("  After region QC: ", ncol(counts_preg), " samples")
print(table(meta_preg$region, meta_preg$species))

# Clean factors
meta_preg$cell_type <- as.factor(make.names(as.character(meta_preg$cell_type)))
meta_preg$species   <- factor(meta_preg$species, levels = SPECIES)
meta_preg$region    <- factor(meta_preg$region)

# Restrict to shared peaks
avail_peaks        <- intersect(shared_peak_ids, rownames(counts_preg))
counts_preg_shared <- counts_preg[avail_peaks, , drop = FALSE]
message("  Shared-peak per-region matrix: ", nrow(counts_preg_shared),
        " peaks x ", ncol(counts_preg_shared), " samples")

saveRDS(list(counts = counts_preg_shared, meta = meta_preg),
        file.path(SAVE_DIR, "pb_data_shared_per_region.rds"))
message("  Saved pb_data_shared_per_region.rds")

# ---------------------------------------------------------------------------
# Union join for per-region data
# ---------------------------------------------------------------------------
message("Building union-join matrix (non-aggregated)...")
union_raw <- merge_pseudobulk(raw$all_counts, raw$all_meta, join_type = "union")

meta_union_preg   <- union_raw$meta[union_raw$meta$age == "Adult", ]
counts_union_preg <- union_raw$counts[, rownames(meta_union_preg), drop = FALSE]

meta_union_preg$total_counts <- colSums(counts_union_preg)
keep_up <- meta_union_preg$total_counts >= MIN_COUNTS_REGION
if ("n_cells" %in% colnames(meta_union_preg)) {
  keep_up <- keep_up & (meta_union_preg$n_cells >= MIN_CELLS_REGION)
}
meta_union_preg   <- meta_union_preg[keep_up, ]
counts_union_preg <- counts_union_preg[, rownames(meta_union_preg), drop = FALSE]

meta_union_preg$cell_type <- as.factor(make.names(as.character(meta_union_preg$cell_type)))
meta_union_preg$species   <- factor(meta_union_preg$species, levels = SPECIES)
meta_union_preg$region    <- factor(meta_union_preg$region)

message("  Union per-region matrix: ", nrow(counts_union_preg),
        " peaks x ", ncol(counts_union_preg), " samples")

saveRDS(list(counts = counts_union_preg, meta = meta_union_preg),
        file.path(SAVE_DIR, "pb_data_union_per_region.rds"))
message("  Saved pb_data_union_per_region.rds")

message("\n=== Per-region data generation complete ===")
message("  pb_data_shared_per_region.rds : ", nrow(counts_preg_shared), " peaks x ",
        ncol(counts_preg_shared), " samples")
message("  pb_data_union_per_region.rds  : ", nrow(counts_union_preg), " peaks x ",
        ncol(counts_union_preg), " samples")
message("  Regions: ", paste(levels(meta_preg$region), collapse = ", "))
