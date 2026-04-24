# =============================================================================
# DESeq2 Pipeline Utility Functions for Cross-Species ATAC-seq
# =============================================================================
# Source this file at the top of analysis notebooks:
#   source("../src/deseq2_utils.R")
#
# Provides: data loading, pseudobulk aggregation, filtering, PCA, DESeq2
#           helpers, volcano plots, heatmaps, motif enrichment (parallel),
#           module identification, BED conversion, GO/GREAT wrappers.
# =============================================================================

# =============================================================================
# SECTION 1: DATA LOADING & PREPROCESSING
# =============================================================================

#' Load pseudobulk counts and metadata from parquet files
#'
#' @param quant_dir Base directory containing per-species subdirectories
#' @param species Character vector of species to load
#' @return List with elements: all_counts, all_meta
load_pseudobulk_data <- function(quant_dir, species) {
  suppressPackageStartupMessages({
    library(arrow)
    library(dplyr)
  })

  all_counts <- list()
  all_meta   <- list()

  for (sp in species) {
    sp_dir      <- file.path(quant_dir, sp)
    counts_file <- file.path(sp_dir, "pseudobulk_counts.parquet")
    meta_file   <- file.path(sp_dir, "pseudobulk_groups.parquet")

    if (!file.exists(counts_file)) {
      message("Skipping ", sp, ": files not found at ", sp_dir)
      next
    }

    counts <- as.data.frame(read_parquet(counts_file))

    # Safety net: deduplicate region_ids if any remain
    if (any(duplicated(counts$region_id))) {
      dup_ids <- unique(counts$region_id[duplicated(counts$region_id)])
      message("  Warning: ", length(dup_ids), " duplicated region IDs in ", sp, " — merging by sum.")
      numeric_mat <- as.matrix(counts[, -which(names(counts) == "region_id")])
      merged_mat  <- rowsum(numeric_mat, group = counts$region_id)
      counts <- as.data.frame(merged_mat)
      counts$region_id <- rownames(counts)
    } else {
      rownames(counts) <- counts$region_id
    }
    counts$region_id <- NULL

    meta <- as.data.frame(read_parquet(meta_file))
    new_labels       <- paste0(sp, "_", meta$label)
    colnames(counts) <- new_labels
    meta$sample_id   <- new_labels
    meta$species     <- sp

    all_counts[[sp]] <- counts
    all_meta[[sp]]   <- meta
    message("Loaded ", sp, ": ", ncol(counts), " pseudobulk samples, ",
            nrow(counts), " peaks")
  }

  return(list(all_counts = all_counts, all_meta = all_meta))
}


#' Merge species count matrices on shared (intersected) peaks
#'
#' @param all_counts Named list of count data.frames (from load_pseudobulk_data)
#' @param all_meta   Named list of metadata data.frames
#' @param join_type  "inner" (shared peaks only) or "union" (all peaks, 0-filled)
#' @return List with elements: counts_merged, meta_merged
merge_pseudobulk <- function(all_counts, all_meta, join_type = "inner") {
  if (join_type == "inner") {
    shared_peaks   <- Reduce(intersect, lapply(all_counts, rownames))
    counts_merged  <- do.call(cbind, unname(lapply(all_counts,
                        function(x) x[shared_peaks, , drop = FALSE])))
    message("Inner join: ", length(shared_peaks), " shared peaks")
  } else {
    union_peaks <- Reduce(union, lapply(all_counts, rownames))
    counts_list <- lapply(all_counts, function(mat) {
      new_mat <- matrix(0, nrow = length(union_peaks), ncol = ncol(mat),
                        dimnames = list(union_peaks, colnames(mat)))
      common <- intersect(rownames(mat), union_peaks)
      new_mat[common, ] <- as.matrix(mat[common, ])
      return(as.data.frame(new_mat))
    })
    counts_merged <- do.call(cbind, unname(counts_list))
    message("Union join: ", length(union_peaks), " total peaks")
  }

  counts_merged[is.na(counts_merged)] <- 0
  meta_merged <- do.call(rbind, unname(all_meta))
  rownames(meta_merged) <- meta_merged$sample_id
  counts_merged <- counts_merged[, rownames(meta_merged)]

  message("Merged matrix: ", nrow(counts_merged), " peaks x ",
          ncol(counts_merged), " samples")
  return(list(counts = counts_merged, meta = meta_merged))
}


#' Aggregate pseudobulk samples by collapsing specified metadata columns
#'
#' @param counts_df  Count matrix (peaks x samples)
#' @param meta_df    Metadata data.frame (rows = samples)
#' @param merge_cols Columns to collapse (set to "Merged")
#' @param sum_cols   Numeric metadata columns to sum (e.g. "n_cells")
#' @return List with elements: counts, meta
aggregate_pseudobulk <- function(counts_df, meta_df,
                                 merge_cols, sum_cols = c("n_cells")) {
  suppressPackageStartupMessages(library(dplyr))

  exclude_cols <- c(merge_cols, sum_cols, "sample_id", "label",
                    "total_counts", "agg_id")
  group_cols   <- setdiff(colnames(meta_df), exclude_cols)

  message("Aggregating across: ", paste(merge_cols, collapse = ", "))
  message("Grouping by: ", paste(group_cols, collapse = ", "))

  # Create unique aggregation ID
  temp_meta       <- meta_df[, group_cols, drop = FALSE]
  meta_df$agg_id  <- apply(temp_meta, 1, paste, collapse = "_")

  # Sum counts using fast C-level rowsum
  counts_t   <- t(counts_df)
  merged_t   <- rowsum(counts_t, group = meta_df$agg_id)
  new_counts <- as.data.frame(t(merged_t))

  # Rebuild metadata
  new_meta_base <- meta_df %>%
    group_by(agg_id) %>%
    slice(1) %>%
    select(-any_of(sum_cols)) %>%
    ungroup()

  actual_sum_cols <- intersect(sum_cols, colnames(meta_df))
  if (length(actual_sum_cols) > 0) {
    new_meta_sums <- meta_df %>%
      group_by(agg_id) %>%
      summarize(across(all_of(actual_sum_cols), sum, na.rm = TRUE),
                .groups = "drop")
    new_meta <- left_join(new_meta_base, new_meta_sums, by = "agg_id")
  } else {
    new_meta <- new_meta_base
  }

  for (col in merge_cols) {
    if (col %in% colnames(new_meta)) {
      new_meta[[col]] <- "Merged"
    }
  }

  new_meta$sample_id <- new_meta$agg_id
  new_meta$label     <- new_meta$agg_id
  new_meta           <- as.data.frame(new_meta)
  rownames(new_meta) <- new_meta$sample_id
  new_counts         <- new_counts[, rownames(new_meta)]

  message("Reduced from ", ncol(counts_df), " to ", ncol(new_counts), " samples.")
  return(list(counts = new_counts, meta = new_meta))
}


#' Filter samples by minimum counts and minimum cells
#'
#' @param counts   Count matrix
#' @param meta     Metadata data.frame
#' @param min_counts Minimum total counts per sample
#' @param min_cells  Minimum cells per pseudobulk
#' @return List with elements: counts, meta
filter_samples <- function(counts, meta,
                           min_counts = 50000, min_cells = 100) {
  meta$total_counts <- colSums(counts)

  if ("n_cells" %in% colnames(meta)) {
    keep <- meta$total_counts >= min_counts & meta$n_cells >= min_cells
    message("Filtering: counts >= ", min_counts, " AND cells >= ", min_cells)
  } else {
    keep <- meta$total_counts >= min_counts
    message("Filtering: counts >= ", min_counts, " (n_cells column not found)")
  }

  meta_f   <- meta[keep, ]
  counts_f <- counts[, rownames(meta_f)]
  message("Kept ", sum(keep), " / ", length(keep), " samples.")
  return(list(counts = counts_f, meta = meta_f))
}


#' Subset to cell types shared across all species
#'
#' @param counts Count matrix
#' @param meta   Metadata with species and cell_type columns
#' @return List with elements: counts, meta
subset_shared_celltypes <- function(counts, meta) {
  ct_per_species <- split(as.character(meta$cell_type), meta$species)
  shared_ct      <- Reduce(intersect, ct_per_species)

  keep         <- meta$cell_type %in% shared_ct
  counts_out   <- counts[, keep]
  meta_out     <- meta[keep, ]
  meta_out$cell_type <- droplevels(meta_out$cell_type)

  message("Shared cell types (", length(shared_ct), "): ",
          paste(shared_ct, collapse = ", "))
  message("Kept ", ncol(counts_out), " samples.")
  return(list(counts = counts_out, meta = meta_out))
}


#' Species-aware peak filtering
#'
#' @param counts  Count matrix
#' @param meta    Metadata with species column
#' @param species Character vector of all species
#' @param min_samples_per_species Minimum samples with signal within a species
#' @param min_species_active Minimum species that must have signal
#' @param max_count_threshold Peak must have >= this max count in some sample
#' @return Logical vector of peaks to keep
species_aware_peak_filter <- function(counts, meta, species,
                                      min_samples_per_species = 2,
                                      min_species_active = 6,
                                      max_count_threshold = 150) {
  n_before <- nrow(counts)

  keep_count <- apply(counts, 1, max) >= max_count_threshold

  active_in_species <- sapply(species, function(sp) {
    sp_cols <- meta$species == sp
    if (sum(sp_cols) > 0) {
      rowSums(counts[, sp_cols, drop = FALSE] >= 10) >= min_samples_per_species
    } else {
      rep(FALSE, nrow(counts))
    }
  })

  keep_species <- rowSums(active_in_species) >= min_species_active
  keep_peaks   <- keep_count & keep_species

  message("Species-aware filtering: kept ", sum(keep_peaks), " / ", n_before,
          " peaks (", round(100 * sum(keep_peaks) / n_before, 1), "%)")
  return(keep_peaks)
}


# =============================================================================
# SECTION 2: PEAK / ANNOTATION UTILITIES
# =============================================================================

#' Swiss-army-knife peak info lookup from master annotation
#'
#' @param peak_ids  Character vector of peak IDs
#' @param species   Species name (e.g. "Human")
#' @param info_type One of "coordinates", "gene", "distance", "status"
#' @param anno_df   Master annotation data.frame
#' @return data.frame with requested information
get_peak_info <- function(peak_ids, species, info_type = "coordinates",
                          anno_df) {
  chr_col   <- paste0(species, "_chr")
  start_col <- paste0(species, "_start")
  end_col   <- paste0(species, "_end")
  gene_col  <- paste0(species, "_gene")
  dist_col  <- paste0(species, "_gene_dist")
  det_col   <- paste0(species, "_det")
  orth_col  <- paste0(species, "_orth")

  sub_anno <- anno_df[anno_df$peak_id %in% peak_ids, ]

  if (info_type == "coordinates") {
    if (!all(c(chr_col, start_col, end_col) %in% colnames(anno_df))) {
      stop("Coordinate columns for ", species, " not found in annotation.")
    }
    res <- sub_anno[!is.na(sub_anno[[chr_col]]) & sub_anno[[chr_col]] != ".", ]
    if (nrow(res) == 0) return(NULL)
    return(data.frame(
      seqnames = res[[chr_col]],
      start    = as.integer(res[[start_col]]),
      end      = as.integer(res[[end_col]]),
      name     = res$peak_id
    ))
  } else if (info_type == "gene") {
    return(data.frame(peak_id = sub_anno$peak_id, gene = sub_anno[[gene_col]]))
  } else if (info_type == "distance") {
    return(data.frame(peak_id = sub_anno$peak_id,
                      distance = sub_anno[[dist_col]]))
  } else if (info_type == "status") {
    return(data.frame(peak_id = sub_anno$peak_id,
                      status = paste0(sub_anno[[det_col]], "/",
                                      sub_anno[[orth_col]])))
  }
}


#' Write peaks to BED file
#'
#' @param peak_ids  Character vector of peak IDs
#' @param species   Species name for coordinate lookup
#' @param anno_df   Master annotation data.frame
#' @param out_file  Output BED file path
#' @return Invisibly returns the BED data.frame (or NULL)
write_peaks_bed <- function(peak_ids, species, anno_df, out_file) {
  options(scipen = 999)
  bed_df <- get_peak_info(peak_ids, species, "coordinates", anno_df)
  if (!is.null(bed_df) && nrow(bed_df) > 0) {
    dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
    write.table(bed_df, out_file, sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = FALSE)
    message("  Saved ", nrow(bed_df), " regions to: ", basename(out_file))
  }
  invisible(bed_df)
}


# =============================================================================
# SECTION 3: DESeq2 HELPERS
# =============================================================================

#' Create a numeric contrast vector for DESeq2 (zero-intercept model)
#'
#' @param sp_pos    Character vector of positive-side species
#' @param sp_neg    Character vector of negative-side species
#' @param res_names Character vector from resultsNames(dds)
#' @return Numeric contrast vector
make_contrast_vector <- function(sp_pos, sp_neg, res_names) {
  vec <- rep(0, length(res_names))
  names(vec) <- res_names

  pos_names <- paste0("species", sp_pos)
  neg_names <- paste0("species", sp_neg)

  valid_pos <- intersect(pos_names, res_names)
  valid_neg <- intersect(neg_names, res_names)

  if (length(valid_pos) > 0) vec[valid_pos] <-  1 / length(valid_pos)
  if (length(valid_neg) > 0) vec[valid_neg] <- -1 / length(valid_neg)

  return(vec)
}


#' Run DESeq2 for each cell type with all-vs-one species contrasts
#'
#' Uses shared peaks (peaks present in all species). Design: ~ 0 + species
#'
#' @param counts_shared Filtered count matrix (shared peaks x samples)
#' @param meta_shared   Metadata with cell_type, species columns
#' @param species       Character vector of species names
#' @param out_dir       Directory to save per-cell-type CSV results
#' @return Named list: res_list[[cell_type]][[species]] = DESeq2 results
run_deseq2_shared_peaks <- function(counts_shared, meta_shared, species,
                                    out_dir,
                                    filter_outliers = TRUE, sd_thresh = 2.5) {
  suppressPackageStartupMessages({
    library(DESeq2)
    library(dplyr)
  })

  summary_dir <- file.path(out_dir, "_summary")
  dir.create(summary_dir, showWarnings = FALSE, recursive = TRUE)

  res_list       <- list()
  outlier_log_df <- data.frame(cell_type = character(), contrast = character(),
                               sample_id = character(), species = character(),
                               log2_lib_size = numeric(),
                               stringsAsFactors = FALSE)
  cell_types <- levels(meta_shared$cell_type)

  for (ct in cell_types) {
    message(sprintf("\n=== [Shared Peaks] %s ===", ct))

    meta_ct   <- meta_shared[meta_shared$cell_type == ct, ]
    counts_ct <- counts_shared[, rownames(meta_ct)]

    present_species <- unique(as.character(meta_ct$species))
    if (length(present_species) < 2) {
      message("  Skipping: fewer than 2 species.")
      next
    }

    # Local peak filter for this cell type
    keep_ct            <- rowSums(counts_ct >= 10) >= 3
    counts_ct_filtered <- counts_ct[keep_ct, ]

    # Per-cell-type outlier filtering (applied once, before fitting the model)
    if (filter_outliers) {
      outlier_flags <- detect_per_contrast_outliers(counts_ct_filtered, meta_ct,
                                                    sd_thresh = sd_thresh)
      removed_ids <- names(outlier_flags)[outlier_flags]
      if (length(removed_ids) > 0) {
        lib_vec        <- log2(colSums(counts_ct_filtered) + 1)
        outlier_log_df <- append_outlier_log(outlier_log_df, ct, "shared_peaks_all_species",
                                             removed_ids, meta_ct, lib_vec)
        keep_samp          <- colnames(counts_ct_filtered)[!outlier_flags]
        meta_ct            <- meta_ct[keep_samp, ]
        counts_ct_filtered <- counts_ct_filtered[, keep_samp, drop = FALSE]
        meta_ct$species    <- droplevels(meta_ct$species)
        message(sprintf("  Outliers removed: %s", paste(removed_ids, collapse = ", ")))
      }
    }

    if (ncol(counts_ct_filtered) < 4 ||
        length(unique(as.character(meta_ct$species))) < 2) {
      message("  Skipping: too few samples after outlier removal.")
      next
    }

    dds_ct <- DESeqDataSetFromMatrix(
      countData = round(counts_ct_filtered),
      colData   = meta_ct,
      design    = ~ 0 + species
    )
    dds_ct <- estimateSizeFactors(dds_ct, type = "poscounts")
    dds_ct <- DESeq(dds_ct, quiet = TRUE)

    res_list[[ct]] <- list()
    ct_dir <- file.path(out_dir, ct, "atac_shared")
    dir.create(ct_dir, showWarnings = FALSE, recursive = TRUE)

    res_names <- resultsNames(dds_ct)

    for (target_sp in present_species) {
      target_term <- paste0("species", target_sp)
      if (!target_term %in% res_names) next

      contrast_vec <- rep(-1 / (length(res_names) - 1), length(res_names))
      target_idx   <- which(res_names == target_term)
      contrast_vec[target_idx] <- 1

      res_sp         <- results(dds_ct, contrast = contrast_vec)
      res_sp_ordered <- res_sp[order(res_sp$padj), ]
      res_list[[ct]][[target_sp]] <- res_sp_ordered

      res_df            <- as.data.frame(res_sp_ordered)
      res_df$peak_id    <- rownames(res_df)
      res_df            <- res_df[, c("peak_id", setdiff(names(res_df), "peak_id"))]
      out_file          <- file.path(ct_dir, paste0(target_sp, "_vs_rest_DA.csv"))
      write.csv(res_df, out_file, row.names = FALSE)

      sig_hits <- sum(res_sp$padj < 0.05 & res_sp$log2FoldChange > 1, na.rm = TRUE)
      message(sprintf("  %-12s specific: %d (padj<0.05, LFC>1)", target_sp, sig_hits))
    }
  }

  saveRDS(res_list, file.path(out_dir, "_summary", "DESeq2_res_list_shared.rds"))
  if (nrow(outlier_log_df) > 0) {
    write.csv(outlier_log_df,
              file.path(out_dir, "_summary", "outlier_shared_DA.csv"),
              row.names = FALSE)
    message("Shared-peaks outlier log saved (", nrow(outlier_log_df), " removals).")
  }
  message("\nCheckpoint: shared-peaks res_list saved.")
  return(res_list)
}


#' Run DESeq2 using only peaks that exist in the contrast species
#'
#' For each cell type and each target species, uses a union-joined count
#' matrix (0-filled for missing peaks) so species-specific insertions are
#' included. Design: ~ 0 + species.
#'
#' @param counts_union   Union-joined count matrix (peaks x samples)
#' @param meta_shared    Metadata with cell_type, species columns
#' @param species        Character vector of species names
#' @param out_dir        Directory to save results
#' @param min_samples    Minimum samples with signal for a peak
#' @return Named list: res_list[[cell_type]][[species]] = DESeq2 results
run_deseq2_contrast_peaks <- function(counts_union, meta_shared, species,
                                      out_dir, min_samples = 3) {
  suppressPackageStartupMessages(library(DESeq2))

  dir.create(file.path(out_dir, "_summary"), showWarnings = FALSE, recursive = TRUE)

  res_list   <- list()
  cell_types <- levels(meta_shared$cell_type)

  message("Union peaks for contrast-specific testing: ", nrow(counts_union))

  # Subset to valid samples (must exist in both counts and metadata)
  valid_samples <- intersect(rownames(meta_shared), colnames(counts_union))
  counts_union  <- counts_union[, valid_samples, drop = FALSE]

  for (ct in cell_types) {
    message(sprintf("\n=== [Contrast Peaks] %s ===", ct))

    meta_ct   <- meta_shared[meta_shared$cell_type == ct, ]
    ct_samples <- intersect(rownames(meta_ct), colnames(counts_union))
    counts_ct <- counts_union[, ct_samples, drop = FALSE]
    meta_ct   <- meta_ct[ct_samples, ]

    present_species <- unique(as.character(meta_ct$species))
    if (length(present_species) < 2) {
      message("  Skipping: fewer than 2 species.")
      next
    }

    # Light filter: peak active in at least min_samples donors
    keep_ct   <- rowSums(counts_ct >= 10) >= min_samples
    counts_ct <- counts_ct[keep_ct, ]

    dds_ct <- DESeqDataSetFromMatrix(
      countData = round(counts_ct),
      colData   = meta_ct,
      design    = ~ 0 + species
    )
    dds_ct <- estimateSizeFactors(dds_ct, type = "poscounts")
    dds_ct <- DESeq(dds_ct, quiet = TRUE)

    res_list[[ct]] <- list()
    ct_dir <- file.path(out_dir, ct, "atac_contrast")
    dir.create(ct_dir, showWarnings = FALSE, recursive = TRUE)

    res_names <- resultsNames(dds_ct)

    for (target_sp in present_species) {
      target_term <- paste0("species", target_sp)
      if (!target_term %in% res_names) next

      contrast_vec <- rep(-1 / (length(res_names) - 1), length(res_names))
      target_idx   <- which(res_names == target_term)
      contrast_vec[target_idx] <- 1

      res_sp         <- results(dds_ct, contrast = contrast_vec)
      res_sp_ordered <- res_sp[order(res_sp$padj), ]
      res_list[[ct]][[target_sp]] <- res_sp_ordered

      res_df         <- as.data.frame(res_sp_ordered)
      res_df$peak_id <- rownames(res_df)
      res_df         <- res_df[, c("peak_id", setdiff(names(res_df), "peak_id"))]
      out_file       <- file.path(ct_dir, paste0(target_sp, "_vs_rest_DA.csv"))
      write.csv(res_df, out_file, row.names = FALSE)

      sig_hits <- sum(res_sp$padj < 0.05 & res_sp$log2FoldChange > 1, na.rm = TRUE)
      message(sprintf("  %-12s specific: %d (padj<0.05, LFC>1)", target_sp, sig_hits))
    }
  }

  saveRDS(res_list, file.path(out_dir, "_summary", "DESeq2_res_list_contrast.rds"))
  message("\nCheckpoint: contrast-peaks res_list saved.")
  return(res_list)
}


# =============================================================================
# SECTION 3b: EVOLUTIONARY BRANCH TESTING
# =============================================================================

#' Build orthology index from master annotation
#'
#' Returns a logical matrix (peak_id x species) indicating whether each peak
#' has physical orthologous sequence in each species.
#'
#' @param anno_df  Master annotation data.frame
#' @param species  Character vector of species names
#' @return Logical matrix (peaks x species)
build_orthology_index <- function(anno_df, species) {
  ortho_mat <- matrix(FALSE, nrow = nrow(anno_df), ncol = length(species),
                      dimnames = list(anno_df$peak_id, species))
  for (sp in species) {
    orth_col <- paste0(sp, "_orth")
    if (orth_col %in% colnames(anno_df)) {
      ortho_mat[, sp] <- as.logical(anno_df[[orth_col]])
      ortho_mat[is.na(ortho_mat[, sp]), sp] <- FALSE
    }
  }
  message("Orthology index: ", nrow(ortho_mat), " peaks x ", ncol(ortho_mat),
          " species")
  return(ortho_mat)
}


#' Define evolutionary branch contrasts (phylogenetic nodes, ILS, pairwise)
#'
#' Returns a named list of contrasts, each with pos/neg species vectors.
#'
#' @return Named list of list(pos = ..., neg = ...)
define_evolutionary_contrasts <- function() {
  contrasts <- list(
    # --- Phylogenetic Branch Nodes (Synapomorphies) ---
    "Node1_Human_vs_Pan"           = list(pos = c("Human"),
                                          neg = c("Chimpanzee", "Bonobo")),
    "Node2_AfricanApes_vs_Gorilla" = list(pos = c("Human", "Chimpanzee", "Bonobo"),
                                          neg = c("Gorilla")),
    "Node3_GreatApes_vs_Macaque"   = list(pos = c("Human", "Chimpanzee", "Bonobo", "Gorilla"),
                                          neg = c("Macaque")),
    "Node4_OldWorld_vs_Marmoset"   = list(pos = c("Human", "Chimpanzee", "Bonobo", "Gorilla", "Macaque"),
                                          neg = c("Marmoset")),

    # --- ILS / Topology Discordance (African Apes) ---
    "ILS_HumanGorilla_vs_Pan"         = list(pos = c("Human", "Gorilla"),
                                             neg = c("Chimpanzee", "Bonobo")),
    "ILS_HumanChimp_vs_GorillaBonobo" = list(pos = c("Human", "Chimpanzee"),
                                             neg = c("Gorilla", "Bonobo")),
    "ILS_HumanBonobo_vs_ChimpGorilla" = list(pos = c("Human", "Bonobo"),
                                             neg = c("Chimpanzee", "Gorilla")),

    # --- Single-Species Divergence ---
    "Div_Human_vs_Apes"        = list(pos = c("Human"),
                                      neg = c("Chimpanzee", "Bonobo", "Gorilla")),
    "Div_Chimp_vs_Apes"        = list(pos = c("Chimpanzee"),
                                      neg = c("Human", "Bonobo", "Gorilla")),
    "Div_Bonobo_vs_Apes"       = list(pos = c("Bonobo"),
                                      neg = c("Human", "Chimpanzee", "Gorilla")),
    "Div_Gorilla_vs_Apes"      = list(pos = c("Gorilla"),
                                      neg = c("Human", "Chimpanzee", "Bonobo")),
    "Div_Macaque_vs_GreatApes" = list(pos = c("Macaque"),
                                      neg = c("Human", "Chimpanzee", "Bonobo", "Gorilla")),

    # --- Pairwise Direct Comparisons ---
    "Pair_Human_vs_Gorilla"  = list(pos = c("Human"), neg = c("Gorilla")),
    "Pair_Human_vs_Chimp"    = list(pos = c("Human"), neg = c("Chimpanzee")),
    "Pair_Human_vs_Bonobo"   = list(pos = c("Human"), neg = c("Bonobo")),
    "Pair_Human_vs_Macaque"  = list(pos = c("Human"), neg = c("Macaque")),
    "Pair_Human_vs_Marmoset" = list(pos = c("Human"), neg = c("Marmoset")),

    # --- Broader Clade Contrasts ---
    # Human vs all non-human primates (broadest human-specific test)
    "Div_Human_vs_AllPrimates" = list(pos = c("Human"),
                                      neg = c("Chimpanzee", "Bonobo", "Gorilla",
                                              "Macaque", "Marmoset")),
    # Apes (Hominoidea) vs monkeys (Cercopithecidae + Cebidae) — major taxonomic split
    "Node_Apes_vs_Monkeys"     = list(pos = c("Human", "Chimpanzee", "Bonobo", "Gorilla"),
                                      neg = c("Macaque", "Marmoset"))
  )
  message("Defined ", length(contrasts), " evolutionary contrasts.")
  return(contrasts)
}


# =============================================================================
# SECTION 3c: PER-CONTRAST OUTLIER DETECTION
# =============================================================================

#' Detect outlier samples within a single contrast's count matrix
#'
#' For each species group in the contrast, flags samples whose log2 library size
#' falls more than sd_thresh SDs below the species mean. Never reduces any
#' species below min_per_species surviving samples.
#'
#' @param counts      Count matrix for this contrast (peaks x samples)
#' @param meta        Metadata for those samples (must have 'species' column)
#' @param sd_thresh   SD threshold below species mean to flag (default 2.5)
#' @param min_per_species Minimum samples to retain per species (default 2)
#' @return Named logical vector: TRUE = outlier. Names = sample IDs.
detect_per_contrast_outliers <- function(counts, meta,
                                         sd_thresh = 2.5,
                                         min_per_species = 2) {
  lib_sizes     <- log2(colSums(counts) + 1)
  outlier_flags <- setNames(rep(FALSE, ncol(counts)), colnames(counts))

  for (sp in unique(as.character(meta$species))) {
    sp_cols <- intersect(rownames(meta)[meta$species == sp], colnames(counts))
    if (length(sp_cols) < 3) next          # can't estimate spread with < 3

    sp_lib  <- lib_sizes[sp_cols]
    sp_sd   <- sd(sp_lib)
    if (sp_sd == 0) next

    low_signal <- sp_lib < (mean(sp_lib) - sd_thresh * sp_sd)

    # Safety: never drop below min_per_species
    if (sum(!low_signal) >= min_per_species) {
      outlier_flags[sp_cols[low_signal]] <- TRUE
    }
  }
  return(outlier_flags)
}


#' Append to an outlier log data frame
#'
#' @param log_df      Existing log data.frame (may be empty)
#' @param cell_type   Cell type label
#' @param contrast    Contrast name
#' @param sample_ids  Character vector of removed sample IDs
#' @param meta        Metadata for species lookup
#' @param lib_sizes   Named numeric vector of log2 library sizes
#' @return Updated log data.frame
append_outlier_log <- function(log_df, cell_type, contrast,
                               sample_ids, meta, lib_sizes) {
  if (length(sample_ids) == 0) return(log_df)
  new_rows <- data.frame(
    cell_type    = cell_type,
    contrast     = contrast,
    sample_id    = sample_ids,
    species      = as.character(meta[sample_ids, "species"]),
    log2_lib_size = round(lib_sizes[sample_ids], 3),
    stringsAsFactors = FALSE
  )
  rbind(log_df, new_rows)
}


#' Save outlier removal summary to CSV and print a compact report
#'
#' @param outlier_log_df Data frame from append_outlier_log() calls
#' @param out_dir        Output directory
save_outlier_summary <- function(outlier_log_df, out_dir) {
  out_file <- file.path(out_dir, "_summary", "outlier_samples_removed.csv")
  dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)

  if (nrow(outlier_log_df) == 0) {
    message("\nOutlier filtering: no samples removed across all contrasts.")
    write.csv(data.frame(cell_type = character(), contrast = character(),
                         sample_id = character(), species = character(),
                         log2_lib_size = numeric()),
              out_file, row.names = FALSE)
    return(invisible(NULL))
  }

  write.csv(outlier_log_df, out_file, row.names = FALSE)

  # Print compact summary grouped by sample
  sample_summary <- outlier_log_df %>%
    group_by(sample_id, species, log2_lib_size) %>%
    summarize(n_contrasts_removed = n(),
              cell_types = paste(unique(cell_type), collapse = "; "),
              .groups = "drop") %>%
    arrange(species, log2_lib_size)

  message("\n=== Outlier samples removed (", nrow(outlier_log_df),
          " sample×contrast removals across ",
          length(unique(outlier_log_df$sample_id)), " unique samples) ===")
  for (i in seq_len(nrow(sample_summary))) {
    r <- sample_summary[i, ]
    message(sprintf("  %-45s [%s]  log2-lib=%.2f  removed from %d contrast(s)",
                    r$sample_id, r$species, r$log2_lib_size,
                    r$n_contrasts_removed))
  }
  message("Full log saved: ", basename(out_file))
  invisible(sample_summary)
}


#' Run DESeq2 along evolutionary branches with orthology-aware peak subsetting
#'
#' For each cell type and each contrast, subsets to peaks where ALL involved
#' species have physical orthologous sequence, then runs DESeq2.
#'
#' @param counts_union   Union-joined count matrix (peaks x samples)
#' @param meta_shared    Metadata with cell_type, species columns
#' @param contrasts      Named list from define_evolutionary_contrasts()
#' @param ortho_mat      Logical matrix from build_orthology_index()
#' @param out_dir        Output directory
#' @param min_samples    Min samples with signal for a peak (default 2)
#' @return Named list: branch_res[[cell_type]][[contrast_name]] = DESeq2 results
run_deseq2_evolutionary <- function(counts_union, meta_shared, contrasts,
                                    ortho_mat, out_dir, min_samples = 2,
                                    filter_outliers = TRUE, sd_thresh = 2.5) {
  suppressPackageStartupMessages({
    library(DESeq2)
    library(dplyr)
  })

  dir.create(file.path(out_dir, "_summary"), showWarnings = FALSE, recursive = TRUE)

  branch_res    <- list()
  outlier_log_df <- data.frame(cell_type = character(), contrast = character(),
                               sample_id = character(), species = character(),
                               log2_lib_size = numeric(),
                               stringsAsFactors = FALSE)
  cell_types <- levels(meta_shared$cell_type)

  valid_samples <- intersect(rownames(meta_shared), colnames(counts_union))
  counts_union  <- counts_union[, valid_samples, drop = FALSE]

  for (ct in cell_types) {
    message(sprintf("\n=== [Evolutionary Branches] %s ===", ct))
    branch_res[[ct]] <- list()

    meta_ct    <- meta_shared[meta_shared$cell_type == ct, ]
    ct_samples <- intersect(rownames(meta_ct), colnames(counts_union))
    meta_ct    <- meta_ct[ct_samples, ]

    for (node_name in names(contrasts)) {
      pos_sp <- contrasts[[node_name]]$pos
      neg_sp <- contrasts[[node_name]]$neg
      req_sp <- c(pos_sp, neg_sp)

      # Check we have samples from all required species
      avail_sp <- unique(as.character(meta_ct$species))
      missing  <- setdiff(req_sp, avail_sp)
      if (length(missing) > 0) next

      # Subset samples to only the species in this contrast
      node_samples <- ct_samples[meta_ct$species %in% req_sp]
      if (length(node_samples) < 4) next

      node_meta   <- meta_ct[node_samples, ]
      node_meta$species <- factor(node_meta$species)

      # Orthology-aware peak subsetting: only peaks present in ALL contrast species
      valid_peaks <- rownames(ortho_mat)[
        rowSums(ortho_mat[, req_sp, drop = FALSE]) == length(req_sp)
      ]
      valid_peaks <- intersect(valid_peaks, rownames(counts_union))

      node_counts <- counts_union[valid_peaks, node_samples, drop = FALSE]

      # Per-contrast outlier filtering
      if (filter_outliers) {
        outlier_flags <- detect_per_contrast_outliers(node_counts, node_meta,
                                                      sd_thresh = sd_thresh)
        removed_ids <- names(outlier_flags)[outlier_flags]
        if (length(removed_ids) > 0) {
          lib_vec <- log2(colSums(node_counts) + 1)
          outlier_log_df <- append_outlier_log(outlier_log_df, ct, node_name,
                                               removed_ids, node_meta, lib_vec)
          keep_samp   <- node_samples[!outlier_flags[node_samples]]
          node_meta   <- node_meta[keep_samp, ]
          node_counts <- node_counts[, keep_samp, drop = FALSE]
          node_meta$species <- droplevels(node_meta$species)
          message(sprintf("    Outliers removed: %s", paste(removed_ids, collapse = ", ")))
        }
        if (length(node_samples[!outlier_flags[node_samples]]) < 4) next
      }

      # Active peak filter
      keep <- rowSums(node_counts >= 10) >= min_samples
      node_counts <- node_counts[keep, , drop = FALSE]

      if (nrow(node_counts) < 50 || length(unique(node_meta$species)) < 2) next

      # Run DESeq2
      tryCatch({
        dds <- DESeqDataSetFromMatrix(
          countData = round(node_counts),
          colData   = node_meta,
          design    = ~ 0 + species
        )
        dds <- estimateSizeFactors(dds, type = "poscounts")
        dds <- DESeq(dds, quiet = TRUE)

        res_names    <- resultsNames(dds)
        contrast_vec <- make_contrast_vector(pos_sp, neg_sp, res_names)

        res            <- results(dds, contrast = contrast_vec, alpha = 0.05)
        res_ordered    <- res[order(res$padj), ]
        branch_res[[ct]][[node_name]] <- res_ordered

        # Save CSV
        ct_dir <- file.path(out_dir, ct, "atac_evolutionary")
        dir.create(ct_dir, showWarnings = FALSE, recursive = TRUE)
        res_df         <- as.data.frame(res_ordered)
        res_df$peak_id <- rownames(res_df)
        res_df <- res_df[, c("peak_id", setdiff(names(res_df), "peak_id"))]
        write.csv(res_df, file.path(ct_dir, paste0(node_name, "_DA.csv")),
                  row.names = FALSE)

        sig_up <- sum(res$padj < 0.05 & res$log2FoldChange > 1, na.rm = TRUE)
        message(sprintf("  %-40s: %d peaks tested, %d sig UP",
                        node_name, nrow(node_counts), sig_up))
      }, error = function(e) {
        message(sprintf("  %-40s: FAILED (%s)", node_name, e$message))
      })
    }
  }

  saveRDS(branch_res, file.path(out_dir, "_summary", "branch_res_list.rds"))
  save_outlier_summary(outlier_log_df, out_dir)
  message("\nCheckpoint: evolutionary branch results saved.")
  return(branch_res)
}


#' Ultra-robust filtering: min(positive donors) > max(negative donors)
#'
#' Applies stringent sample-level filter on DESeq2-significant peaks.
#' A peak survives only if its LOWEST log2CPM in the positive-side donors
#' exceeds its HIGHEST log2CPM in the negative-side donors.
#'
#' @param branch_res  Output from run_deseq2_evolutionary()
#' @param counts_union Union count matrix
#' @param meta_shared  Metadata
#' @param contrasts    Named list from define_evolutionary_contrasts()
#' @param out_dir      Output directory
#' @param padj_thresh  Adjusted p-value threshold (default 0.05)
#' @param lfc_thresh   Log2FC threshold (default 1)
#' @param min_logcpm   Minimum log2CPM in positive donors (default 1)
#' @return Named list: robust_peaks[[ct]][[node]] = character vector of peak IDs
ultra_robust_filter <- function(branch_res, counts_union, meta_shared,
                                contrasts, out_dir,
                                padj_thresh = 0.05, lfc_thresh = 1,
                                min_logcpm = 1) {
  suppressPackageStartupMessages(library(matrixStats))

  dir.create(file.path(out_dir, "_summary"), showWarnings = FALSE, recursive = TRUE)

  valid_samples <- intersect(rownames(meta_shared), colnames(counts_union))
  counts_union  <- counts_union[, valid_samples, drop = FALSE]

  robust_peaks <- list()

  for (ct in names(branch_res)) {
    message("\nUltra-robust filtering: ", ct)
    robust_peaks[[ct]] <- list()

    meta_ct    <- meta_shared[meta_shared$cell_type == ct, ]
    ct_samples <- intersect(rownames(meta_ct), colnames(counts_union))
    meta_ct    <- meta_ct[ct_samples, ]

    # Build log2CPM using full library sizes for this cell type
    ct_counts  <- counts_union[, ct_samples, drop = FALSE]
    lib_sizes  <- colSums(ct_counts)
    cpm_mat    <- t(t(ct_counts) / lib_sizes) * 1e6
    logcpm_mat <- log2(cpm_mat + 1)

    for (node_name in names(branch_res[[ct]])) {
      res <- branch_res[[ct]][[node_name]]
      if (is.null(res)) next

      # Get significant UP peaks
      res_df <- as.data.frame(res)
      sig_up <- rownames(res_df)[!is.na(res_df$padj) &
                                   res_df$padj < padj_thresh &
                                   res_df$log2FoldChange > lfc_thresh]
      sig_up <- intersect(sig_up, rownames(logcpm_mat))
      if (length(sig_up) == 0) next

      pos_sp <- contrasts[[node_name]]$pos
      neg_sp <- contrasts[[node_name]]$neg

      pos_samples <- ct_samples[meta_ct$species %in% pos_sp]
      neg_samples <- ct_samples[meta_ct$species %in% neg_sp]

      if (length(pos_samples) < 1 || length(neg_samples) < 1) next

      pos_mat <- as.matrix(logcpm_mat[sig_up, pos_samples, drop = FALSE])
      neg_mat <- as.matrix(logcpm_mat[sig_up, neg_samples, drop = FALSE])

      min_pos <- rowMins(pos_mat)
      max_neg <- rowMaxs(neg_mat)

      is_robust <- (min_pos > max_neg) & (min_pos > min_logcpm)
      is_robust[is.na(is_robust)] <- FALSE

      surviving <- sig_up[is_robust]
      robust_peaks[[ct]][[node_name]] <- surviving

      message(sprintf("  [%s] DESeq2 UP: %d → Ultra-Robust: %d",
                      node_name, length(sig_up), length(surviving)))

      # Save BED + CSV
      if (length(surviving) > 0) {
        ct_node_dir <- file.path(out_dir, ct, "atac_ultra_robust")
        dir.create(ct_node_dir, showWarnings = FALSE, recursive = TRUE)

        out_df <- data.frame(
          peak_id        = surviving,
          min_pos_logCPM = min_pos[is_robust],
          max_neg_logCPM = max_neg[is_robust]
        )
        write.csv(out_df,
                  file.path(ct_node_dir, paste0(node_name, "_DA.csv")),
                  row.names = FALSE)
      }
    }
  }

  saveRDS(robust_peaks, file.path(out_dir, "_summary", "ultra_robust_peaks_list.rds"))
  message("\nUltra-robust filtering complete. Saved checkpoint.")
  return(robust_peaks)
}


# =============================================================================
# SECTION 3d: PER-REGION EVOLUTIONARY BRANCH ANALYSIS
# =============================================================================

#' Run evolutionary branch DESeq2 separately within each anatomical region
#'
#' Uses the same orthology-aware peak subsetting and outlier filtering as
#' run_deseq2_evolutionary(), but restricts samples to one anatomical region
#' at a time (e.g. "Duodenum", "Colon"). Only adult samples are used.
#' A lower per-sample cell threshold is applied since per-region pseudobulks
#' have fewer cells than region-merged ones.
#'
#' Input data must NOT have been region-aggregated (region column preserved).
#'
#' @param counts_union_regions Union count matrix — samples NOT region-aggregated
#' @param meta_regions         Metadata with 'region', 'age', 'cell_type', 'species'
#' @param contrasts            Named list from define_evolutionary_contrasts()
#' @param ortho_mat            Logical orthology matrix from build_orthology_index()
#' @param out_dir              Output base directory
#' @param target_regions       Regions to process (default c("Duodenum", "Colon"))
#' @param min_samples          Min samples with signal per peak (default 2)
#' @param min_cells_region     Min cells per pseudobulk for this analysis (default 50)
#' @param min_counts_region    Min total counts per pseudobulk (default 30000)
#' @param filter_outliers      Apply per-contrast outlier filtering (default TRUE)
#' @param sd_thresh            SD threshold for outlier detection (default 2.5)
#' @return Nested list: region_res[[region]][[cell_type]][[contrast]] = DESeq2 results
run_deseq2_per_region <- function(counts_union_regions, meta_regions,
                                  contrasts, ortho_mat, out_dir,
                                  target_regions    = c("Duodenum", "Colon"),
                                  min_samples       = 2,
                                  min_cells_region  = 50,
                                  min_counts_region = 30000,
                                  filter_outliers   = TRUE,
                                  sd_thresh         = 2.5) {
  suppressPackageStartupMessages({
    library(DESeq2)
    library(dplyr)
  })

  region_res      <- list()
  all_outlier_log <- data.frame(cell_type = character(), contrast = character(),
                                sample_id = character(), species = character(),
                                log2_lib_size = numeric(), region = character(),
                                stringsAsFactors = FALSE)

  # Enforce: adult samples only
  meta_regions <- meta_regions[meta_regions$age == "Adult", ]
  valid_cols   <- intersect(rownames(meta_regions), colnames(counts_union_regions))
  counts_union_regions <- counts_union_regions[, valid_cols, drop = FALSE]
  meta_regions         <- meta_regions[valid_cols, ]

  message("Adult samples available: ", nrow(meta_regions))
  message("Regions in data: ", paste(unique(meta_regions$region), collapse = ", "))

  for (region in target_regions) {
    message(sprintf("\n\n========== REGION: %s ==========", region))

    region_res[[region]] <- list()

    # ---- Subset to this region ----
    reg_rows <- rownames(meta_regions)[meta_regions$region == region]
    if (length(reg_rows) == 0) {
      message("  No samples found for region: ", region, " — skipping.")
      next
    }

    meta_reg   <- meta_regions[reg_rows, ]
    counts_reg <- counts_union_regions[, reg_rows, drop = FALSE]

    # ---- QC filters for region-level pseudobulks ----
    meta_reg$total_counts <- colSums(counts_reg)
    keep_qc <- meta_reg$total_counts >= min_counts_region
    if ("n_cells" %in% colnames(meta_reg)) {
      keep_qc <- keep_qc & (meta_reg$n_cells >= min_cells_region)
    }
    meta_reg   <- meta_reg[keep_qc, ]
    counts_reg <- counts_reg[, rownames(meta_reg), drop = FALSE]
    message(sprintf("  After QC: %d samples (counts >= %d, cells >= %d)",
                    nrow(meta_reg), min_counts_region, min_cells_region))

    # ---- Cell types present in this region ----
    cell_types_reg <- sort(unique(as.character(meta_reg$cell_type)))

    for (ct in cell_types_reg) {
      meta_ct    <- meta_reg[meta_reg$cell_type == ct, ]
      ct_samples <- rownames(meta_ct)
      counts_ct  <- counts_reg[, ct_samples, drop = FALSE]

      region_res[[region]][[ct]] <- list()

      for (node_name in names(contrasts)) {
        pos_sp <- contrasts[[node_name]]$pos
        neg_sp <- contrasts[[node_name]]$neg
        req_sp <- c(pos_sp, neg_sp)

        avail_sp <- unique(as.character(meta_ct$species))
        if (length(setdiff(req_sp, avail_sp)) > 0) next

        node_samples <- ct_samples[meta_ct$species %in% req_sp]
        if (length(node_samples) < 4) next

        node_meta <- meta_ct[node_samples, ]
        node_meta$species <- factor(node_meta$species)

        # Orthology-aware peak subsetting
        valid_peaks <- rownames(ortho_mat)[
          rowSums(ortho_mat[, req_sp, drop = FALSE]) == length(req_sp)
        ]
        valid_peaks <- intersect(valid_peaks, rownames(counts_ct))
        node_counts <- counts_ct[valid_peaks, node_samples, drop = FALSE]

        # Per-contrast outlier filtering
        if (filter_outliers) {
          outlier_flags <- detect_per_contrast_outliers(node_counts, node_meta,
                                                        sd_thresh = sd_thresh)
          removed_ids <- names(outlier_flags)[outlier_flags]
          if (length(removed_ids) > 0) {
            lib_vec  <- log2(colSums(node_counts) + 1)
            new_rows <- append_outlier_log(data.frame(), ct, node_name,
                                           removed_ids, node_meta, lib_vec)
            new_rows$region    <- region
            all_outlier_log    <- rbind(all_outlier_log, new_rows)
            keep_samp   <- node_samples[!outlier_flags[node_samples]]
            node_meta   <- node_meta[keep_samp, ]
            node_counts <- node_counts[, keep_samp, drop = FALSE]
            node_meta$species <- droplevels(node_meta$species)
          }
          if (ncol(node_counts) < 4) next
        }

        # Active peak filter
        keep <- rowSums(node_counts >= 10) >= min_samples
        node_counts <- node_counts[keep, , drop = FALSE]
        if (nrow(node_counts) < 50) next

        tryCatch({
          dds <- DESeqDataSetFromMatrix(
            countData = round(node_counts),
            colData   = node_meta,
            design    = ~ 0 + species
          )
          dds <- estimateSizeFactors(dds, type = "poscounts")
          dds <- DESeq(dds, quiet = TRUE)

          res_names    <- resultsNames(dds)
          contrast_vec <- make_contrast_vector(pos_sp, neg_sp, res_names)
          res          <- results(dds, contrast = contrast_vec, alpha = 0.05)
          res_ordered  <- res[order(res$padj), ]

          region_res[[region]][[ct]][[node_name]] <- res_ordered

          # Save CSV
          ct_dir <- file.path(out_dir, ct, "atac_region")
          dir.create(ct_dir, showWarnings = FALSE, recursive = TRUE)
          res_df         <- as.data.frame(res_ordered)
          res_df$peak_id <- rownames(res_df)
          res_df <- res_df[, c("peak_id", setdiff(names(res_df), "peak_id"))]
          write.csv(res_df, file.path(ct_dir, paste0(region, "_", node_name, "_DA.csv")),
                    row.names = FALSE)

          sig_up <- sum(res$padj < 0.05 & res$log2FoldChange > 1, na.rm = TRUE)
          message(sprintf("  [%s | %s] %-40s: %d peaks, %d sig UP",
                          region, ct, node_name, nrow(node_counts), sig_up))
        }, error = function(e) {
          message(sprintf("  [%s | %s] %-40s: FAILED (%s)",
                          region, ct, node_name, e$message))
        })
      }
    }
  }

  # ---- Save outputs ----
  saveRDS(region_res, file.path(out_dir, "_summary", "region_res_list.rds"))

  if (nrow(all_outlier_log) > 0) {
    write.csv(all_outlier_log,
              file.path(out_dir, "_summary", "outlier_region_DA.csv"),
              row.names = FALSE)
    message("\nRegion outlier log: ", nrow(all_outlier_log), " removals saved.")
  }

  message("\nPer-region DESeq2 complete.")
  return(region_res)
}


#' Build a region-comparison summary table for one contrast and cell type
#'
#' For a given contrast, reports which peaks are significant in Duodenum only,
#' Colon only, or both, with their log2FC values in each region.
#'
#' @param region_res     Output from run_deseq2_per_region()
#' @param cell_type      Cell type to query
#' @param contrast       Contrast name
#' @param padj_thresh    Significance threshold (default 0.05)
#' @param lfc_thresh     Log2FC threshold (default 1)
#' @return data.frame with columns: peak_id, sig_duodenum, sig_colon,
#'         lfc_duodenum, lfc_colon, category
compare_regions <- function(region_res, cell_type, contrast,
                            padj_thresh = 0.05, lfc_thresh = 1) {
  regions <- names(region_res)
  res_list <- lapply(regions, function(rg) {
    r <- region_res[[rg]][[cell_type]][[contrast]]
    if (is.null(r)) return(NULL)
    df <- as.data.frame(r)
    df$peak_id <- rownames(df)
    df
  })
  names(res_list) <- regions

  res_list <- Filter(Negate(is.null), res_list)
  if (length(res_list) < 2) {
    message("Need results from at least 2 regions.")
    return(invisible(NULL))
  }

  all_peaks <- Reduce(union, lapply(res_list, function(x) x$peak_id))
  out_df    <- data.frame(peak_id = all_peaks, stringsAsFactors = FALSE)

  for (rg in names(res_list)) {
    df <- res_list[[rg]]
    rownames(df) <- df$peak_id
    lfc_col  <- paste0("lfc_", tolower(rg))
    padj_col <- paste0("padj_", tolower(rg))
    sig_col  <- paste0("sig_", tolower(rg))

    out_df[[lfc_col]]  <- df[out_df$peak_id, "log2FoldChange"]
    out_df[[padj_col]] <- df[out_df$peak_id, "padj"]
    out_df[[sig_col]]  <- !is.na(out_df[[padj_col]]) &
                          out_df[[padj_col]] < padj_thresh &
                          !is.na(out_df[[lfc_col]]) &
                          out_df[[lfc_col]] > lfc_thresh
  }

  rg_names  <- names(res_list)
  sig_cols  <- paste0("sig_", tolower(rg_names))
  sig_mat   <- as.matrix(out_df[, sig_cols, drop = FALSE])
  out_df$n_regions_sig <- rowSums(sig_mat, na.rm = TRUE)
  out_df$category <- apply(sig_mat, 1, function(x) {
    which_sig <- rg_names[which(x)]
    if (length(which_sig) == 0) "Neither"
    else if (length(which_sig) == length(rg_names)) "Both"
    else paste0(which_sig, "_only")
  })

  out_df <- out_df[out_df$n_regions_sig > 0, ]
  out_df <- out_df[order(-out_df$n_regions_sig), ]
  return(out_df)
}


#' Plot a region-comparison heatmap (signed p-value per region) for one cell type
#'
#' Shows peaks significant in either region, colored by signed -log10(padj).
#' Columns are region × contrast combinations; rows are peaks.
#'
#' @param region_res   Output from run_deseq2_per_region()
#' @param cell_type    Cell type to plot
#' @param out_file     Output PDF path
#' @param padj_thresh  Significance threshold
#' @param lfc_thresh   Log2FC threshold
plot_region_comparison_heatmap <- function(region_res, cell_type, out_file,
                                           padj_thresh = 0.05, lfc_thresh = 1,
                                           max_peaks = 2000) {
  suppressPackageStartupMessages({
    library(ComplexHeatmap)
    library(circlize)
  })

  regions    <- names(region_res)
  contrasts  <- unique(unlist(lapply(region_res, function(rg) names(rg[[cell_type]]))))
  contrasts  <- contrasts[!is.null(contrasts)]

  # Build signed-p matrix: rows = peaks, columns = region_contrast
  sig_peaks_all <- c()
  col_data      <- list()

  for (rg in regions) {
    for (cn in contrasts) {
      r <- region_res[[rg]][[cell_type]][[cn]]
      if (is.null(r)) next
      df <- as.data.frame(r)
      sig_idx <- which(!is.na(df$padj) & df$padj < padj_thresh &
                         df$log2FoldChange > lfc_thresh)
      if (length(sig_idx) > 0)
        sig_peaks_all <- c(sig_peaks_all, rownames(df)[sig_idx])
      padj_safe <- df$padj; padj_safe[is.na(padj_safe)] <- 1
      padj_safe[padj_safe == 0] <- 1e-300
      col_vals <- df$log2FoldChange * -log10(padj_safe)
      col_data[[paste0(rg, "_", cn)]] <- setNames(col_vals, rownames(df))
    }
  }

  sig_peaks_all <- unique(sig_peaks_all)
  if (length(sig_peaks_all) == 0) {
    message("No significant peaks for region comparison: ", cell_type)
    return(invisible(NULL))
  }

  # Subsample to max_peaks by highest absolute signed p-value sum across columns
  if (length(sig_peaks_all) > max_peaks) {
    message("  Subsampling ", length(sig_peaks_all), " → ", max_peaks,
            " peaks for region comparison heatmap")
    all_vals <- lapply(col_data, function(v) abs(v[sig_peaks_all]))
    score_mat <- do.call(cbind, all_vals)
    row_scores <- rowSums(score_mat, na.rm = TRUE)
    sig_peaks_all <- sig_peaks_all[order(row_scores, decreasing = TRUE)[seq_len(max_peaks)]]
  }

  mat <- matrix(0, nrow = length(sig_peaks_all), ncol = length(col_data),
                dimnames = list(sig_peaks_all, names(col_data)))
  for (cn in names(col_data)) {
    common <- intersect(sig_peaks_all, names(col_data[[cn]]))
    mat[common, cn] <- col_data[[cn]][common]
  }

  cap_val <- quantile(abs(mat), 0.99, na.rm = TRUE)
  mat[mat >  cap_val] <-  cap_val
  mat[mat < -cap_val] <- -cap_val

  col_fun <- colorRamp2(c(-cap_val, 0, cap_val),
                        c("#377eb8", "white", "#e41a1c"))

  ht <- Heatmap(mat,
                name             = "Signed\n-log10(padj)",
                col              = col_fun,
                cluster_rows     = TRUE,
                cluster_columns  = FALSE,
                show_row_names   = FALSE,
                column_names_rot = 45,
                column_names_gp  = gpar(fontsize = 8),
                column_title     = paste("Region Comparison:", cell_type),
                column_title_gp  = gpar(fontsize = 13, fontface = "bold"),
                border           = TRUE)

  dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
  pdf(out_file, width = max(10, length(col_data) * 0.4 + 4),
      height = max(8, length(sig_peaks_all) * 0.01 + 4))
  draw(ht, merge_legend = TRUE)
  dev.off()
  message("Region comparison heatmap saved: ", basename(out_file))
}


# =============================================================================
# SECTION 4: PCA
# =============================================================================

#' Run global VST and generate PCA plots for multiple colorings
#'
#' @param counts    Count matrix
#' @param meta      Metadata data.frame
#' @param out_dir   Output directory for plots
#' @param prefix    Filename prefix (e.g. "Shared_Peaks" or "Enterocytes")
#' @param subset_ct Optional: character vector of cell types to plot (NULL = all)
#' @return VST DESeqTransform object (for downstream use)
run_pca_suite <- function(counts, meta, out_dir, prefix = "Global",
                          subset_ct = NULL) {
  suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
  })

  plots_dir <- file.path(out_dir, "plots")
  dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

  # Run VST — dynamically build design formula based on available factors
  design_terms <- c()
  if (length(unique(meta$species)) > 1)   design_terms <- c(design_terms, "species")
  if (length(unique(meta$cell_type)) > 1) design_terms <- c(design_terms, "cell_type")
  if (length(design_terms) == 0) design_terms <- "1"  # intercept-only fallback
  design_formula <- as.formula(paste("~", paste(design_terms, collapse = " + ")))
  message("  PCA design: ", deparse(design_formula))

  dds <- DESeqDataSetFromMatrix(
    countData = round(counts),
    colData   = meta,
    design    = design_formula
  )
  dds <- estimateSizeFactors(dds, type = "poscounts")
  vsd <- vst(dds, blind = FALSE)

  # Extract PCA coordinates
  intgroup <- intersect(c("cell_type", "species", "donor", "region", "age"),
                        colnames(meta))
  pca_data   <- plotPCA(vsd, intgroup = intgroup, returnData = TRUE)
  percentVar <- round(100 * attr(pca_data, "percentVar"))

  # Optionally filter to a subset of cell types
  if (!is.null(subset_ct)) {
    pca_data <- pca_data[pca_data$cell_type %in% subset_ct, ]
  }

  # Internal plot builder
  .build_pca <- function(df, color_var, shape_var = NULL, title,
                         add_labels = FALSE) {
    p <- ggplot(df, aes_string(x = "PC1", y = "PC2",
                               color = color_var, shape = shape_var)) +
      geom_point(size = 4, alpha = 0.85) +
      xlab(paste0("PC1: ", percentVar[1], "% variance")) +
      ylab(paste0("PC2: ", percentVar[2], "% variance")) +
      theme_bw() +
      ggtitle(title) +
      theme(legend.position = "right",
            legend.title = element_text(face = "bold"),
            plot.title = element_text(face = "bold", size = 14))

    if (add_labels && "donor" %in% colnames(df)) {
      p <- p + ggrepel::geom_text_repel(aes(label = donor), size = 3,
                                         show.legend = FALSE,
                                         max.overlaps = 20)
    }

    if (!is.null(shape_var)) {
      n_shapes <- length(unique(df[[shape_var]]))
      p <- p + scale_shape_manual(
        values = c(15, 16, 17, 18, 3, 4, 8, 9)[1:n_shapes])
    }
    return(p)
  }

  # Generate standard PCA panels
  color_vars <- intersect(c("species", "cell_type", "region", "age", "donor"),
                          colnames(pca_data))

  for (cv in color_vars) {
    title <- paste(prefix, "PCA:", tools::toTitleCase(cv))
    p     <- .build_pca(pca_data, color_var = cv, title = title)
    fname <- file.path(plots_dir,
                       paste0(prefix, "_PCA_", cv, ".pdf"))
    ggsave(fname, p, width = 10, height = 7)
    message("  Saved: ", basename(fname))
  }

  # Combined species + cell_type
  if (all(c("species", "cell_type") %in% colnames(pca_data))) {
    p_comb <- .build_pca(pca_data, color_var = "species",
                         shape_var = "cell_type",
                         title = paste(prefix, "PCA: Species & Cell Type"),
                         add_labels = TRUE)
    ggsave(file.path(plots_dir, paste0(prefix, "_PCA_Combined.pdf")),
           p_comb, width = 12, height = 8)
    message("  Saved: ", prefix, "_PCA_Combined.pdf")
  }

  message("PCA suite complete for: ", prefix)
  return(vsd)
}


# =============================================================================
# SECTION 5: VOLCANO PLOTS
# =============================================================================

#' Generate a volcano plot with labeled top peaks
#'
#' @param res_df       DESeq2 results data.frame (must have log2FoldChange, padj)
#' @param title        Plot title
#' @param out_file     Output PDF path
#' @param n_label      Number of top peaks to label on each side
#' @param lfc_thresh   Log2 fold-change threshold for significance
#' @param padj_thresh  Adjusted p-value threshold
#' @return ggplot object (invisibly)
plot_volcano <- function(res_df, title, out_file,
                         n_label = 10, lfc_thresh = 1, padj_thresh = 0.05) {
  suppressPackageStartupMessages({
    library(ggplot2)
    library(ggrepel)
  })

  df <- as.data.frame(res_df)
  df$peak_id <- rownames(df)
  df <- df[!is.na(df$padj), ]

  # Classify significance

  df$Significance <- "Not Significant"
  df$Significance[df$padj < padj_thresh & df$log2FoldChange >  lfc_thresh] <- "Up"
  df$Significance[df$padj < padj_thresh & df$log2FoldChange < -lfc_thresh] <- "Down"
  df$Significance <- factor(df$Significance,
                            levels = c("Up", "Down", "Not Significant"))

  # Identify top peaks to label (by combined score: |LFC| * -log10(padj))
  df$rank_score <- abs(df$log2FoldChange) * -log10(df$padj)

  up_df   <- df[df$Significance == "Up", ]
  down_df <- df[df$Significance == "Down", ]

  top_up   <- head(up_df[order(up_df$rank_score, decreasing = TRUE), ], n_label)
  top_down <- head(down_df[order(down_df$rank_score, decreasing = TRUE), ], n_label)

  df$label <- ""
  df$label[df$peak_id %in% top_up$peak_id]   <- df$peak_id[df$peak_id %in% top_up$peak_id]
  df$label[df$peak_id %in% top_down$peak_id] <- df$peak_id[df$peak_id %in% top_down$peak_id]

  volcano_colors <- c("Up" = "#e41a1c", "Down" = "#377eb8",
                       "Not Significant" = "grey80")

  p <- ggplot(df, aes(x = log2FoldChange, y = -log10(padj),
                      color = Significance)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = volcano_colors) +
    geom_vline(xintercept = c(-lfc_thresh, lfc_thresh),
               linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(padj_thresh),
               linetype = "dashed", color = "black") +
    geom_text_repel(
      data = df[df$label != "", ],
      aes(label = label),
      size = 2.5, max.overlaps = 30,
      fontface = "italic",
      segment.color = "grey40",
      segment.alpha = 0.6,
      show.legend = FALSE
    ) +
    theme_bw() +
    labs(title = title,
         subtitle = sprintf("Up: %d | Down: %d",
                            sum(df$Significance == "Up"),
                            sum(df$Significance == "Down")),
         x = "Log2 Fold Change",
         y = "-Log10(Adjusted P-Value)") +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          plot.title = element_text(face = "bold", size = 14))

  dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
  ggsave(out_file, p, width = 8, height = 8)
  message("  Volcano saved: ", basename(out_file))
  invisible(p)
}


# =============================================================================
# SECTION 5b: ACCESSIBILITY DOTPLOT
# =============================================================================

#' Plot per-species accessibility for a set of peaks as a species-level dotplot
#'
#' For each peak, shows mean log2CPM per species (dot size = n donors, color =
#' mean log2CPM). Peaks are clustered by their species accessibility profile.
#' Useful for inspecting top hits from differential analysis.
#'
#' @param peak_ids     Character vector of peak IDs to plot
#' @param counts       Count matrix (peaks x samples)
#' @param meta         Metadata with 'species' and 'donor' columns
#' @param out_file     Output PDF path
#' @param title        Plot title
#' @param species_order Species order for x-axis
#' @param max_peaks    Cap for number of peaks to show (default 60)
plot_accessibility_dotplot <- function(peak_ids, counts, meta, out_file,
                                       title = "Accessibility per Species",
                                       species_order = c("Human", "Chimpanzee",
                                                          "Bonobo", "Gorilla",
                                                          "Macaque", "Marmoset"),
                                       max_peaks = 60) {
  suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
  })

  peak_ids <- intersect(peak_ids, rownames(counts))
  if (length(peak_ids) == 0) {
    message("No peaks found in counts matrix for dotplot.")
    return(invisible(NULL))
  }
  if (length(peak_ids) > max_peaks) {
    peak_ids <- peak_ids[1:max_peaks]
    message("  Dotplot: capping at ", max_peaks, " peaks")
  }

  # Compute log2CPM
  lib_sizes  <- colSums(counts)
  cpm_mat    <- t(t(counts[peak_ids, , drop = FALSE]) / lib_sizes) * 1e6
  logcpm_mat <- log2(cpm_mat + 1)

  # Build long data frame
  long_df <- as.data.frame(logcpm_mat) %>%
    tibble::rownames_to_column("peak_id") %>%
    tidyr::pivot_longer(-peak_id, names_to = "sample_id", values_to = "logcpm") %>%
    dplyr::left_join(data.frame(sample_id = rownames(meta),
                                species   = as.character(meta$species),
                                stringsAsFactors = FALSE),
                     by = "sample_id")

  # Summarise per peak × species
  sum_df <- long_df %>%
    group_by(peak_id, species) %>%
    summarize(mean_logcpm = mean(logcpm, na.rm = TRUE),
              n_donors    = n(),
              .groups = "drop")

  # Order species
  present_sp  <- intersect(species_order, unique(sum_df$species))
  sum_df$species <- factor(sum_df$species, levels = present_sp)

  # Cluster peaks by their mean-logcpm vector
  wide_mat <- sum_df %>%
    pivot_wider(id_cols = peak_id, names_from = species,
                values_from = mean_logcpm, values_fill = 0) %>%
    tibble::column_to_rownames("peak_id")
  if (nrow(wide_mat) > 1) {
    peak_order <- rownames(wide_mat)[
      hclust(dist(wide_mat))$order
    ]
  } else {
    peak_order <- rownames(wide_mat)
  }
  sum_df$peak_id <- factor(sum_df$peak_id, levels = peak_order)

  p <- ggplot(sum_df, aes(x = species, y = peak_id,
                          size = n_donors, color = mean_logcpm)) +
    geom_point(alpha = 0.85) +
    scale_color_gradient2(low  = "#2166ac", mid = "#f7f7f7", high = "#d73027",
                          midpoint = median(sum_df$mean_logcpm, na.rm = TRUE),
                          name = "Mean\nlog2CPM") +
    scale_size_continuous(range = c(1.5, 6), name = "N donors") +
    theme_bw() +
    theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 11),
          axis.text.y  = element_text(size = 5),
          panel.grid.major = element_line(color = "grey92"),
          plot.title   = element_text(face = "bold", size = 13)) +
    labs(title = title, x = NULL, y = "Peak ID")

  dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
  dynamic_h <- max(6, length(peak_order) * 0.12 + 2)
  ggsave(out_file, p, width = 8, height = dynamic_h, limitsize = FALSE)
  message("  Accessibility dotplot saved: ", basename(out_file))
  invisible(p)
}


#' Strip-plot of per-donor log2CPM for a single peak, grouped by species
#'
#' Shows raw data points (one per donor) plus mean ± SD. Useful for inspecting
#' individual top hits and confirming they pass the ultra-robust filter visually.
#'
#' @param peak_id   Single peak ID
#' @param counts    Count matrix
#' @param meta      Metadata with 'species' and 'donor'
#' @param out_file  Output PDF path (or NULL to return plot)
#' @param title     Plot title (default = peak_id)
#' @param species_order Species order for x-axis
plot_peak_strip <- function(peak_id, counts, meta, out_file = NULL,
                             title = NULL,
                             species_order = c("Human", "Chimpanzee", "Bonobo",
                                               "Gorilla", "Macaque", "Marmoset")) {
  suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
  })

  if (!peak_id %in% rownames(counts)) {
    message("Peak not found: ", peak_id)
    return(invisible(NULL))
  }

  lib_sizes <- colSums(counts)
  cpm       <- counts[peak_id, ] / lib_sizes * 1e6
  logcpm    <- log2(cpm + 1)

  df <- data.frame(sample_id = names(logcpm),
                   logcpm    = as.numeric(logcpm),
                   species   = as.character(meta[names(logcpm), "species"]),
                   stringsAsFactors = FALSE)
  present_sp <- intersect(species_order, unique(df$species))
  df$species <- factor(df$species, levels = present_sp)

  # Summarise for mean ± SD overlay
  sum_df <- df %>%
    group_by(species) %>%
    summarize(mean_val = mean(logcpm), sd_val = sd(logcpm), .groups = "drop")

  p <- ggplot(df, aes(x = species, y = logcpm, color = species)) +
    geom_jitter(width = 0.15, size = 2.5, alpha = 0.8) +
    geom_crossbar(data = sum_df,
                  aes(x = species, y = mean_val,
                      ymin = mean_val - sd_val, ymax = mean_val + sd_val),
                  width = 0.4, color = "black", fill = NA, fatten = 2) +
    theme_bw() +
    theme(legend.position  = "none",
          axis.text.x      = element_text(angle = 45, hjust = 1, size = 11),
          plot.title       = element_text(face = "bold", size = 12)) +
    labs(title = if (is.null(title)) peak_id else title,
         x = NULL, y = "log2(CPM + 1)")

  if (!is.null(out_file)) {
    dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
    ggsave(out_file, p, width = 6, height = 4)
    message("  Strip plot saved: ", basename(out_file))
  }
  invisible(p)
}


# =============================================================================
# SECTION 6: HEATMAPS
# =============================================================================

#' Species-specific marker heatmap for a single cell type
#'
#' @param res_ct     Named list of DESeq2 results for one cell type
#' @param vsd        VST object from global normalization
#' @param meta       Metadata data.frame
#' @param cell_type  Cell type to plot
#' @param out_file   Output PDF path
#' @param top_n      Number of top peaks per species
plot_species_heatmap <- function(res_ct, vsd, meta, cell_type, out_file,
                                 top_n = 25) {
  suppressPackageStartupMessages({
    library(pheatmap)
    library(viridis)
  })

  top_peaks_list <- list()
  for (sp in names(res_ct)) {
    res_sp  <- as.data.frame(res_ct[[sp]])
    # Filter significant UP peaks
    sig_res <- res_sp[!is.na(res_sp$padj) &
                        res_sp$padj < 0.05 &
                        res_sp$log2FoldChange > 1, ]
    if (nrow(sig_res) == 0) next

    # Rank by combined score: |LFC| * -log10(padj)
    sig_res$rank_score <- abs(sig_res$log2FoldChange) *
                          -log10(pmax(sig_res$padj, 1e-300))
    sig_res <- sig_res[order(sig_res$rank_score, decreasing = TRUE), ]
    top_peaks_list[[sp]] <- rownames(head(sig_res, top_n))

    message(sprintf("  %s: %d sig peaks, top %d selected",
                    sp, nrow(sig_res), length(top_peaks_list[[sp]])))
  }

  top_regions <- unique(unlist(top_peaks_list))
  message("  Total unique top regions: ", length(top_regions))

  if (length(top_regions) < 3) {
    message("  Skipping heatmap: too few significant markers for ", cell_type)
    return(invisible(NULL))
  }

  meta_ct <- meta[meta$cell_type == cell_type, ]
  meta_ct <- meta_ct[order(meta_ct$species), ]
  samples <- rownames(meta_ct)

  # Intersect with available VST features
  top_regions <- intersect(top_regions, rownames(vsd))
  message("  After VST intersect: ", length(top_regions), " regions")
  if (length(top_regions) < 3) {
    message("  Skipping: peaks not in VST object for ", cell_type)
    return(invisible(NULL))
  }

  mat        <- SummarizedExperiment::assay(vsd)[top_regions, samples]
  mat_scaled <- t(scale(t(mat)))

  # Replace NA/NaN from scaling (zero-variance rows)
  mat_scaled[is.nan(mat_scaled)] <- 0
  mat_scaled[is.na(mat_scaled)]  <- 0

  anno_col <- data.frame(Species = meta_ct$species, row.names = samples)

  dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
  pdf(out_file, width = 9, height = 7)
  pheatmap(mat_scaled,
           annotation_col  = anno_col,
           cluster_cols    = FALSE,
           cluster_rows    = TRUE,
           show_rownames   = FALSE,
           show_colnames   = FALSE,
           color           = viridis(100),
           main            = paste("Top Species Markers in", cell_type))
  dev.off()
  message("  Heatmap saved: ", basename(out_file))
}


#' Cross-cell-type signed p-value heatmap (for one target species)
#'
#' @param res_list  Nested list: res_list[[ct]][[species]]
#' @param target_sp Target species (default "Human")
#' @param out_file  Output PDF path
#' @return Matrix used for plotting (for downstream module extraction)
plot_cross_celltype_heatmap <- function(res_list, target_sp = "Human",
                                        out_file) {
  suppressPackageStartupMessages(library(pheatmap))

  signed_p_list  <- list()
  sig_regions_all <- c()

  for (ct in names(res_list)) {
    if (target_sp %in% names(res_list[[ct]])) {
      res_df <- as.data.frame(res_list[[ct]][[target_sp]])

      sig_idx <- which(!is.na(res_df$padj) & res_df$padj < 0.05 &
                         abs(res_df$log2FoldChange) > 1)
      sig_regions_all <- c(sig_regions_all, rownames(res_df)[sig_idx])

      padj_safe <- res_df$padj
      padj_safe[is.na(padj_safe)] <- 1
      padj_safe[padj_safe == 0]   <- 1e-300

      signed_pval <- res_df$log2FoldChange * -log10(padj_safe)

      df <- data.frame(region = rownames(res_df), metric = signed_pval)
      colnames(df)[2] <- ct
      signed_p_list[[ct]] <- df
    }
  }

  sig_regions_all <- unique(sig_regions_all)
  message("Regions significant in >= 1 cell type: ", length(sig_regions_all))

  if (length(sig_regions_all) == 0) {
    message("No significant regions found. Skipping heatmap.")
    return(invisible(NULL))
  }

  merged_df          <- Reduce(function(x, y) merge(x, y, by = "region", all = TRUE),
                               signed_p_list)
  rownames(merged_df) <- merged_df$region
  merged_df$region    <- NULL

  plot_mat <- as.matrix(merged_df[sig_regions_all, ])
  plot_mat[is.na(plot_mat)] <- 0

  # Cap at 99th percentile
  cap_val <- quantile(abs(plot_mat), 0.99, na.rm = TRUE)
  plot_mat[plot_mat >  cap_val] <-  cap_val
  plot_mat[plot_mat < -cap_val] <- -cap_val

  bwr_palette <- colorRampPalette(c("#377eb8", "white", "#e41a1c"))(100)
  max_abs     <- max(abs(plot_mat))
  breaks      <- seq(-max_abs, max_abs, length.out = 101)

  dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
  pdf(out_file, width = 10, height = 8)
  pheatmap(plot_mat,
           color         = bwr_palette,
           breaks        = breaks,
           cluster_cols  = TRUE,
           cluster_rows  = TRUE,
           show_rownames = FALSE,
           show_colnames = TRUE,
           main          = paste0(target_sp, "-Specific Regions Across Cell Types\n",
                                  "(Log2FC x -log10(padj))"))
  dev.off()

  message("Cross-cell-type heatmap saved: ", basename(out_file))
  return(plot_mat)
}


# =============================================================================
# SECTION 7: K-MEANS MODULE EXTRACTION
# =============================================================================

#' Extract K-means modules from a signed-pvalue matrix
#'
#' @param plot_mat     Numeric matrix (regions x cell_types)
#' @param num_modules  Number of K-means clusters
#' @param min_regions  Minimum regions per module to keep
#' @param out_dir      Directory to save module CSVs
#' @param seed         Random seed for reproducibility
#' @return List: km_res (kmeans object), master_df (assignments),
#'         keep_modules (passing module numbers)
extract_kmeans_modules <- function(plot_mat, num_modules = 40,
                                   min_regions = 5, out_dir, seed = 42) {
  suppressPackageStartupMessages({
    library(ComplexHeatmap)
    library(circlize)
  })

  cluster_dir <- file.path(out_dir, "_summary", "cross_celltype_modules")
  dir.create(cluster_dir, showWarnings = FALSE, recursive = TRUE)

  set.seed(seed)
  km_res       <- kmeans(plot_mat, centers = num_modules, iter.max = 100,
                         nstart = 25)
  module_sizes <- table(km_res$cluster)
  keep_modules <- as.numeric(names(module_sizes[module_sizes >= min_regions]))

  message(sprintf("K-means: keeping %d / %d modules (>= %d regions)",
                  length(keep_modules), num_modules, min_regions))

  # Build ordered peak list with behavior annotation
  row_order           <- c()
  ordered_assignments <- c()
  behavior_list       <- list()

  for (mod_num in keep_modules) {
    mod_regions  <- names(km_res$cluster)[km_res$cluster == mod_num]
    center_vals  <- km_res$centers[mod_num, ]
    top_ct       <- names(center_vals)[which.max(center_vals)]
    max_val      <- max(center_vals)
    min_val      <- min(center_vals)

    if (min_val > (max_val * 0.5) && max_val > 0) {
      behavior <- "Universal (Up)"
    } else if (max_val < 0) {
      behavior <- "Universal (Down)"
    } else {
      behavior <- paste("Strong in", top_ct)
    }
    behavior_list[[as.character(mod_num)]] <- behavior

    # Sort internally by max absolute signal
    mod_mat        <- plot_mat[mod_regions, , drop = FALSE]
    internal_order <- order(apply(abs(mod_mat), 1, max), decreasing = TRUE)
    sorted_regions <- mod_regions[internal_order]

    row_order           <- c(row_order, sorted_regions)
    ordered_assignments <- c(ordered_assignments,
                             rep(mod_num, length(sorted_regions)))

    # Save per-module CSV
    write.csv(data.frame(region_id = sorted_regions),
              file.path(cluster_dir, paste0("Module_", mod_num, "_regions.csv")),
              row.names = FALSE)

    message(sprintf("  Module %d: %-25s (%d regions)",
                    mod_num, behavior, length(mod_regions)))
  }

  # Master mapping
  master_df <- data.frame(
    region_id = row_order,
    module    = paste0("Module_", ordered_assignments),
    behavior  = sapply(as.character(ordered_assignments),
                       function(x) behavior_list[[x]])
  )
  write.csv(master_df,
            file.path(cluster_dir, "All_Regions_Module_Assignments.csv"),
            row.names = FALSE)

  # Draw ComplexHeatmap
  plot_mat_ordered <- plot_mat[row_order, , drop = FALSE]
  split_factor     <- factor(paste0("Module ", ordered_assignments),
                             levels = paste0("Module ", keep_modules))

  max_abs  <- max(abs(plot_mat_ordered))
  col_fun  <- colorRamp2(c(-max_abs, 0, max_abs),
                          c("#377eb8", "white", "#e41a1c"))
  lbl_size <- ifelse(length(keep_modules) > 20, 8, 12)

  ht <- Heatmap(plot_mat_ordered,
                name             = "Signed_Pval",
                col              = col_fun,
                row_split        = split_factor,
                row_title_rot    = 0,
                row_title_gp     = gpar(fontsize = lbl_size, fontface = "bold"),
                row_gap          = unit(2, "mm"),
                cluster_rows     = FALSE,
                show_row_names   = FALSE,
                cluster_columns  = TRUE,
                column_names_rot = 45,
                column_names_gp  = gpar(fontsize = 11),
                border           = TRUE,
                column_title     = paste("K-Means Modules (k =",
                                         length(keep_modules), ")"),
                column_title_gp  = gpar(fontsize = 14, fontface = "bold"),
                heatmap_legend_param = list(title = "Signed\n-log10(padj)"))

  heatmap_file <- file.path(out_dir, "plots",
                            "Cross_CellType_ComplexHeatmap_Modules.pdf")
  dynamic_height <- max(9, length(keep_modules) * 0.4)

  pdf(heatmap_file, width = 11, height = dynamic_height)
  draw(ht, merge_legend = TRUE)
  dev.off()

  message("Module ComplexHeatmap saved: ", basename(heatmap_file))
  return(list(km_res       = km_res,
              master_df    = master_df,
              keep_modules = keep_modules))
}


# =============================================================================
# SECTION 8: ENRICHMENT (GO / GREAT)
# =============================================================================

#' Run GO enrichment for a list of gene symbols
#'
#' @param genes     Character vector of gene symbols
#' @param label     Label for this analysis (used in filenames)
#' @param out_dir   Directory to save results
#' @param universe  Character vector of background gene symbols (all tested genes).
#'                  If NULL, defaults to all annotated genes in org.Hs.eg.db.
#' @return enrichResult object (or NULL)
run_go_enrichment <- function(genes, label, out_dir, universe = NULL) {
  suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Hs.eg.db)
  })

  enrich_dir <- file.path(out_dir, "enrichment")
  dir.create(enrich_dir, showWarnings = FALSE, recursive = TRUE)

  genes <- unique(genes[genes != "." & !is.na(genes)])
  if (length(genes) < 5) {
    message("  Too few genes for ", label, " (", length(genes), ")")
    return(invisible(NULL))
  }

  universe <- if (!is.null(universe)) unique(universe[universe != "." & !is.na(universe)]) else NULL

  ego <- enrichGO(gene          = genes,
                  universe      = universe,
                  OrgDb         = org.Hs.eg.db,
                  keyType       = "SYMBOL",
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.2)

  if (!is.null(ego) && nrow(ego) > 0) {
    write.csv(as.data.frame(ego),
              file.path(enrich_dir, paste0(label, "_GO_BP.csv")))
    p <- dotplot(ego, showCategory = 15) + ggplot2::ggtitle(paste(label, "GO"))
    ggplot2::ggsave(file.path(enrich_dir, paste0(label, "_GO_dotplot.pdf")),
                    p, width = 8, height = 7)
    message("  GO enrichment saved for: ", label)
  } else {
    message("  No significant GO terms for: ", label)
  }
  return(invisible(ego))
}


#' Run disease and pathway enrichment (KEGG, Disease Ontology, DisGeNET, Cancer genes)
#'
#' @param genes       Character vector of gene symbols
#' @param label       Label for this analysis (used in filenames)
#' @param out_dir     Directory to save results
#' @param padj_cutoff FDR threshold (default 0.05)
#' @param universe    Character vector of background gene symbols (all tested genes).
#'                    If NULL, defaults to all annotated genes in org.Hs.eg.db.
#' @return Named list of enrichResult objects (kegg, do, dgn, ncg)
run_disease_enrichment <- function(genes, label, out_dir, padj_cutoff = 0.05,
                                   universe = NULL) {
  suppressPackageStartupMessages({
    library(clusterProfiler)
    library(DOSE)
    library(org.Hs.eg.db)
    library(enrichplot)
    library(ggplot2)
  })

  enrich_dir <- file.path(out_dir, "enrichment")
  dir.create(enrich_dir, showWarnings = FALSE, recursive = TRUE)

  genes <- unique(genes[genes != "." & !is.na(genes) & genes != ""])
  if (length(genes) < 5) {
    message("  Too few genes for disease enrichment: ", label, " (", length(genes), ")")
    return(invisible(NULL))
  }

  # Convert SYMBOL → ENTREZID (required for KEGG, DOSE, DisGeNET)
  id_map <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db, drop = TRUE)
  entrez <- unique(id_map$ENTREZID)
  message("  ", label, ": ", length(genes), " symbols → ", length(entrez), " ENTREZ IDs")

  # Convert universe SYMBOL → ENTREZID if provided
  universe_entrez <- NULL
  if (!is.null(universe)) {
    universe <- unique(universe[universe != "." & !is.na(universe) & universe != ""])
    uni_map  <- tryCatch(
      bitr(universe, fromType = "SYMBOL", toType = "ENTREZID",
           OrgDb = org.Hs.eg.db, drop = TRUE),
      error = function(e) NULL
    )
    if (!is.null(uni_map)) universe_entrez <- unique(uni_map$ENTREZID)
    message("  Universe: ", length(universe), " symbols → ",
            length(universe_entrez), " ENTREZ IDs")
  }

  results <- list()

  # Helper: save result and dotplot
  save_enrich <- function(er, tag) {
    if (!is.null(er) && nrow(er) > 0) {
      er_readable <- tryCatch(setReadable(er, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
                              error = function(e) er)
      write.csv(as.data.frame(er_readable),
                file.path(enrich_dir, paste0(label, "_", tag, ".csv")),
                row.names = FALSE)
      p <- dotplot(er_readable, showCategory = 15) +
        ggtitle(paste(label, gsub("_", " ", tag))) +
        theme(axis.text.y = element_text(size = 8))
      ggsave(file.path(enrich_dir, paste0(label, "_", tag, "_dotplot.pdf")),
             p, width = 9, height = 7)
      message("    [OK] ", tag, ": ", nrow(er), " terms")
      return(er_readable)
    } else {
      message("    [--] ", tag, ": no significant terms")
      return(NULL)
    }
  }

  # 1. KEGG pathways
  results$kegg <- tryCatch({
    er <- enrichKEGG(gene          = entrez,
                     universe      = universe_entrez,
                     organism      = "hsa",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = padj_cutoff)
    save_enrich(er, "KEGG_Pathways")
  }, error = function(e) { message("    [!] KEGG failed: ", e$message); NULL })

  # 2. Disease Ontology (DO)
  results$do <- tryCatch({
    er <- enrichDO(gene          = entrez,
                   universe      = universe_entrez,
                   ont           = "DO",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = padj_cutoff,
                   readable      = TRUE)
    save_enrich(er, "Disease_Ontology")
  }, error = function(e) { message("    [!] DO failed: ", e$message); NULL })

  # 3. DisGeNET (curated disease–gene associations)
  results$dgn <- tryCatch({
    er <- enrichDGN(gene          = entrez,
                    universe      = universe_entrez,
                    pAdjustMethod = "BH",
                    pvalueCutoff  = padj_cutoff,
                    readable      = TRUE)
    save_enrich(er, "DisGeNET")
  }, error = function(e) { message("    [!] DisGeNET failed: ", e$message); NULL })

  # 4. Network of Cancer Genes (NCG)
  results$ncg <- tryCatch({
    er <- enrichNCG(gene          = entrez,
                    universe      = universe_entrez,
                    pAdjustMethod = "BH",
                    pvalueCutoff  = padj_cutoff,
                    readable      = TRUE)
    save_enrich(er, "Cancer_Genes")
  }, error = function(e) { message("    [!] NCG failed: ", e$message); NULL })

  # Combined summary table: top 5 terms per database
  summary_rows <- list()
  db_labels <- c(kegg = "KEGG", do = "DiseaseOntology", dgn = "DisGeNET", ncg = "CancerGenes")
  for (db in names(db_labels)) {
    er <- results[[db]]
    if (!is.null(er) && nrow(er) > 0) {
      df <- as.data.frame(er)[1:min(5, nrow(er)), c("Description", "GeneRatio", "p.adjust")]
      df$Database <- db_labels[[db]]
      df$Label    <- label
      summary_rows[[db]] <- df
    }
  }
  if (length(summary_rows) > 0) {
    summary_df <- do.call(rbind, summary_rows)
    write.csv(summary_df,
              file.path(enrich_dir, paste0(label, "_top_disease_terms.csv")),
              row.names = FALSE)
  }

  return(invisible(results))
}


#' Run GO + disease enrichment together (convenience wrapper)
#'
#' @param peaks          Character vector of significant peak IDs
#' @param species        Species for gene lookup (default "Human")
#' @param label          Label prefix for filenames
#' @param out_dir        Output directory
#' @param anno_df        Master annotation data.frame
#' @param ct             Cell type name (used to route output to {out_dir}/{ct}/enrichment/)
#' @param universe_peaks Character vector of ALL tested peak IDs for this contrast
#'                       (used as enrichment background). If NULL, defaults to all
#'                       annotated genes in org.Hs.eg.db — too broad; always pass this.
run_full_enrichment <- function(peaks, species = "Human", label, out_dir, anno_df,
                                ct = NULL, universe_peaks = NULL) {
  if (length(peaks) < 5) {
    message("  Skipping enrichment for ", label, ": too few peaks (", length(peaks), ")")
    return(invisible(NULL))
  }
  genes <- get_peak_info(peaks, species, "gene", anno_df)$gene
  genes <- unique(genes[genes != "." & !is.na(genes) & genes != ""])

  # Build universe from all tested peaks
  universe_genes <- NULL
  if (!is.null(universe_peaks) && length(universe_peaks) > 0) {
    ug <- get_peak_info(universe_peaks, species, "gene", anno_df)$gene
    universe_genes <- unique(ug[ug != "." & !is.na(ug) & ug != ""])
  }

  ct_out_dir <- if (!is.null(ct)) file.path(out_dir, ct) else out_dir
  message("\n--- Enrichment: ", label, " (", length(peaks), " peaks, ",
          length(genes), " genes; universe: ",
          if (is.null(universe_genes)) "org.Hs.eg.db default" else paste0(length(universe_genes), " genes"),
          ") ---")
  go_res      <- run_go_enrichment(genes, label, ct_out_dir, universe = universe_genes)
  disease_res <- run_disease_enrichment(genes, label, ct_out_dir, universe = universe_genes)
  invisible(list(go = go_res, disease = disease_res))
}


#' Run GREAT analysis for a set of peak coordinates
#'
#' @param peak_ids  Character vector of peak IDs
#' @param species   Species for coordinate lookup
#' @param anno_df   Master annotation data.frame
#' @param label     Label for filenames
#' @param out_dir   Directory to save results
#' @param genome    Genome assembly (default "hg38")
run_great_analysis <- function(peak_ids, species, anno_df, label, out_dir,
                               genome = "hg38") {
  suppressPackageStartupMessages(library(rGREAT))

  enrich_dir <- file.path(out_dir, "enrichment")
  dir.create(enrich_dir, showWarnings = FALSE, recursive = TRUE)

  coords <- get_peak_info(peak_ids, species, "coordinates", anno_df)
  if (is.null(coords) || nrow(coords) < 5) {
    message("  Too few peaks for GREAT: ", label)
    return(invisible(NULL))
  }

  tryCatch({
    job       <- submitGreatJob(coords, species = genome)
    res_great <- getEnrichmentTables(job)

    if ("GO Biological Process" %in% names(res_great)) {
      write.csv(res_great[["GO Biological Process"]],
                file.path(enrich_dir, paste0(label, "_GREAT_GO_BP.csv")))
      message("  GREAT results saved for: ", label)
    }
  }, error = function(e) {
    message("  GREAT failed for ", label, ": ", e$message)
  })
}


# =============================================================================
# SECTION 9: MOTIF ENRICHMENT (PARALLEL)
# =============================================================================

#' Run monaLisa binned motif enrichment with parallelization
#'
#' @param peak_list_named Named list of peak ID vectors (one per group)
#' @param anno_df         Master annotation data.frame
#' @param species         Species for coordinate lookup (default "Human")
#' @param n_cores         Number of parallel cores (default 10)
#' @param genome          BSgenome object or string (default hg38)
#' @param out_rds         Path to save the SummarizedExperiment RDS
#' @return SummarizedExperiment with enrichment results
run_parallel_motif_enrichment <- function(peak_list_named, anno_df,
                                          species    = "Human",
                                          n_cores    = 4L,
                                          genome     = NULL,
                                          out_rds    = NULL,
                                          batch_size = 25L) {
    suppressPackageStartupMessages({
      library(monaLisa)
      library(JASPAR2022)
      library(TFBSTools)
      library(BSgenome.Hsapiens.UCSC.hg38)
      library(Biostrings)
      library(BiocParallel)
      library(SummarizedExperiment)
    })

    if (is.null(genome)) genome <- BSgenome.Hsapiens.UCSC.hg38

    # Load JASPAR PWMs once
    pwms <- getMatrixSet(JASPAR2022,
                         opts = list(tax_group = "vertebrates",
                                     matrixtype = "PWM"))

    # Build full mapping and fetch sequences once
    mapping_df <- data.frame(
      peak_id = unlist(peak_list_named, use.names = FALSE),
      group   = rep(names(peak_list_named), lengths(peak_list_named))
    )

    coords     <- get_peak_info(mapping_df$peak_id, species, "coordinates", anno_df)
    mapping_df <- mapping_df[mapping_df$peak_id %in% coords$name, ]
    rownames(coords) <- coords$name
    coords     <- coords[mapping_df$peak_id, ]

    gr           <- GenomicRanges::makeGRangesFromDataFrame(coords)
    seqs_all     <- Biostrings::getSeq(genome, gr)
    group_factor <- factor(mapping_df$group, levels = names(peak_list_named))

    n_groups <- length(levels(group_factor))
    message("Total groups: ", n_groups, " — processing in batches of ", batch_size)
    message("Peak memory per batch: ~",
            round(length(seqs_all) / n_groups * batch_size * 400 / 1e6, 1),
            " MB sequences")

    # Use SnowParam for PWM-level parallelism within each batch.
    # Keep workers low (2-4): SnowParam serialises ALL seqs to each worker,
    # so more workers = proportionally more memory.
    n_workers <- min(as.integer(n_cores), 4L)
    bp <- if (n_workers > 1L) {
      SnowParam(workers = n_workers, type = "SOCK", progressbar = FALSE)
    } else {
      SerialParam()
    }
    message("Parallel backend: ", class(bp)[1], " (", n_workers, " workers)")

    # Split groups into batches and run each batch separately
    group_names <- levels(group_factor)
    batch_idx   <- split(seq_along(group_names),
                         ceiling(seq_along(group_names) / batch_size))

    se_batches <- vector("list", length(batch_idx))

    for (bi in seq_along(batch_idx)) {
      batch_groups <- group_names[batch_idx[[bi]]]
      keep_rows    <- as.character(group_factor) %in% batch_groups

      batch_seqs   <- seqs_all[keep_rows]
      batch_factor <- droplevels(group_factor[keep_rows])

      message("  Batch ", bi, "/", length(batch_idx),
              ": ", length(batch_groups), " groups, ",
              length(batch_seqs), " sequences")

      se_batches[[bi]] <- tryCatch(
        calcBinnedMotifEnrR(seqs      = batch_seqs,
                            bins      = batch_factor,
                            pwmL      = pwms,
                            min.score = 10,
                            BPPARAM   = bp,
                            verbose   = FALSE),
        error = function(e) {
          warning("Batch ", bi, " failed (", conditionMessage(e),
                  ") — retrying serially")
          calcBinnedMotifEnrR(seqs      = batch_seqs,
                              bins      = batch_factor,
                              pwmL      = pwms,
                              min.score = 10,
                              BPPARAM   = SerialParam(),
                              verbose   = FALSE)
        }
      )
      gc()  # release memory between batches
    }

    # Combine batch SE objects by column-binding assays
    message("Combining ", length(se_batches), " batch results...")
    se_enrich <- do.call(cbind, se_batches)

    if (!is.null(out_rds)) {
      dir.create(dirname(out_rds), showWarnings = FALSE, recursive = TRUE)
      saveRDS(se_enrich, out_rds)
      message("Saved SE object: ", basename(out_rds))
    }

    return(se_enrich)
}


#' Plot motif enrichment heatmap from SummarizedExperiment
#'
#' @param se_obj    SummarizedExperiment from monaLisa
#' @param title     Plot title
#' @param out_file  Output PDF path
#' @param top_per_col Number of top motifs to keep per column
plot_motif_heatmap <- function(se_obj, title, out_file, top_per_col = 3) {
  suppressPackageStartupMessages({
    library(ComplexHeatmap)
    library(circlize)
    library(SummarizedExperiment)
  })

  mat          <- assay(se_obj, "negLog10Padj")
  rownames(mat) <- gsub(".*_", "", rownames(mat))

  top_indices <- unique(as.vector(
    apply(mat, 2, function(x) order(x, decreasing = TRUE)[1:top_per_col])
  ))
  plot_sub <- mat[top_indices, , drop = FALSE]
  plot_sub[plot_sub > 20] <- 20

  col_fun <- colorRamp2(c(0, 2, 20), c("white", "orange", "red"))

  ht <- Heatmap(plot_sub,
                name             = "Enrichment",
                col              = col_fun,
                column_title     = title,
                cluster_rows     = TRUE,
                cluster_columns  = TRUE,
                show_row_names   = TRUE,
                border           = TRUE,
                column_names_rot = 45,
                heatmap_legend_param = list(title = "-log10(p-adj)"))

  dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
  pdf(out_file, width = 12, height = 14)
  draw(ht, merge_legend = TRUE)
  dev.off()
  message("  Motif heatmap saved: ", basename(out_file))
}


# =============================================================================
# SECTION 10: QC PLOTS
# =============================================================================

#' Plot library size distribution per species
#'
#' @param meta       Metadata with total_counts and species columns
#' @param out_file   Output PDF path
#' @param min_counts Threshold line to draw
#' @param min_cells  Cell threshold (for subtitle)
plot_count_distribution <- function(meta, out_file,
                                    min_counts = 50000, min_cells = 100) {
  suppressPackageStartupMessages(library(ggplot2))

  p <- ggplot(meta, aes(x = total_counts, fill = species)) +
    geom_histogram(bins = 30, color = "black", alpha = 0.8) +
    geom_vline(xintercept = min_counts, color = "red",
               linetype = "dashed", linewidth = 1) +
    facet_wrap(~ species, scales = "free_y") +
    theme_bw() +
    scale_x_log10(labels = scales::comma) +
    labs(title    = "Total Counts per Pseudobulk Sample",
         subtitle = paste("Thresholds: Counts >=", min_counts,
                          "| Cells >=", min_cells),
         x = "Total Counts (log10 scale)",
         y = "Number of Samples") +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))

  dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
  ggsave(out_file, p, width = 10, height = 6)
  message("  Count distribution plot saved.")
}


# =============================================================================
# SECTION 11: RNA PSEUDOBULK DESeq2
# =============================================================================
# Mirrors the ATAC pseudobulk pipeline but loads from parquet files produced
# by scripts/create_rna_pseudobulk.py (generated from the merged h5ad).
#
# Pseudobulk unit: Individual × cell_type × Region  (non-aggregated)
#                  Individual × cell_type            (region-aggregated, default)
# Design: ~ 0 + species (same no-intercept model as ATAC)
# ==============================================================================

#' Load RNA pseudobulk data from parquet files generated by create_rna_pseudobulk.py
#'
#' @param pb_dir   Root directory with per-species subdirs (each has
#'                 pseudobulk_counts.parquet + pseudobulk_meta.parquet)
#' @param species  Character vector of species to load
#' @return list(all_counts = list(species -> matrix), all_meta = list(species -> df))
load_rna_pseudobulk <- function(pb_dir, species) {
  suppressPackageStartupMessages(library(arrow))

  all_counts <- list()
  all_meta   <- list()

  for (sp in species) {
    sp_dir     <- file.path(pb_dir, sp)
    counts_f   <- file.path(sp_dir, "pseudobulk_counts.parquet")
    meta_f     <- file.path(sp_dir, "pseudobulk_meta.parquet")

    if (!file.exists(counts_f) || !file.exists(meta_f)) {
      message("Skipping ", sp, ": files not found at ", sp_dir)
      next
    }

    counts_mat <- as.data.frame(read_parquet(counts_f))
    meta_df    <- as.data.frame(read_parquet(meta_f))

    # counts: rows=genes, cols=sample_ids (from parquet index)
    # Restore row names from the first column if it was saved as data column
    if (!is.null(counts_mat[["__index_level_0__"]])) {
      rownames(counts_mat) <- counts_mat[["__index_level_0__"]]
      counts_mat[["__index_level_0__"]] <- NULL
    }
    counts_mat <- as.matrix(counts_mat)
    storage.mode(counts_mat) <- "integer"

    # Meta: index is sample_id
    if (!is.null(meta_df[["sample_id"]])) {
      rownames(meta_df) <- meta_df[["sample_id"]]
    }

    # Align columns
    shared_samples <- intersect(colnames(counts_mat), rownames(meta_df))
    counts_mat <- counts_mat[, shared_samples, drop = FALSE]
    meta_df    <- meta_df[shared_samples, , drop = FALSE]

    all_counts[[sp]] <- counts_mat
    all_meta[[sp]]   <- meta_df
    message("Loaded ", sp, ": ", nrow(counts_mat), " genes x ", ncol(counts_mat), " samples")
  }
  return(list(all_counts = all_counts, all_meta = all_meta))
}


#' Merge RNA pseudobulk matrices across species (inner = shared genes)
#'
#' @param all_counts Named list of count matrices (genes × samples)
#' @param all_meta   Named list of metadata data.frames
#' @return list(counts = merged_matrix, meta = merged_meta)
merge_rna_pseudobulk <- function(all_counts, all_meta, join_type = "inner") {
  if (join_type == "inner") {
    shared_genes <- Reduce(intersect, lapply(all_counts, rownames))
    message("Inner join: ", length(shared_genes), " shared genes")
  } else {
    shared_genes <- Reduce(union, lapply(all_counts, rownames))
    message("Union join: ", length(shared_genes), " total genes")
  }

  mats <- lapply(all_counts, function(m) {
    missing <- setdiff(shared_genes, rownames(m))
    if (length(missing) > 0) {
      pad <- matrix(0L, nrow = length(missing), ncol = ncol(m),
                    dimnames = list(missing, colnames(m)))
      m <- rbind(m[intersect(shared_genes, rownames(m)), , drop = FALSE], pad)
    }
    m[shared_genes, , drop = FALSE]
  })

  merged_counts <- do.call(cbind, unname(mats))
  # unname() prevents do.call(rbind/cbind) from prefixing dimension names with list element names
  merged_meta   <- do.call(rbind, unname(all_meta))
  # Restore rownames from sample_id column if needed
  if (!is.null(merged_meta[["sample_id"]]) &&
      !identical(rownames(merged_meta), as.character(merged_meta[["sample_id"]]))) {
    rownames(merged_meta) <- as.character(merged_meta[["sample_id"]])
  }
  merged_meta$species <- factor(merged_meta$species)

  message("Merged RNA matrix: ", nrow(merged_counts), " genes x ", ncol(merged_counts), " samples")
  return(list(counts = merged_counts, meta = merged_meta))
}


#' Aggregate RNA pseudobulk across regions (sum counts, merge metadata)
#'
#' Collapses Individual × cell_type × Region → Individual × cell_type
#' by summing counts and summing n_cells.
#'
#' @param counts  Genes × samples count matrix
#' @param meta    Samples × metadata data.frame (must have donor, cell_type, species)
#' @return list(counts, meta)
aggregate_rna_pseudobulk <- function(counts, meta) {
  meta$agg_key <- paste(meta$donor, meta$cell_type, meta$species, sep = "__")
  keys <- unique(meta$agg_key)

  new_counts <- matrix(0L, nrow = nrow(counts), ncol = length(keys),
                       dimnames = list(rownames(counts), keys))
  n_keys     <- length(keys)
  new_meta   <- data.frame(
    donor     = rep(NA_character_, n_keys),
    cell_type = rep(NA_character_, n_keys),
    species   = rep(NA_character_, n_keys),
    n_cells   = rep(0L,            n_keys),
    n_counts  = rep(0L,            n_keys),
    row.names = keys,
    stringsAsFactors = FALSE
  )

  for (k in keys) {
    cols <- rownames(meta)[meta$agg_key == k]
    new_counts[, k] <- rowSums(counts[, cols, drop = FALSE])
    m0 <- meta[cols[1], , drop = FALSE]
    new_meta[k, "donor"]     <- as.character(m0$donor)
    new_meta[k, "cell_type"] <- as.character(m0$cell_type)
    new_meta[k, "species"]   <- as.character(m0$species)
    new_meta[k, "n_cells"]   <- sum(meta[cols, "n_cells"], na.rm = TRUE)
    new_meta[k, "n_counts"]  <- sum(new_counts[, k])
  }
  storage.mode(new_counts) <- "integer"
  return(list(counts = new_counts, meta = new_meta))
}


#' Run DESeq2 on RNA pseudobulk for evolutionary contrasts
#'
#' Replicates run_deseq2_evolutionary() but for RNA expression.
#' Uses the same no-intercept design, per-contrast orthology subsetting is skipped
#' (all genes are shared by construction from the inner-join merge).
#'
#' @param counts_rna    Genes × samples integer count matrix
#' @param meta_rna      Sample metadata (must have species, cell_type columns)
#' @param evo_contrasts Named list of pos/neg species vectors (from define_evolutionary_contrasts())
#' @param out_dir       Base output directory
#' @param min_samples   Minimum samples per species group (default 2)
#' @param filter_outliers  Apply per-contrast outlier detection (default TRUE)
#' @return Nested list: cell_type → contrast → DESeqResults
run_deseq2_rna_evolutionary <- function(counts_rna, meta_rna, evo_contrasts, out_dir,
                                        min_samples = 2, filter_outliers = TRUE,
                                        sd_thresh = 2.5) {
  suppressPackageStartupMessages(library(DESeq2))

  dir.create(file.path(out_dir, "_summary"), showWarnings = FALSE, recursive = TRUE)

  meta_rna$species   <- as.character(meta_rna$species)
  meta_rna$cell_type <- make.names(as.character(meta_rna$cell_type))
  cell_types <- unique(meta_rna$cell_type)

  results <- list()

  for (ct in cell_types) {
    ct_idx    <- rownames(meta_rna)[meta_rna$cell_type == ct]
    ct_counts <- counts_rna[, ct_idx, drop = FALSE]
    ct_meta   <- meta_rna[ct_idx, , drop = FALSE]

    message("\n=== [RNA] ", ct, " ===")

    # Remove all-zero genes for this cell type
    keep_genes <- rowSums(ct_counts) > 0
    ct_counts  <- ct_counts[keep_genes, , drop = FALSE]

    results[[ct]] <- list()
    outlier_log   <- data.frame()
    ct_out_dir    <- file.path(out_dir, ct, "rna_evolutionary")
    dir.create(ct_out_dir, showWarnings = FALSE, recursive = TRUE)

    for (contrast_name in names(evo_contrasts)) {
      pos_sp <- evo_contrasts[[contrast_name]]$pos
      neg_sp <- evo_contrasts[[contrast_name]]$neg
      all_sp <- c(pos_sp, neg_sp)

      ct_sp_idx  <- rownames(ct_meta)[ct_meta$species %in% all_sp]
      if (length(ct_sp_idx) == 0) next

      sub_counts <- ct_counts[, ct_sp_idx, drop = FALSE]
      sub_meta   <- ct_meta[ct_sp_idx, , drop = FALSE]

      # Per-contrast outlier removal
      if (filter_outliers) {
        flags <- detect_per_contrast_outliers(sub_counts, sub_meta,
                                              sd_thresh = sd_thresh)
        if (any(flags)) {
          removed <- names(flags)[flags]
          message("  [outlier] ", contrast_name, ": removing ", paste(removed, collapse=", "))
          outlier_log <- append_outlier_log(outlier_log, ct, contrast_name, removed,
                                            colSums(sub_counts)[removed])
          keep_s      <- !flags
          sub_counts  <- sub_counts[, keep_s, drop = FALSE]
          sub_meta    <- sub_meta[keep_s, , drop = FALSE]
        }
      }

      # Check min samples per group
      pos_n <- sum(sub_meta$species %in% pos_sp)
      neg_n <- sum(sub_meta$species %in% neg_sp)
      if (pos_n < min_samples || neg_n < min_samples) next

      # Filter low-count genes for this contrast
      keep_g    <- rowSums(sub_counts >= 5) >= min_samples
      sub_counts <- sub_counts[keep_g, , drop = FALSE]
      if (nrow(sub_counts) < 100) next

      sub_meta$species_group <- factor(
        ifelse(sub_meta$species %in% pos_sp, "pos", "neg"),
        levels = c("neg", "pos")
      )

      tryCatch({
        dds <- DESeqDataSetFromMatrix(
          countData = sub_counts,
          colData   = sub_meta,
          design    = ~ species_group
        )
        sizeFactors(dds) <- estimateSizeFactors(dds, type = "poscounts")$sizeFactor
        dds <- DESeq(dds, fitType = "local", quiet = TRUE)
        res <- results(dds, contrast = c("species_group", "pos", "neg"),
                       independentFiltering = TRUE)
        res$gene <- rownames(res)

        # Save CSV
        res_df <- as.data.frame(res)
        res_df <- res_df[order(res_df$padj, na.last = TRUE), ]
        write.csv(res_df, file.path(ct_out_dir, paste0(contrast_name, "_DE.csv")))

        n_up   <- sum(!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange > 1, na.rm=TRUE)
        n_down <- sum(!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange < -1, na.rm=TRUE)
        message("  ", sprintf("%-40s", contrast_name), ": ",
                n_up, " UP, ", n_down, " DOWN  (padj<0.05, |LFC|>1)")
        results[[ct]][[contrast_name]] <- res

      }, error = function(e) {
        message("  [!] DESeq2 failed for ", ct, " / ", contrast_name, ": ", e$message)
      })
    }

    # Save outlier log
    if (nrow(outlier_log) > 0) {
      save_outlier_summary(outlier_log, ct_out_dir, prefix = "RNA")
    }
  }

  # Save full results list
  saveRDS(results, file.path(out_dir, "_summary", "RNA_DE_res_list.rds"))
  message("\nRNA evolutionary DE complete. Saved RNA_DE_res_list.rds")
  return(invisible(results))
}


#' Run DESeq2 RNA DE for shared-peak-style all-vs-one contrasts (species_X vs rest)
#'
#' @param counts_rna  Genes × samples count matrix
#' @param meta_rna    Sample metadata with species, cell_type
#' @param species     Character vector of species
#' @param out_dir     Output directory
#' @return Nested list: cell_type → species → DESeqResults
run_deseq2_rna_species <- function(counts_rna, meta_rna, species, out_dir,
                                   min_samples = 2, filter_outliers = TRUE,
                                   sd_thresh = 2.5) {
  suppressPackageStartupMessages(library(DESeq2))

  dir.create(file.path(out_dir, "_summary"), showWarnings = FALSE, recursive = TRUE)

  meta_rna$species   <- factor(as.character(meta_rna$species), levels = species)
  meta_rna$cell_type <- make.names(as.character(meta_rna$cell_type))
  cell_types <- unique(meta_rna$cell_type)
  results    <- list()

  for (ct in cell_types) {
    ct_idx    <- rownames(meta_rna)[meta_rna$cell_type == ct]
    ct_counts <- counts_rna[, ct_idx, drop = FALSE]
    ct_meta   <- meta_rna[ct_idx, , drop = FALSE]
    message("\n=== [RNA species] ", ct, " (", ncol(ct_counts), " samples) ===")

    keep_g    <- rowSums(ct_counts) > 0
    ct_counts <- ct_counts[keep_g, , drop = FALSE]

    results[[ct]] <- list()
    ct_out        <- file.path(out_dir, ct, "rna_species")
    dir.create(ct_out, showWarnings = FALSE, recursive = TRUE)

    for (sp in species) {
      if (sum(ct_meta$species == sp) < min_samples) next
      if (sum(ct_meta$species != sp) < min_samples) next

      sub_meta             <- ct_meta
      sub_meta$is_focal    <- factor(ifelse(sub_meta$species == sp, sp, "Other"),
                                     levels = c("Other", sp))

      tryCatch({
        dds <- DESeqDataSetFromMatrix(
          countData = ct_counts,
          colData   = sub_meta,
          design    = ~ is_focal
        )
        sizeFactors(dds) <- estimateSizeFactors(dds, type = "poscounts")$sizeFactor
        dds <- DESeq(dds, fitType = "local", quiet = TRUE)
        res <- results(dds, contrast = c("is_focal", sp, "Other"),
                       independentFiltering = TRUE)
        res$gene <- rownames(res)

        res_df <- as.data.frame(res)
        res_df <- res_df[order(res_df$padj, na.last = TRUE), ]
        write.csv(res_df, file.path(ct_out, paste0(sp, "_vs_rest_DE.csv")))

        n_up <- sum(!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange > 1, na.rm=TRUE)
        message("  ", sprintf("%-12s", sp), ": ", n_up, " UP (padj<0.05, LFC>1)")
        results[[ct]][[sp]] <- res
      }, error = function(e) {
        message("  [!] ", sp, " failed: ", e$message)
      })
    }
  }

  saveRDS(results, file.path(out_dir, "_summary", "RNA_DE_species_res_list.rds"))
  message("RNA species DE complete.")
  return(invisible(results))
}


#' Ultra-robust filter for RNA DE genes
#'
#' Applies the same min(positive donors) > max(negative donors) criterion as
#' ultra_robust_filter() for ATAC peaks, but on RNA pseudobulk log2CPM.
#' A gene is "robust UP" if its lowest log2CPM across all positive-species donors
#' exceeds the highest log2CPM across all negative-species donors.
#' A gene is "robust DOWN" if its highest log2CPM across positive donors is below
#' the lowest log2CPM across all negative donors.
#'
#' @param rna_res    Nested list: cell_type → contrast → DESeqResults
#' @param counts_rna Genes × samples raw count matrix (all CTs combined)
#' @param meta_rna   Sample metadata with cell_type and species columns
#' @param contrasts  Named list from define_evolutionary_contrasts()
#' @param out_dir    Root output directory (cell-type-first: {out_dir}/{CT}/rna_robust/)
#' @param padj_thresh Adjusted p-value threshold (default 0.05)
#' @param lfc_thresh  |log2FC| threshold (default 1)
#' @param min_logcpm  Minimum log2CPM required in the positive group (default 0)
#' @return Nested list: cell_type → contrast → list(up=genes, down=genes)
ultra_robust_filter_rna <- function(rna_res, counts_rna, meta_rna, contrasts, out_dir,
                                    padj_thresh = 0.05, lfc_thresh = 1, min_logcpm = 0) {
  suppressPackageStartupMessages(library(matrixStats))

  valid_samples <- intersect(rownames(meta_rna), colnames(counts_rna))
  counts_rna    <- counts_rna[, valid_samples, drop = FALSE]

  robust_genes <- list()

  for (ct in names(rna_res)) {
    message("\nRobust RNA DE filter: ", ct)
    robust_genes[[ct]] <- list()

    meta_ct    <- meta_rna[meta_rna$cell_type == ct, , drop = FALSE]
    ct_samples <- intersect(rownames(meta_ct), colnames(counts_rna))
    if (length(ct_samples) < 2) next
    meta_ct   <- meta_ct[ct_samples, , drop = FALSE]
    ct_counts <- counts_rna[, ct_samples, drop = FALSE]

    lib_sizes  <- colSums(ct_counts)
    cpm_mat    <- t(t(ct_counts) / lib_sizes) * 1e6
    logcpm_mat <- log2(cpm_mat + 1)

    for (contrast_name in names(rna_res[[ct]])) {
      if (!contrast_name %in% names(contrasts)) next
      res <- rna_res[[ct]][[contrast_name]]
      if (is.null(res)) next

      res_df  <- as.data.frame(res)
      pos_sp  <- contrasts[[contrast_name]]$pos
      neg_sp  <- contrasts[[contrast_name]]$neg

      pos_samp <- ct_samples[meta_ct$species %in% pos_sp]
      neg_samp <- ct_samples[meta_ct$species %in% neg_sp]
      if (length(pos_samp) < 1 || length(neg_samp) < 1) next

      # --- Robust UP ---
      sig_up <- rownames(res_df)[!is.na(res_df$padj) &
                                   res_df$padj < padj_thresh &
                                   res_df$log2FoldChange > lfc_thresh]
      sig_up <- intersect(sig_up, rownames(logcpm_mat))

      up_robust <- character(0)
      if (length(sig_up) > 0) {
        min_pos  <- rowMins(as.matrix(logcpm_mat[sig_up, pos_samp, drop = FALSE]))
        max_neg  <- rowMaxs(as.matrix(logcpm_mat[sig_up, neg_samp, drop = FALSE]))
        keep     <- (min_pos > max_neg) & (min_pos > min_logcpm)
        keep[is.na(keep)] <- FALSE
        up_robust <- sig_up[keep]
      }

      # --- Robust DOWN ---
      sig_dn <- rownames(res_df)[!is.na(res_df$padj) &
                                   res_df$padj < padj_thresh &
                                   res_df$log2FoldChange < -lfc_thresh]
      sig_dn <- intersect(sig_dn, rownames(logcpm_mat))

      dn_robust <- character(0)
      if (length(sig_dn) > 0) {
        max_pos  <- rowMaxs(as.matrix(logcpm_mat[sig_dn, pos_samp, drop = FALSE]))
        min_neg  <- rowMins(as.matrix(logcpm_mat[sig_dn, neg_samp, drop = FALSE]))
        keep     <- (max_pos < min_neg) & (min_neg > min_logcpm)
        keep[is.na(keep)] <- FALSE
        dn_robust <- sig_dn[keep]
      }

      message(sprintf("  [%s]  UP: %d sig -> %d robust  |  DOWN: %d sig -> %d robust",
                      contrast_name,
                      length(sig_up), length(up_robust),
                      length(sig_dn), length(dn_robust)))

      robust_genes[[ct]][[contrast_name]] <- list(up = up_robust, down = dn_robust)

      n_robust <- length(up_robust) + length(dn_robust)
      if (n_robust > 0) {
        ct_dir <- file.path(out_dir, ct, "rna_robust")
        dir.create(ct_dir, showWarnings = FALSE, recursive = TRUE)

        parts <- list()
        if (length(up_robust) > 0) {
          df <- res_df[up_robust, , drop = FALSE]
          df$gene <- rownames(df); df$direction <- "UP"
          parts[[1]] <- df
        }
        if (length(dn_robust) > 0) {
          df <- res_df[dn_robust, , drop = FALSE]
          df$gene <- rownames(df); df$direction <- "DOWN"
          parts[[2]] <- df
        }
        out_df <- do.call(rbind, parts)
        out_df <- out_df[order(out_df$padj, na.last = TRUE), ]
        write.csv(out_df,
                  file.path(ct_dir, paste0(contrast_name, "_robust_DE.csv")),
                  row.names = FALSE)
      }
    }
  }

  summary_dir <- file.path(out_dir, "_summary")
  dir.create(summary_dir, showWarnings = FALSE, recursive = TRUE)
  saveRDS(robust_genes, file.path(summary_dir, "RNA_DE_robust_list.rds"))
  message("\nRobust RNA DE complete. Saved RNA_DE_robust_list.rds")
  return(invisible(robust_genes))
}
