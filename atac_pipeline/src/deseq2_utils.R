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
                                    out_dir) {
  suppressPackageStartupMessages(library(DESeq2))

  res_dir <- file.path(out_dir, "differential_results", "shared_peaks")
  dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

  res_list   <- list()
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
    keep_ct           <- rowSums(counts_ct >= 10) >= 3
    counts_ct_filtered <- counts_ct[keep_ct, ]

    dds_ct <- DESeqDataSetFromMatrix(
      countData = round(counts_ct_filtered),
      colData   = meta_ct,
      design    = ~ 0 + species
    )
    dds_ct <- estimateSizeFactors(dds_ct, type = "poscounts")
    dds_ct <- DESeq(dds_ct, quiet = TRUE)

    res_list[[ct]] <- list()
    ct_dir <- file.path(res_dir, ct)
    dir.create(ct_dir, showWarnings = FALSE)

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
      out_file          <- file.path(ct_dir, paste0(target_sp, "_vs_All_Others.csv"))
      write.csv(res_df, out_file, row.names = FALSE)

      sig_hits <- sum(res_sp$padj < 0.05 & res_sp$log2FoldChange > 1, na.rm = TRUE)
      message(sprintf("  %-12s specific: %d (padj<0.05, LFC>1)", target_sp, sig_hits))
    }
  }

  saveRDS(res_list, file.path(res_dir, "DESeq2_res_list_shared.rds"))
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

  res_dir <- file.path(out_dir, "differential_results", "contrast_peaks")
  dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

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
    ct_dir <- file.path(res_dir, ct)
    dir.create(ct_dir, showWarnings = FALSE)

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
      out_file       <- file.path(ct_dir, paste0(target_sp, "_vs_All_Others.csv"))
      write.csv(res_df, out_file, row.names = FALSE)

      sig_hits <- sum(res_sp$padj < 0.05 & res_sp$log2FoldChange > 1, na.rm = TRUE)
      message(sprintf("  %-12s specific: %d (padj<0.05, LFC>1)", target_sp, sig_hits))
    }
  }

  saveRDS(res_list, file.path(res_dir, "DESeq2_res_list_contrast.rds"))
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
    "Pair_Human_vs_Marmoset" = list(pos = c("Human"), neg = c("Marmoset"))
  )
  message("Defined ", length(contrasts), " evolutionary contrasts.")
  return(contrasts)
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
                                    ortho_mat, out_dir, min_samples = 2) {
  suppressPackageStartupMessages(library(DESeq2))

  branch_dir <- file.path(out_dir, "differential_results", "evolutionary_branches")
  dir.create(branch_dir, showWarnings = FALSE, recursive = TRUE)

  branch_res <- list()
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
        ct_dir <- file.path(branch_dir, ct)
        dir.create(ct_dir, showWarnings = FALSE)
        res_df         <- as.data.frame(res_ordered)
        res_df$peak_id <- rownames(res_df)
        res_df <- res_df[, c("peak_id", setdiff(names(res_df), "peak_id"))]
        write.csv(res_df, file.path(ct_dir, paste0(node_name, ".csv")),
                  row.names = FALSE)

        sig_up <- sum(res$padj < 0.05 & res$log2FoldChange > 1, na.rm = TRUE)
        message(sprintf("  %-35s: %d peaks tested, %d sig UP",
                        node_name, nrow(node_counts), sig_up))
      }, error = function(e) {
        message(sprintf("  %-35s: FAILED (%s)", node_name, e$message))
      })
    }
  }

  saveRDS(branch_res, file.path(branch_dir, "branch_res_list.rds"))
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

  robust_dir <- file.path(out_dir, "differential_results", "ultra_robust_branches")
  dir.create(robust_dir, showWarnings = FALSE, recursive = TRUE)

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
        ct_node_dir <- file.path(robust_dir, ct)
        dir.create(ct_node_dir, showWarnings = FALSE)

        out_df <- data.frame(
          peak_id        = surviving,
          min_pos_logCPM = min_pos[is_robust],
          max_neg_logCPM = max_neg[is_robust]
        )
        write.csv(out_df,
                  file.path(ct_node_dir, paste0(node_name, "_UltraRobust.csv")),
                  row.names = FALSE)
      }
    }
  }

  saveRDS(robust_peaks, file.path(robust_dir, "ultra_robust_peaks_list.rds"))
  message("\nUltra-robust filtering complete. Saved checkpoint.")
  return(robust_peaks)
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

  cluster_dir <- file.path(out_dir, "differential_results", "cross_celltype_modules")
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
#' @return enrichResult object (or NULL)
run_go_enrichment <- function(genes, label, out_dir) {
  suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Hs.eg.db)
  })

  enrich_dir <- file.path(out_dir, "differential_results", "enrichment")
  dir.create(enrich_dir, showWarnings = FALSE, recursive = TRUE)

  genes <- unique(genes[genes != "." & !is.na(genes)])
  if (length(genes) < 5) {
    message("  Too few genes for ", label, " (", length(genes), ")")
    return(invisible(NULL))
  }

  ego <- enrichGO(gene          = genes,
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

  enrich_dir <- file.path(out_dir, "differential_results", "enrichment")
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
