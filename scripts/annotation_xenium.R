#!/usr/bin/env Rscript
# ==============================================================================
# annotate_xenium.R
#
# Consensus cell type annotation of Proseg-segmented Xenium spatial data
# using a pre-built scRNA-seq reference (from build_reference.R).
#
# Key design: the full-genome reference (~52K genes) is used to annotate
# ~193K reference cells at build time. But SingleR and RCTD are trained
# here on only the genes present in the Xenium panel (~500-5000 genes).
# This ensures all marker genes selected by SingleR's DE step are actually
# measurable in the spatial data, producing cleaner Spearman correlations
# than zero-padding absent genes would.
#
# Pipeline:
#   1. Auto-discover proseg regions
#   2. Scan panel genes (union across all regions)
#   3. Load full consensus reference, subset to panel genes
#   4. Train SingleR + build RCTD components on panel-gene subset
#   5. For each Xenium region:
#      a. Load proseg counts + cell metadata
#      b. Run SingleR (per-cell correlation-based)
#      c. Run RCTD doublet mode (probabilistic deconvolution)
#      d. Consensus scoring -- label only cells where BOTH methods agree
#         and BOTH pass their individual confidence thresholds
#      e. Tag unlabeled cells with reason codes
#      f. Generate diagnostics
#
# Inputs:
#   From build_reference.R output directory:
#     - consensus_reference.rds    Full integrated Seurat object
#
#   Proseg output per region (r_seurat/ subdirectory):
#     - counts.mtx.gz
#     - cell-metadata.csv.gz
#     - gene-metadata.csv.gz
#
# Outputs per region (in output_root/<region>/):
#   - <region>_annotated.rds            Seurat object with all annotations
#   - all_cells_annotation.csv.gz       Full per-cell table
#   - unlabeled_cells.csv.gz            Unlabeled cells + reason codes
#   - singler_raw_results.rds           Raw SingleR output for reanalysis
#   - singler_vs_rctd_confusion.csv     Agreement matrix
#   - disagreement_matrix.csv           Label pairs where methods disagree
#   - spatial_consensus.pdf             Spatial map of consensus labels
#   - spatial_unlabeled_reasons.pdf     Spatial map of reason codes
#   - singler_delta_distribution.pdf    Per-label delta.next violins
#   - confidence_scatter.pdf            2D confidence space (delta vs weight)
#
# Panel-specific cached models (in output_root/):
#   - singler_trained_panel.rds         SingleR model trained on panel genes
#   - rctd_reference_panel.rds          RCTD components on panel genes
#   - panel_genes.txt                   Gene list used for training
#   - global_summary.csv                Per-region annotation summary
#
# Dependencies:
#   install.packages(c("Seurat", "SeuratObject", "Matrix", "data.table",
#                      "ggplot2", "patchwork", "dplyr"))
#   BiocManager::install(c("SingleR", "SingleCellExperiment", "spacexr",
#                          "BiocParallel", "SummarizedExperiment"))
#   devtools::install_github("immunogenomics/presto")  # optional
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(SingleR)
  library(SingleCellExperiment)
  library(spacexr)
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(BiocParallel)
})

# Source shared helpers
source("./scripts/helper_functions.R")

# ==============================================================================
# CONFIGURATION
# ==============================================================================

config <- list(

  # --- Reference (from build_reference.R) ---
  reference_dir = "./reference/",

  # --- Proseg results root ---
  # Expected structure: proseg_root/<dataset>/r_seurat/{counts,cell-metadata,gene-metadata}
  proseg_root = "./proseg_results",
  output_root = "./data",

  # --- Xenium QC ---
  xen_min_features = 5,
  xen_min_counts   = 10,

  # --- Coordinate column names in proseg cell-metadata.csv.gz ---
  coord_x = "x",
  coord_y = "y",

  # --- SingleR parameters ---
  singler_quantile  = 0.8,
  singler_fine_tune = TRUE,
  singler_de_method = "wilcox",

  # --- RCTD parameters ---
  rctd_mode         = "doublet",
  rctd_max_cores    = 8,
  rctd_gene_cutoff  = 0.000125,
  rctd_fc_cutoff    = 0.5,
  rctd_fc_cutoff_reg = 0.75,
  rctd_UMI_min      = 0,
  rctd_counts_MIN   = 0,
  rctd_CELL_MIN_INSTANCE = 10,

  # --- Consensus thresholds ---
  singler_delta_floor = 0.05,
  singler_delta_pctile_threshold = 0.95,
  rctd_weight_threshold = 0.70,
  rctd_confident_classes = c("singlet", "doublet_certain"),

  # --- Panel-specific model caching ---
  # If TRUE and cached models exist with matching gene list, skip retraining.
  use_cache = TRUE
)

dir.create(config$output_root, showWarnings = FALSE, recursive = TRUE)


# ==============================================================================
# AUTO-DISCOVER PROSEG REGIONS
# ==============================================================================

log_msg("================================================================")
log_msg(" Discovering proseg regions in: ", config$proseg_root)
log_msg("================================================================\n")

required_files <- c("counts.mtx.gz", "cell-metadata.csv.gz", "gene-metadata.csv.gz")

r_seurat_dirs <- list.dirs(config$proseg_root, recursive = TRUE, full.names = TRUE)
r_seurat_dirs <- r_seurat_dirs[basename(r_seurat_dirs) == "r_seurat"]

xenium_regions <- list()
for (rdir in r_seurat_dirs) {
  present <- required_files %in% list.files(rdir)
  if (all(present)) {
    region_name <- basename(dirname(rdir))
    xenium_regions[[region_name]] <- rdir
    log_msg("  Found: ", region_name, " -> ", rdir)
  } else {
    missing <- required_files[!present]
    log_msg("  [SKIP] ", rdir, " -- missing: ", paste(missing, collapse = ", "))
  }
}

if (length(xenium_regions) == 0) {
  stop("[ERROR] No valid proseg r_seurat directories found under: ",
       config$proseg_root)
}

log_msg("\nDiscovered ", length(xenium_regions), " regions: ",
        paste(names(xenium_regions), collapse = ", "), "\n")


# ==============================================================================
# DISCOVER PANEL GENES (union across all regions)
# ==============================================================================

log_msg("================================================================")
log_msg(" Discovering Xenium panel genes")
log_msg("================================================================\n")

panel_genes <- character(0)
for (nm in names(xenium_regions)) {
  gm <- fread(file.path(xenium_regions[[nm]], "gene-metadata.csv.gz"))
  region_genes <- gm$gene
  log_msg("  ", nm, ": ", length(region_genes), " genes")
  panel_genes <- union(panel_genes, region_genes)
}

panel_genes <- sort(unique(panel_genes))
log_msg("\nPanel gene union: ", length(panel_genes), " genes across ",
        length(xenium_regions), " regions\n")


# ==============================================================================
# LOAD REFERENCE + PANEL-SPECIFIC TRAINING
# ==============================================================================

log_msg("================================================================")
log_msg(" Building panel-specific reference models")
log_msg("================================================================\n")

# Cache paths
cache_singler <- file.path(config$output_root, "singler_trained_panel.rds")
cache_rctd    <- file.path(config$output_root, "rctd_reference_panel.rds")
cache_genes   <- file.path(config$output_root, "panel_genes.txt")

# Check if cached models exist and match current panel
use_cached <- FALSE
if (config$use_cache &&
    file.exists(cache_singler) &&
    file.exists(cache_rctd) &&
    file.exists(cache_genes)) {

  cached_genes <- readLines(cache_genes)
  if (identical(sort(cached_genes), sort(panel_genes))) {
    log_msg("[CACHE] Panel genes match cached models. Loading from cache.")
    use_cached <- TRUE
  } else {
    log_msg("[CACHE] Panel genes changed (",
            length(cached_genes), " cached vs ", length(panel_genes),
            " current). Retraining.")
  }
}

if (use_cached) {

  singler_trained <- readRDS(cache_singler)
  log_msg("[CACHE] SingleR panel model loaded.")

  rctd_ref <- readRDS(cache_rctd)
  log_msg("[CACHE] RCTD panel components loaded.")
  log_msg("[CACHE]   ", ncol(rctd_ref$counts), " cells, ",
          nrow(rctd_ref$counts), " panel genes, ",
          length(levels(rctd_ref$cell_types)), " types\n")

} else {

  # Load full consensus reference
  ref_path <- file.path(config$reference_dir, "consensus_reference.rds")
  log_msg("[REF] Loading full reference: ", ref_path)
  ref_obj <- readRDS(ref_path)
  log_msg("[REF] Full reference: ", ncol(ref_obj), " cells, ",
          nrow(ref_obj), " genes, ",
          length(unique(ref_obj$consensus_label)), " types")

  # Intersect reference genes with panel genes
  shared_genes <- intersect(rownames(ref_obj), panel_genes)
  panel_only   <- setdiff(panel_genes, rownames(ref_obj))

  log_msg("[REF] Panel genes in reference: ", length(shared_genes), " / ",
          length(panel_genes))
  if (length(panel_only) > 0) {
    log_msg("[REF] Panel genes NOT in reference (", length(panel_only),
            "): ", paste(head(panel_only, 20), collapse = ", "),
            if (length(panel_only) > 20) " ..." else "")
  }

  if (length(shared_genes) < 50) {
    stop("[ERROR] Only ", length(shared_genes),
         " panel genes found in reference. Check gene name format.")
  }

  # Subset reference to panel genes
  log_msg("[REF] Subsetting reference to ", length(shared_genes), " panel genes...")
  ref_panel <- subset(ref_obj, features = shared_genes)

  # Ensure counts layer exists and is clean
  ref_panel <- ensure_clean_assay(ref_panel, "ref_panel")

  # Log-normalize for SingleR training
  log_msg("[REF] Normalizing panel-subset reference...")
  ref_panel <- NormalizeData(ref_panel, verbose = FALSE)

  # Free full reference
  rm(ref_obj); gc(verbose = FALSE)

  # --- Train SingleR on panel genes ---
  log_msg("[SINGLER] Training on ", length(shared_genes),
          " panel genes, ", ncol(ref_panel), " cells...")

  ref_sce <- as.SingleCellExperiment(ref_panel)
  if (!"logcounts" %in% assayNames(ref_sce)) {
    logcounts(ref_sce) <- GetAssayData(ref_panel, assay = "RNA", layer = "data")
  }

  singler_trained <- trainSingleR(
    ref       = ref_sce,
    labels    = ref_sce$consensus_label,
    de.method = config$singler_de_method,
    BPPARAM   = MulticoreParam(config$rctd_max_cores)
  )

  log_msg("[SINGLER] Training complete.")
  saveRDS(singler_trained, cache_singler)
  log_msg("[CACHE] Saved: ", cache_singler)

  rm(ref_sce); gc(verbose = FALSE)

  # --- Build RCTD components on panel genes ---
  log_msg("[RCTD] Building panel-gene reference components...")

  ref_counts <- GetAssayData(ref_panel, assay = "RNA", layer = "counts")
  ref_types  <- setNames(factor(ref_panel$consensus_label), colnames(ref_panel))
  ref_numi   <- setNames(colSums(ref_counts), colnames(ref_panel))

  rctd_ref <- list(
    counts     = ref_counts,
    cell_types = ref_types,
    nUMI       = ref_numi
  )

  log_msg("[RCTD] Panel reference: ", ncol(ref_counts), " cells, ",
          nrow(ref_counts), " genes, ",
          length(levels(ref_types)), " types")

  saveRDS(rctd_ref, cache_rctd)
  log_msg("[CACHE] Saved: ", cache_rctd)

  # Save gene list for cache validation
  writeLines(panel_genes, cache_genes)

  rm(ref_panel, ref_counts); gc(verbose = FALSE)
}

# Sanity check: warn if RCTD CELL_MIN_INSTANCE might drop types
ref_type_counts <- table(rctd_ref$cell_types)
rctd_would_drop <- names(ref_type_counts[ref_type_counts < config$rctd_CELL_MIN_INSTANCE])
if (length(rctd_would_drop) > 0) {
  log_msg("[WARN] RCTD CELL_MIN_INSTANCE = ", config$rctd_CELL_MIN_INSTANCE,
          " will drop ", length(rctd_would_drop), " types:")
  log_msg("       ", paste(rctd_would_drop, collapse = ", "))
  log_msg("       SingleR can still assign these -- systematic disagreement.\n")
}


# ==============================================================================
# HELPER: resolve coordinate column names
# ==============================================================================

resolve_coord_col <- function(meta_colnames, preferred, alternatives) {
  if (preferred %in% meta_colnames) return(preferred)
  for (alt in alternatives) {
    if (alt %in% meta_colnames) return(alt)
  }
  return(preferred)
}


# ==============================================================================
# PER-REGION ANNOTATION
# ==============================================================================

annotate_region <- function(region_name, region_path, config,
                            singler_trained, rctd_ref, panel_genes) {

  log_msg("\n========================================")
  log_msg(" REGION: ", region_name)
  log_msg("========================================\n")

  region_outdir <- file.path(config$output_root, region_name)
  dir.create(region_outdir, showWarnings = FALSE, recursive = TRUE)


  # ============================================================================
  # LOAD PROSEG DATA
  # ============================================================================

  log_msg("[LOAD] Reading proseg outputs from: ", region_path)

  counts    <- readMM(file.path(region_path, "counts.mtx.gz"))
  cell_meta <- fread(file.path(region_path, "cell-metadata.csv.gz"))
  gene_meta <- fread(file.path(region_path, "gene-metadata.csv.gz"))

  # Proseg outputs cells x genes; transpose to genes x cells
  if (nrow(counts) == nrow(cell_meta) && ncol(counts) == nrow(gene_meta)) {
    counts <- t(counts)
  }

  rownames(counts) <- gene_meta$gene
  cell_ids <- paste0(region_name, "_", seq_len(ncol(counts)))
  colnames(counts) <- cell_ids

  meta_df <- as.data.frame(cell_meta)
  rownames(meta_df) <- cell_ids

  xen <- CreateSeuratObject(
    counts    = counts,
    meta.data = meta_df,
    project   = region_name
  )

  log_msg("[LOAD] Raw: ", ncol(xen), " cells, ", nrow(xen), " genes")

  # QC
  xen <- subset(xen,
    nFeature_RNA >= config$xen_min_features &
    nCount_RNA   >= config$xen_min_counts
  )
  log_msg("[QC]   Post-filter: ", ncol(xen), " cells")

  # Normalize
  xen <- NormalizeData(xen, verbose = FALSE)


  # ============================================================================
  # SINGLER
  # ============================================================================

  log_msg("[SINGLER] Classifying ", ncol(xen), " cells...")

  xen_sce <- as.SingleCellExperiment(xen)
  if (!"logcounts" %in% assayNames(xen_sce)) {
    logcounts(xen_sce) <- GetAssayData(xen, layer = "data")
  }

  # Subset test SCE to the panel genes used during training.
  # Since SingleR was trained on panel genes, test must have the same genes.
  # This region may have a subset of the union (if panels differ across regions).
  # Genes in the trained model but not in this region get zero-padded
  # (minimal impact since training selected markers from the panel union,
  # and most/all should be present if panels are consistent).
  test_genes  <- rownames(xen_sce)
  train_genes <- rownames(singler_trained$ref)

  shared <- intersect(train_genes, test_genes)
  missing_in_test <- setdiff(train_genes, test_genes)

  log_msg("[SINGLER] Trained on: ", length(train_genes),
          " | This region: ", length(test_genes),
          " | Shared: ", length(shared))

  if (length(missing_in_test) > 0) {
    log_msg("[SINGLER] Padding ", length(missing_in_test),
            " genes absent from this region (panel union vs region subset)")
    zero_mat <- Matrix(0,
                       nrow = length(missing_in_test),
                       ncol = ncol(xen_sce),
                       sparse = TRUE,
                       dimnames = list(missing_in_test, colnames(xen_sce)))
    logcounts(xen_sce) <- rbind(logcounts(xen_sce), zero_mat)[train_genes, ]
    if ("counts" %in% assayNames(xen_sce)) {
      zero_counts <- Matrix(0, nrow = length(missing_in_test), ncol = ncol(xen_sce),
                            sparse = TRUE, dimnames = list(missing_in_test, colnames(xen_sce)))
      counts(xen_sce) <- rbind(counts(xen_sce), zero_counts)[train_genes, ]
    }
  } else {
    logcounts(xen_sce) <- logcounts(xen_sce)[train_genes, ]
    if ("counts" %in% assayNames(xen_sce)) {
      counts(xen_sce) <- counts(xen_sce)[train_genes, ]
    }
  }

  singler_results <- classifySingleR(
    test      = xen_sce,
    trained   = singler_trained,
    quantile  = config$singler_quantile,
    fine.tune = config$singler_fine_tune,
    BPPARAM   = MulticoreParam(config$rctd_max_cores)
  )

  # Store results
  xen$singler_label      <- singler_results$labels
  xen$singler_pruned     <- singler_results$pruned.labels
  xen$singler_delta_next <- singler_results$delta.next
  xen$singler_max_score  <- apply(singler_results$scores, 1, max)

  # Per-label percentile rank of delta.next
  xen$singler_delta_pctile <- ave(
    xen$singler_delta_next,
    xen$singler_label,
    FUN = function(x) rank(x) / length(x)
  )

  n_pruned <- sum(is.na(xen$singler_pruned))
  log_msg("[SINGLER] Complete. ", n_pruned, " cells pruned by outlier detection.")

  saveRDS(singler_results, file.path(region_outdir, "singler_raw_results.rds"))
  rm(xen_sce, singler_results); gc(verbose = FALSE)


  # ============================================================================
  # RCTD
  # ============================================================================

  log_msg("[RCTD] Running doublet-mode deconvolution...")

  common_genes <- intersect(rownames(xen), rownames(rctd_ref$counts))
  log_msg("[RCTD] Gene intersection: ", length(common_genes), " genes")

  spatial_counts <- GetAssayData(xen, layer = "counts")[common_genes, ]

  # Resolve coordinate columns
  cx <- resolve_coord_col(colnames(xen@meta.data), config$coord_x,
                          c("x_centroid", "X", "centroid_x"))
  cy <- resolve_coord_col(colnames(xen@meta.data), config$coord_y,
                          c("y_centroid", "Y", "centroid_y"))

  if (!(cx %in% colnames(xen@meta.data))) {
    stop("[ERROR] Cannot find x-coordinate column.")
  }
  if (!(cy %in% colnames(xen@meta.data))) {
    stop("[ERROR] Cannot find y-coordinate column.")
  }

  log_msg("[RCTD] Coordinate columns: ", cx, ", ", cy)

  coords <- data.frame(
    x = xen@meta.data[[cx]],
    y = xen@meta.data[[cy]],
    row.names = colnames(xen)
  )

  ref_rctd_obj <- Reference(
    counts     = rctd_ref$counts[common_genes, ],
    cell_types = rctd_ref$cell_types,
    nUMI       = rctd_ref$nUMI
  )

  query_rctd <- SpatialRNA(
    coords = coords,
    counts = spatial_counts,
    nUMI   = setNames(colSums(spatial_counts), colnames(spatial_counts))
  )

  rctd_obj <- create.RCTD(
    spatialRNA        = query_rctd,
    reference         = ref_rctd_obj,
    max_cores         = config$rctd_max_cores,
    gene_cutoff       = config$rctd_gene_cutoff,
    fc_cutoff         = config$rctd_fc_cutoff,
    fc_cutoff_reg     = config$rctd_fc_cutoff_reg,
    UMI_min           = config$rctd_UMI_min,
    counts_MIN        = config$rctd_counts_MIN,
    CELL_MIN_INSTANCE = config$rctd_CELL_MIN_INSTANCE
  )

  rctd_obj <- run.RCTD(rctd_obj, doublet_mode = config$rctd_mode)

  # Extract results
  rctd_df      <- rctd_obj@results$results_df
  rctd_weights <- rctd_obj@results$weights
  rctd_cells   <- rownames(rctd_df)

  # Initialize RCTD columns
  xen$rctd_spot_class   <- NA_character_
  xen$rctd_first_type   <- NA_character_
  xen$rctd_second_type  <- NA_character_
  xen$rctd_first_weight <- NA_real_

  matching  <- intersect(colnames(xen), rctd_cells)
  match_idx <- match(matching, colnames(xen))

  xen$rctd_spot_class[match_idx]  <- as.character(rctd_df[matching, "spot_class"])
  xen$rctd_first_type[match_idx]  <- as.character(rctd_df[matching, "first_type"])
  xen$rctd_second_type[match_idx] <- as.character(rctd_df[matching, "second_type"])

  # Vectorized first_type weight extraction
  if (!is.null(rctd_weights) && nrow(rctd_weights) > 0) {
    weight_cells <- intersect(matching, rownames(rctd_weights))
    if (length(weight_cells) > 0) {
      first_types <- as.character(rctd_df[weight_cells, "first_type"])
      valid <- first_types %in% colnames(rctd_weights)
      if (any(valid)) {
        vc <- weight_cells[valid]
        ft <- first_types[valid]
        weight_vals <- rctd_weights[cbind(vc, ft)]
        xen$rctd_first_weight[match(vc, colnames(xen))] <- weight_vals
      }
    }
  }

  n_rctd_dropped <- ncol(xen) - length(matching)
  if (n_rctd_dropped > 0) {
    log_msg("[RCTD] Dropped ", n_rctd_dropped, " cells (below UMI/gene thresholds)")
  }

  log_msg("[RCTD] Complete. Spot class distribution:")
  rctd_class_tbl <- table(xen$rctd_spot_class, useNA = "ifany")
  for (cl in names(rctd_class_tbl)) {
    log_msg("  ", cl, ": ", rctd_class_tbl[cl])
  }

  rm(rctd_obj, ref_rctd_obj, query_rctd, rctd_weights, rctd_df)
  gc(verbose = FALSE)


  # ============================================================================
  # CONSENSUS SCORING
  # ============================================================================

  log_msg("[CONSENSUS] Applying dual-method confidence filter...")

  # --- SingleR confidence ---
  # Pass if: not pruned, delta >= floor, within-label percentile above threshold.
  # With threshold = 0.95: pass = percentile >= 0.05 (removes bottom 5% per label)
  singler_pass <- (
    !is.na(xen$singler_pruned) &
    !is.na(xen$singler_delta_next) &
    xen$singler_delta_next >= config$singler_delta_floor &
    xen$singler_delta_pctile >= (1 - config$singler_delta_pctile_threshold)
  )
  singler_pass[is.na(singler_pass)] <- FALSE

  # --- RCTD confidence ---
  # Pass if: singlet/doublet_certain AND first_type weight >= threshold
  rctd_pass <- (
    xen$rctd_spot_class %in% config$rctd_confident_classes &
    !is.na(xen$rctd_first_weight) &
    xen$rctd_first_weight >= config$rctd_weight_threshold
  )
  rctd_pass[is.na(rctd_pass)] <- FALSE

  # --- Label agreement ---
  labels_agree <- (!is.na(xen$singler_label) &
                   !is.na(xen$rctd_first_type) &
                   xen$singler_label == xen$rctd_first_type)
  labels_agree[is.na(labels_agree)] <- FALSE

  # --- Consensus: agree AND both confident ---
  consensus_pass <- labels_agree & singler_pass & rctd_pass
  xen$consensus_label <- ifelse(consensus_pass, xen$singler_label, NA_character_)

  # --- Reason codes ---
  xen$unlabeled_reason <- NA_character_
  idx <- !consensus_pass
  if (any(idx)) {
    xen$unlabeled_reason[idx] <- case_when(
      is.na(xen$rctd_first_type[idx])
        ~ "rctd_dropped",
      labels_agree[idx] & !singler_pass[idx] & !rctd_pass[idx]
        ~ "agree_both_low_confidence",
      labels_agree[idx] & !singler_pass[idx] & rctd_pass[idx]
        ~ "agree_singler_low_confidence",
      labels_agree[idx] & singler_pass[idx] & !rctd_pass[idx]
        ~ "agree_rctd_low_confidence",
      !labels_agree[idx] & singler_pass[idx] & rctd_pass[idx]
        ~ "disagree_both_confident",
      !labels_agree[idx] & !singler_pass[idx] & !rctd_pass[idx]
        ~ "disagree_both_low_confidence",
      !labels_agree[idx] & singler_pass[idx] & !rctd_pass[idx]
        ~ "disagree_only_singler_confident",
      !labels_agree[idx] & !singler_pass[idx] & rctd_pass[idx]
        ~ "disagree_only_rctd_confident",
      TRUE ~ "other"
    )
  }

  # --- Summary ---
  n_total     <- ncol(xen)
  n_consensus <- sum(!is.na(xen$consensus_label))
  pct         <- round(100 * n_consensus / n_total, 1)

  log_msg("[CONSENSUS] ", region_name, " results:")
  log_msg("  Total cells:       ", n_total)
  log_msg("  Consensus labeled: ", n_consensus, " (", pct, "%)")
  log_msg("  Unlabeled:         ", n_total - n_consensus)

  reason_tbl <- table(xen$unlabeled_reason, useNA = "ifany")
  log_msg("  Reason breakdown:")
  for (r in names(reason_tbl)) {
    log_msg("    ", r, ": ", reason_tbl[r])
  }

  consensus_tbl <- sort(table(xen$consensus_label), decreasing = TRUE)
  log_msg("  Consensus label composition:")
  for (ct in names(consensus_tbl)) {
    log_msg("    ", ct, ": ", consensus_tbl[ct])
  }


  # ============================================================================
  # EXPORTS
  # ============================================================================

  log_msg("[EXPORT] Writing tables...")

  export_cols <- intersect(
    c("consensus_label", "unlabeled_reason",
      "singler_label", "singler_pruned", "singler_delta_next",
      "singler_delta_pctile", "singler_max_score",
      "rctd_spot_class", "rctd_first_type", "rctd_second_type",
      "rctd_first_weight",
      "nCount_RNA", "nFeature_RNA",
      cx, cy),
    colnames(xen@meta.data)
  )

  fwrite(xen@meta.data[, export_cols, drop = FALSE],
         file.path(region_outdir, "all_cells_annotation.csv.gz"))

  unlabeled_mask <- is.na(xen$consensus_label)
  if (any(unlabeled_mask)) {
    fwrite(xen@meta.data[unlabeled_mask, export_cols, drop = FALSE],
           file.path(region_outdir, "unlabeled_cells.csv.gz"))
  }

  confusion <- table(
    SingleR = xen$singler_label,
    RCTD    = xen$rctd_first_type,
    useNA   = "ifany"
  )
  fwrite(as.data.frame.matrix(confusion),
         file.path(region_outdir, "singler_vs_rctd_confusion.csv"),
         row.names = TRUE)

  disagree_mask <- !labels_agree & !is.na(xen$singler_label) &
                   !is.na(xen$rctd_first_type)
  if (any(disagree_mask)) {
    dtab <- table(
      SingleR = xen$singler_label[disagree_mask],
      RCTD    = xen$rctd_first_type[disagree_mask]
    )
    fwrite(as.data.frame.matrix(dtab),
           file.path(region_outdir, "disagreement_matrix.csv"),
           row.names = TRUE)
  }


  # ============================================================================
  # DIAGNOSTIC PLOTS
  # ============================================================================

  log_msg("[DIAG] Generating plots...")

  has_coords <- (cx %in% colnames(xen@meta.data) &&
                 cy %in% colnames(xen@meta.data))

  if (has_coords) {
    plot_df <- data.frame(
      x      = xen@meta.data[[cx]],
      y      = xen@meta.data[[cy]],
      label  = ifelse(is.na(xen$consensus_label), "Unlabeled",
                      xen$consensus_label),
      reason = xen$unlabeled_reason,
      stringsAsFactors = FALSE
    )

    p1 <- ggplot(plot_df, aes(x, y, color = label)) +
      geom_point(size = 0.1, alpha = 0.5) +
      coord_fixed() + theme_minimal(base_size = 10) +
      ggtitle(paste(region_name, "-- Consensus labels")) +
      theme(legend.key.size = unit(0.3, "cm"))
    ggsave(file.path(region_outdir, "spatial_consensus.pdf"),
           p1, width = 14, height = 10)

    unlabeled_df <- plot_df[plot_df$label == "Unlabeled", ]
    if (nrow(unlabeled_df) > 0) {
      p2 <- ggplot(unlabeled_df, aes(x, y, color = reason)) +
        geom_point(size = 0.2, alpha = 0.6) +
        coord_fixed() + theme_minimal(base_size = 10) +
        ggtitle(paste(region_name, "-- Unlabeled: reason codes"))
      ggsave(file.path(region_outdir, "spatial_unlabeled_reasons.pdf"),
             p2, width = 14, height = 10)
    }
  }

  p3 <- ggplot(xen@meta.data,
    aes(x = singler_label, y = singler_delta_next)) +
    geom_violin(fill = "lightblue", alpha = 0.4, scale = "width") +
    geom_jitter(
      aes(color = ifelse(is.na(consensus_label), "Unlabeled", "Labeled")),
      size = 0.05, alpha = 0.2, width = 0.2
    ) +
    scale_color_manual(values = c(Labeled = "grey30", Unlabeled = "red")) +
    geom_hline(yintercept = config$singler_delta_floor,
               lty = 2, color = "darkred", linewidth = 0.3) +
    coord_flip() + theme_minimal(base_size = 9) +
    ggtitle(paste(region_name, "-- SingleR delta.next by label")) +
    labs(color = "Consensus", y = "delta.next", x = NULL)
  ggsave(file.path(region_outdir, "singler_delta_distribution.pdf"),
         p3, width = 10, height = max(6, 0.4 * length(unique(xen$singler_label))))

  p4 <- ggplot(xen@meta.data,
    aes(x = singler_delta_next, y = rctd_first_weight)) +
    geom_point(
      aes(color = ifelse(is.na(consensus_label), "Unlabeled", "Labeled")),
      size = 0.15, alpha = 0.2
    ) +
    scale_color_manual(values = c(Labeled = "steelblue", Unlabeled = "red")) +
    geom_hline(yintercept = config$rctd_weight_threshold,
               lty = 2, color = "grey40") +
    geom_vline(xintercept = config$singler_delta_floor,
               lty = 2, color = "grey40") +
    annotate("rect",
             xmin = config$singler_delta_floor, xmax = Inf,
             ymin = config$rctd_weight_threshold, ymax = Inf,
             fill = "steelblue", alpha = 0.05) +
    theme_minimal(base_size = 10) +
    ggtitle(paste(region_name, "-- Confidence space")) +
    labs(x = "SingleR delta.next", y = "RCTD first_type weight",
         color = "Consensus")
  ggsave(file.path(region_outdir, "confidence_scatter.pdf"),
         p4, width = 8, height = 7)

  # Save annotated object
  saveRDS(xen, file.path(region_outdir, paste0(region_name, "_annotated.rds")))
  log_msg("[DONE] ", region_name, " complete.\n")

  return(xen)
}


# ==============================================================================
# RUN ALL REGIONS
# ==============================================================================

log_msg("================================================================")
log_msg(" Annotating ", length(xenium_regions), " Xenium regions")
log_msg("================================================================")

annotated <- mapply(
  annotate_region,
  region_name = names(xenium_regions),
  region_path = xenium_regions,
  MoreArgs = list(
    config          = config,
    singler_trained = singler_trained,
    rctd_ref        = rctd_ref,
    panel_genes     = panel_genes
  ),
  SIMPLIFY = FALSE
)


# ==============================================================================
# GLOBAL SUMMARY
# ==============================================================================

log_msg("================================================================")
log_msg(" GLOBAL SUMMARY")
log_msg("================================================================\n")

summary_df <- do.call(rbind, lapply(names(annotated), function(nm) {
  obj <- annotated[[nm]]
  data.frame(
    region        = nm,
    total_cells   = ncol(obj),
    consensus     = sum(!is.na(obj$consensus_label)),
    unlabeled     = sum(is.na(obj$consensus_label)),
    pct_consensus = round(100 * sum(!is.na(obj$consensus_label)) / ncol(obj), 1),
    n_types       = length(unique(na.omit(obj$consensus_label))),
    n_disagree_both_confident = sum(
      obj$unlabeled_reason == "disagree_both_confident", na.rm = TRUE
    ),
    stringsAsFactors = FALSE
  )
}))

log_msg("Per-region summary:")
for (i in seq_len(nrow(summary_df))) {
  r <- summary_df[i, ]
  log_msg("  ", r$region, ": ", r$consensus, "/", r$total_cells,
          " (", r$pct_consensus, "%) consensus, ",
          r$n_types, " types, ",
          r$n_disagree_both_confident, " high-confidence disagreements")
}

fwrite(summary_df, file.path(config$output_root, "global_summary.csv"))

log_msg("\nOutputs in: ", normalizePath(config$output_root))
log_msg("Done.")