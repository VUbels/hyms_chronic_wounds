#!/usr/bin/env Rscript
# ==============================================================================
# build_reference.R
#
# Build a consensus scRNA-seq / snRNA-seq reference from a mixture of:
#   (A) Pre-processed Seurat .rds objects  (type = "rds")
#   (B) Raw CellRanger/SpaceRanger h5 outputs (type = "h5")
#
# Usage:
#   Rscript build_reference.R                    # use label_mapping.csv if present
#   Rscript build_reference.R --consensus        # auto-harmonize (Jaccard >= 0.4)
#   Rscript build_reference.R --consensus 0.5    # custom threshold
#   Rscript build_reference.R --no-consensus     # keep all labels as-is
#
# Outputs (in output_dir):
#   consensus_reference.rds          Integrated Seurat object
#   singler_trained.rds              Pre-trained SingleR model
#   rctd_reference_components.rds    Counts + cell_types + nUMI for spacexr
#   reference_umap.pdf               Integration QC plots
#   label_mapping_template.csv       Edit -> label_mapping.csv to remap labels
#   reference_type_counts.csv        Cells per type per dataset
#   reference_build_log.txt          Provenance log
#
# Dependencies:
#   install.packages(c("Seurat","SeuratObject","Matrix","data.table",
#                      "ggplot2","patchwork","harmony","hdf5r"))
#   BiocManager::install(c("SingleR","SingleCellExperiment","scDblFinder",
#                          "scuttle","BiocParallel"))
#   devtools::install_github("immunogenomics/presto")  # optional, faster DE
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(SingleR)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(scuttle)
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(harmony)
  library(BiocParallel)
})

if (requireNamespace("presto", quietly = TRUE)) {
  cat("presto detected -- wilcox DE will use fast C++ implementation\n")
} else {
  cat("NOTE: install 'presto' for faster SingleR training\n")
  cat("  devtools::install_github('immunogenomics/presto')\n")
}

# ==============================================================================
# COMMAND-LINE ARGUMENTS
# ==============================================================================
#
# Usage:
#   Rscript build_reference.R                     # default: use label_mapping.csv if present
#   Rscript build_reference.R --consensus          # auto-harmonize labels via Jaccard (J>=0.4)
#   Rscript build_reference.R --consensus 0.5      # custom Jaccard threshold
#   Rscript build_reference.R --no-consensus       # keep all labels as-is (ignore mapping file)
#

cli_args <- commandArgs(trailingOnly = TRUE)

# Defaults
consensus_mode     <- NULL    # NULL = use label_mapping.csv if exists, else as-is
jaccard_threshold  <- 0.4

if (length(cli_args) > 0) {
  i <- 1
  while (i <= length(cli_args)) {
    arg <- cli_args[i]
    if (arg == "--consensus") {
      consensus_mode <- TRUE
      # Check if next arg is a number (optional threshold)
      if (i + 1 <= length(cli_args) && !grepl("^--", cli_args[i + 1])) {
        jaccard_threshold <- as.numeric(cli_args[i + 1])
        if (is.na(jaccard_threshold) || jaccard_threshold < 0 || jaccard_threshold > 1) {
          stop("[ERROR] --consensus threshold must be between 0 and 1")
        }
        i <- i + 1
      }
    } else if (arg == "--no-consensus") {
      consensus_mode <- FALSE
    } else {
      cat("Unknown argument:", arg, "\n")
      cat("Usage: Rscript build_reference.R [--consensus [threshold]] [--no-consensus]\n")
      quit(status = 1)
    }
    i <- i + 1
  }
}

if (isTRUE(consensus_mode)) {
  cat("Mode: CONSENSUS (Jaccard threshold =", jaccard_threshold, ")\n")
} else if (isFALSE(consensus_mode)) {
  cat("Mode: NO-CONSENSUS (all labels kept as-is)\n")
} else {
  cat("Mode: DEFAULT (use label_mapping.csv if present)\n")
}

# ==============================================================================
# CONFIGURATION
# ==============================================================================

config <- list(

  inputs = list(
    HSCA_ref = list(
      type           = "rds",
      path           = "./reference/HSCA_slim.rds",
      annotation_col = "inherited_celltype_lvl_5_extended",
      batch_col      = "sample",
      skip_qc        = TRUE,
      skip_doublets  = TRUE
    ),
    CW_ref = list(
      type           = "rds",
      path           = "./reference/CW_slim.rds",
      annotation_col = "CellType",
      batch_col      = "samples",
      skip_qc        = TRUE,
      skip_doublets  = TRUE
    )
  ),

  min_genes          = 200,
  max_genes          = 8000,
  max_mt_pct         = 30,
  min_cells_per_type = 10,

  integration_method  = "harmony",
  integration_dims    = 1:30,
  n_variable_features = 3000,

  # Normalization: "sctransform" or "lognormalize"
  # sctransform:  better variance stabilization, HIGH memory (~2x object size)
  # lognormalize: lightweight, sufficient for reference building + SingleR/RCTD
  normalization_method = "lognormalize",

  transfer_dims       = 1:30,
  transfer_k_filter   = 200,

  singler_de_method   = "wilcox",
  n_cores             = 8,

  # --- Consensus mode settings ---
  # anchor_source: the input whose labels are the "gold standard".
  # In consensus mode, non-anchor inputs have their labels mapped to
  # the anchor vocabulary where Jaccard >= threshold.
  # Labels below threshold are kept as-is from the original dataset.
  anchor_source = "HSCA_ref",

  output_dir = "./reference/"
)


# ==============================================================================
# INITIALISE
# ==============================================================================

dir.create(config$output_dir, showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(config$output_dir, "reference_build_log.txt")
log_con  <- file(log_file, open = "wt")

# Source shared helpers (log_msg, qc_filter, remove_doublets, ensure_clean_assay,
# load_h5_to_seurat, sample_name_from_h5, get_top_markers)
source("./scripts/helper_functions.R")


# ==============================================================================
# STEP 1 -- LOAD ALL INPUTS
# ==============================================================================

log_msg("========================================================")
log_msg("   build_reference.R")
log_msg("========================================================\n")
log_msg("=== STEP 1: Loading inputs ===\n")

loaded <- list()
annotated_names <- c()

for (nm in names(config$inputs)) {

  inp <- config$inputs[[nm]]
  log_msg("--- ", nm, " (type = ", inp$type, ") ---")

  # ................................................................
  # TYPE: RDS
  # ................................................................
  if (inp$type == "rds") {

    stopifnot("path must exist" = file.exists(inp$path))

    obj <- readRDS(inp$path)
    log_msg(nm, " | Loaded RDS: ", ncol(obj), " cells, ", nrow(obj), " genes")

    obj <- ensure_clean_assay(obj, nm)

    # Annotation
    if (!is.null(inp$annotation_col) &&
        inp$annotation_col %in% colnames(obj@meta.data)) {
      obj$consensus_label <- obj[[inp$annotation_col]][, 1]
      annotated_names <- c(annotated_names, nm)
      log_msg(nm, " | Annotation column '", inp$annotation_col, "': ",
              length(unique(obj$consensus_label)), " types")
    } else if (!is.null(inp$annotation_col)) {
      warning(nm, ": annotation_col '", inp$annotation_col,
              "' not found. Will attempt label transfer.")
      obj$consensus_label <- NA_character_
    } else {
      obj$consensus_label <- NA_character_
    }

    # Batch
    if (!is.null(inp$batch_col) &&
        inp$batch_col %in% colnames(obj@meta.data)) {
      obj$batch <- as.character(obj[[inp$batch_col]][, 1])
    } else {
      obj$batch <- nm
    }

    obj$ref_source <- nm

    # QC
    if (!isTRUE(inp$skip_qc)) {
      obj <- qc_filter(obj, config, nm)
    } else {
      log_msg(nm, " | QC skipped (skip_qc = TRUE)")
      if (!"percent.mt" %in% colnames(obj@meta.data)) {
        obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-|^mt-")
      }
    }

    # Doublets
    if (!isTRUE(inp$skip_doublets)) {
      obj <- remove_doublets(obj, "batch", nm)
    } else {
      log_msg(nm, " | Doublet removal skipped (skip_doublets = TRUE)")
    }

    loaded[[nm]] <- obj

  # ................................................................
  # TYPE: H5
  # ................................................................
  } else if (inp$type == "h5") {

    stopifnot("h5_dir must exist" = dir.exists(inp$h5_dir))

    h5_files <- list.files(inp$h5_dir,
                           pattern = "filtered_feature_bc_matrix\\.h5$",
                           full.names = TRUE, recursive = FALSE)

    if (length(h5_files) == 0) {
      stop("[ERROR] No *filtered_feature_bc_matrix.h5 found in: ", inp$h5_dir)
    }

    log_msg(nm, " | Found ", length(h5_files), " h5 files in ", inp$h5_dir)

    # Auto-detect spatial
    spatial_zips <- list.files(inp$h5_dir, pattern = "spatial_images\\.zip$",
                               full.names = FALSE)
    is_spatial <- isTRUE(inp$is_spatial)
    if (is.null(inp$is_spatial) && length(spatial_zips) > 0) {
      is_spatial <- TRUE
    }

    if (is_spatial) {
      log_msg(nm, " | *** SPATIAL DATA DETECTED ***")
      log_msg(nm, " |     Visium spots -- each barcode treated as one observation.")
    }

    # Load each h5
    h5_objects <- list()
    for (h5 in h5_files) {
      sname <- sample_name_from_h5(basename(h5))
      log_msg(nm, " | Loading: ", basename(h5), " -> '", sname, "'")
      tryCatch({
        h5_obj <- load_h5_to_seurat(h5, sname)
        h5_objects[[sname]] <- h5_obj
        log_msg(nm, " |   ", ncol(h5_obj), " barcodes, ", nrow(h5_obj), " genes")
      }, error = function(e) {
        log_msg(nm, " |   [WARN] Failed: ", conditionMessage(e))
      })
    }

    if (length(h5_objects) == 0) {
      stop("[ERROR] No h5 files loaded for ", nm)
    }

    if (length(h5_objects) == 1) {
      obj <- h5_objects[[1]]
    } else {
      obj <- merge(h5_objects[[1]], y = h5_objects[-1],
                   add.cell.ids = names(h5_objects))
    }

    obj <- ensure_clean_assay(obj, nm)

    obj$ref_source <- nm
    obj$batch      <- obj$sample_name
    obj$is_spatial <- is_spatial

    log_msg(nm, " | Merged: ", ncol(obj), " barcodes, ", nrow(obj), " genes")

    # External metadata
    if (!is.null(inp$metadata_file) && file.exists(inp$metadata_file)) {
      log_msg(nm, " | Reading metadata: ", inp$metadata_file)
      ext_meta <- fread(inp$metadata_file)
      log_msg(nm, " |   ", nrow(ext_meta), " rows, ", ncol(ext_meta), " columns")

      barcode_cols <- intersect(tolower(colnames(ext_meta)),
                                c("barcode", "cell_id", "cellid", "cell_barcode"))

      if (length(barcode_cols) > 0) {
        bc_col <- colnames(ext_meta)[tolower(colnames(ext_meta)) == barcode_cols[1]]
        cell_barcodes <- sub("^.*_", "", colnames(obj))
        idx <- match(cell_barcodes, ext_meta[[bc_col]])
        matched <- sum(!is.na(idx))
        log_msg(nm, " |   Matched ", matched, " / ", ncol(obj), " barcodes")
        if (matched > 0) {
          for (col in setdiff(colnames(ext_meta), bc_col)) {
            obj[[col]] <- ext_meta[[col]][idx]
          }
        }
      } else {
        log_msg(nm, " |   No barcode column -- sample-level metadata only")
      }
    }

    # Annotation
    if (!is.null(inp$annotation_col) &&
        inp$annotation_col %in% colnames(obj@meta.data)) {
      obj$consensus_label <- obj[[inp$annotation_col]][, 1]
      annotated_names <- c(annotated_names, nm)
      log_msg(nm, " | Annotations from '", inp$annotation_col, "': ",
              length(unique(na.omit(obj$consensus_label))), " types")
    } else {
      obj$consensus_label <- NA_character_
      log_msg(nm, " | No annotations. Will label-transfer after integration.")
    }

    obj <- qc_filter(obj, config, nm)

    if (is_spatial) {
      log_msg(nm, " | Doublet removal skipped (spatial data)")
    } else {
      obj <- remove_doublets(obj, "batch", nm)
    }

    loaded[[nm]] <- obj

  } else {
    stop("[ERROR] Unknown input type '", inp$type, "' for ", nm)
  }

  log_msg(nm, " | Final: ", ncol(loaded[[nm]]), " cells\n")
}


# ==============================================================================
# STEP 2 -- LABEL HARMONIZATION
# ==============================================================================

log_msg("=== STEP 2: Label harmonization ===\n")

needs_transfer <- setdiff(names(loaded), annotated_names)

if (length(needs_transfer) > 0 && length(annotated_names) == 0) {
  stop("[ERROR] No annotated inputs. At least one must have annotation_col.")
}

# Generate label template from annotated inputs (always, for reference)
annotated_labels <- do.call(rbind, lapply(annotated_names, function(nm) {
  data.frame(
    source          = nm,
    original_label  = sort(unique(loaded[[nm]]$consensus_label)),
    consensus_label = NA_character_,
    stringsAsFactors = FALSE
  )
}))

template_file  <- file.path(config$output_dir, "label_mapping_template.csv")
completed_file <- file.path(config$output_dir, "label_mapping.csv")
fwrite(annotated_labels, template_file)
log_msg("Label template written: ", template_file)


# --------------------------------------------------------------------------
# MODE A: --consensus (auto Jaccard-based harmonization)
# --------------------------------------------------------------------------
if (isTRUE(consensus_mode)) {

  log_msg("CONSENSUS MODE: auto-harmonizing via Jaccard (threshold = ",
          jaccard_threshold, ")")
  log_msg("  Anchor source: ", config$anchor_source)

  if (!(config$anchor_source %in% annotated_names)) {
    stop("[ERROR] anchor_source '", config$anchor_source,
         "' not found in annotated inputs: ",
         paste(annotated_names, collapse = ", "))
  }

  non_anchor <- setdiff(annotated_names, config$anchor_source)

  if (length(non_anchor) == 0) {
    log_msg("  Only one annotated input -- nothing to harmonize.")
  } else {

    # Normalize each input for marker detection
    log_msg("  Normalizing inputs for marker detection...")
    loaded_norm <- list()
    for (nm in annotated_names) {
      obj_tmp <- NormalizeData(loaded[[nm]], verbose = FALSE)
      loaded_norm[[nm]] <- obj_tmp
    }

    # Get markers for anchor
    log_msg("  Finding markers for anchor: ", config$anchor_source)
    anchor_markers <- get_top_markers(
      loaded_norm[[config$anchor_source]],
      "consensus_label", n_top = 200
    )
    log_msg("    ", length(anchor_markers), " types, median ",
            median(sapply(anchor_markers, length)), " markers/type")

    # For each non-anchor input, compute Jaccard and auto-map
    auto_mapping <- data.frame(
      source          = character(),
      original_label  = character(),
      consensus_label = character(),
      jaccard         = numeric(),
      anchor_match    = character(),
      stringsAsFactors = FALSE
    )

    for (nm in non_anchor) {
      log_msg("\n  Computing Jaccard: ", nm, " vs ", config$anchor_source)

      query_markers <- get_top_markers(
        loaded_norm[[nm]], "consensus_label", n_top = 200
      )
      log_msg("    ", length(query_markers), " types")

      jac_mat <- compute_jaccard(query_markers, anchor_markers)
      log_msg("    Jaccard range: [", round(min(jac_mat), 3), ", ",
              round(max(jac_mat), 3), "]")

      # For each query type, find best anchor match
      for (qtype in rownames(jac_mat)) {
        scores <- jac_mat[qtype, ]
        best_idx   <- which.max(scores)
        best_name  <- names(scores)[best_idx]
        best_score <- scores[best_idx]

        if (best_score >= jaccard_threshold) {
          new_label <- best_name
          flag <- "MAPPED"
        } else {
          new_label <- qtype    # keep original
          flag <- "KEPT"
        }

        auto_mapping <- rbind(auto_mapping, data.frame(
          source          = nm,
          original_label  = qtype,
          consensus_label = new_label,
          jaccard         = round(best_score, 4),
          anchor_match    = best_name,
          stringsAsFactors = FALSE
        ))

        log_msg(sprintf("    [%s] %-20s -> %-30s (J=%.3f, best=%s)",
                        flag, qtype, new_label, best_score, best_name))
      }
    }

    rm(loaded_norm); gc(verbose = FALSE)

    # Save auto-generated mapping
    auto_map_file <- file.path(config$output_dir, "label_mapping_auto.csv")
    fwrite(auto_mapping, auto_map_file)
    log_msg("\n  Auto mapping saved: ", auto_map_file)

    # Also write as label_mapping.csv for reproducibility on re-runs
    map_for_csv <- auto_mapping[, c("source", "original_label", "consensus_label")]
    fwrite(map_for_csv, completed_file)
    log_msg("  Written to: ", completed_file,
            " (will be used on default-mode re-runs)")

    # Apply the mapping
    for (nm in non_anchor) {
      obj <- loaded[[nm]]
      lm  <- auto_mapping[auto_mapping$source == nm, ]

      obj$consensus_label <- lm$consensus_label[
        match(obj$consensus_label, lm$original_label)
      ]

      n_before <- ncol(obj)
      obj <- subset(obj, !is.na(consensus_label))
      if (ncol(obj) < n_before) {
        log_msg(nm, " | Excluded ", n_before - ncol(obj), " cells mapped to NA")
      }

      n_mapped <- sum(lm$original_label != lm$consensus_label)
      n_kept   <- sum(lm$original_label == lm$consensus_label)
      log_msg(nm, " | Applied: ", n_mapped, " types mapped to anchor, ",
              n_kept, " kept original. ",
              length(unique(obj$consensus_label)), " consensus types, ",
              ncol(obj), " cells")

      loaded[[nm]] <- obj
    }

    # Anchor keeps its labels as-is
    log_msg(config$anchor_source, " | Anchor -- labels unchanged (",
            length(unique(loaded[[config$anchor_source]]$consensus_label)),
            " types)")
  }


# --------------------------------------------------------------------------
# MODE B: --no-consensus (all labels as-is, ignore mapping file)
# --------------------------------------------------------------------------
} else if (isFALSE(consensus_mode)) {

  log_msg("NO-CONSENSUS MODE: all labels kept as-is.")
  log_msg("  label_mapping.csv will be ignored even if present.")


# --------------------------------------------------------------------------
# MODE C: default (use label_mapping.csv if present, else as-is)
# --------------------------------------------------------------------------
} else {

  if (file.exists(completed_file)) {

    log_msg("Found label_mapping.csv -- applying harmonization")
    label_map <- fread(completed_file)

    mapped_sources <- unique(label_map$source)
    log_msg("  Sources in mapping file: ", paste(mapped_sources, collapse = ", "))

    for (nm in annotated_names) {
      obj <- loaded[[nm]]

      if (nm %in% mapped_sources) {
        lm <- label_map[source == nm]

        missing <- setdiff(unique(obj$consensus_label), lm$original_label)
        if (length(missing) > 0) {
          stop("[ERROR] ", nm, " has labels not in mapping: ",
               paste(missing, collapse = ", "))
        }

        obj$consensus_label <- lm$consensus_label[
          match(obj$consensus_label, lm$original_label)
        ]

        n_before <- ncol(obj)
        obj <- subset(obj, !is.na(consensus_label))
        if (ncol(obj) < n_before) {
          log_msg(nm, " | Excluded ", n_before - ncol(obj), " cells mapped to NA")
        }

        log_msg(nm, " | Mapped: ", length(unique(obj$consensus_label)),
                " consensus types, ", ncol(obj), " cells")

      } else {
        log_msg(nm, " | No mapping entries -- using original labels as-is (",
                length(unique(obj$consensus_label)), " types)")
      }

      loaded[[nm]] <- obj
    }

  } else {
    log_msg("No label_mapping.csv found. Using original labels as-is.")
    log_msg("  To remap, edit label_mapping_template.csv -> label_mapping.csv",
            " and re-run.")
    log_msg("  Or run with --consensus to auto-harmonize via Jaccard.")
  }
}

# Remove rare cell types
for (nm in annotated_names) {
  obj <- loaded[[nm]]
  type_counts <- table(obj$consensus_label)
  rare <- names(type_counts[type_counts < config$min_cells_per_type])
  if (length(rare) > 0) {
    log_msg(nm, " | Dropping rare types (< ", config$min_cells_per_type,
            " cells): ", paste(rare, collapse = ", "))
    obj <- subset(obj, !(consensus_label %in% rare))
    loaded[[nm]] <- obj
  }
}

for (nm in annotated_names) {
  log_msg(nm, " | ", length(unique(loaded[[nm]]$consensus_label)),
          " consensus types, ", ncol(loaded[[nm]]), " cells")
}


# ==============================================================================
# STEP 3 -- MERGE
# ==============================================================================

log_msg("\n=== STEP 3: Merge and integrate ===\n")

all_objects <- unname(loaded)

if (length(all_objects) == 1) {
  ref_merged <- all_objects[[1]]
} else {
  ref_merged <- merge(all_objects[[1]],
                      y = all_objects[-1],
                      add.cell.ids = names(loaded))
}

# Seurat v5 merge() auto-splits layers -- rebuild assay from counts only
ref_merged[["RNA"]] <- JoinLayers(ref_merged[["RNA"]])

all_layers <- Layers(ref_merged[["RNA"]])
non_counts <- all_layers[!grepl("^counts$", all_layers)]
if (length(non_counts) > 0) {
  log_msg("Rebuilding RNA assay (stripping: ",
          paste(non_counts, collapse = ", "), ")")
  counts_mat <- LayerData(ref_merged[["RNA"]], layer = "counts")
  new_assay  <- CreateAssay5Object(counts = counts_mat)
  ref_merged[["RNA"]] <- new_assay
}
log_msg("RNA layers after cleanup: ",
        paste(Layers(ref_merged[["RNA"]]), collapse = ", "))

# Free memory from input objects before normalization
rm(all_objects, loaded)
gc(verbose = FALSE)

log_msg("Merged: ", ncol(ref_merged), " cells, ", nrow(ref_merged), " genes")
log_msg("Sources: ", paste(unique(ref_merged$ref_source), collapse = ", "))


# ==============================================================================
# STEP 4 -- NORMALIZE + INTEGRATE
# ==============================================================================

log_msg("\nSplitting layers by ref_source for per-batch normalization...")
ref_merged[["RNA"]] <- split(ref_merged[["RNA"]], f = ref_merged$ref_source)

log_msg("Normalization method: ", config$normalization_method)

if (config$normalization_method == "sctransform") {

  # SCTransform: memory-heavy. Force sequential to avoid duplicating the
  # object across future workers. conserve.memory processes genes in chunks.
  options(future.globals.maxSize = 16 * 1024^3)
  future::plan("sequential")

  log_msg("Running SCTransform (sequential, conserve.memory = TRUE)...")
  ref_merged <- SCTransform(ref_merged, verbose = FALSE,
                            conserve.memory = TRUE)

} else {

  # LogNormalize: lightweight, sufficient for Harmony integration and
  # SingleR / RCTD reference building.
  log_msg("Running NormalizeData + FindVariableFeatures + ScaleData...")
  ref_merged <- NormalizeData(ref_merged, verbose = FALSE)
  ref_merged <- FindVariableFeatures(ref_merged,
                                      nfeatures = config$n_variable_features,
                                      verbose = FALSE)
  ref_merged <- ScaleData(ref_merged, verbose = FALSE)

}

log_msg("Running PCA...")
ref_merged <- RunPCA(ref_merged, npcs = 50, verbose = FALSE)

log_msg("Integrating with: ", config$integration_method)

if (config$integration_method == "harmony") {
  ref_merged <- IntegrateLayers(
    ref_merged,
    method         = HarmonyIntegration,
    orig.reduction = "pca",
    new.reduction  = "harmony",
    verbose        = FALSE
  )
  red_use <- "harmony"
} else if (config$integration_method == "rpca") {
  ref_merged <- IntegrateLayers(
    ref_merged,
    method         = RPCAIntegration,
    orig.reduction = "pca",
    new.reduction  = "integrated.rpca",
    verbose        = FALSE
  )
  red_use <- "integrated.rpca"
} else {
  ref_merged <- IntegrateLayers(
    ref_merged,
    method         = CCAIntegration,
    orig.reduction = "pca",
    new.reduction  = "integrated.cca",
    verbose        = FALSE
  )
  red_use <- "integrated.cca"
}

ref_merged[["RNA"]] <- JoinLayers(ref_merged[["RNA"]])

ref_merged <- RunUMAP(ref_merged, reduction = red_use,
                      dims = config$integration_dims, verbose = FALSE)
ref_merged <- FindNeighbors(ref_merged, reduction = red_use,
                            dims = config$integration_dims, verbose = FALSE)

log_msg("Integration complete. Reduction: ", red_use)


# ==============================================================================
# STEP 5 -- LABEL TRANSFER (for unannotated inputs)
# ==============================================================================

if (length(needs_transfer) > 0) {

  log_msg("\n=== STEP 5: Label transfer for unannotated inputs ===\n")

  has_label <- !is.na(ref_merged$consensus_label) &
               ref_merged$consensus_label != ""

  n_unlabeled <- sum(!has_label)
  log_msg("Cells needing transfer: ", n_unlabeled, " / ", ncol(ref_merged))

  if (n_unlabeled > 0 && sum(has_label) > 0) {

    ref_sub   <- subset(ref_merged, cells = colnames(ref_merged)[has_label])
    query_sub <- subset(ref_merged, cells = colnames(ref_merged)[!has_label])

    log_msg("Reference subset: ", ncol(ref_sub), " cells")
    log_msg("Query subset:     ", ncol(query_sub), " cells")

    anchors <- FindTransferAnchors(
      reference           = ref_sub,
      query               = query_sub,
      reference.reduction = red_use,
      dims                = config$transfer_dims,
      k.filter            = config$transfer_k_filter,
      verbose             = FALSE
    )

    predictions <- TransferData(
      anchorset = anchors,
      refdata   = ref_sub$consensus_label,
      dims      = config$transfer_dims,
      verbose   = FALSE
    )

    transfer_labels <- predictions$predicted.id
    transfer_scores <- predictions$prediction.score.max
    names(transfer_labels) <- colnames(query_sub)
    names(transfer_scores) <- colnames(query_sub)

    ref_merged$consensus_label[!has_label] <-
      transfer_labels[colnames(ref_merged)[!has_label]]

    ref_merged$transfer_score <- NA_real_
    ref_merged$transfer_score[!has_label] <-
      transfer_scores[colnames(ref_merged)[!has_label]]
    ref_merged$label_source <- ifelse(has_label, "original", "transfer")

    log_msg("Transfer complete.")
    log_msg("  Min:    ", round(min(transfer_scores, na.rm = TRUE), 3))
    log_msg("  Median: ", round(median(transfer_scores, na.rm = TRUE), 3))
    log_msg("  Mean:   ", round(mean(transfer_scores, na.rm = TRUE), 3))

    low_conf <- sum(transfer_scores < 0.5, na.rm = TRUE)
    if (low_conf > 0) {
      log_msg("  WARNING: ", low_conf, " cells with transfer score < 0.5")
    }

  } else if (sum(has_label) == 0) {
    stop("[ERROR] No cells have annotations.")
  } else {
    log_msg("All cells already labeled. Skipping transfer.")
  }

} else {
  log_msg("\n=== STEP 5: Label transfer -- not needed (all inputs annotated) ===")
  ref_merged$label_source   <- "original"
  ref_merged$transfer_score <- NA_real_
}


# ==============================================================================
# STEP 6 -- POST-TRANSFER QC
# ==============================================================================

log_msg("\n=== STEP 6: Post-transfer QC ===\n")

type_counts <- table(ref_merged$consensus_label)
rare <- names(type_counts[type_counts < config$min_cells_per_type])
if (length(rare) > 0) {
  log_msg("Dropping rare types (< ", config$min_cells_per_type,
          " cells): ", paste(rare, collapse = ", "))
  ref_merged <- subset(ref_merged, !(consensus_label %in% rare))
}

log_msg("Final: ", ncol(ref_merged), " cells, ",
        length(unique(ref_merged$consensus_label)), " types")


# ==============================================================================
# STEP 7 -- EXPORT
# ==============================================================================

log_msg("\n=== STEP 7: Exporting ===\n")

# 7a. Seurat object
ref_rds <- file.path(config$output_dir, "consensus_reference.rds")
saveRDS(ref_merged, ref_rds)
log_msg("Seurat reference        -> ", ref_rds)

# 7b. SingleR
log_msg("Training SingleR model (", config$singler_de_method, " DE)...")

ref_sce <- as.SingleCellExperiment(ref_merged)
if (!"logcounts" %in% assayNames(ref_sce)) {
  # SCTransform stores corrected counts, not log-normalized.
  # LogNormalize stores in the 'data' layer which should map to logcounts.
  # If logcounts is missing regardless of method, generate it.
  ref_tmp <- NormalizeData(ref_merged, assay = "RNA", verbose = FALSE)
  logcounts(ref_sce) <- GetAssayData(ref_tmp, assay = "RNA", layer = "data")
  rm(ref_tmp); gc(verbose = FALSE)
}

singler_trained <- trainSingleR(
  ref       = ref_sce,
  labels    = ref_sce$consensus_label,
  de.method = config$singler_de_method,
  BPPARAM   = MulticoreParam(config$n_cores)
)

singler_rds <- file.path(config$output_dir, "singler_trained.rds")
saveRDS(singler_trained, singler_rds)
log_msg("SingleR trained model   -> ", singler_rds)

# 7c. RCTD
ref_counts <- GetAssayData(ref_merged, assay = "RNA", layer = "counts")
ref_types  <- setNames(factor(ref_merged$consensus_label), colnames(ref_merged))
ref_numi   <- setNames(colSums(ref_counts), colnames(ref_merged))

rctd_components <- list(
  counts     = ref_counts,
  cell_types = ref_types,
  nUMI       = ref_numi
)

rctd_rds <- file.path(config$output_dir, "rctd_reference_components.rds")
saveRDS(rctd_components, rctd_rds)
log_msg("RCTD reference          -> ", rctd_rds)


# ==============================================================================
# STEP 8 -- QC PLOTS
# ==============================================================================

log_msg("\n=== STEP 8: QC plots ===\n")

p1 <- DimPlot(ref_merged, group.by = "ref_source", raster = FALSE) +
  ggtitle("By input source")

p2 <- DimPlot(ref_merged, group.by = "consensus_label", raster = FALSE,
              label = TRUE, repel = TRUE, label.size = 2.5) +
  ggtitle("Consensus labels") + NoLegend()

p3 <- DimPlot(ref_merged, group.by = "label_source", raster = FALSE) +
  ggtitle("Label origin (original vs transfer)")

pdf_out <- file.path(config$output_dir, "reference_umap.pdf")
ggsave(pdf_out, (p1 | p2) / p3, width = 20, height = 14)
log_msg("UMAP plots              -> ", pdf_out)

type_summary <- as.data.frame(table(
  source = ref_merged$ref_source,
  label  = ref_merged$consensus_label
))
type_summary <- type_summary[type_summary$Freq > 0, ]
type_summary <- type_summary[order(type_summary$label, -type_summary$Freq), ]
fwrite(type_summary, file.path(config$output_dir, "reference_type_counts.csv"))

if (any(ref_merged$label_source == "transfer", na.rm = TRUE)) {
  tf_cells <- subset(ref_merged, label_source == "transfer")
  p4 <- VlnPlot(tf_cells, features = "transfer_score",
                 group.by = "consensus_label", pt.size = 0) +
    ggtitle("Label transfer confidence") +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
    NoLegend()
  ggsave(file.path(config$output_dir, "transfer_confidence.pdf"),
         p4, width = 14, height = 6)
  log_msg("Transfer confidence     -> transfer_confidence.pdf")
}


# ==============================================================================
# SUMMARY
# ==============================================================================

log_msg("\n========================================================")
log_msg("   REFERENCE BUILD COMPLETE")
log_msg("========================================================")
log_msg("Total cells:         ", ncol(ref_merged))
log_msg("Total genes:         ", nrow(ref_merged))
log_msg("Input sources:       ", paste(unique(ref_merged$ref_source), collapse = ", "))
log_msg("Consensus types:     ", length(unique(ref_merged$consensus_label)))
log_msg("Original labels:     ", sum(ref_merged$label_source == "original", na.rm = TRUE))
log_msg("Transferred labels:  ", sum(ref_merged$label_source == "transfer", na.rm = TRUE))
log_msg("Integration method:  ", config$integration_method)
log_msg("Normalization:       ", config$normalization_method)
if (isTRUE(consensus_mode)) {
  log_msg("Label mode:          consensus (Jaccard >= ", jaccard_threshold, ")")
  log_msg("Anchor source:       ", config$anchor_source)
} else if (isFALSE(consensus_mode)) {
  log_msg("Label mode:          no-consensus (all labels as-is)")
} else {
  log_msg("Label mode:          default (label_mapping.csv)")
}
log_msg("Output directory:    ", normalizePath(config$output_dir))
log_msg("")
log_msg("Files:")
log_msg("  consensus_reference.rds          Integrated Seurat v5 object")
log_msg("  singler_trained.rds              Pre-trained SingleR model")
log_msg("  rctd_reference_components.rds    counts + cell_types + nUMI")
log_msg("  reference_umap.pdf               Integration QC (3 panels)")
log_msg("  reference_type_counts.csv        Per-source cell counts")
log_msg("  label_mapping_template.csv       Edit -> label_mapping.csv")
log_msg("  reference_build_log.txt          This log")

if (isTRUE(consensus_mode)) {
  log_msg("  label_mapping_auto.csv           Auto Jaccard mapping + scores")
  log_msg("  label_mapping.csv                Generated for default-mode re-runs")
}

if (any(ref_merged$label_source == "transfer", na.rm = TRUE)) {
  log_msg("  transfer_confidence.pdf          Per-type transfer scores")
}

log_msg("========================================================")

close(log_con)

cat("\nReference build complete. Outputs in:",
    normalizePath(config$output_dir), "\n")