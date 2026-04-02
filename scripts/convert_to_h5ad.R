# ==============================================================================
# convert_rds_to_h5ad.R
#
# Source this file and call:
#   convert_annotated_rds(base_dir = "./annotated_data")
#
# Requires: Seurat, SeuratObject, scCustomize, reticulate
# Assumes reticulate is already configured via setup_py_env() or equivalent.
# ==============================================================================

#' Convert annotated Seurat v5 Xenium RDS files to h5ad
#'
#' Recursively finds *_annotated.rds files under base_dir, converts each to
#' h5ad using scCustomize::as.anndata, and writes the output alongside the
#' source rds.
#'
#' @param base_dir Character. Root directory to search. Default "./annotated_data".
#' @param overwrite Logical. Re-convert even if h5ad already exists. Default FALSE.
#' @param main_layer Character. Which Seurat layer maps to adata.X. Default "data".
#' @param other_layers Character vector. Additional layers to transfer. Default "counts".
#' @param default_fov Character. FOV to set as default. Default "proseg".
#'
#' @return Data frame summarising conversion results (path, n_cells, n_genes, status).
convert_annotated_rds <- function(
    base_dir     = "./annotated_data",
    overwrite    = FALSE,
    main_layer   = "counts",
    other_layers = NULL,
    default_fov  = "proseg"
) {
  
  require(Seurat)
  require(SeuratObject)
  require(scCustomize)
  
  stopifnot(dir.exists(base_dir))
  
  rds_files <- list.files(
    path       = base_dir,
    pattern    = "_annotated\\.rds$",
    recursive  = TRUE,
    full.names = TRUE
  )
  
  if (length(rds_files) == 0) {
    stop("No *_annotated.rds files found under: ", base_dir)
  }
  
  message(sprintf("Found %d annotated RDS files", length(rds_files)))
  
  results <- data.frame(
    rds_path  = character(),
    h5ad_path = character(),
    n_cells   = integer(),
    n_genes   = integer(),
    status    = character(),
    stringsAsFactors = FALSE
  )
  
  for (rds_path in rds_files) {
    
    sample_name <- sub("_annotated\\.rds$", "", basename(rds_path))
    out_dir     <- dirname(rds_path)
    h5ad_name   <- paste0(sample_name, "_annotated.h5ad")
    h5ad_path   <- file.path(out_dir, h5ad_name)
    
    # Skip if exists and not overwriting
    if (file.exists(h5ad_path) && !overwrite) {
      message(sprintf("[SKIP] %s", h5ad_path))
      results <- rbind(results, data.frame(
        rds_path = rds_path, h5ad_path = h5ad_path,
        n_cells = NA_integer_, n_genes = NA_integer_,
        status = "skipped", stringsAsFactors = FALSE
      ))
      next
    }
    
    message(sprintf("[LOAD] %s", rds_path))
    obj <- readRDS(rds_path)
    
    # Set default FOV to proseg
    if (default_fov %in% Images(obj)) {
      DefaultFOV(obj) <- default_fov
      message(sprintf("  Default FOV -> '%s'", default_fov))
    }
    
    # Join layers (no-op if already joined, safety for split-layer objects)
    tryCatch({
      obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
    }, error = function(e) {
      message(sprintf("  JoinLayers note: %s", conditionMessage(e)))
    })
    
    # Convert
    message(sprintf("[CONVERT] -> %s", h5ad_path))
    tryCatch({
      as.anndata(
        x                = obj,
        file_path        = out_dir,
        file_name        = h5ad_name,
        assay            = "RNA",
        main_layer       = main_layer,
        other_layers     = other_layers,
        transer_dimreduc = TRUE,
        verbose          = TRUE
      )
      
      results <- rbind(results, data.frame(
        rds_path = rds_path, h5ad_path = h5ad_path,
        n_cells = ncol(obj), n_genes = nrow(obj),
        status = "converted", stringsAsFactors = FALSE
      ))
      message(sprintf("[DONE] %d cells, %d genes", ncol(obj), nrow(obj)))
      
    }, error = function(e) {
      warning(sprintf("[FAIL] %s: %s", sample_name, conditionMessage(e)))
      results <<- rbind(results, data.frame(
        rds_path = rds_path, h5ad_path = h5ad_path,
        n_cells = ncol(obj), n_genes = nrow(obj),
        status = paste0("error: ", conditionMessage(e)),
        stringsAsFactors = FALSE
      ))
    })
    
    rm(obj)
    gc(verbose = FALSE)
  }
  
  message(sprintf(
    "\nDone. %d converted, %d skipped, %d failed.",
    sum(results$status == "converted"),
    sum(results$status == "skipped"),
    sum(grepl("^error", results$status))
  ))
  
  return(results)
}