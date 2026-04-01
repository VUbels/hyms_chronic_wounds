################################################################################
# helper_functions_proseg.R
#
# Proseg + Xenium spatial helper functions for Seurat v5.
#
# WORKFLOW:
#   1. Run proseg -> proseg-output.zarr
#   2. Python: proseg_to_seurat.py -> produces:
#        proseg-anndata.h5ad           (counts + metadata + spatial coords)
#        proseg-seg-vertices.csv.gz    (polygon vertices as x,y,cell table)
#        proseg-centroids.csv.gz       (cell centroids)
#   3. R: load_proseg_h5ad()           (Seurat object with counts + metadata)
#   4. R: attach_proseg_fov()          (native FOV with segmentation + centroids)
#   5. R: attach_xenium_fov()          (Xenium boundaries as second FOV)
#
# After this, ImageDimPlot, ImageFeaturePlot, Crop, DefaultBoundary all work.
#
# REQUIREMENTS:
#   - anndataR + rhdf5 (for h5ad reading)
#   - SeuratObject >= 5.0.0, Seurat >= 5.0.0
#   - arrow, dplyr (for Xenium parquet files)
################################################################################

library(data.table)
library(Seurat)
library(SeuratObject)
library(anndataR)
library(rhdf5)
library(arrow)
library(dplyr)
library(Matrix)
library(spacexr)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(SingleR)
library(ggplot2)
library(sf)
library(patchwork)
library(viridis)
library(BiocParallel) 

################################################################################
# STEP 1: Load proseg h5ad -> Seurat object (counts + metadata, no FOV yet)
################################################################################

load_proseg_h5ad <- function(h5ad_path) {

  if (!file.exists(h5ad_path)) stop("h5ad not found: ", h5ad_path)
  
  message("Reading h5ad: ", h5ad_path)
  adata <- read_h5ad(h5ad_path)
  seu <- adata$as_Seurat()
  
  # Fix layer naming: anndataR stores as "X", Seurat expects "counts"
  available_layers <- Layers(seu[["RNA"]])
  if ("X" %in% available_layers && !"counts" %in% available_layers) {
    seu[["RNA"]]$counts <- seu[["RNA"]]$X
    seu[["RNA"]]$X <- NULL
    message("Renamed layer 'X' -> 'counts'")
  }
  
  # Store spatial coords in metadata if available
  spatial_coords <- adata$obsm[["spatial"]]
  if (!is.null(spatial_coords)) {
    seu$centroid_x <- spatial_coords[, 1]
    seu$centroid_y <- spatial_coords[, 2]
  }
  
  message("Seurat object: ", ncol(seu), " cells, ", nrow(seu), " genes")
  return(seu)
}


################################################################################
# STEP 2: Attach proseg segmentation as native FOV
#
# Reads the vertex CSV produced by proseg_to_seurat.py.
# Format: x, y, cell (one row per vertex, cell name repeated).
# This is exactly what CreateSegmentation.data.frame() expects.
#
# Also reads centroids CSV or computes from metadata.
################################################################################

attach_proseg_fov <- function(seu,
                              vertices_csv_path,
                              centroids_csv_path = NULL,
                              fov_name = "proseg") {
  require(SeuratObject)
  require(Seurat)
  require(data.table)
  
  if (!file.exists(vertices_csv_path)) {
    stop("Vertex CSV not found: ", vertices_csv_path)
  }
  
  # ---- Load vertices ----
  message("Reading polygon vertices: ", vertices_csv_path)
  vtx <- fread(vertices_csv_path)
  
  if (!all(c("x", "y", "cell") %in% colnames(vtx))) {
    stop("Vertex CSV must have columns: x, y, cell")
  }
  
  # Filter to cells in Seurat object
  seurat_cells <- colnames(seu)
  vtx <- vtx[vtx$cell %in% seurat_cells, ]
  n_cells_with_seg <- length(unique(vtx$cell))
  message("  ", nrow(vtx), " vertices for ", n_cells_with_seg, " cells")
  
  # ---- Build Segmentation ----
  seg <- CreateSegmentation(as.data.frame(vtx))
  
  # ---- Build Centroids ----
  if (!is.null(centroids_csv_path) && file.exists(centroids_csv_path)) {
    message("Reading centroids: ", centroids_csv_path)
    cents_df <- fread(centroids_csv_path)
    cents_df <- as.data.frame(cents_df[cents_df$cell %in% seurat_cells, ])
  } else if (all(c("centroid_x", "centroid_y") %in% colnames(seu@meta.data))) {
    cents_df <- data.frame(
      x    = seu$centroid_x,
      y    = seu$centroid_y,
      cell = seurat_cells
    )
  } else {
    stop("No centroids available. Provide centroids_csv_path or ensure ",
         "centroid_x/centroid_y are in metadata (from load_proseg_h5ad).")
  }
  
  cents <- CreateCentroids(cents_df)
  
  # ---- Build FOV ----
  fov <- CreateFOV(
    coords = list("segmentation" = seg, "centroids" = cents),
    type   = c("segmentation", "centroids"),
    assay  = DefaultAssay(seu),
    key    = paste0(fov_name, "_")
  )
  
  # Subset to shared cells
  shared <- intersect(Cells(fov[["segmentation"]]), seurat_cells)
  fov <- subset(fov, cells = shared)
  
  seu[[fov_name]] <- fov
  message("Attached FOV '", fov_name, "' with segmentation + centroids (",
          length(shared), " cells).")
  
  return(seu)
}


################################################################################
# STEP 3: Attach Xenium boundaries as a native FOV
#
# Reads cell_boundaries.parquet from original Xenium output.
# Converts to vertex data.frame and creates a proper FOV.
#
# Cell name mapping: Xenium cell_id -> proseg original_cell_id -> colnames(seu)
################################################################################

attach_xenium_fov <- function(seu,
                              xenium_dir,
                              fov_name = "xenium",
                              boundary_type = "cell") {
  require(arrow)
  require(dplyr)
  require(SeuratObject)
  require(Seurat)
  
  parquet_file <- switch(boundary_type,
                         "cell"    = file.path(xenium_dir, "cell_boundaries.parquet"),
                         "nucleus" = file.path(xenium_dir, "nucleus_boundaries.parquet"),
                         stop("boundary_type must be 'cell' or 'nucleus'")
  )
  if (!file.exists(parquet_file)) stop("Not found: ", parquet_file)
  
  message("Reading Xenium ", boundary_type, " boundaries...")
  bdf <- read_parquet(parquet_file)
  
  required <- c("cell_id", "vertex_x", "vertex_y")
  if (!all(required %in% colnames(bdf))) {
    stop("Parquet must contain: ", paste(required, collapse = ", "))
  }
  
  if ("label_id" %in% colnames(bdf)) {
    bdf <- bdf %>% arrange(cell_id, label_id)
  }
  
  # ---- Map Xenium cell_id -> Seurat colnames ----
  seurat_names <- colnames(seu)
  xenium_ids <- as.character(bdf$cell_id)
  unique_xenium <- unique(xenium_ids)
  
  # Try direct match
  direct <- sum(unique_xenium %in% seurat_names)
  
  if (direct > 0) {
    message("  Direct match: ", direct, " / ", length(unique_xenium), " Xenium cells")
    id_map <- setNames(unique_xenium, unique_xenium)
  } else {
    # Map via original_cell_id in metadata
    meta <- seu@meta.data
    ocid_col <- intersect(c("original_cell_id", "cell_id"), colnames(meta))
    if (length(ocid_col) > 0) {
      id_map <- setNames(seurat_names, as.character(meta[[ocid_col[1]]]))
      n_mapped <- sum(!is.na(id_map[unique_xenium]))
      message("  Mapped via ", ocid_col[1], ": ", n_mapped,
              " / ", length(unique_xenium), " Xenium cells")
    } else {
      stop("No cell name mapping found. Need 'original_cell_id' in metadata.")
    }
  }
  
  # Apply mapping to all vertex rows
  bdf$cell <- as.character(id_map[xenium_ids])
  bdf <- bdf[!is.na(bdf$cell), ]
  
  # ---- Build vertex data.frame for CreateSegmentation ----
  vtx_df <- data.frame(
    x    = bdf$vertex_x,
    y    = bdf$vertex_y,
    cell = bdf$cell
  )
  
  n_cells <- length(unique(vtx_df$cell))
  message("  Building segmentation from ", nrow(vtx_df),
          " vertices (", n_cells, " cells)")
  
  seg <- CreateSegmentation(vtx_df)
  
  # ---- Build Centroids ----
  cells_parquet <- file.path(xenium_dir, "cells.parquet")
  if (file.exists(cells_parquet)) {
    cmeta <- read_parquet(cells_parquet)
    if (all(c("cell_id", "x_centroid", "y_centroid") %in% colnames(cmeta))) {
      cmeta$cell <- as.character(id_map[as.character(cmeta$cell_id)])
      cmeta <- cmeta[!is.na(cmeta$cell) & cmeta$cell %in% seurat_names, ]
      cents_df <- data.frame(
        x    = cmeta$x_centroid,
        y    = cmeta$y_centroid,
        cell = cmeta$cell
      )
    } else {
      cents_df <- NULL
    }
  } else {
    cents_df <- NULL
  }
  
  # Fallback: compute centroids from vertices
  if (is.null(cents_df) || nrow(cents_df) == 0) {
    message("  Computing centroids from vertex means...")
    cents_df <- vtx_df %>%
      group_by(cell) %>%
      summarise(x = mean(x), y = mean(y), .groups = "drop") %>%
      as.data.frame()
  }
  
  cents <- CreateCentroids(cents_df)
  
  # ---- Build FOV ----
  fov <- CreateFOV(
    coords = list("segmentation" = seg, "centroids" = cents),
    type   = c("segmentation", "centroids"),
    assay  = DefaultAssay(seu),
    key    = paste0(fov_name, "_")
  )
  
  shared <- intersect(Cells(fov[["segmentation"]]), seurat_names)
  fov <- subset(fov, cells = shared)
  
  seu[[fov_name]] <- fov
  message("Attached FOV '", fov_name, "' with segmentation + centroids (",
          length(shared), " cells).")
  
  return(seu)
}


################################################################################
# STEP 3b: Attach Xenium nuclei as additional boundary within existing FOV
################################################################################

attach_xenium_nuclei <- function(seu,
                                 xenium_dir,
                                 fov_name = "xenium",
                                 boundary_name = "nuclei") {
  require(arrow)
  require(dplyr)
  require(SeuratObject)
  
  parquet_file <- file.path(xenium_dir, "nucleus_boundaries.parquet")
  if (!file.exists(parquet_file)) stop("Not found: ", parquet_file)
  if (!fov_name %in% Images(seu)) {
    stop("FOV '", fov_name, "' not found. Run attach_xenium_fov() first.")
  }
  
  message("Reading Xenium nucleus boundaries...")
  bdf <- read_parquet(parquet_file)
  
  # ---- Same cell name mapping as attach_xenium_fov ----
  seurat_names <- colnames(seu)
  xenium_ids <- as.character(bdf$cell_id)
  unique_xenium <- unique(xenium_ids)
  
  direct <- sum(unique_xenium %in% seurat_names)
  if (direct > 0) {
    id_map <- setNames(unique_xenium, unique_xenium)
  } else {
    meta <- seu@meta.data
    ocid_col <- intersect(c("original_cell_id", "cell_id"), colnames(meta))
    if (length(ocid_col) > 0) {
      id_map <- setNames(seurat_names, as.character(meta[[ocid_col[1]]]))
    } else {
      stop("No cell name mapping found.")
    }
  }
  
  bdf$cell <- as.character(id_map[as.character(bdf$cell_id)])
  bdf <- bdf[!is.na(bdf$cell) & bdf$cell %in% seurat_names, ]
  
  vtx_df <- data.frame(x = bdf$vertex_x, y = bdf$vertex_y, cell = bdf$cell)
  
  nuclei_seg <- CreateSegmentation(vtx_df)
  seu[[fov_name]][[boundary_name]] <- nuclei_seg
  
  n_cells <- length(unique(vtx_df$cell))
  message("Added '", boundary_name, "' boundary to FOV '", fov_name,
          "' (", n_cells, " cells).")
  return(seu)
}


################################################################################
# CONVENIENCE: Full proseg load in one call - For segmentation comparison
################################################################################

load_proseg_full <- function(h5ad_path,
                             vertices_csv_path,
                             centroids_csv_path = NULL,
                             xenium_dir = NULL,
                             proseg_fov_name = "proseg",
                             xenium_fov_name = "xenium") {
  
  # Load counts + metadata
  seu <- load_proseg_h5ad(h5ad_path)
  
  # Attach proseg segmentation FOV
  seu <- attach_proseg_fov(
    seu, vertices_csv_path,
    centroids_csv_path = centroids_csv_path,
    fov_name = proseg_fov_name
  )
  
  # Optionally attach Xenium boundaries
  if (!is.null(xenium_dir)) {
    seu <- attach_xenium_fov(seu, xenium_dir, fov_name = xenium_fov_name)
    seu <- attach_xenium_nuclei(seu, xenium_dir, fov_name = xenium_fov_name)
  }
  
  return(seu)
}

################################################################################
# LOAD AND OUTPUT A PREFORMATTED SEURAT.RDS AFTER PROSEG COMPUTATION
################################################################################

build_proseg_seurat <- function(proseg_dir,
                                xenium_dir,
                                resolution = 0.7,
                                npcs       = 30,
                                umap_dims  = 1:20,
                                overwrite  = FALSE) {
  
  # Discover all zarr outputs
  zarr_hits <- list.files(
    proseg_dir,
    pattern      = "^proseg-output\\.zarr$",
    recursive    = TRUE,
    full.names   = TRUE,
    include.dirs = TRUE
  )
  
  # zarr sits directly in <proseg_dir>/<sample>/proseg-output.zarr
  sample_dirs <- unique(dirname(zarr_hits))
  
  if (length(sample_dirs) == 0) {
    stop("No proseg-output.zarr found under: ", proseg_dir)
  }
  
  cat("Found", length(sample_dirs), "sample(s)\n")
  
  for (sample_dir in sample_dirs) {
    
    sample_name <- basename(sample_dir)
    rds_path    <- file.path(sample_dir, paste0(sample_name, "_proseg_seurat.rds"))
    
    cat("\n========================================\n")
    cat("Sample:", sample_name, "\n")
    cat("========================================\n")
    
    if (!overwrite && file.exists(rds_path)) {
      cat("  Skipping — RDS already exists:", rds_path, "\n")
      next
    }
    
    zarr_path      <- file.path(sample_dir, "proseg-output.zarr")
    xenium_sample  <- file.path(xenium_dir, sample_name)
    
    # --- Validate ---
    if (!dir.exists(xenium_sample)) {
      warning(sample_name, ": xenium dir not found at ", xenium_sample, " — skipping.")
      next
    }
    
    # --- Python: zarr -> h5ad + vertex CSV ---
    cat("  Converting zarr to h5ad...\n")
    zconv$convert(
      zarr_path  = zarr_path,
      output_dir = sample_dir
    )
    
    h5ad_path      <- file.path(sample_dir, "proseg-anndata.h5ad")
    vertices_path  <- file.path(sample_dir, "proseg-seg-vertices.csv.gz")
    centroids_path <- file.path(sample_dir, "proseg-centroids.csv.gz")
    
    missing <- c(h5ad_path, vertices_path, centroids_path)
    missing <- missing[!file.exists(missing)]
    if (length(missing) > 0) {
      warning(sample_name, ": conversion produced missing files — skipping.\n",
              paste("  ", missing, collapse = "\n"))
      next
    }
    
    # --- Build Seurat with dual FOVs ---
    cat("  Building Seurat object...\n")
    obj <- load_proseg_full(
      h5ad_path          = h5ad_path,
      vertices_csv_path  = vertices_path,
      centroids_csv_path = centroids_path,
      xenium_dir         = xenium_sample,
      proseg_fov_name    = "proseg",
      xenium_fov_name    = "xenium"
    )
    
    obj <- subset(obj, subset = nCount_RNA > 0)
    
    cat("  Cells after subset:", ncol(obj), "\n")
    cat("  FOV proseg cells:",  length(Cells(obj[["proseg"]])), "\n")
    cat("  FOV xenium cells:",  length(Cells(obj[["xenium"]])), "\n")
    
    proseg_fov_cells <- Cells(obj[["proseg"]], boundary = "centroids")
    xenium_fov_cells <- Cells(obj[["xenium"]], boundary = "centroids")
    cat("  proseg FOV cells in Seurat:",
        sum(proseg_fov_cells %in% colnames(obj)), "/",
        length(proseg_fov_cells), "\n")
    cat("  xenium FOV cells in Seurat:",
        sum(xenium_fov_cells %in% colnames(obj)), "/",
        length(xenium_fov_cells), "\n")
    
    # --- Standard Seurat v5 workflow ---
    cat("  Running SCTransform...\n")
    obj <- SCTransform(obj, assay = "RNA", return.only.var.genes = FALSE)
    
    cat("  Running PCA...\n")
    obj <- RunPCA(obj, npcs = npcs, features = rownames(obj))
    
    cat("  Running UMAP...\n")
    obj <- RunUMAP(obj, dims = umap_dims)
    
    cat("  Clustering...\n")
    obj <- FindNeighbors(obj, reduction = "pca", dims = 1:npcs)
    obj <- FindClusters(obj, resolution = resolution)
    
    # --- Save ---
    saveRDS(obj, file = rds_path)
    cat("  Saved:", rds_path, "\n")
    
    rm(obj)
    gc()
  }
}

###############################################################################
# ANNOTATE PROSEG-SEURAT OBJECT BASED ON REFERENCE MODEL
###############################################################################

annotate_proseg_seurat <- function(proseg_dir,
                                   reference_dir,
                                   # Method selection: "both", "singler", "rctd"
                                   annotation_method  = "both",
                                   # QC
                                   xen_min_features = 5,
                                   xen_min_counts   = 10,
                                   coord_x          = "x",
                                   coord_y          = "y",
                                   # SingleR
                                   singler_quantile  = 0.8,
                                   singler_fine_tune = FALSE,
                                   singler_de_method = "wilcox",
                                   # RCTD
                                   rctd_mode              = "doublet",
                                   rctd_max_cores         = 8,
                                   rctd_gene_cutoff       = 0.000125,
                                   rctd_fc_cutoff         = 0.5,
                                   rctd_fc_cutoff_reg     = 0.75,
                                   rctd_UMI_min           = 0,
                                   rctd_counts_MIN        = 0,
                                   rctd_CELL_MIN_INSTANCE = 10,
                                   # Consensus
                                   singler_delta_floor            = 0.05,
                                   singler_delta_pctile_threshold = 0.95,
                                   rctd_weight_threshold          = 0.70,
                                   rctd_confident_classes         = c("singlet", "doublet_certain"),
                                   # Caching
                                   use_cache = TRUE,
                                   overwrite = FALSE) {
  
  # Validate annotation_method
  annotation_method <- match.arg(annotation_method, c("both", "singler", "rctd"))
  run_singler <- annotation_method %in% c("both", "singler")
  run_rctd    <- annotation_method %in% c("both", "rctd")
  
  config <- as.list(environment())
  config <- config[!names(config) %in% c("proseg_dir", "reference_dir", "overwrite")]
  
  ts <- function() format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  
  cat(ts(), "Annotation method:", annotation_method, "\n")
  cat(ts(), "  Run SingleR:", run_singler, "| Run RCTD:", run_rctd, "\n\n")
  
  # ============================================================================
  # DISCOVER BUILT SEURAT OBJECTS
  # ============================================================================
  
  cat(ts(), "================================================================\n")
  cat(ts(), " Discovering built Seurat objects in:", proseg_dir, "\n")
  cat(ts(), "================================================================\n\n")
  
  rds_hits <- list.files(
    proseg_dir,
    pattern    = "_proseg_seurat\\.rds$",
    recursive  = TRUE,
    full.names = TRUE
  )
  
  if (length(rds_hits) == 0) {
    stop("[ERROR] No *_proseg_seurat.rds files found under: ", proseg_dir,
         "\n  Run build_proseg_seurat() first.")
  }
  
  regions <- list()
  for (rds_path in rds_hits) {
    sample_dir  <- dirname(rds_path)
    region_name <- basename(sample_dir)
    regions[[region_name]] <- list(
      dir      = sample_dir,
      rds_path = rds_path
    )
    cat(ts(), "  Found:", region_name, "->", rds_path, "\n")
  }
  
  cat(ts(), "\nDiscovered", length(regions), "regions:",
      paste(names(regions), collapse = ", "), "\n\n")
  
  # ============================================================================
  # PANEL GENES (union across regions)
  # ============================================================================
  
  cat(ts(), "================================================================\n")
  cat(ts(), " Discovering Xenium panel genes\n")
  cat(ts(), "================================================================\n\n")
  
  panel_genes <- character(0)
  for (nm in names(regions)) {
    gm_path <- file.path(regions[[nm]]$dir, "gene-metadata.csv.gz")
    if (file.exists(gm_path)) {
      gm <- fread(gm_path)
      region_genes <- gm$gene
    } else {
      cat(ts(), " ", nm, ": gene-metadata.csv.gz not found, reading RDS header...\n")
      tmp <- readRDS(regions[[nm]]$rds_path)
      region_genes <- rownames(tmp)
      rm(tmp); gc(verbose = FALSE)
    }
    cat(ts(), " ", nm, ":", length(region_genes), "genes\n")
    panel_genes <- union(panel_genes, region_genes)
  }
  
  panel_genes <- sort(unique(panel_genes))
  cat(ts(), "\nPanel gene union:", length(panel_genes), "genes\n\n")
  
  # ============================================================================
  # REFERENCE + PANEL-SPECIFIC TRAINING
  # ============================================================================
  
  cat(ts(), "================================================================\n")
  cat(ts(), " Building panel-specific reference models\n")
  cat(ts(), "================================================================\n\n")
  
  cache_singler <- file.path(proseg_dir, "singler_trained_panel.rds")
  cache_rctd    <- file.path(proseg_dir, "rctd_reference_panel.rds")
  cache_genes   <- file.path(proseg_dir, "panel_genes.txt")
  
  # Determine what's cached
  genes_match <- FALSE
  if (config$use_cache && file.exists(cache_genes)) {
    cached_genes <- readLines(cache_genes)
    genes_match  <- identical(sort(cached_genes), sort(panel_genes))
    if (!genes_match) {
      cat(ts(), "[CACHE] Panel genes changed. Retraining needed models.\n")
    }
  }
  
  singler_cached <- genes_match && file.exists(cache_singler)
  rctd_cached    <- genes_match && file.exists(cache_rctd)
  
  # Load or train SingleR
  singler_trained <- NULL
  if (run_singler) {
    if (singler_cached) {
      singler_trained <- readRDS(cache_singler)
      cat(ts(), "[CACHE] SingleR panel model loaded.\n")
    }
  }
  
  # Load or prepare RCTD ref
  rctd_ref <- NULL
  if (run_rctd) {
    if (rctd_cached) {
      rctd_ref <- readRDS(cache_rctd)
      cat(ts(), "[CACHE] RCTD panel components loaded.\n")
    }
  }
  
  # If anything still needs training, load the reference
  need_singler_train <- run_singler && is.null(singler_trained)
  need_rctd_train    <- run_rctd && is.null(rctd_ref)
  
  if (need_singler_train || need_rctd_train) {
    ref_path <- file.path(reference_dir, "consensus_reference.rds")
    cat(ts(), "[REF] Loading full reference:", ref_path, "\n")
    ref_obj <- readRDS(ref_path)
    cat(ts(), "[REF] Full reference:", ncol(ref_obj), "cells,",
        nrow(ref_obj), "genes,",
        length(unique(ref_obj$consensus_label)), "types\n")
    
    shared_genes <- intersect(rownames(ref_obj), panel_genes)
    panel_only   <- setdiff(panel_genes, rownames(ref_obj))
    
    cat(ts(), "[REF] Panel genes in reference:", length(shared_genes), "/",
        length(panel_genes), "\n")
    if (length(panel_only) > 0) {
      cat(ts(), "[REF] Panel genes NOT in reference (", length(panel_only),
          "):", paste(head(panel_only, 20), collapse = ", "),
          if (length(panel_only) > 20) " ..." else "", "\n")
    }
    if (length(shared_genes) < 50) {
      stop("[ERROR] Only ", length(shared_genes),
           " panel genes found in reference.")
    }
    
    cat(ts(), "[REF] Subsetting reference to", length(shared_genes), "panel genes...\n")
    ref_panel <- subset(ref_obj, features = shared_genes)
    ref_panel <- ensure_clean_assay(ref_panel, "ref_panel")
    ref_panel <- NormalizeData(ref_panel, verbose = FALSE)
    rm(ref_obj); gc(verbose = FALSE)
    
    if (need_singler_train) {
      cat(ts(), "[SINGLER] Training on", length(shared_genes), "panel genes...\n")
      ref_counts_mat <- GetAssayData(ref_panel, assay = "RNA", layer = "counts")
      ref_data_mat   <- GetAssayData(ref_panel, assay = "RNA", layer = "data")
      ref_sce <- SingleCellExperiment(
        assays  = list(counts = ref_counts_mat, logcounts = ref_data_mat),
        colData = ref_panel@meta.data
      )
      rm(ref_counts_mat, ref_data_mat)
      
      singler_trained <- trainSingleR(
        ref       = ref_sce,
        labels    = ref_sce$consensus_label,
        de.method = config$singler_de_method,
        BPPARAM   = MulticoreParam(config$rctd_max_cores)
      )
      cat(ts(), "[SINGLER] Training complete.\n")
      saveRDS(singler_trained, cache_singler)
      rm(ref_sce); gc(verbose = FALSE)
    }
    
    if (need_rctd_train) {
      cat(ts(), "[RCTD] Building panel-gene reference components...\n")
      ref_counts <- GetAssayData(ref_panel, assay = "RNA", layer = "counts")
      ref_types  <- setNames(factor(ref_panel$consensus_label), colnames(ref_panel))
      ref_numi   <- setNames(colSums(ref_counts), colnames(ref_panel))
      
      rctd_ref <- list(
        counts     = ref_counts,
        cell_types = ref_types,
        nUMI       = ref_numi
      )
      
      cat(ts(), "[RCTD] Panel reference:", ncol(ref_counts), "cells,",
          nrow(ref_counts), "genes,",
          length(levels(ref_types)), "types\n")
      saveRDS(rctd_ref, cache_rctd)
      rm(ref_counts); gc(verbose = FALSE)
    }
    
    writeLines(panel_genes, cache_genes)
    rm(ref_panel); gc(verbose = FALSE)
  }
  
  # RCTD drop warning
  if (run_rctd) {
    ref_type_counts <- table(rctd_ref$cell_types)
    rctd_would_drop <- names(ref_type_counts[ref_type_counts < config$rctd_CELL_MIN_INSTANCE])
    if (length(rctd_would_drop) > 0) {
      cat(ts(), "[WARN] RCTD CELL_MIN_INSTANCE =", config$rctd_CELL_MIN_INSTANCE,
          "will drop", length(rctd_would_drop), "types:",
          paste(rctd_would_drop, collapse = ", "), "\n")
    }
  }
  
  # ============================================================================
  # HELPERS
  # ============================================================================
  
  resolve_coord_col <- function(meta_colnames, preferred, alternatives) {
    if (preferred %in% meta_colnames) return(preferred)
    for (alt in alternatives) {
      if (alt %in% meta_colnames) return(alt)
    }
    return(preferred)
  }
  
  # ============================================================================
  # PER-REGION ANNOTATION
  # ============================================================================
  
  annotate_region <- function(region_name, region_info) {
    
    cat(ts(), "\n========================================\n")
    cat(ts(), " REGION:", region_name, "\n")
    cat(ts(), "========================================\n\n")
    
    region_outdir <- file.path(region_info$dir, "annotation")
    dir.create(region_outdir, showWarnings = FALSE, recursive = TRUE)
    
    annotated_rds <- region_info$rds_path
    
    # --- Load ---
    if (!overwrite) {
      tmp <- readRDS(annotated_rds)
      
      if ("SCT" %in% names(tmp@assays)) {
        tmp@assays[["SCT"]] <- NULL
        tmp@active.assay <- "RNA"
      }
      
      if ("consensus_label" %in% colnames(tmp@meta.data)) {
        cat(ts(), "  Skipping -- already annotated:", annotated_rds, "\n")
        
        region_summary <- data.frame(
          region        = region_name,
          total_cells   = ncol(tmp),
          consensus     = sum(!is.na(tmp$consensus_label)),
          unlabeled     = sum(is.na(tmp$consensus_label)),
          pct_consensus = round(100 * sum(!is.na(tmp$consensus_label)) / ncol(tmp), 1),
          n_types       = length(unique(na.omit(tmp$consensus_label))),
          n_disagree_both_confident = if (annotation_method == "both") {
            sum(tmp$unlabeled_reason == "disagree_both_confident", na.rm = TRUE)
          } else { NA_integer_ },
          stringsAsFactors = FALSE
        )
        
        rm(tmp); gc(verbose = FALSE)
        return(region_summary)
      }
      xen <- tmp
      rm(tmp)
    } else {
      xen <- readRDS(annotated_rds)
      if ("SCT" %in% names(xen@assays)) {
        xen@assays[["SCT"]] <- NULL
        xen@active.assay <- "RNA"
      }
    }
    
    cat(ts(), "[LOAD] Loaded:", ncol(xen), "cells,", nrow(xen), "genes\n")
    cat(ts(), "[LOAD] Assays:", paste(names(xen@assays), collapse = ", "), "\n")
    cat(ts(), "[LOAD] Default assay:", xen@active.assay, "\n")
    
    # --- Normalize RNA ---
    
    available_layers <- Layers(xen[["RNA"]])
    if (!"data" %in% available_layers) {
      cat(ts(), "[NORM] No data layer found. Running NormalizeData...\n")
      xen <- NormalizeData(xen, assay = "RNA", verbose = FALSE)
    } else {
      rna_counts <- GetAssayData(xen, assay = "RNA", layer = "counts")
      rna_data   <- GetAssayData(xen, assay = "RNA", layer = "data")
      if (identical(rna_counts[1:min(100, nrow(rna_counts)), 1:min(10, ncol(rna_counts))],
                    rna_data[1:min(100, nrow(rna_counts)), 1:min(10, ncol(rna_counts))])) {
        cat(ts(), "[NORM] Data layer matches counts. Running NormalizeData...\n")
        xen <- NormalizeData(xen, assay = "RNA", verbose = FALSE)
      } else {
        cat(ts(), "[NORM] RNA data layer already log-normalized.\n")
      }
      rm(rna_counts, rna_data); gc(verbose = FALSE)
    }
    
    # ========================================================================
    # SingleR
    # ========================================================================
    
    if (run_singler) {
      cat(ts(), "[SINGLER] Classifying", ncol(xen), "cells...\n")
      
      counts_mat <- GetAssayData(xen, assay = "RNA", layer = "counts")
      data_mat   <- GetAssayData(xen, assay = "RNA", layer = "data")
      
      xen_sce <- SingleCellExperiment(
        assays = list(counts = counts_mat, logcounts = data_mat),
        colData = xen@meta.data
      )
      
      test_genes  <- rownames(xen_sce)
      train_genes <- rownames(singler_trained$ref)
      shared <- intersect(train_genes, test_genes)
      missing_in_test <- setdiff(train_genes, test_genes)
      
      cat(ts(), "[SINGLER] Trained on:", length(train_genes),
          "| This region:", length(test_genes),
          "| Shared:", length(shared), "\n")
      
      # ---- Gene alignment: rebuild SCE from reordered matrices ----
      lc_mat <- logcounts(xen_sce)
      ct_mat <- if ("counts" %in% assayNames(xen_sce)) counts(xen_sce) else NULL
      xen_coldata <- colData(xen_sce)
      
      if (length(missing_in_test) > 0) {
        cat(ts(), "[SINGLER] Padding", length(missing_in_test), "absent genes\n")
        zero_mat <- Matrix(0, nrow = length(missing_in_test), ncol = ncol(xen_sce),
                           sparse = TRUE,
                           dimnames = list(missing_in_test, colnames(xen_sce)))
        lc_mat <- rbind(lc_mat, zero_mat)[train_genes, ]
        if (!is.null(ct_mat)) {
          ct_mat <- rbind(ct_mat, zero_mat)[train_genes, ]
        }
      } else {
        lc_mat <- lc_mat[train_genes, ]
        if (!is.null(ct_mat)) {
          ct_mat <- ct_mat[train_genes, ]
        }
      }
      
      assay_list <- list(logcounts = lc_mat)
      if (!is.null(ct_mat)) assay_list$counts <- ct_mat
      
      xen_sce <- SingleCellExperiment(
        assays  = assay_list,
        colData = xen_coldata
      )
      rm(lc_mat, ct_mat, xen_coldata); gc(verbose = FALSE)
      # ---- End gene alignment ----
      
      singler_results <- classifySingleR(
        test      = xen_sce,
        trained   = singler_trained,
        quantile  = config$singler_quantile,
        fine.tune = config$singler_fine_tune,
        BPPARAM   = MulticoreParam(config$rctd_max_cores)
      )
      
      xen$singler_label      <- singler_results$labels
      xen$singler_pruned     <- singler_results$pruned.labels
      xen$singler_delta_next <- singler_results$delta.next
      xen$singler_max_score  <- apply(singler_results$scores, 1, max)
      xen$singler_delta_pctile <- ave(
        xen$singler_delta_next, xen$singler_label,
        FUN = function(x) rank(x) / length(x)
      )
      
      cat(ts(), "[SINGLER] Complete.", sum(is.na(xen$singler_pruned)), "cells pruned.\n")
      saveRDS(singler_results, file.path(region_outdir, "singler_raw_results.rds"))
      rm(xen_sce, singler_results, counts_mat, data_mat); gc(verbose = FALSE)
    }
    
    # ========================================================================
    # RCTD
    # ========================================================================
    
    if (run_rctd) {
      cat(ts(), "[RCTD] Running", config$rctd_mode, "mode deconvolution...\n")
      
      common_genes <- intersect(rownames(xen), rownames(rctd_ref$counts))
      cat(ts(), "[RCTD] Gene intersection:", length(common_genes), "genes\n")
      
      spatial_counts <- GetAssayData(xen, assay = "RNA", layer = "counts")[common_genes, ]
      
      cx <- resolve_coord_col(colnames(xen@meta.data), config$coord_x,
                              c("x_centroid", "X", "centroid_x"))
      cy <- resolve_coord_col(colnames(xen@meta.data), config$coord_y,
                              c("y_centroid", "Y", "centroid_y"))
      
      if (!(cx %in% colnames(xen@meta.data))) stop("[ERROR] Cannot find x-coordinate column.")
      if (!(cy %in% colnames(xen@meta.data))) stop("[ERROR] Cannot find y-coordinate column.")
      
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
      
      xen$rctd_spot_class   <- NA_character_
      xen$rctd_first_type   <- NA_character_
      xen$rctd_second_type  <- NA_character_
      xen$rctd_first_weight <- NA_real_
      
      if (config$rctd_mode == "full") {
        rctd_weights <- rctd_obj@results$weights
        rctd_cells   <- rownames(rctd_weights)
        matching     <- intersect(colnames(xen), rctd_cells)
        match_idx    <- match(matching, colnames(xen))
        
        row_sums <- rowSums(rctd_weights[matching, , drop = FALSE])
        norm_weights <- rctd_weights[matching, , drop = FALSE] / row_sums
        
        top_type   <- colnames(norm_weights)[apply(norm_weights, 1, which.max)]
        top_weight <- apply(norm_weights, 1, max)
        
        xen$rctd_first_type[match_idx]   <- top_type
        xen$rctd_first_weight[match_idx] <- top_weight
        xen$rctd_spot_class[match_idx]   <- "full"
        
        rm(rctd_weights, norm_weights)
        
      } else if (config$rctd_mode == "doublet") {
        rctd_df      <- rctd_obj@results$results_df
        rctd_weights <- rctd_obj@results$weights
        rctd_cells   <- rownames(rctd_df)
        matching     <- intersect(colnames(xen), rctd_cells)
        match_idx    <- match(matching, colnames(xen))
        
        xen$rctd_first_type[match_idx]   <- as.character(rctd_df[matching, "first_type"])
        xen$rctd_spot_class[match_idx]   <- as.character(rctd_df[matching, "spot_class"])
        xen$rctd_second_type[match_idx]  <- as.character(rctd_df[matching, "second_type"])
        
        if (!is.null(rctd_weights) && nrow(rctd_weights) > 0) {
          weight_cells <- intersect(matching, rownames(rctd_weights))
          if (length(weight_cells) > 0) {
            first_types <- as.character(rctd_df[weight_cells, "first_type"])
            valid <- first_types %in% colnames(rctd_weights)
            if (any(valid)) {
              vc <- weight_cells[valid]
              ft <- first_types[valid]
              xen$rctd_first_weight[match(vc, colnames(xen))] <- rctd_weights[cbind(vc, ft)]
            }
          }
        }
        rm(rctd_df, rctd_weights)
        
      } else {
        # Multi mode
        rctd_weights <- rctd_obj@results$weights
        rctd_cells   <- rownames(rctd_weights)
        matching     <- intersect(colnames(xen), rctd_cells)
        match_idx    <- match(matching, colnames(xen))
        
        row_sums <- rowSums(rctd_weights[matching, , drop = FALSE])
        norm_weights <- rctd_weights[matching, , drop = FALSE] / row_sums
        
        top_type   <- colnames(norm_weights)[apply(norm_weights, 1, which.max)]
        top_weight <- apply(norm_weights, 1, max)
        
        xen$rctd_first_type[match_idx]   <- top_type
        xen$rctd_first_weight[match_idx] <- top_weight
        xen$rctd_spot_class[match_idx]   <- "multi"
        
        rm(rctd_weights, norm_weights)
      }
      
      n_rctd_dropped <- ncol(xen) - length(matching)
      if (n_rctd_dropped > 0) {
        cat(ts(), "[RCTD] Dropped", n_rctd_dropped, "cells (below UMI/gene thresholds)\n")
      }
      
      cat(ts(), "[RCTD] Complete.", length(matching), "cells assigned.\n")
      cat(ts(), "[RCTD] Mode:", config$rctd_mode, "\n")
      
      rm(rctd_obj, ref_rctd_obj, query_rctd, spatial_counts)
      gc(verbose = FALSE)
    }
    
    # ========================================================================
    # CONSENSUS / LABELING
    # ========================================================================
    
    if (annotation_method == "both") {
      # --- Dual-method consensus ---
      cat(ts(), "[CONSENSUS] Applying dual-method confidence filter...\n")
      
      singler_pass <- (
        !is.na(xen$singler_pruned) &
          !is.na(xen$singler_delta_next) &
          xen$singler_delta_next >= config$singler_delta_floor &
          xen$singler_delta_pctile >= (1 - config$singler_delta_pctile_threshold)
      )
      singler_pass[is.na(singler_pass)] <- FALSE
      
      if (config$rctd_mode == "doublet") {
        rctd_pass <- (
          xen$rctd_spot_class %in% config$rctd_confident_classes &
            !is.na(xen$rctd_first_weight) &
            xen$rctd_first_weight >= config$rctd_weight_threshold
        )
      } else {
        rctd_pass <- (
          !is.na(xen$rctd_first_weight) &
            xen$rctd_first_weight >= config$rctd_weight_threshold
        )
      }
      rctd_pass[is.na(rctd_pass)] <- FALSE
      
      labels_agree <- (!is.na(xen$singler_label) &
                         !is.na(xen$rctd_first_type) &
                         xen$singler_label == xen$rctd_first_type)
      labels_agree[is.na(labels_agree)] <- FALSE
      
      consensus_pass <- labels_agree & singler_pass & rctd_pass
      xen$consensus_label <- ifelse(consensus_pass, xen$singler_label, NA_character_)
      
      idx <- !consensus_pass
      xen$unlabeled_reason <- NA_character_
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
      
    } else if (annotation_method == "singler") {
      # --- SingleR only: use pruned labels with delta threshold ---
      cat(ts(), "[LABEL] Applying SingleR confidence filter...\n")
      
      singler_pass <- (
        !is.na(xen$singler_pruned) &
          !is.na(xen$singler_delta_next) &
          xen$singler_delta_next >= config$singler_delta_floor &
          xen$singler_delta_pctile >= (1 - config$singler_delta_pctile_threshold)
      )
      singler_pass[is.na(singler_pass)] <- FALSE
      
      xen$consensus_label <- ifelse(singler_pass, xen$singler_label, NA_character_)
      
      xen$unlabeled_reason <- NA_character_
      idx <- !singler_pass
      if (any(idx)) {
        xen$unlabeled_reason[idx] <- case_when(
          is.na(xen$singler_pruned[idx])
          ~ "singler_pruned",
          !is.na(xen$singler_delta_next[idx]) &
            xen$singler_delta_next[idx] < config$singler_delta_floor
          ~ "singler_low_delta",
          !is.na(xen$singler_delta_pctile[idx]) &
            xen$singler_delta_pctile[idx] < (1 - config$singler_delta_pctile_threshold)
          ~ "singler_low_delta_pctile",
          TRUE ~ "other"
        )
      }
      
    } else {
      # --- RCTD only ---
      cat(ts(), "[LABEL] Applying RCTD confidence filter...\n")
      
      if (config$rctd_mode == "doublet") {
        rctd_pass <- (
          xen$rctd_spot_class %in% config$rctd_confident_classes &
            !is.na(xen$rctd_first_weight) &
            xen$rctd_first_weight >= config$rctd_weight_threshold
        )
      } else {
        rctd_pass <- (
          !is.na(xen$rctd_first_weight) &
            xen$rctd_first_weight >= config$rctd_weight_threshold
        )
      }
      rctd_pass[is.na(rctd_pass)] <- FALSE
      
      xen$consensus_label <- ifelse(rctd_pass, xen$rctd_first_type, NA_character_)
      
      xen$unlabeled_reason <- NA_character_
      idx <- !rctd_pass
      if (any(idx)) {
        xen$unlabeled_reason[idx] <- case_when(
          is.na(xen$rctd_first_type[idx])
          ~ "rctd_dropped",
          !is.na(xen$rctd_first_weight[idx]) &
            xen$rctd_first_weight[idx] < config$rctd_weight_threshold
          ~ "rctd_low_weight",
          config$rctd_mode == "doublet" &
            !(xen$rctd_spot_class[idx] %in% config$rctd_confident_classes)
          ~ "rctd_uncertain_spot_class",
          TRUE ~ "other"
        )
      }
    }
    
    n_total     <- ncol(xen)
    n_consensus <- sum(!is.na(xen$consensus_label))
    pct         <- round(100 * n_consensus / n_total, 1)
    cat(ts(), "[LABEL]", region_name, ":", n_consensus, "/", n_total,
        "(", pct, "%) labeled\n")
    
    reason_tbl <- table(xen$unlabeled_reason, useNA = "ifany")
    cat(ts(), "  Reason breakdown:\n")
    for (r in names(reason_tbl)) {
      cat(ts(), "    ", r, ":", reason_tbl[r], "\n")
    }
    
    consensus_tbl <- sort(table(xen$consensus_label), decreasing = TRUE)
    cat(ts(), "  Label composition:\n")
    for (ct in names(consensus_tbl)) {
      cat(ts(), "    ", ct, ":", consensus_tbl[ct], "\n")
    }
    
    # --- Exports ---
    cat(ts(), "[EXPORT] Writing tables...\n")
    
    cx <- resolve_coord_col(colnames(xen@meta.data), config$coord_x,
                            c("x_centroid", "X", "centroid_x"))
    cy <- resolve_coord_col(colnames(xen@meta.data), config$coord_y,
                            c("y_centroid", "Y", "centroid_y"))
    
    export_cols <- intersect(
      c("consensus_label", "unlabeled_reason",
        "singler_label", "singler_pruned", "singler_delta_next",
        "singler_delta_pctile", "singler_max_score",
        "rctd_spot_class", "rctd_first_type", "rctd_second_type",
        "rctd_first_weight",
        "nCount_RNA", "nFeature_RNA", cx, cy),
      colnames(xen@meta.data)
    )
    
    fwrite(xen@meta.data[, export_cols, drop = FALSE],
           file.path(region_outdir, "all_cells_annotation.csv.gz"))
    
    unlabeled_mask <- is.na(xen$consensus_label)
    if (any(unlabeled_mask)) {
      fwrite(xen@meta.data[unlabeled_mask, export_cols, drop = FALSE],
             file.path(region_outdir, "unlabeled_cells.csv.gz"))
    }
    
    # Confusion matrix only if both methods ran
    if (annotation_method == "both") {
      labels_agree <- (!is.na(xen$singler_label) &
                         !is.na(xen$rctd_first_type) &
                         xen$singler_label == xen$rctd_first_type)
      labels_agree[is.na(labels_agree)] <- FALSE
      
      confusion <- table(SingleR = xen$singler_label,
                         RCTD = xen$rctd_first_type, useNA = "ifany")
      fwrite(as.data.frame.matrix(confusion),
             file.path(region_outdir, "singler_vs_rctd_confusion.csv"),
             row.names = TRUE)
      
      disagree_mask <- !labels_agree & !is.na(xen$singler_label) &
        !is.na(xen$rctd_first_type)
      if (any(disagree_mask)) {
        dtab <- table(SingleR = xen$singler_label[disagree_mask],
                      RCTD = xen$rctd_first_type[disagree_mask])
        fwrite(as.data.frame.matrix(dtab),
               file.path(region_outdir, "disagreement_matrix.csv"),
               row.names = TRUE)
      }
    }
    
    # --- Diagnostic plots ---
    cat(ts(), "[DIAG] Generating plots...\n")
    
    has_coords <- cx %in% colnames(xen@meta.data) && cy %in% colnames(xen@meta.data)
    
    if (has_coords) {
      plot_df <- data.frame(
        x = xen@meta.data[[cx]], y = xen@meta.data[[cy]],
        label = ifelse(is.na(xen$consensus_label), "Unlabeled", xen$consensus_label),
        reason = xen$unlabeled_reason, stringsAsFactors = FALSE
      )
      
      p1 <- ggplot(plot_df, aes(x, y, color = label)) +
        geom_point(size = 0.1, alpha = 0.5) +
        coord_fixed() + theme_minimal(base_size = 10) +
        ggtitle(paste(region_name, "--", annotation_method, "labels")) +
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
    
    # SingleR delta plot (only if SingleR ran)
    if (run_singler) {
      p3 <- ggplot(xen@meta.data, aes(x = singler_label, y = singler_delta_next)) +
        geom_violin(fill = "lightblue", alpha = 0.4, scale = "width") +
        geom_jitter(aes(color = ifelse(is.na(consensus_label), "Unlabeled", "Labeled")),
                    size = 0.05, alpha = 0.2, width = 0.2) +
        scale_color_manual(values = c(Labeled = "grey30", Unlabeled = "red")) +
        geom_hline(yintercept = config$singler_delta_floor, lty = 2, color = "darkred") +
        coord_flip() + theme_minimal(base_size = 9) +
        ggtitle(paste(region_name, "-- SingleR delta.next by label"))
      ggsave(file.path(region_outdir, "singler_delta_distribution.pdf"),
             p3, width = 10, height = max(6, 0.4 * length(unique(xen$singler_label))))
    }
    
    # Confidence scatter (only if both ran)
    if (annotation_method == "both") {
      p4 <- ggplot(xen@meta.data, aes(x = singler_delta_next, y = rctd_first_weight)) +
        geom_point(aes(color = ifelse(is.na(consensus_label), "Unlabeled", "Labeled")),
                   size = 0.15, alpha = 0.2) +
        scale_color_manual(values = c(Labeled = "steelblue", Unlabeled = "red")) +
        geom_hline(yintercept = config$rctd_weight_threshold, lty = 2, color = "grey40") +
        geom_vline(xintercept = config$singler_delta_floor, lty = 2, color = "grey40") +
        annotate("rect",
                 xmin = config$singler_delta_floor, xmax = Inf,
                 ymin = config$rctd_weight_threshold, ymax = Inf,
                 fill = "steelblue", alpha = 0.05) +
        theme_minimal(base_size = 10) +
        ggtitle(paste(region_name, "-- Confidence space"))
      ggsave(file.path(region_outdir, "confidence_scatter.pdf"),
             p4, width = 8, height = 7)
    }
    
    # --- Save annotated object ---
    saveRDS(xen, annotated_rds)
    cat(ts(), "[DONE]", region_name, "->", annotated_rds, "\n")
    
    # --- Build summary, then free object ---
    region_summary <- data.frame(
      region        = region_name,
      total_cells   = ncol(xen),
      consensus     = sum(!is.na(xen$consensus_label)),
      unlabeled     = sum(is.na(xen$consensus_label)),
      pct_consensus = pct,
      n_types       = length(unique(na.omit(xen$consensus_label))),
      n_disagree_both_confident = if (annotation_method == "both") {
        sum(xen$unlabeled_reason == "disagree_both_confident", na.rm = TRUE)
      } else { NA_integer_ },
      stringsAsFactors = FALSE
    )
    
    rm(xen); gc(verbose = FALSE)
    cat("\n")
    
    return(region_summary)
  }
  
  # ============================================================================
  # RUN ALL REGIONS
  # ============================================================================
  
  cat(ts(), "================================================================\n")
  cat(ts(), " Annotating", length(regions), "Xenium regions\n")
  cat(ts(), "================================================================\n")
  
  summaries <- mapply(
    annotate_region,
    region_name = names(regions),
    region_info = regions,
    SIMPLIFY = FALSE
  )
  
  summaries <- Filter(Negate(is.null), summaries)
  
  # ============================================================================
  # GLOBAL SUMMARY
  # ============================================================================
  
  if (length(summaries) > 0) {
    cat(ts(), "================================================================\n")
    cat(ts(), " GLOBAL SUMMARY\n")
    cat(ts(), "================================================================\n\n")
    
    summary_df <- do.call(rbind, summaries)
    fwrite(summary_df, file.path(proseg_dir, "annotation_summary.csv"))
    
    for (i in seq_len(nrow(summary_df))) {
      r <- summary_df[i, ]
      cat(ts(), " ", r$region, ":", r$consensus, "/", r$total_cells,
          "(", r$pct_consensus, "%) labeled,",
          r$n_types, "types")
      if (annotation_method == "both") {
        cat(",", r$n_disagree_both_confident, "high-confidence disagreements")
      }
      cat("\n")
    }
  }
  
  cat(ts(), "\nAnnotation complete. Method:", annotation_method, "\n")
  invisible(summary_df)
}

########################################################
# Ensure RNA counts layer is a clean sparse matrix
########################################################

ensure_clean_assay <- function(obj, label = "object") {
  
  library(Matrix)
  
  counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
  
  if (!inherits(counts, "dgCMatrix")) {
    cat("[", label, "] Converting counts to dgCMatrix\n")
    counts <- as(counts, "dgCMatrix")
  }
  
  if (any(is.na(counts@x))) {
    n_na <- sum(is.na(counts@x))
    cat("[", label, "] Replacing", n_na, "NA values with 0\n")
    counts@x[is.na(counts@x)] <- 0
  }
  
  if (any(counts@x < 0)) {
    n_neg <- sum(counts@x < 0)
    cat("[", label, "] Flooring", n_neg, "negative values to 0\n")
    counts@x[counts@x < 0] <- 0
  }
  
  counts <- drop0(counts)
  
  obj[["RNA"]]$counts <- counts
  return(obj)
}

################################################################################
# LOGGING MESSAGE FUNCTION
################################################################################

log_msg <- function(...) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), paste0(...), "\n")
}

################################################################################
# PLOTTING
#
# Centroid-level: use native ImageDimPlot (works via FOV centroids)
# Polygon-level: use ggplot + geom_sf (works via @misc sf objects)
################################################################################

# Plot polygons from @misc with optional metadata fill
plot_polygons <- function(seu, polygon_slot, fill_by = NULL, line_col = "grey30",
                          line_width = 0.1) {
  require(ggplot2)
  require(sf)
  
  poly_sf <- seu@misc[[polygon_slot]]
  if (is.null(poly_sf)) stop("No polygons in @misc$", polygon_slot)
  
  if (!is.null(fill_by)) {
    if (!fill_by %in% colnames(seu@meta.data)) {
      stop("Metadata column '", fill_by, "' not found.")
    }
    meta_df <- data.frame(
      cell = rownames(seu@meta.data),
      fill_var = seu@meta.data[[fill_by]]
    )
    poly_sf <- merge(poly_sf, meta_df, by = "cell", all.x = TRUE, sort = FALSE)
    poly_sf$fill_var <- as.factor(poly_sf$fill_var)
    
    p <- ggplot(poly_sf) +
      geom_sf(aes(fill = fill_var), color = line_col, linewidth = line_width) +
      theme_void() + coord_sf() +
      labs(fill = fill_by)
  } else {
    p <- ggplot(poly_sf) +
      geom_sf(fill = NA, color = line_col, linewidth = line_width) +
      theme_void() + coord_sf()
  }
  
  return(p)
}


# Plot feature expression on polygons
plot_spatial_feature <- function(seu, polygon_slot, feature,
                                 palette = viridis::viridis(100),
                                 line_col = "grey50", line_width = 0.05) {
  require(ggplot2)
  require(sf)
  
  if (!feature %in% rownames(seu)) stop("Feature '", feature, "' not found.")
  
  poly_sf <- seu@misc[[polygon_slot]]
  if (is.null(poly_sf)) stop("No polygons in @misc$", polygon_slot)
  
  expr <- FetchData(seu, vars = feature)
  expr$cell <- rownames(expr)
  poly_sf <- merge(poly_sf, expr, by = "cell", all.x = TRUE, sort = FALSE)
  
  ggplot(poly_sf) +
    geom_sf(aes(fill = .data[[feature]]), color = line_col, linewidth = line_width) +
    scale_fill_gradientn(colors = palette, na.value = "#eeeeee", name = feature) +
    theme_void() + coord_sf()
}


# Side-by-side segmentation comparison
compare_segmentations <- function(seu,
                                  slot1 = "proseg_polygons",
                                  slot2 = "xenium_polygons",
                                  fill_by = NULL,
                                  labels = c("ProSeg", "Xenium")) {
  require(patchwork)
  
  p1 <- plot_polygons(seu, slot1, fill_by = fill_by) + ggtitle(labels[1])
  p2 <- plot_polygons(seu, slot2, fill_by = fill_by) + ggtitle(labels[2])
  
  p1 + p2
}


# Zoomed comparison
compare_segmentations_zoomed <- function(seu,
                                         slot1 = "proseg_polygons",
                                         slot2 = "xenium_polygons",
                                         x_range, y_range,
                                         fill_by = NULL,
                                         labels = c("ProSeg", "Xenium")) {
  require(patchwork)
  require(sf)
  
  crop_and_plot <- function(slot, label) {
    poly_sf <- seu@misc[[slot]]
    poly_sf <- st_crop(poly_sf,
                       xmin = x_range[1], xmax = x_range[2],
                       ymin = y_range[1], ymax = y_range[2])
    
    if (!is.null(fill_by)) {
      meta_df <- data.frame(
        cell = rownames(seu@meta.data),
        fill_var = as.factor(seu@meta.data[[fill_by]])
      )
      poly_sf <- merge(poly_sf, meta_df, by = "cell", all.x = TRUE, sort = FALSE)
      ggplot(poly_sf) +
        geom_sf(aes(fill = fill_var), color = "grey30", linewidth = 0.2) +
        theme_void() + coord_sf() + ggtitle(label) + labs(fill = fill_by)
    } else {
      ggplot(poly_sf) +
        geom_sf(fill = NA, color = "grey30", linewidth = 0.2) +
        theme_void() + coord_sf() + ggtitle(label)
    }
  }
  
  p1 <- crop_and_plot(slot1, labels[1])
  p2 <- crop_and_plot(slot2, labels[2])
  p1 + p2
}


################################################################################
# NUCLEI-PER-CELL COMPUTATION
################################################################################

compute_nuclei_per_cell <- function(seu,
                                    cell_slot,
                                    nuclei_slot) {
  require(sf)
  
  cells_sf  <- seu@misc[[cell_slot]]
  nuclei_sf <- seu@misc[[nuclei_slot]]
  
  if (is.null(cells_sf))  stop("No polygons in @misc$", cell_slot)
  if (is.null(nuclei_sf)) stop("No polygons in @misc$", nuclei_slot)
  
  cells_sf  <- st_make_valid(cells_sf)
  nuclei_sf <- st_make_valid(nuclei_sf)
  
  message("Computing nuclei per cell...")
  overlaps <- st_within(nuclei_sf, cells_sf, sparse = TRUE)
  
  cell_hits <- sapply(overlaps, function(hit) {
    if (length(hit) == 0) return(NA_integer_)
    hit[1]
  })
  
  counts <- table(na.omit(cell_hits))
  data.frame(
    cell     = cells_sf$cell[as.integer(names(counts))],
    n_nuclei = as.integer(counts)
  )
}


compare_nuclei_segmentation <- function(seu,
                                        cell_slot1, cell_slot2,
                                        nuclei_slot,
                                        labels = c("Segmentation 1", "Segmentation 2")) {
  require(patchwork)
  require(ggplot2)
  
  counts1 <- compute_nuclei_per_cell(seu, cell_slot1, nuclei_slot)
  counts2 <- compute_nuclei_per_cell(seu, cell_slot2, nuclei_slot)
  
  plot_dist <- function(counts_df, title) {
    ggplot(counts_df, aes(x = n_nuclei)) +
      geom_histogram(binwidth = 1, fill = "steelblue", color = "white") +
      theme_minimal() + ggtitle(title) +
      labs(x = "Nuclei per cell", y = "Count")
  }
  
  p1 <- plot_dist(counts1, labels[1])
  p2 <- plot_dist(counts2, labels[2])
  p1 + p2
}


################################################################################
# TWO-PANEL SEGMENTATION WITH NUCLEI OVERLAY
#
# Replaces plot_two_segmentation_panels from the original code.
################################################################################

plot_two_segmentation_panels <- function(seu,
                                         cell_seg_slot1,
                                         cell_seg_slot2,
                                         nuclei_slot,
                                         labels = c("ProSeg", "Xenium"),
                                         nuclei_fill = "#4a90e2",
                                         nuclei_border = "#d9d9d9",
                                         cell_border_color = "black",
                                         cell_fill1 = "#fdbf6f",
                                         cell_fill2 = "#b2df8a",
                                         zoom_bbox = NULL) {
  require(ggplot2)
  require(sf)
  require(patchwork)
  
  cells1 <- seu@misc[[cell_seg_slot1]]
  cells2 <- seu@misc[[cell_seg_slot2]]
  nuclei <- seu@misc[[nuclei_slot]]
  
  if (!is.null(zoom_bbox)) {
    cells1 <- st_crop(cells1, xmin = zoom_bbox[1], xmax = zoom_bbox[2],
                      ymin = zoom_bbox[3], ymax = zoom_bbox[4])
    cells2 <- st_crop(cells2, xmin = zoom_bbox[1], xmax = zoom_bbox[2],
                      ymin = zoom_bbox[3], ymax = zoom_bbox[4])
    nuclei <- st_crop(nuclei, xmin = zoom_bbox[1], xmax = zoom_bbox[2],
                      ymin = zoom_bbox[3], ymax = zoom_bbox[4])
  }
  
  plot_panel <- function(cells, cell_fill, label) {
    ggplot() +
      geom_sf(data = cells, fill = cell_fill,
              color = cell_border_color, alpha = 0.4) +
      geom_sf(data = nuclei, fill = nuclei_fill,
              color = nuclei_border, alpha = 0.6) +
      theme_minimal() + ggtitle(label) + coord_sf(expand = FALSE)
  }
  
  p1 <- plot_panel(cells1, cell_fill1, labels[1])
  p2 <- plot_panel(cells2, cell_fill2, labels[2])
  p1 + p2
}


###############################################################################
# Plot heatmap (pheatmap or ggplot fallback)
###############################################################################

plot_jaccard_heatmap <- function(jaccard_mat, nm_row, nm_col,
                                 pair_tag, output_dir, n_top_markers,
                                 filtered = TRUE) {

  if (filtered) {
    max_per_col <- apply(jaccard_mat, 2, max)
    relevant_cols <- names(max_per_col[max_per_col >= 0.02])
    if (length(relevant_cols) < 3) {
      relevant_cols <- names(sort(max_per_col, decreasing = TRUE))[
        1:min(30, length(max_per_col))]
    }
    plot_mat <- jaccard_mat[, relevant_cols, drop = FALSE]
    suffix   <- ""
  } else {
    plot_mat <- jaccard_mat
    suffix   <- "_full"
  }

  fname <- file.path(output_dir,
                      paste0("label_jaccard_heatmap", suffix, "_", pair_tag, ".pdf"))

  if (requireNamespace("pheatmap", quietly = TRUE)) {
    library(pheatmap)

    pdf(fname,
        width  = max(14, ncol(plot_mat) * 0.25),
        height = max(8, nrow(plot_mat) * 0.5))

    pheatmap(plot_mat,
             color           = colorRampPalette(c("white", "cornflowerblue",
                                                  "darkblue"))(100),
             cluster_rows    = TRUE,
             cluster_cols    = TRUE,
             clustering_method = "ward.D2",
             display_numbers = filtered,
             number_format   = "%.2f",
             number_color    = if (filtered) ifelse(plot_mat > 0.15, "white", "grey30") else "grey30",
             fontsize_number = 7,
             fontsize_row    = 10,
             fontsize_col    = if (filtered) 8 else 6,
             angle_col       = 45,
             main            = paste0("Jaccard: ", nm_row, " (rows) vs ",
                                      nm_col, " (cols)",
                                      if (filtered) paste0("\nTop ", n_top_markers, " markers") else ""),
             border_color    = NA)
    dev.off()

  } else {
    plot_df <- as.data.frame(as.table(plot_mat))
    colnames(plot_df) <- c("type_row", "type_col", "Jaccard")

    p <- ggplot(plot_df, aes(x = type_col, y = type_row, fill = Jaccard)) +
      geom_tile(color = "grey90") +
      geom_text(aes(label = ifelse(Jaccard >= 0.05,
                                    sprintf("%.2f", Jaccard), "")),
                size = 2.5) +
      scale_fill_gradient(low = "white", high = "darkblue") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
            axis.text.y = element_text(size = 10)) +
      labs(title = paste0("Jaccard: ", nm_row, " vs ", nm_col),
           x = nm_col, y = nm_row)

    ggsave(fname, p,
           width = max(14, ncol(plot_mat) * 0.3), height = 10)
  }

  cat("  Saved:", basename(fname), "\n")
}

################################################################################
# HELPER: resolve coordinate column names
################################################################################

resolve_coord_col <- function(meta_colnames, preferred, alternatives) {
  if (preferred %in% meta_colnames) return(preferred)
  for (alt in alternatives) {
    if (alt %in% meta_colnames) return(alt)
  }
  return(preferred)  # fall through, will be caught later
}
