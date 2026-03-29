#!/usr/bin/env Rscript
################################################################################
# helper_functions.R
#
# Shared helper functions for the reference building pipeline.
# Sourced by: build_reference.R, label_jaccard_heatmap.R, subcluster_labels.R
#
# Expects the following to be defined before sourcing:
#   - log_con   : file connection for logging (build_reference.R sets this up)
#   - config    : list with QC thresholds (min_genes, max_genes, max_mt_pct)
#
# If log_con does not exist, log_msg falls back to cat() only.
################################################################################

################################################################################
# LOGGING
################################################################################

log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)
  cat(msg, "\n")
  if (exists("log_con") && isOpen(log_con)) {
    writeLines(msg, log_con)
    flush(log_con)
  }
}


################################################################################
# Load a single h5 file into a Seurat object
#
# Handles both CellRanger and SpaceRanger h5 outputs.
# If the h5 contains multiple modalities (e.g. Gene Expression + Antibody
# Capture), only Gene Expression is retained.
################################################################################

load_h5_to_seurat <- function(h5_path, sample_name) {

  counts <- Read10X_h5(h5_path)

  if (is.list(counts)) {
    if ("Gene Expression" %in% names(counts)) {
      counts <- counts[["Gene Expression"]]
    } else {
      counts <- counts[[1]]
    }
  }

  obj <- CreateSeuratObject(
    counts    = counts,
    project   = sample_name,
    min.cells = 3
  )

  obj$sample_name <- sample_name
  obj$orig_h5     <- basename(h5_path)

  return(obj)
}


################################################################################
# Extract a human-readable sample name from h5 file names
#   e.g. "GSM7717016_Skin_filtered_feature_bc_matrix.h5" -> "GSM7717016_Skin"
################################################################################

sample_name_from_h5 <- function(h5_basename) {
  sub("_filtered_feature_bc_matrix\\.h5$", "", h5_basename, ignore.case = TRUE)
}


# ==============================================================================
# Standard QC filter: gene count, UMI count, mito percentage
# ==============================================================================

qc_filter <- function(obj, config, label) {

  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-|^mt-")
  n_before <- ncol(obj)

  obj <- subset(obj,
    nFeature_RNA >= config$min_genes &
    nFeature_RNA <= config$max_genes &
    percent.mt   <= config$max_mt_pct
  )

  log_msg(label, " | QC: ", n_before, " -> ", ncol(obj), " cells")
  return(obj)
}


################################################################################
# Doublet removal via scDblFinder
#
# Falls back to no-batch mode if any batch has < 2 cells.
################################################################################

remove_doublets <- function(obj, batch_field = "sample_name", label = "") {

  sce <- as.SingleCellExperiment(obj)
  set.seed(42)

  batch_vec <- NULL
  if (batch_field %in% colnames(obj@meta.data)) {
    tbl <- table(obj[[batch_field]][, 1])
    if (all(tbl >= 2)) {
      batch_vec <- batch_field
    }
  }

  sce <- scDblFinder(sce, samples = batch_vec)
  obj$scDblFinder_class <- sce$scDblFinder.class

  n_dbl <- sum(obj$scDblFinder_class == "doublet")
  obj   <- subset(obj, scDblFinder_class == "singlet")
  log_msg(label, " | Doublets removed: ", n_dbl,
          " | Remaining: ", ncol(obj))

  return(obj)
}

################################################################################
# Import proseg output into single seurat RDS file
################################################################################

load_proseg <- function(proseg_dir) {
  require(Seurat)
  require(Matrix)
  require(data.table)
  require(dplyr)
  require(sf)
  
  # File paths
  counts_path       <- file.path(proseg_dir, "counts.mtx.gz")
  cell_meta_path    <- file.path(proseg_dir, "cell-metadata.csv.gz")
  gene_meta_path    <- file.path(proseg_dir, "gene-metadata.csv.gz")
  transcript_path   <- file.path(proseg_dir, "transcript-metadata.csv.gz")
  cell_poly_path    <- file.path(proseg_dir, "cell-polygons.geojson.gz")
  
  for (f in c(counts_path, cell_meta_path, gene_meta_path, transcript_path, cell_poly_path)) {
    if (!file.exists(f)) stop("Missing required file: ", f)
  }
  
  # Load counts
  counts <- readMM(gzfile(counts_path))
  counts <- t(counts)  # genes x cells
  
  # Load metadata
  genes <- fread(cmd = paste("zcat", gene_meta_path))
  transcripts <- fread(cmd = paste("zcat", transcript_path))
  cells <- fread(cmd = paste("zcat", cell_meta_path))
  
  # Map gene names
  if ("gene_id" %in% colnames(transcripts) && "gene_name" %in% colnames(transcripts)) {
    gene_map <- transcripts %>% select(gene_id, gene_name)
    rownames(counts) <- make.unique(gene_map$gene_name[match(genes$gene_id, gene_map$gene_id)])
  } else if ("gene" %in% colnames(genes)) {
    rownames(counts) <- make.unique(genes$gene)
  } else stop("Gene/transcript metadata must contain gene_name or gene column")
  
  # Assign cell names
  colnames(counts) <- make.unique(cells$original_cell_id)
  rownames(cells) <- colnames(counts)
  
  # Create Seurat object
  seu <- CreateSeuratObject(
    counts = counts,
    meta.data = as.data.frame(cells)
  )
  
  # Add centroid coordinates if available
  coord_candidates <- list(
    c("x_centroid", "y_centroid"),
    c("x", "y"),
    c("center_x", "center_y")
  )
  for (coords in coord_candidates) {
    if (all(coords %in% colnames(cells))) {
      seu[["x"]] <- cells[[coords[1]]]
      seu[["y"]] <- cells[[coords[2]]]
      break
    }
  }
  
  # Load polygons
  tmp_poly <- tempfile(fileext = ".geojson")
  gz_con <- gzfile(cell_poly_path, "rb")
  writeLines(readLines(gz_con, warn = FALSE), con = tmp_poly)
  close(gz_con)
  cell_polygons <- st_read(tmp_poly, quiet = TRUE)
  unlink(tmp_poly)
  
  # Map polygon cell column to original_cell_id
  if (!"cell_id" %in% colnames(cell_polygons)) {
    if ("cell" %in% colnames(cell_polygons)) {
      # polygon 0-based → original_cell_id
      cell_polygons$cell_id <- cells$original_cell_id[cell_polygons$cell + 1]
    } else stop("Polygons must have 'cell' or 'cell_id' column")
  }
  
  seu@misc$proseg_polygons <- cell_polygons
  message("Stored ", nrow(cell_polygons), " polygons in @misc$proseg_polygons")
  
  message("Seurat object created successfully.")
  return(seu)
}

############################################################################
# PLOT PROSEG CENTROIDS
############################################################################

plot_proseg_centroids <- function(seu, color_by = NULL, point_size = 1) {
  # Dependencies
  require(ggplot2)
  
  # Extract meta.data
  meta <- seu@meta.data
  
  # Check that centroid columns exist
  if (!all(c("centroid_x", "centroid_y") %in% colnames(meta))) {
    stop("Seurat object must contain 'x' and 'y' centroid columns in meta.data")
  }
  
  # Base plot
  p <- ggplot(meta, aes(x = centroid_x, y = centroid_y)) +
    geom_point(aes_string(color = color_by), size = point_size) +
    coord_fixed() +   # ensures aspect ratio is 1:1
    theme_minimal() +
    labs(x = "Centroid X", y = "Centroid Y") +
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10))
  
  return(p)
}

###############################################################################
# STORE PROSEG POLYGONS
###############################################################################

store_proseg_polygons <- function(seu, geojson_gz_file, cell_col = "cell") {
  require(geojsonsf)
  require(sf)
  
  # Read GeoJSON from gzipped file
  poly_sf <- geojson_sf(gzfile(geojson_gz_file))
  
  # Ensure the polygon cell column exists
  if (!cell_col %in% colnames(poly_sf)) {
    stop(sprintf("GeoJSON must have a '%s' property for each polygon", cell_col))
  }
  
  # Map polygon indices to original_cell_id using Seurat metadata
  meta_df <- seu@meta.data
  meta_df$cell_index <- as.character(seq_len(nrow(meta_df)) - 1)  # 0-based index in polygon file
  poly_sf$cell_id <- meta_df$original_cell_id[match(as.character(poly_sf[[cell_col]]),
                                                    meta_df$cell_index)]
  
  # Store the polygon sf object
  seu@misc$proseg_polygons <- poly_sf
  
  return(seu)
}


###############################################################################
# STORE XENIUM POLYGONS (from parquet vertex table)
###############################################################################
store_xenium_polygons <- function(seu, parquet_file) {
  require(arrow)
  require(dplyr)
  require(sf)
  
  message("Reading Xenium cell boundaries (parquet)...")
  cell_boundaries <- read_parquet(parquet_file)
  
  # Basic validation
  required_cols <- c("cell_id", "vertex_x", "vertex_y")
  if (!all(required_cols %in% colnames(cell_boundaries))) {
    stop("Parquet file must contain: cell_id, vertex_x, vertex_y")
  }
  
  message("Constructing polygons...")
  
  # Ensure proper ordering (important for valid polygons)
  cell_boundaries <- cell_boundaries %>%
    arrange(cell_id, label_id)
  
  # Build polygons per cell
  polygons_sf <- cell_boundaries %>%
    group_by(cell_id) %>%
    summarise(
      geometry = list(st_polygon(list(as.matrix(cbind(vertex_x, vertex_y)))))
    ) %>%
    st_as_sf()
  
  # Ensure valid geometries (important!)
  polygons_sf <- st_make_valid(polygons_sf)
  
  # Store in Seurat object
  seu@misc$xenium_polygons <- polygons_sf
  
  message("Stored Xenium polygons in seu@misc$xenium_polygons")
  
  return(seu)
}

###############################################################################
# GENERIC POLYGON PLOTTER
###############################################################################

plot_polygons <- function(seu, polygon_slot, fill_by = NULL, line_col = "black") {
  require(ggplot2)
  require(sf)
  
  poly_sf <- seu@misc[[polygon_slot]]
  
  if (is.null(poly_sf)) {
    stop(sprintf("No polygons found in seu@misc$%s", polygon_slot))
  }
  
  if (!is.null(fill_by)) {
    if (!fill_by %in% colnames(seu@meta.data)) {
      stop(sprintf("Metadata column '%s' not found in Seurat object", fill_by))
    }
    
    meta_df <- seu@meta.data
    meta_df$cell_id <- rownames(meta_df)
    
    # Ensure polygon object has cell_id
    if (!"cell_id" %in% colnames(poly_sf)) {
      stop(sprintf("Polygon slot '%s' must contain a 'cell_id' column", polygon_slot))
    }
    
    poly_sf <- merge(poly_sf, meta_df[, c("cell_id", fill_by)],
                     by = "cell_id", all.x = TRUE, sort = FALSE)
    
    # Force discrete coloring (DimPlot behaviour)
    poly_sf[[fill_by]] <- as.factor(poly_sf[[fill_by]])
  }
  
  ggplot(poly_sf) +
    geom_sf(aes_string(fill = fill_by), color = line_col, size = 0.2, show.legend = TRUE) +
    theme_void() +
    coord_sf()
}

###############################################################################
# SIDE-BY-SIDE SEGMENTATION COMPARISON
###############################################################################
compare_segmentations <- function(seu, fill_by = "cluster") {
  require(patchwork)
  
  p1 <- plot_polygons(seu, "proseg_polygons", fill_by = fill_by) +
    ggtitle("ProSeg Segmentation")
  
  p2 <- plot_polygons(seu, "xenium_polygons", fill_by = fill_by) +
    ggtitle("Xenium Segmentation")
  
  p1 + p2
}

###############################################################################
# STORE XENIUM NUCLEI POLYGONS
###############################################################################
store_xenium_nuclei <- function(seu, parquet_file) {
  require(arrow)
  require(dplyr)
  require(sf)
  
  message("Reading Xenium nucleus boundaries...")
  nuclei_df <- read_parquet(parquet_file)
  
  required_cols <- c("cell_id", "vertex_x", "vertex_y")
  if (!all(required_cols %in% colnames(nuclei_df))) {
    stop("Parquet must contain: cell_id, vertex_x, vertex_y")
  }
  
  nuclei_df <- nuclei_df %>%
    arrange(cell_id)
  
  # Function to safely build polygon
  build_polygon <- function(df) {
    coords <- as.matrix(df[, c("vertex_x", "vertex_y")])
    
    # Need at least 3 points
    if (nrow(coords) < 3) return(NULL)
    
    # Close polygon if needed
    if (!all(coords[1, ] == coords[nrow(coords), ])) {
      coords <- rbind(coords, coords[1, ])
    }
    
    # Try building polygon safely
    tryCatch(
      st_polygon(list(coords)),
      error = function(e) NULL
    )
  }
  
  nuclei_list <- nuclei_df %>%
    group_split(cell_id)
  
  geometries <- lapply(nuclei_list, build_polygon)
  cell_ids   <- sapply(nuclei_list, function(x) x$cell_id[1])
  
  # Remove failed geometries
  valid_idx <- !sapply(geometries, is.null)
  
  nuclei_sf <- st_sf(
    cell_id = cell_ids[valid_idx],
    geometry = st_sfc(geometries[valid_idx])
  )
  
  nuclei_sf <- st_make_valid(nuclei_sf)
  
  seu@misc$xenium_nuclei <- nuclei_sf
  
  message("Stored ", nrow(nuclei_sf), " nucleus polygons")
  
  return(seu)
}
###############################################################################
# GENERIC FUNCTION: COMPUTE NUCLEI PER CELL
###############################################################################

compute_nuclei_per_cell <- function(seu, cell_segmentation_slot, nuclei_segmentation_slot) {
  require(sf)
  require(dplyr)
  
  # Retrieve polygons
  cells  <- seu@misc[[cell_segmentation_slot]]
  nuclei <- seu@misc[[nuclei_segmentation_slot]]
  
  if (!inherits(cells, "sf") || !inherits(nuclei, "sf")) {
    stop("Both cell and nuclei objects must be sf objects stored in misc slot")
  }
  
  # Assign CRS if missing (assume longlat degrees)
  if (is.na(st_crs(cells))) st_crs(cells) <- 4326
  if (is.na(st_crs(nuclei))) st_crs(nuclei) <- 4326
  
  # Transform nuclei to match cells if needed
  if (st_crs(nuclei) != st_crs(cells)) {
    nuclei <- st_transform(nuclei, st_crs(cells))
  }
  
  # Make polygons valid robustly
  cells  <- st_make_valid(cells)
  nuclei <- st_make_valid(nuclei)
  
  # Additional fix: buffer by zero to resolve self-intersections (works in most cases)
  suppressWarnings({
    cells  <- st_buffer(cells, 0)
    nuclei <- st_buffer(nuclei, 0)
  })
  
  message("Computing nuclei per cell...")
  
  # Compute which nuclei fall within which cells
  overlaps <- st_within(nuclei, cells, sparse = TRUE)
  
  # Flatten mapping
  cell_hits <- sapply(overlaps, function(hit) {
    if (length(hit) == 0) return(NA_integer_)
    hit[1]
  })
  
  # Count nuclei per cell
  counts <- table(na.omit(cell_hits))
  
  # Map back to cell IDs
  cell_ids <- rownames(seu@meta.data)[as.integer(names(counts))]
  
  data.frame(
    cell_id = cell_ids,
    n_nuclei = as.integer(counts)
  )
}

###############################################################################
# GENERIC FUNCTION: COMPARE NUCLEI DISTRIBUTION BETWEEN TWO SEGMENTATIONS
###############################################################################
compare_nuclei_segmentation <- function(seu, 
                                        cell_segmentation_slot1, 
                                        cell_segmentation_slot2, 
                                        nuclei_segmentation_slot,
                                        labels = c("Segmentation 1", "Segmentation 2")) {
  require(patchwork)
  
  counts1 <- compute_nuclei_per_cell(seu, cell_segmentation_slot1, nuclei_segmentation_slot)
  counts2 <- compute_nuclei_per_cell(seu, cell_segmentation_slot2, nuclei_segmentation_slot)
  
  p1 <- plot_nuclei_distribution(counts1, title = labels[1])
  p2 <- plot_nuclei_distribution(counts2, title = labels[2])
  
  p1 + p2
}

###############################################################################
# PLOT NUCLEI DISTRIBUTION BETWEEN TWO SEGMENTATIONS
###############################################################################

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
  require(gridExtra)
  
  # Extract sf objects
  cells1 <- seu@misc[[cell_seg_slot1]]
  cells2 <- seu@misc[[cell_seg_slot2]]
  nuclei <- seu@misc[[nuclei_slot]]
  
  # Ensure same CRS
  if (is.na(st_crs(cells1))) st_crs(cells1) <- 3857
  if (is.na(st_crs(cells2))) st_crs(cells2) <- 3857
  if (is.na(st_crs(nuclei))) st_crs(nuclei) <- 3857
  cells2 <- st_transform(cells2, st_crs(cells1))
  nuclei <- st_transform(nuclei, st_crs(cells1))
  
  # Apply zoom if provided
  if (!is.null(zoom_bbox)) {
    xmin <- zoom_bbox[1]; xmax <- zoom_bbox[2]; ymin <- zoom_bbox[3]; ymax <- zoom_bbox[4]
    cells1 <- st_crop(cells1, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
    cells2 <- st_crop(cells2, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
    nuclei <- st_crop(nuclei, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
  }
  
  # Base plotting function
  plot_panel <- function(cells, label) {
    ggplot() +
      geom_sf(data = cells, aes(), fill = ifelse(label == labels[1], cell_fill1, cell_fill2),
              color = cell_border_color, alpha = 0.4) +
      geom_sf(data = nuclei, aes(), fill = nuclei_fill, color = nuclei_border, alpha = 0.6) +
      theme_minimal() +
      ggtitle(label) +
      coord_sf(expand = FALSE)
  }
  
  p1 <- plot_panel(cells1, labels[1])
  p2 <- plot_panel(cells2, labels[2])
  
  gridExtra::grid.arrange(p1, p2, ncol = 2)
}


################################################################################
# ZOOMED IN SEGMENTATION COMPARISON
################################################################################

scout_and_plot_comparison <- function(seu,
                                      cell_seg_slot1,
                                      cell_seg_slot2,
                                      nuclei_slot,
                                      labels = c("ProSeg", "Xenium"),
                                      nuclei_fill = "#4a90e2",
                                      nuclei_border = "#d9d9d9",
                                      cell_border_color = "black",
                                      cell_fill1 = "#fdbf6f",
                                      cell_fill2 = "#b2df8a") {
  require(sf)
  require(leaflet)
  require(leaflet.extras)
  
  # Extract sf objects
  cells1 <- seu@misc[[cell_seg_slot1]]
  cells2 <- seu@misc[[cell_seg_slot2]]
  nuclei <- seu@misc[[nuclei_slot]]
  
  # Ensure CRS
  default_crs <- 4326
  if (is.na(st_crs(cells1))) st_crs(cells1) <- default_crs
  if (is.na(st_crs(cells2))) st_crs(cells2) <- default_crs
  if (is.na(st_crs(nuclei))) st_crs(nuclei) <- default_crs
  cells2 <- st_transform(cells2, st_crs(cells1))
  nuclei <- st_transform(nuclei, st_crs(cells1))
  
  # Open in browser
  options(viewer = NULL)
  
  # Launch interactive map for scouting
  message("Interactive scout: hover to see coordinates. Pan/zoom as needed.")
  leaflet() %>%
    addTiles() %>%
    addPolygons(data = cells1, color = cell_fill1, weight = 1, fillOpacity = 0.3) %>%
    addPolygons(data = cells2, color = cell_fill2, weight = 1, fillOpacity = 0.3) %>%
    addPolygons(data = nuclei, color = nuclei_fill, weight = 0.5, fillOpacity = 0.6) %>%
    addMouseCoordinates()
  
  # Ask user for zoom bounding box
  cat("Enter zoom coordinates as numeric xmin, xmax, ymin, ymax (separated by commas): ")
  coords <- scan(what = numeric(), sep = ",", nmax = 4)
  if (length(coords) != 4) stop("You must enter exactly 4 numeric values: xmin, xmax, ymin, ymax")
  
  # Call plotting function
  plot_two_segmentation_panels(
    seu = seu,
    cell_seg_slot1 = cell_seg_slot1,
    cell_seg_slot2 = cell_seg_slot2,
    nuclei_slot = nuclei_slot,
    labels = labels,
    nuclei_fill = nuclei_fill,
    nuclei_border = nuclei_border,
    cell_border_color = cell_border_color,
    cell_fill1 = cell_fill1,
    cell_fill2 = cell_fill2,
    zoom_bbox = coords
  )
}

################################################################################
# FEATURE PLOT ON CUSTOM SEGMENTATION
################################################################################

plot_spatial_feature <- function(seu,
                                         cell_seg_slot = "proseg_polygons",
                                         feature,
                                         nuclei_slot = NULL,
                                         nuclei_fill = "#4a90e2",
                                         nuclei_border = "#d9d9d9",
                                         cell_fill = "#f2f2f2",
                                         cell_border_color = "black",
                                         feature_palette = viridis::viridis(100)) {
  require(ggplot2)
  require(sf)
  require(dplyr)
  
  # Validate feature
  if (!feature %in% rownames(seu)) {
    stop("Feature '", feature, "' not found in Seurat object assay")
  }
  
  # Get cell polygons
  if (!cell_seg_slot %in% names(seu@misc)) {
    stop("Cell segmentation slot '", cell_seg_slot, "' not found in @misc")
  }
  cells_sf <- seu@misc[[cell_seg_slot]]
  
  if (!"cell_id" %in% colnames(cells_sf)) {
    stop("Cell polygons must have a 'cell_id' column")
  }
  
  # Pull feature expression
  expr <- FetchData(seu, vars = feature)
  expr$cell_id <- rownames(expr)
  
  # Join feature expression to polygons
  cells_sf <- left_join(cells_sf, expr, by = "cell_id")
  
  # Start ggplot
  p <- ggplot()
  
  # Plot cells with feature fill
  p <- p + geom_sf(data = cells_sf, aes(fill = .data[[feature]]), color = cell_border_color, size = 0.2)
  
  # Optionally add nuclei
  if (!is.null(nuclei_slot) && nuclei_slot %in% names(seu@misc)) {
    nuclei_sf <- seu@misc[[nuclei_slot]]
    p <- p + geom_sf(data = nuclei_sf, fill = nuclei_fill, color = nuclei_border, size = 0.1, alpha = 0.6)
  }
  
  # Feature palette
  p <- p + scale_fill_gradientn(colors = feature_palette, na.value = "#eeeeee", name = feature)
  
  p <- p + theme_void() + theme(legend.position = "right")
  
  return(p)
}


# ==============================================================================
# Fast marker detection (presto if available, else FindAllMarkers)
#
# Returns a named list: cell_type -> character vector of gene names
# ==============================================================================

get_top_markers <- function(obj, label_col, n_top = 200,
                            min_pct = 0.1, logfc_thresh = 0.25) {

  Idents(obj) <- label_col

  if (requireNamespace("presto", quietly = TRUE)) {
    cat("  Using presto for marker detection...\n")

    expr_mat <- GetAssayData(obj, layer = "data")
    labels   <- obj@meta.data[[label_col]]

    res <- presto::wilcoxauc(expr_mat, labels)
    res <- data.table::as.data.table(res)

    res <- res[logFC > 0 & padj < 0.05]
    res <- res[order(group, -auc)]
    markers <- res[, head(.SD, n_top), by = group]

    marker_list <- split(markers$feature, markers$group)

  } else {
    cat("  presto not available. Using FindAllMarkers (slow)...\n")

    all_markers <- FindAllMarkers(obj,
                                  only.pos  = TRUE,
                                  min.pct   = min_pct,
                                  logfc.threshold = logfc_thresh,
                                  test.use  = "wilcox",
                                  verbose   = FALSE)
    all_markers <- data.table::as.data.table(all_markers)
    all_markers <- all_markers[order(cluster, -avg_log2FC)]
    top <- all_markers[, head(.SD, n_top), by = cluster]

    marker_list <- split(top$gene, top$cluster)
  }

  return(marker_list)
}


# ==============================================================================
# Compute Jaccard matrix between two marker lists
# ==============================================================================

compute_jaccard <- function(markers_row, markers_col) {
  types_row <- sort(names(markers_row))
  types_col <- sort(names(markers_col))

  mat <- matrix(0, nrow = length(types_row), ncol = length(types_col),
                dimnames = list(types_row, types_col))

  for (tr in types_row) {
    for (tc in types_col) {
      g_row <- markers_row[[tr]]
      g_col <- markers_col[[tc]]
      inter <- length(intersect(g_row, g_col))
      uni   <- length(union(g_row, g_col))
      mat[tr, tc] <- if (uni > 0) inter / uni else 0
    }
  }
  return(mat)
}


# ==============================================================================
# Plot heatmap (pheatmap or ggplot fallback)
# ==============================================================================

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

# ==============================================================================
# HELPER: resolve coordinate column names
# ==============================================================================

resolve_coord_col <- function(meta_colnames, preferred, alternatives) {
  if (preferred %in% meta_colnames) return(preferred)
  for (alt in alternatives) {
    if (alt %in% meta_colnames) return(alt)
  }
  return(preferred)  # fall through, will be caught later
}
