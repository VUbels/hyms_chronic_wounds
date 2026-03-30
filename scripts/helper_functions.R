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


################################################################################
# STEP 1: Load proseg h5ad -> Seurat object (counts + metadata, no FOV yet)
################################################################################

load_proseg_h5ad <- function(h5ad_path) {
  require(anndataR)
  require(SeuratObject)
  require(Seurat)
  
  if (!file.exists(h5ad_path)) stop("h5ad not found: ", h5ad_path)
  
  message("Reading h5ad: ", h5ad_path)
  adata <- read_h5ad(h5ad_path)
  seu <- adata$as_Seurat()
  
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
# CONVENIENCE: Full proseg load in one call
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
