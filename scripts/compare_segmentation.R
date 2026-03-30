source("./scripts/helper_functions.R")
source("./scripts/setup_py_environment.R")

# ---- Python: proseg zarr -> h5ad + vertex CSV ----
setup_py_env("hyms_chronic_wounds", "/home/uvictor/miniconda3/condabin/conda")

reticulate::py_run_string("
import sys
sys.path.append('./scripts')
")

zconv <- reticulate::import("zarr_h5ad_conversion", convert = TRUE)

zconv$convert(
  zarr_path  = "./proseg_results/HS_Region_1/python/proseg-output.zarr",
  output_dir = "./proseg_results/HS_Region_1/seurat_ready/"
)


# ---- R: Build Seurat with native FOVs ----
seurat_proseg <- load_proseg_full(
  h5ad_path         = "./proseg_results/HS_Region_1/seurat_ready/proseg-anndata.h5ad",
  vertices_csv_path = "./proseg_results/HS_Region_1/seurat_ready/proseg-seg-vertices.csv.gz",
  centroids_csv_path = "./proseg_results/HS_Region_1/seurat_ready/proseg-centroids.csv.gz",
  xenium_dir        = "/mnt/d/HYMS/chronic_wounds/HS_Region_1",
  proseg_fov_name   = "proseg",
  xenium_fov_name   = "xenium"
)

# ---- Native Seurat functions now work ----
ImageDimPlot(seurat_proseg, fov = "proseg")
ImageDimPlot(seurat_proseg, fov = "xenium")
#ImageFeaturePlot(seurat_proseg, features = "SERPINB7", fov = "proseg")

# Zoomed with polygon boundaries visible
DefaultBoundary(seurat_proseg[["proseg"]]) <- "segmentation"
crop <- Crop(seurat_proseg[["proseg"]], x = c(1000, 1500), y = c(1500, 2000))
seurat_proseg[["proseg_zoom"]] <- crop
ImageDimPlot(seurat_proseg, fov = "proseg_zoom")

# Compare segmentations
DefaultBoundary(seurat_proseg[["xenium"]]) <- "segmentation"
ImageDimPlot(seurat_proseg, fov = "xenium")

# Switch to nuclei view
DefaultBoundary(seurat_proseg[["xenium"]]) <- "nuclei"
ImageDimPlot(seurat_proseg, fov = "xenium")

#########################
# DIAGNOSTICS
#########################

# # 1. What does the proseg GeoJSON sf look like after reading?
# tmp_poly <- tempfile(fileext = ".geojson")
# gz_con <- gzfile("./proseg_results/HS_Region_1/r_seurat/cell-polygons.geojson.gz", "rb")
# writeLines(readLines(gz_con, warn = FALSE), con = tmp_poly)
# close(gz_con)
# cell_polygons <- sf::st_read(tmp_poly, quiet = TRUE)
# 
# cat("--- GeoJSON sf object ---\n")
# print(str(cell_polygons[1:5, ]))
# cat("\nGeometry types:", unique(as.character(sf::st_geometry_type(cell_polygons))), "\n")
# cat("Columns:", colnames(cell_polygons), "\n")
# cat("nrow:", nrow(cell_polygons), "\n")
# cat("First 5 'cell' values:", head(cell_polygons$cell, 5), "\n")
# 
# # 2. What does cell-metadata.csv.gz look like?
# library(data.table)
# cells <- fread(cmd = "zcat ./proseg_results/HS_Region_1/r_seurat/cell-metadata.csv.gz")
# cat("\n--- cell-metadata.csv.gz ---\n")
# cat("Columns:", colnames(cells), "\n")
# cat("nrow:", nrow(cells), "\n")
# cat("First 5 original_cell_id:", head(cells$original_cell_id, 5), "\n")
# 
# # 3. After casting MULTIPOLYGON -> POLYGON, how many rows?
# poly_cast <- sf::st_cast(cell_polygons, "POLYGON", warn = FALSE)
# cat("\n--- After st_cast to POLYGON ---\n")
# cat("nrow before:", nrow(cell_polygons), " nrow after:", nrow(poly_cast), "\n")
# cat("Geometry types:", unique(as.character(sf::st_geometry_type(poly_cast))), "\n")
# 
# # 4. What does the Seurat object look like after CreateSeuratObject?
# # (just the colnames)
# cat("\n--- Seurat colnames ---\n")
# # If you already ran load_proseg partially and it failed, 
# # just recreate the counts + colnames part:
# library(Matrix)
# counts <- readMM(gzfile("./proseg_results/HS_Region_1/r_seurat/counts.mtx.gz"))
# counts <- t(counts)
# colnames(counts) <- make.unique(as.character(cells$original_cell_id))
# cat("ncol(counts):", ncol(counts), "\n")
# cat("First 5 colnames:", head(colnames(counts), 5), "\n")

# 1. What does the proseg centroid range look like?
cat("--- Centroids from metadata ---\n")
cat("centroid_x range:", range(seurat_proseg$centroid_x), "\n")
cat("centroid_y range:", range(seurat_proseg$centroid_y), "\n")

# 2. What does the FOV think the coordinates are?
fov_coords <- GetTissueCoordinates(seurat_proseg[["proseg"]])
cat("\n--- FOV coordinates ---\n")
cat("x range:", range(fov_coords$x), "\n")
cat("y range:", range(fov_coords$y), "\n")
cat("nrow:", nrow(fov_coords), "\n")
head(fov_coords)

# 3. What does your original plotting approach look like for comparison?
# (using the raw sf from the GeoJSON, not through Seurat)
library(ggplot2)
library(sf)
tmp_poly <- tempfile(fileext = ".geojson")
gz_con <- gzfile("./proseg_results/HS_Region_1/r_seurat/cell-polygons.geojson.gz", "rb")
writeLines(readLines(gz_con, warn = FALSE), con = tmp_poly)
close(gz_con)
raw_polys <- st_read(tmp_poly, quiet = TRUE)
p_raw <- ggplot(raw_polys) + geom_sf(fill = NA, color = "black", linewidth = 0.1) + 
  theme_void() + ggtitle("Raw GeoJSON polygons")

p_seurat <- ImageDimPlot(seurat_proseg, fov = "proseg") + ggtitle("ImageDimPlot")

p_raw | p_seurat