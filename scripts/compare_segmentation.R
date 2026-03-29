source("./scripts/helper_functions.R")

seurat_proseg <- load_proseg(proseg_dir = "./proseg_results/HS_Region_1/r_seurat/")

seurat_proseg <- store_xenium_polygons(seurat_proseg, "/mnt/d/hyms/chronic_wounds/hs_region_1/cell_boundaries.parquet")
seurat_proseg <- store_xenium_nuclei(seurat_proseg, "/mnt/d/hyms/chronic_wounds/hs_region_1/nucleus_boundaries.parquet")

plot_polygons(seurat_proseg, polygon_slot = "proseg_polygons")

plot_two_segmentation_panels(
  seu = seurat_proseg,
  cell_seg_slot1 = "proseg_polygons",
  cell_seg_slot2 = "xenium_polygons",
  nuclei_slot = "xenium_nuclei",
  labels = c("ProSeg", "Xenium")
)

scout_and_plot_comparison(
  seu = seurat_proseg,
  cell_seg_slot1 = "proseg_polygons",
  cell_seg_slot2 = "xenium_polygons",
  nuclei_slot = "xenium_nuclei"
)

compare_nuclei_segmentation(
  seu = seurat_proseg,
  cell_segmentation_slot1 = "proseg_polygons",
  cell_segmentation_slot2 = "xenium_polygons",
  nuclei_segmentation_slot = "xenium_nuclei",
  labels = c("ProSeg", "Xenium")
)

plot_spatial_feature(
  seu = seurat_proseg,
  cell_seg_slot = "proseg_polygons",
  feature = "SERPINB7",
  nuclei_slot = "xenium_nuclei"
)

# Replace with your actual parquet file path
parquet_file <- "/mnt/d/HYMS/chronic_wounds/HS_Region_1/cell_boundaries.parquet"

cell_boundaries <- read_parquet(parquet_file)

# Check the data
head(cell_boundaries)