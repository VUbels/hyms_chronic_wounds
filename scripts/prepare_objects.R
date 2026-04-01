source("./scripts/helper_functions.R")
source("./scripts/setup_py_environment.R")

setup_py_env("hyms_chronic_wounds", "/home/uvictor/miniconda3/condabin/conda")

reticulate::py_run_string("
import sys
sys.path.append('./scripts')
")
zconv <- reticulate::import("zarr_h5ad_conversion", convert = TRUE)

build_proseg_seurat(
  proseg_dir = "./proseg_results",
  xenium_dir = "/mnt/d/HYMS/chronic_wounds"
)

annotate_proseg_seurat (proseg_dir = "./proseg_results",
                        reference_dir = "./reference",
                        # Method selection: "both", "singler", "rctd"
                        annotation_method  = "singler",
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
                        overwrite = TRUE)
