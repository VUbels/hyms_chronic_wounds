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

annotate_proseg_seurat(
  proseg_dir    = "./proseg_results",
  reference_dir = "./reference"
)