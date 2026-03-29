#!/usr/bin/env Rscript
# ==============================================================================
# prep_CW.R
#
# Build a single Seurat RDS from the GSE265972 (Liu, Landén et al.) dataset.
#
# Input files:
#   - GSE265972_processed.tar.gz  Contains per-sample .zip files, each holding
#                                  CellRanger output (barcodes/features/matrix)
#   - CW_anno_metadata.txt        Post-QC metadata with cell type annotations
#                                  for 48,346 cells across 9 samples
#
# Pipeline:
#   1. Extract tar.gz -> per-sample .zip files
#   2. Extract each .zip -> find barcodes.tsv.gz / features.tsv.gz / matrix.mtx.gz
#   3. Load each sample with Read10X, prefix barcodes to match metadata
#   4. Merge all samples into one Seurat object
#   5. Filter to cells present in metadata (already QC'd by original authors)
#   6. Attach metadata (cell types, condition, QC stats)
#   7. Strip to counts + metadata only
#   8. Save slim RDS
#
# Output:
#   - CW_slim.rds           Ready for build_reference.R
#
# Dependencies:
#   library(Seurat), library(SeuratObject), library(data.table)
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(data.table)
  library(Matrix)
})

# --- CONFIG ---
tar_path      <- "/mnt/d/scRNA_datasets/HYMS_reference_dataset/GSE265972_processed.tar.gz"
metadata_path <- "/mnt/d/scRNA_datasets/HYMS_reference_dataset/CW_anno_metadata.txt"
output_path   <- "./reference/CW_slim.rds"
scratch_dir   <- "./tmp_GSE265972"    # temporary extraction directory
seed          <- 123
# --------------

dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)


# ==============================================================================
# 1. Read metadata
# ==============================================================================

cat("=== Step 1: Reading metadata ===\n")

meta <- fread(metadata_path, sep = "\t", header = TRUE)
cat("  Metadata rows:", nrow(meta), "\n")
cat("  Columns:", paste(colnames(meta), collapse = ", "), "\n")
cat("  Samples:", paste(sort(unique(meta$samples)), collapse = ", "), "\n")
cat("  Conditions:", paste(sort(unique(meta$Condition)), collapse = ", "), "\n")
cat("  Cell types:", length(unique(meta$CellType)), "\n\n")

# The barcodeindex column is the cell identifier: e.g. NS23_AAACCCAAGGGCTTCC-1
# This is what we need to match after loading the 10X data
expected_samples <- sort(unique(meta$samples))
cat("  Expected samples from metadata:\n")
for (s in expected_samples) {
  cat("    ", s, ":", sum(meta$samples == s), "cells\n")
}
cat("\n")


# ==============================================================================
# 2. Extract tar.gz -> per-sample .zip files
# ==============================================================================

cat("=== Step 2: Extracting tar.gz ===\n")

dir.create(scratch_dir, showWarnings = FALSE, recursive = TRUE)

cat("  Extracting:", tar_path, "\n")
untar(tar_path, exdir = scratch_dir)

# List extracted contents
extracted <- list.files(scratch_dir, recursive = FALSE)
cat("  Extracted files:\n")
for (f in extracted) {
  cat("    ", f, "\n")
}
cat("\n")

# Find .zip files
zip_files <- list.files(scratch_dir, pattern = "\\.zip$",
                        full.names = TRUE, recursive = FALSE)
cat("  Found", length(zip_files), "zip files\n\n")


# ==============================================================================
# 3. Extract each .zip and load 10X data
# ==============================================================================

cat("=== Step 3: Loading per-sample 10X data ===\n\n")

sample_objects <- list()

for (zf in zip_files) {
  
  # Derive sample name from zip filename: NS23.zip -> NS23
  sample_name <- sub("\\.zip$", "", basename(zf))
  cat("--- ", sample_name, " ---\n")
  
  # Create per-sample extraction directory
  sample_dir <- file.path(scratch_dir, sample_name)
  dir.create(sample_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Extract zip
  unzip(zf, exdir = sample_dir)
  
  # Recursively find the 10X matrix files
  # CellRanger outputs: barcodes.tsv.gz, features.tsv.gz (or genes.tsv.gz),
  #                     matrix.mtx.gz
  # They may be nested inside a subdirectory (e.g. filtered_feature_bc_matrix/)
  mtx_files <- list.files(sample_dir, pattern = "matrix\\.mtx(\\.gz)?$",
                          recursive = TRUE, full.names = TRUE)
  
  if (length(mtx_files) == 0) {
    cat("  [WARN] No matrix.mtx found for ", sample_name, " — skipping\n\n")
    next
  }
  
  # Use the directory containing the matrix file as the Read10X input
  tenx_dir <- dirname(mtx_files[1])
  cat("  10X dir: ", tenx_dir, "\n")
  
  # List contents for verification
  tenx_contents <- list.files(tenx_dir)
  cat("  Contents:", paste(tenx_contents, collapse = ", "), "\n")
  
  # Load with Read10X
  tryCatch({
    counts <- Read10X(data.dir = tenx_dir)
    
    # Handle multi-modal returns
    if (is.list(counts) && !is(counts, "dgCMatrix")) {
      if ("Gene Expression" %in% names(counts)) {
        counts <- counts[["Gene Expression"]]
      } else {
        counts <- counts[[1]]
      }
    }
    
    cat("  Raw barcodes:", ncol(counts), "  Genes:", nrow(counts), "\n")
    
    # Prefix barcodes with sample name to match metadata format
    # Metadata format: NS23_AAACCCAAGGGCTTCC-1
    colnames(counts) <- paste0(sample_name, "_", colnames(counts))
    
    # Check overlap with metadata
    meta_barcodes <- meta$barcodeindex[meta$samples == sample_name]
    overlap <- sum(colnames(counts) %in% meta_barcodes)
    cat("  Overlap with metadata:", overlap, "/", length(meta_barcodes),
        "expected cells\n")
    
    # Create Seurat object
    obj <- CreateSeuratObject(counts = counts, project = sample_name,
                              min.cells = 0, min.features = 0)
    obj$sample_name <- sample_name
    
    sample_objects[[sample_name]] <- obj
    cat("  Loaded successfully\n\n")
    
  }, error = function(e) {
    cat("  [ERROR] Failed to load: ", conditionMessage(e), "\n\n")
  })
}

if (length(sample_objects) == 0) {
  stop("[ERROR] No samples loaded successfully.")
}

cat("Successfully loaded", length(sample_objects), "samples\n\n")


# ==============================================================================
# 4. Merge all samples
# ==============================================================================

cat("=== Step 4: Merging samples ===\n")

if (length(sample_objects) == 1) {
  merged <- sample_objects[[1]]
} else {
  merged <- merge(sample_objects[[1]],
                  y = sample_objects[-1])
}

cat("  Merged object:", ncol(merged), "barcodes,", nrow(merged), "genes\n\n")

# Free memory
rm(sample_objects)
gc(verbose = FALSE)


# ==============================================================================
# 5. Filter to metadata barcodes (authors' QC)
# ==============================================================================

cat("=== Step 5: Filtering to QC'd cells from metadata ===\n")

# Barcodes in metadata
meta_bcs <- meta$barcodeindex

# Barcodes in merged object
merged_bcs <- colnames(merged)

in_both    <- intersect(merged_bcs, meta_bcs)
in_meta_only    <- setdiff(meta_bcs, merged_bcs)
in_merged_only  <- length(merged_bcs) - length(in_both)

cat("  Cells in merged object:", length(merged_bcs), "\n")
cat("  Cells in metadata:     ", length(meta_bcs), "\n")
cat("  Intersection:          ", length(in_both), "\n")
cat("  In metadata only:      ", length(in_meta_only), "\n")
cat("  In merged only:        ", in_merged_only, "\n")

if (length(in_both) == 0) {
  # Debug: show barcode format comparison
  cat("\n  [DEBUG] First 3 merged barcodes:\n")
  cat("    ", paste(head(merged_bcs, 3), collapse = "\n     "), "\n")
  cat("  [DEBUG] First 3 metadata barcodes:\n")
  cat("    ", paste(head(meta_bcs, 3), collapse = "\n     "), "\n")
  stop("[ERROR] No barcode overlap. Check prefix format.")
}

merged <- subset(merged, cells = in_both)
cat("  After filter:", ncol(merged), "cells\n\n")


# ==============================================================================
# 6. Attach metadata
# ==============================================================================

cat("=== Step 6: Attaching metadata ===\n")

# Convert metadata to data.frame indexed by barcodeindex
meta_df <- as.data.frame(meta)
rownames(meta_df) <- meta_df$barcodeindex

# Align to current cell order
meta_df <- meta_df[colnames(merged), , drop = FALSE]

# Attach all metadata columns to the Seurat object
for (col in colnames(meta_df)) {
  if (col == "barcodeindex") next   # already used as rownames
  merged[[col]] <- meta_df[[col]]
}

cat("  Attached", ncol(meta_df) - 1, "metadata columns\n")
cat("  CellType values:", length(unique(merged$CellType)), "types\n")
cat("  Condition values:", paste(unique(merged$Condition), collapse = ", "), "\n\n")


# ==============================================================================
# 7. Set ident and strip non-essential data
# ==============================================================================

cat("=== Step 7: Finalising object ===\n")

# Set CellType as active ident
Idents(merged) <- "CellType"
cat("  Active ident set to: CellType\n")

# Coerce to v5 assay
if (!inherits(merged[["RNA"]], "Assay5")) {
  merged[["RNA"]] <- as(merged[["RNA"]], "Assay5")
  cat("  Coerced RNA to Assay5\n")
}

# Remove all assays except RNA
for (a in setdiff(Assays(merged), "RNA")) {
  merged[[a]] <- NULL
  cat("  Removed assay:", a, "\n")
}

# Remove any reductions
for (r in Reductions(merged)) {
  merged[[r]] <- NULL
  cat("  Removed reduction:", r, "\n")
}

# Remove graphs
if (length(names(merged@graphs)) > 0) {
  for (g in names(merged@graphs)) {
    merged@graphs[[g]] <- NULL
    cat("  Removed graph:", g, "\n")
  }
}

# Drop data/scale.data layers (build_reference.R re-normalizes)
existing_layers <- Layers(merged[["RNA"]])
cat("  RNA layers present:", paste(existing_layers, collapse = ", "), "\n")

for (lyr in intersect(existing_layers, c("data", "scale.data"))) {
  merged[["RNA"]][[lyr]] <- NULL
  cat("  Removed layer:", lyr, "\n")
}

cat("  RNA layers kept:", paste(Layers(merged[["RNA"]]), collapse = ", "), "\n\n")


# ==============================================================================
# 8. Save
# ==============================================================================

cat("=== Step 8: Saving ===\n\n")

cat("Final object:\n")
cat("  Cells:", ncol(merged), "\n")
cat("  Genes:", nrow(merged), "\n")
cat("  Types:", length(unique(merged$CellType)), "\n")
cat("  Samples:", length(unique(merged$samples)), "\n")
cat("  Memory:", format(object.size(merged), units = "GB"), "\n")

cat("\nSaving to:", output_path, "\n")
saveRDS(merged, output_path)

file_size <- file.info(output_path)$size / 1e9
cat("File size on disk:", round(file_size, 2), "GB\n")

# Per-type summary
cat("\nPer-type cell counts:\n")
final_counts <- sort(table(merged$CellType), decreasing = TRUE)
for (i in seq_along(final_counts)) {
  cat(sprintf("  %-20s %5d\n", names(final_counts)[i], final_counts[i]))
}

# Per-condition summary
cat("\nPer-condition cell counts:\n")
cond_counts <- sort(table(merged$Condition), decreasing = TRUE)
for (i in seq_along(cond_counts)) {
  cat(sprintf("  %-10s %5d\n", names(cond_counts)[i], cond_counts[i]))
}

# Clean up scratch directory
cat("\nCleaning up scratch directory:", scratch_dir, "\n")
unlink(scratch_dir, recursive = TRUE)

cat("\nDone. Use this as input to build_reference.R:\n")
cat('  GSE265972 = list(\n')
cat('    type           = "rds",\n')
cat(paste0('    path           = "', output_path, '",\n'))
cat('    annotation_col = "CellType",\n')
cat('    batch_col      = "samples",\n')
cat('    skip_qc        = TRUE,       # already filtered by original authors\n')
cat('    skip_doublets  = TRUE        # already removed by original authors\n')
cat('  )\n')