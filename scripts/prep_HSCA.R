#!/usr/bin/env Rscript
# ==============================================================================
# prep_hsca.R
#
# Slim down the HSCA extended reference object for use in build_reference.R
#
#   1. Clean inherited_celltype_lvl_5_extended: strip _1, _2, _3, _4 suffixes
#   2. Set cleaned column as active ident
#   3. Keep only first 25 metadata columns + inherited_celltype_lvl_5_extended
#   4. Downsample to max 2000 cells per label (keep all if < 2000)
#   5. Strip non-essential slots (reductions, graphs, non-RNA assays)
#   6. Save slim RDS
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})

# --- CONFIG ---
input_path  <- "/mnt/d/scRNA_datasets/HYMS_reference_dataset/HSCA_extended_ref.rds"
output_path <- "./reference/HSCA_slim.rds"
max_cells_per_type <- 2000
seed <- 123
# --------------

cat("Loading HSCA object...\n")
obj <- readRDS(input_path)
cat("Loaded:", ncol(obj), "cells,", nrow(obj), "genes\n")
cat("Memory estimate:", format(object.size(obj), units = "GB"), "\n\n")

# ==========================================================================
# 1. Clean cell type labels: strip trailing _1 through _4
# ==========================================================================

raw_labels <- obj$inherited_celltype_lvl_5_extended

# Remove suffixes _1, _2, _3, _4 at end of string only
cleaned_labels <- sub("_[1-4]$", "", raw_labels)

cat("Label cleaning:\n")
cat("  Unique labels before:", length(unique(raw_labels)), "\n")
cat("  Unique labels after: ", length(unique(cleaned_labels)), "\n")

# Show which labels were collapsed
changed <- raw_labels != cleaned_labels
if (any(changed)) {
  collapse_map <- unique(data.frame(
    original = raw_labels[changed],
    cleaned  = cleaned_labels[changed],
    stringsAsFactors = FALSE
  ))
  collapse_map <- collapse_map[order(collapse_map$cleaned, collapse_map$original), ]
  cat("  Collapsed labels:\n")
  for (i in seq_len(nrow(collapse_map))) {
    cat("    ", collapse_map$original[i], " -> ", collapse_map$cleaned[i], "\n")
  }
}

obj$inherited_celltype_lvl_5_extended <- cleaned_labels

# ==========================================================================
# 2. Set as active ident
# ==========================================================================

Idents(obj) <- "inherited_celltype_lvl_5_extended"
cat("\nActive ident set to: inherited_celltype_lvl_5_extended\n")

# ==========================================================================
# 3. Trim metadata: keep first 25 columns + the label column
# ==========================================================================

meta <- obj@meta.data
all_cols <- colnames(meta)

cat("\nMetadata columns (", length(all_cols), " total):\n")

# First 25
keep_cols <- all_cols[1:min(25, length(all_cols))]

# Ensure the label column is included even if not in first 25
label_col <- "inherited_celltype_lvl_5_extended"
if (!(label_col %in% keep_cols)) {
  keep_cols <- c(keep_cols, label_col)
}

drop_cols <- setdiff(all_cols, keep_cols)

cat("  Keeping:", length(keep_cols), "columns\n")
cat("  Dropping:", length(drop_cols), "columns\n")

if (length(drop_cols) > 0) {
  # Print first few dropped for confirmation
  n_show <- min(10, length(drop_cols))
  cat("  Dropped (first ", n_show, "): ",
      paste(drop_cols[1:n_show], collapse = ", "),
      if (length(drop_cols) > n_show) paste0(", ... +", length(drop_cols) - n_show, " more"),
      "\n")
}

obj@meta.data <- meta[, keep_cols, drop = FALSE]

# ==========================================================================
# 4. Downsample: max 2000 cells per label
# ==========================================================================

cat("\nDownsampling to max", max_cells_per_type, "cells per type...\n")

set.seed(seed)

type_counts <- table(obj$inherited_celltype_lvl_5_extended)
cat("  Types:", length(type_counts), "\n")
cat("  Types already <=", max_cells_per_type, ":",
    sum(type_counts <= max_cells_per_type), "\n")
cat("  Types to downsample:",
    sum(type_counts > max_cells_per_type), "\n")

# Stratified sampling
cells_keep <- c()
for (ct in names(type_counts)) {
  ct_cells <- colnames(obj)[obj$inherited_celltype_lvl_5_extended == ct]
  if (length(ct_cells) > max_cells_per_type) {
    ct_cells <- sample(ct_cells, max_cells_per_type)
  }
  cells_keep <- c(cells_keep, ct_cells)
}

n_before <- ncol(obj)
obj <- subset(obj, cells = cells_keep)

cat("  Cells:", n_before, "->", ncol(obj), "\n")
cat("  Reduction:", round((1 - ncol(obj) / n_before) * 100, 1), "%\n")

# ==========================================================================
# 5. Strip non-essential slots
# ==========================================================================

cat("\nStripping non-essential data...\n")

# Remove all assays except RNA
assay_names <- Assays(obj)
for (a in setdiff(assay_names, "RNA")) {
  obj[[a]] <- NULL
  cat("  Removed assay:", a, "\n")
}

# Remove all reductions
red_names <- Reductions(obj)
for (r in red_names) {
  obj[[r]] <- NULL
  cat("  Removed reduction:", r, "\n")
}

# Remove graphs
graph_names <- names(obj@graphs)
if (length(graph_names) > 0) {
  for (g in graph_names) {
    obj@graphs[[g]] <- NULL
    cat("  Removed graph:", g, "\n")
  }
}

# Remove neighbors
if (length(obj@neighbors) > 0) {
  obj@neighbors <- list()
  cat("  Removed neighbor data\n")
}

# Coerce to v5 assay if needed
if (!inherits(obj[["RNA"]], "Assay5")) {
  obj[["RNA"]] <- as(obj[["RNA"]], "Assay5")
  cat("  Coerced RNA to Assay5\n")
}

# Drop data/scale.data layers if present (build_reference.R re-normalizes)
existing_layers <- Layers(obj[["RNA"]])
cat("  RNA layers present:", paste(existing_layers, collapse = ", "), "\n")

for (lyr in intersect(existing_layers, c("data", "scale.data"))) {
  obj[["RNA"]][[lyr]] <- NULL
  cat("  Removed layer:", lyr, "\n")
}

cat("  RNA layers kept:", paste(Layers(obj[["RNA"]]), collapse = ", "), "\n")

# ==========================================================================
# 6. Save
# ==========================================================================

cat("\nFinal object:\n")
cat("  Cells:", ncol(obj), "\n")
cat("  Genes:", nrow(obj), "\n")
cat("  Types:", length(unique(obj$inherited_celltype_lvl_5_extended)), "\n")
cat("  Memory:", format(object.size(obj), units = "GB"), "\n")

cat("\nSaving to:", output_path, "\n")
saveRDS(obj, output_path)

file_size <- file.info(output_path)$size / 1e9
cat("File size on disk:", round(file_size, 2), "GB\n")

# Print per-type summary
cat("\nPer-type cell counts:\n")
final_counts <- sort(table(obj$inherited_celltype_lvl_5_extended), decreasing = TRUE)
for (i in seq_along(final_counts)) {
  cat(sprintf("  %-40s %5d\n", names(final_counts)[i], final_counts[i]))
}

cat("\nDone. Use this as input to build_reference.R:\n")
cat('  HSCA = list(\n')
cat('    type           = "rds",\n')
cat(paste0('    path           = "', output_path, '",\n'))
cat('    annotation_col = "inherited_celltype_lvl_5_extended",\n')
cat('    batch_col      = "sample",   # <-- verify column name\n')
cat('    skip_qc        = TRUE,\n')
cat('    skip_doublets  = TRUE\n')
cat('  )\n')