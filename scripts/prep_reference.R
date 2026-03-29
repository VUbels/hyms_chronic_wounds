#!/usr/bin/env Rscript
# ==============================================================================
# prep_reference.R
#
# Compute Jaccard similarity between cell type labels across reference
# datasets based on shared marker genes.  Supports splitting inputs by
# a metadata column (e.g. condition) to compare healthy vs disease subsets
# independently against a reference atlas.
#
# Uses the same config$inputs structure as build_reference.R, extended
# with optional subset_col / subset_values fields.
#
# Approach:
#   1. Load all annotated inputs; optionally split by condition
#   2. Normalize each independently
#   3. Find top marker genes per cell type per subset
#   4. Compute pairwise Jaccard for every cross-dataset label pair
#   5. Plot heatmaps; produce delta-Jaccard comparison when subsets exist
#   6. Suggest label mappings
#
# Output (in output_dir):
#   - label_jaccard_heatmap_<A>_vs_<B>.pdf         Filtered heatmap per pair
#   - label_jaccard_heatmap_full_<A>_vs_<B>.pdf     Full heatmap per pair
#   - label_jaccard_matrix_<A>_vs_<B>.csv            Raw Jaccard matrix per pair
#   - label_jaccard_delta_<subset1>_vs_<subset2>.pdf Delta heatmap (if subsets)
#   - suggested_label_mapping.csv                    Best matches all pairs
#
# Dependencies:
#   install.packages(c("Seurat","SeuratObject","Matrix","data.table",
#                      "ggplot2","pheatmap"))
#   devtools::install_github("immunogenomics/presto")  # optional, faster DE
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(ggplot2)
})

source("./scripts/helper_functions.R")

# ==============================================================================
# CONFIGURATION
# ==============================================================================

config <- list(

  inputs = list(

    # ---- Input 1: pre-processed atlas (reference baseline) ----
    HSCA_ref = list(
      type           = "rds",
      path           = "./reference/HSCA_slim.rds",
      annotation_col = "inherited_celltype_lvl_5_extended",
      batch_col      = "sample",
      skip_qc        = TRUE,
      skip_doublets  = TRUE
      # No subset_col -> used as a single block
    ),

    # ---- Input 2: chronic wound data, SPLIT by condition ----
    # subset_col:    metadata column to split on
    # subset_values: named list mapping {display_name -> column_value}
    #                Each subset becomes its own entry in the Jaccard analysis.
    #                If omitted, the input is used as a single block.
    CW_ref = list(
      type           = "rds",
      path           = "./reference/CW_slim.rds",
      annotation_col = "CellType",
      batch_col      = "samples",
      skip_qc        = TRUE,
      skip_doublets  = TRUE,
      subset_col     = "Condition",
      subset_values  = list(
        CW_healthy = "NS",     # normal skin
        CW_VU      = "VU"      # venous ulcer
      )
    )
    # ---- Add more inputs here ----
  ),

  # --- Jaccard-specific settings ---
  n_top_markers  = 200,
  jaccard_thresh = 0.10,

  # --- Output ---
  output_dir = "./reference/"
)


# ==============================================================================
# INITIALISE
# ==============================================================================

dir.create(config$output_dir, showWarnings = FALSE, recursive = TRUE)

cat("========================================================\n")
cat("   Label Jaccard Similarity Analysis\n")
cat("   (with condition-based subsetting)\n")
cat("========================================================\n\n")


# ==============================================================================
# HELPER: fast marker detection
# ==============================================================================

get_top_markers <- function(obj, label_col, n_top = 200) {

  Idents(obj) <- label_col

  if (requireNamespace("presto", quietly = TRUE)) {
    cat("  Using presto for marker detection...\n")

    expr_mat <- GetAssayData(obj, layer = "data")
    labels   <- obj@meta.data[[label_col]]

    res <- presto::wilcoxauc(expr_mat, labels)
    res <- as.data.table(res)

    res <- res[logFC > 0 & padj < 0.05]
    res <- res[order(group, -auc)]
    markers <- res[, head(.SD, n_top), by = group]

    marker_list <- split(markers$feature, markers$group)

  } else {
    cat("  presto not available. Using FindAllMarkers (this may be slow)...\n")

    all_markers <- FindAllMarkers(obj,
                                  only.pos  = TRUE,
                                  min.pct   = 0.1,
                                  logfc.threshold = 0.25,
                                  test.use  = "wilcox",
                                  verbose   = FALSE)
    all_markers <- as.data.table(all_markers)
    all_markers <- all_markers[order(cluster, -avg_log2FC)]
    top <- all_markers[, head(.SD, n_top), by = cluster]

    marker_list <- split(top$gene, top$cluster)
  }

  return(marker_list)
}

# ==============================================================================
# 1. Load inputs and expand subsets
# ==============================================================================

cat("=== Step 1: Loading inputs ===\n\n")

# datasets: flat list of name -> list(obj, label_col, markers, parent, subset_value)
datasets <- list()

for (nm in names(config$inputs)) {

  inp <- config$inputs[[nm]]

  if (is.null(inp$annotation_col)) {
    cat(nm, "| annotation_col is NULL — skipping\n\n")
    next
  }

  if (inp$type == "rds") {
    cat(nm, "| Loading:", inp$path, "\n")
    obj <- readRDS(inp$path)
  } else {
    cat(nm, "| type '", inp$type, "' — run prep script first.\n\n")
    next
  }

  # Coerce + join
  if (!inherits(obj[["RNA"]], "Assay5")) {
    obj[["RNA"]] <- as(obj[["RNA"]], "Assay5")
  }
  if (length(Layers(obj[["RNA"]])) > 1) {
    obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
  }

  # Validate annotation column
  if (!(inp$annotation_col %in% colnames(obj@meta.data))) {
    cat(nm, "| [ERROR] annotation_col '", inp$annotation_col,
        "' not found. Skipping.\n\n")
    next
  }

  cat(nm, "|", ncol(obj), "cells,", nrow(obj), "genes,",
      length(unique(obj[[inp$annotation_col]][, 1])), "types\n")

  # ------------------------------------------------------------------
  # Subset expansion
  # ------------------------------------------------------------------
  if (!is.null(inp$subset_col) && !is.null(inp$subset_values)) {

    if (!(inp$subset_col %in% colnames(obj@meta.data))) {
      cat(nm, "| [ERROR] subset_col '", inp$subset_col,
          "' not found. Loading as single block.\n")
      inp$subset_col <- NULL
    }
  }

  if (!is.null(inp$subset_col) && !is.null(inp$subset_values)) {

    cat(nm, "| Splitting by '", inp$subset_col, "':\n")
    condition_vec <- obj[[inp$subset_col]][, 1]

    for (sub_name in names(inp$subset_values)) {
      sub_val   <- inp$subset_values[[sub_name]]
      sub_cells <- which(condition_vec == sub_val)

      if (length(sub_cells) == 0) {
        cat("  ", sub_name, " (", sub_val, "): 0 cells — skipping\n")
        next
      }

      sub_obj <- subset(obj, cells = colnames(obj)[sub_cells])
      n_types <- length(unique(sub_obj[[inp$annotation_col]][, 1]))

      cat("  ", sub_name, " (", sub_val, "):", ncol(sub_obj), "cells,",
          n_types, "types\n")

      # Normalize
      sub_obj <- NormalizeData(sub_obj, verbose = FALSE)

      datasets[[sub_name]] <- list(
        obj          = sub_obj,
        label_col    = inp$annotation_col,
        parent       = nm,
        subset_value = sub_val
      )
    }

    # Also keep the combined (unsplit) version for reference
    cat("  ", nm, " (combined):", ncol(obj), "cells\n")
    obj <- NormalizeData(obj, verbose = FALSE)
    datasets[[nm]] <- list(
      obj          = obj,
      label_col    = inp$annotation_col,
      parent       = nm,
      subset_value = "all"
    )

  } else {
    # No subsetting — single block
    obj <- NormalizeData(obj, verbose = FALSE)
    datasets[[nm]] <- list(
      obj          = obj,
      label_col    = inp$annotation_col,
      parent       = nm,
      subset_value = "all"
    )
  }

  cat(nm, "| Ready\n\n")
}

if (length(datasets) < 2) {
  stop("[ERROR] Need at least 2 entries for comparison. Found: ",
       length(datasets))
}

cat("Loaded", length(datasets), "dataset entries:",
    paste(names(datasets), collapse = ", "), "\n\n")


# ==============================================================================
# 2. Find markers for each entry
# ==============================================================================

cat("=== Step 2: Finding markers ===\n\n")

for (nm in names(datasets)) {
  ds <- datasets[[nm]]
  n_types <- length(unique(ds$obj[[ds$label_col]][, 1]))
  cat(nm, " (", n_types, " types):\n", sep = "")

  markers <- get_top_markers(ds$obj, ds$label_col, config$n_top_markers)
  datasets[[nm]]$markers <- markers

  cat("  Got markers for", length(markers), "types\n")
  cat("  Median markers per type:", median(sapply(markers, length)), "\n\n")
}


# ==============================================================================
# 3. Pairwise Jaccard: each subset vs HSCA_ref
# ==============================================================================

cat("=== Step 3: Pairwise Jaccard comparisons ===\n\n")

all_suggestions <- data.frame(
  source_row    = character(),
  source_col    = character(),
  type_row      = character(),
  best_col      = character(),
  jaccard       = numeric(),
  second_col    = character(),
  jaccard_2nd   = numeric(),
  stringsAsFactors = FALSE
)

# Store Jaccard matrices keyed by pair name for delta computation later
jaccard_store <- list()

ds_names <- names(datasets)

# Compare every entry against every other (excluding self)
# But skip comparisons between subsets of the same parent (e.g. CW_healthy vs CW_VU)
# unless you want those too — uncomment the skip below to enable
done_pairs <- c()

for (i in seq_along(ds_names)) {
  for (j in seq_along(ds_names)) {
    if (i == j) next

    nm_col <- ds_names[i]   # columns (reference side)
    nm_row <- ds_names[j]   # rows (query side)

    # Skip if this is a subset compared against its own combined parent
    if (nm_col == datasets[[nm_row]]$parent) next
    if (nm_row == datasets[[nm_col]]$parent) next

    # Skip duplicate pairs (A vs B == B vs A for our purposes)
    pair_key <- paste(sort(c(nm_col, nm_row)), collapse = "___")
    if (pair_key %in% done_pairs) next
    done_pairs <- c(done_pairs, pair_key)

    # Skip subsets from same parent against each other
    # (CW_healthy vs CW_VU is interesting but not the primary comparison)
    # Remove the next 2 lines if you want those comparisons too:
    # if (datasets[[nm_col]]$parent == datasets[[nm_row]]$parent) next

    markers_col <- datasets[[nm_col]]$markers
    markers_row <- datasets[[nm_row]]$markers

    cat("--- ", nm_row, " (rows) vs ", nm_col, " (columns) ---\n")
    cat("  ", length(names(markers_row)), " x ", length(names(markers_col)), "\n")

    # Compute
    jaccard_mat <- compute_jaccard(markers_row, markers_col)
    pair_tag    <- paste0(nm_col, "_vs_", nm_row)
    jaccard_store[[pair_tag]] <- jaccard_mat

    cat("  Range: [", round(min(jaccard_mat), 4), ",",
        round(max(jaccard_mat), 4), "]\n")

    # Save matrix
    jac_df <- as.data.frame(jaccard_mat)
    jac_df <- cbind(type = rownames(jac_df), jac_df)
    fwrite(jac_df, file.path(config$output_dir,
                              paste0("label_jaccard_matrix_", pair_tag, ".csv")))

    # Suggestions
    cat("\n  Suggested mappings (Jaccard >=", config$jaccard_thresh, "):\n")

    for (tr in rownames(jaccard_mat)) {
      scores <- jaccard_mat[tr, ]
      top2   <- sort(scores, decreasing = TRUE)[1:min(2, length(scores))]

      best_name  <- names(top2)[1]
      best_score <- top2[1]
      second_name  <- if (length(top2) >= 2) names(top2)[2] else NA_character_
      second_score <- if (length(top2) >= 2) top2[2] else NA_real_

      all_suggestions <- rbind(all_suggestions, data.frame(
        source_row  = nm_row,
        source_col  = nm_col,
        type_row    = tr,
        best_col    = best_name,
        jaccard     = round(best_score, 4),
        second_col  = second_name,
        jaccard_2nd = round(second_score, 4),
        stringsAsFactors = FALSE
      ))

      flag <- if (best_score >= config$jaccard_thresh) "***" else "   "
      cat(sprintf("  %s %-20s -> %-35s (J=%.3f)  | 2nd: %-35s (J=%.3f)\n",
                  flag, tr, best_name, best_score,
                  ifelse(is.na(second_name), "-", second_name),
                  ifelse(is.na(second_score), 0, second_score)))
    }

    # Heatmaps
    cat("\n")
    plot_jaccard_heatmap(jaccard_mat, nm_row, nm_col, pair_tag,
                         config$output_dir, config$n_top_markers, filtered = TRUE)
    plot_jaccard_heatmap(jaccard_mat, nm_row, nm_col, pair_tag,
                         config$output_dir, config$n_top_markers, filtered = FALSE)

    cat("\n")
  }
}


# ==============================================================================
# 4. Delta-Jaccard: compare healthy vs disease divergence from HSCA
# ==============================================================================

cat("=== Step 4: Delta-Jaccard (condition comparison) ===\n\n")

# Find inputs that were subset and identify their paired comparisons against
# the same reference.  For each pair of subsets from the same parent, compute
# delta = J(healthy, HSCA) - J(disease, HSCA) per shared cell type.
# Positive delta = healthy is more similar to HSCA than disease.

# Identify subsets grouped by parent
parents <- unique(sapply(datasets, function(d) d$parent))

for (par in parents) {
  subset_entries <- names(datasets)[sapply(datasets, function(d) {
    d$parent == par && d$subset_value != "all"
  })]

  if (length(subset_entries) < 2) next

  cat("Parent '", par, "' has subsets: ", paste(subset_entries, collapse = ", "), "\n")

  # Find which reference(s) both subsets were compared against
  ref_entries <- setdiff(names(datasets), names(datasets)[sapply(datasets, function(d) {
    d$parent == par
  })])

  for (ref_nm in ref_entries) {
    # Look up Jaccard matrices for each subset vs this reference
    subset_mats <- list()
    for (sub_nm in subset_entries) {
      # Try both orderings of the pair tag
      tag1 <- paste0(ref_nm, "_vs_", sub_nm)
      tag2 <- paste0(sub_nm, "_vs_", ref_nm)
      if (tag1 %in% names(jaccard_store)) {
        subset_mats[[sub_nm]] <- jaccard_store[[tag1]]
      } else if (tag2 %in% names(jaccard_store)) {
        # Transpose: rows and columns are swapped
        subset_mats[[sub_nm]] <- t(jaccard_store[[tag2]])
      }
    }

    if (length(subset_mats) < 2) next

    cat("  Comparing subsets against reference '", ref_nm, "'\n")

    # For each pair of subsets, compute delta on shared row types
    for (k in seq_len(length(subset_entries) - 1)) {
      for (l in (k + 1):length(subset_entries)) {
        sub_a <- subset_entries[k]
        sub_b <- subset_entries[l]

        mat_a <- subset_mats[[sub_a]]
        mat_b <- subset_mats[[sub_b]]

        # Shared cell types (rows) and shared HSCA types (columns)
        shared_rows <- intersect(rownames(mat_a), rownames(mat_b))
        shared_cols <- intersect(colnames(mat_a), colnames(mat_b))

        if (length(shared_rows) == 0 || length(shared_cols) == 0) {
          cat("    No shared types between ", sub_a, " and ", sub_b, "\n")
          next
        }

        mat_a_sub <- mat_a[shared_rows, shared_cols, drop = FALSE]
        mat_b_sub <- mat_b[shared_rows, shared_cols, drop = FALSE]

        # Delta: positive means sub_a (first subset) is MORE similar to HSCA
        delta <- mat_a_sub - mat_b_sub

        cat("    Delta (", sub_a, " - ", sub_b, "):\n")
        cat("      Shared types:", length(shared_rows), "rows x",
            length(shared_cols), "cols\n")
        cat("      Range: [", round(min(delta), 4), ",",
            round(max(delta), 4), "]\n")

        # Per-row summary: for each CW type, what's its best-match delta?
        cat("\n    Per-type divergence (best HSCA match):\n")
        cat(sprintf("    %-20s  J(%s)   J(%s)   Delta   Best HSCA match\n",
                    "CellType", sub_a, sub_b))
        cat("    ", paste(rep("-", 85), collapse = ""), "\n")

        for (tr in shared_rows) {
          best_a <- max(mat_a_sub[tr, ])
          best_b <- max(mat_b_sub[tr, ])
          best_col_a <- colnames(mat_a_sub)[which.max(mat_a_sub[tr, ])]
          best_col_b <- colnames(mat_b_sub)[which.max(mat_b_sub[tr, ])]
          d <- best_a - best_b

          # Flag types where disease diverges substantially
          flag <- ""
          if (d > 0.05)  flag <- "  <-- healthy closer"
          if (d < -0.05) flag <- "  <-- disease closer"

          cat(sprintf("    %-20s  %.3f   %.3f   %+.3f  %s / %s%s\n",
                      tr, best_a, best_b, d, best_col_a, best_col_b, flag))
        }

        # Save delta matrix
        delta_tag <- paste0("delta_", sub_a, "_vs_", sub_b, "_ref_", ref_nm)
        delta_df  <- as.data.frame(delta)
        delta_df  <- cbind(type = rownames(delta_df), delta_df)
        fwrite(delta_df, file.path(config$output_dir,
                                    paste0("label_jaccard_", delta_tag, ".csv")))

        # Delta heatmap
        cat("\n    Generating delta heatmap...\n")

        # Filter to relevant columns
        max_per_col <- apply(abs(delta), 2, max)
        relevant_cols <- names(max_per_col[max_per_col >= 0.01])
        if (length(relevant_cols) < 3) {
          relevant_cols <- names(sort(max_per_col, decreasing = TRUE))[
            1:min(30, length(max_per_col))]
        }
        delta_plot <- delta[, relevant_cols, drop = FALSE]

        fname <- file.path(config$output_dir,
                            paste0("label_jaccard_", delta_tag, ".pdf"))

        if (requireNamespace("pheatmap", quietly = TRUE)) {
          library(pheatmap)

          # Diverging color scale: blue = sub_a closer, red = sub_b closer
          max_abs <- max(abs(delta_plot))
          breaks  <- seq(-max_abs, max_abs, length.out = 101)

          pdf(fname,
              width  = max(14, ncol(delta_plot) * 0.25),
              height = max(8, nrow(delta_plot) * 0.5))

          pheatmap(delta_plot,
                   color           = colorRampPalette(c("firebrick3", "white",
                                                        "steelblue3"))(100),
                   breaks          = breaks,
                   cluster_rows    = TRUE,
                   cluster_cols    = TRUE,
                   clustering_method = "ward.D2",
                   display_numbers = TRUE,
                   number_format   = "%+.2f",
                   number_color    = "black",
                   fontsize_number = 7,
                   fontsize_row    = 10,
                   fontsize_col    = 8,
                   angle_col       = 45,
                   main            = paste0("Delta Jaccard: ", sub_a,
                                            " - ", sub_b,
                                            "\nvs ", ref_nm,
                                            " (blue = ", sub_a,
                                            " closer, red = ", sub_b,
                                            " closer)"),
                   border_color    = NA)
          dev.off()
        } else {
          plot_df <- as.data.frame(as.table(delta_plot))
          colnames(plot_df) <- c("type_row", "type_col", "Delta")

          p <- ggplot(plot_df, aes(x = type_col, y = type_row, fill = Delta)) +
            geom_tile(color = "grey90") +
            geom_text(aes(label = sprintf("%+.2f", Delta)), size = 2.5) +
            scale_fill_gradient2(low = "firebrick3", mid = "white",
                                 high = "steelblue3", midpoint = 0) +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
                  axis.text.y = element_text(size = 10)) +
            labs(title = paste0("Delta: ", sub_a, " - ", sub_b, " vs ", ref_nm),
                 x = ref_nm, y = "CellType")

          ggsave(fname, p,
                 width = max(14, length(relevant_cols) * 0.3), height = 10)
        }

        cat("    Saved:", basename(fname), "\n\n")
      }
    }
  }
}


# ==============================================================================
# 5. Save combined suggestions
# ==============================================================================

fwrite(all_suggestions, file.path(config$output_dir,
                                   "suggested_label_mapping.csv"))

cat("========================================================\n")
cat("  Outputs in:", config$output_dir, "\n")
cat("\n  Per-pair files:\n")
cat("    label_jaccard_matrix_*.csv       Raw Jaccard matrices\n")
cat("    label_jaccard_heatmap_*.pdf      Filtered heatmaps\n")
cat("    label_jaccard_heatmap_full_*.pdf Full heatmaps\n")
cat("\n  Condition comparison:\n")
cat("    label_jaccard_delta_*.csv        Delta matrices\n")
cat("    label_jaccard_delta_*.pdf        Diverging heatmaps\n")
cat("\n  suggested_label_mapping.csv        Best matches all pairs\n")
cat("\n  Interpretation of delta heatmap:\n")
cat("    Blue  = healthy (NS) is MORE similar to HSCA\n")
cat("    Red   = disease (VU) is MORE similar to HSCA\n")
cat("    White = no difference between conditions\n")
cat("    Large red values = VU has NOVEL cell states not in healthy\n")
cat("    Large blue values = VU has DIVERGED from the shared healthy state\n")
cat("========================================================\n")