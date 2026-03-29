# Quick gene name diagnostic
suppressPackageStartupMessages({
  library(data.table)
})

# Panel genes from first region
gm <- fread("./proseg_results/XE765/r_seurat/gene-metadata.csv.gz")
panel <- gm$gene

cat("=== PANEL GENES (first 30) ===\n")
cat(paste(head(panel, 30), collapse = ", "), "\n\n")

# Reference genes
ref <- readRDS("./reference/consensus_reference.rds")
ref_genes <- rownames(ref)

cat("=== REFERENCE GENES (first 30) ===\n")
cat(paste(head(ref_genes, 30), collapse = ", "), "\n\n")

# Direct overlap
shared <- intersect(panel, ref_genes)
cat("Direct overlap:", length(shared), "/", length(panel), "\n\n")

# Case-insensitive overlap
shared_ci <- intersect(toupper(panel), toupper(ref_genes))
cat("Case-insensitive overlap:", length(shared_ci), "/", length(panel), "\n\n")

# Show examples of near-matches
cat("=== NEAR MATCHES (first 20 panel genes) ===\n")
for (g in head(panel, 20)) {
  # Check exact
  exact <- g %in% ref_genes
  # Check uppercase
  upper_match <- toupper(g) %in% toupper(ref_genes)
  if (upper_match && !exact) {
    ref_version <- ref_genes[toupper(ref_genes) == toupper(g)][1]
    cat(sprintf("  %-15s -> %-15s (case mismatch)\n", g, ref_version))
  } else if (exact) {
    cat(sprintf("  %-15s -> EXACT MATCH\n", g))
  } else {
    cat(sprintf("  %-15s -> NOT FOUND\n", g))
  }
}

rm(ref); gc(verbose = FALSE)