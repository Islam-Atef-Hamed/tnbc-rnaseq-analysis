# =======================================================================
# Script 04: Prepare Gene List for STRING Network Analysis
# Project : TNBC vs Normal RNA-seq Analysis
# Author  : Islam Atef Hamed
# Date    : 2025
# Input   : results/tables/DEGs_Upregulated_TNBC.csv (from script 02)
# Output  : results/tables/Top500_Upregulated_for_STRING.txt
# Usage   : Upload the output .txt file to https://string-db.org
#           Paste all gene names -> Organism: Homo sapiens
# =======================================================================

library(dplyr)

# ── 0. Create output directory ───────────────────────────────────────────
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# ── 1. Load Upregulated DEGs ─────────────────────────────────────────────
message("Loading upregulated DEGs...")
sig_up <- read.csv("results/tables/DEGs_Upregulated_TNBC.csv")
message(paste("Total upregulated DEGs loaded:", nrow(sig_up)))

# ── 2. Select Top Genes for STRING ───────────────────────────────────────
# STRING can handle ~2000 proteins per query, but networks become hard to
# interpret beyond 500 nodes. We use the top 500 genes ranked by padj.
# Change head(500) to head(1000) if you want a broader network.
top_genes_for_string <- sig_up %>%
  arrange(padj) %>%
  head(500) %>%
  pull(gene_name)

message(paste("Genes selected for STRING:", length(top_genes_for_string)))
message("Top 10 genes:")
print(head(top_genes_for_string, 10))

# ── 3. Save Gene List ────────────────────────────────────────────────────
# One gene symbol per line — ready to paste into STRING search box
output_path <- "results/tables/Top500_Upregulated_for_STRING.txt"
writeLines(top_genes_for_string, output_path)
message(paste("Saved:", output_path))

message("=================================================")
message("Next step: Upload to https://string-db.org")
message("  -> Search -> Multiple proteins")
message("  -> Organism: Homo sapiens")
message("  -> Minimum interaction score: 0.4 (medium confidence)")
message("=================================================")
