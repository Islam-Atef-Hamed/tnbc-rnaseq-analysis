# =======================================================================
# Script 03: GO and KEGG Enrichment Analysis
# Project : TNBC vs Normal RNA-seq Analysis
# Author  : Islam Atef Hamed
# Date    : 2025
# Input   : DESeq2_Object.rds (from script 02)
# Output  : GO/KEGG enrichment tables and plots
# =======================================================================

library(DESeq2)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# ── 0. Create output directories if they don't exist ─────────────────────
dir.create("results/tables",  recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# ── 1. Load DESeq2 Results ───────────────────────────────────────────────
message("=== [1/6] Loading DESeq2 object ===")
dds <- readRDS("DESeq2_Object.rds")
res <- results(dds, contrast = c("Group", "TNBC", "Normal"))

# Filter upregulated DEGs (padj < 0.01, log2FC > 2)
sig_up     <- subset(res, padj < 0.01 & log2FoldChange > 2 & !is.na(padj))
sig_up_ids <- rownames(sig_up)
message(paste("Upregulated DEGs:", length(sig_up_ids)))

# ── 2. Build Background Gene Universe (critical for unbiased KEGG) ───────
message("=== [2/6] Building background gene universe ===")

# Using all genes tested in DESeq2 as background (not the whole genome).
# This prevents inflation of large, non-cancer KEGG pathways such as
# "Neuroactive ligand-receptor interaction" or "Alcoholism".
all_gene_ids       <- rownames(results(dds))
all_gene_ids_clean <- sub("\\..*$", "", all_gene_ids)  # strip Ensembl version

all_entrez <- bitr(all_gene_ids_clean,
                   fromType = "ENSEMBL",
                   toType   = "ENTREZID",
                   OrgDb    = org.Hs.eg.db)
message(paste("Background universe size:", nrow(all_entrez), "genes"))

# ── 3. Convert Upregulated Gene IDs ──────────────────────────────────────
message("=== [3/6] Converting upregulated gene IDs ===")

# Strip Ensembl version numbers (e.g. ENSG00000141510.18 -> ENSG00000141510)
sig_up_ids_clean <- sub("\\..*$", "", sig_up_ids)

# Ensembl -> Gene Symbol
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys      = sig_up_ids_clean,
                       column    = "SYMBOL",
                       keytype   = "ENSEMBL",
                       multiVals = "first")
gene_symbols <- unique(gene_symbols[!is.na(gene_symbols)])
message(paste("Gene symbols mapped:", length(gene_symbols)))

# Gene Symbol -> Entrez ID (required for KEGG)
gene_entrez <- bitr(gene_symbols,
                    fromType = "SYMBOL",
                    toType   = "ENTREZID",
                    OrgDb    = org.Hs.eg.db)
message(paste("Entrez IDs obtained:", nrow(gene_entrez)))

unmapped <- setdiff(gene_symbols, gene_entrez$SYMBOL)
if (length(unmapped) > 0)
  message(paste("Unmapped genes:", length(unmapped)))

# ── 4. GO Enrichment ─────────────────────────────────────────────────────
message("=== [4/6] GO Enrichment Analysis ===")

run_go <- function(ont_type) {
  enrichGO(gene          = gene_entrez$ENTREZID,
           OrgDb         = org.Hs.eg.db,
           ont           = ont_type,
           pAdjustMethod = "BH",
           pvalueCutoff  = 0.05,
           qvalueCutoff  = 0.05,
           readable      = TRUE)
}

ego_bp <- run_go("BP")
ego_mf <- run_go("MF")
ego_cc <- run_go("CC")

write.csv(as.data.frame(ego_bp), "results/tables/GO_BP_Enrichment.csv")
write.csv(as.data.frame(ego_mf), "results/tables/GO_MF_Enrichment.csv")
write.csv(as.data.frame(ego_cc), "results/tables/GO_CC_Enrichment.csv")

message(paste("GO-BP terms:", nrow(ego_bp)))
message(paste("GO-MF terms:", nrow(ego_mf)))
message(paste("GO-CC terms:", nrow(ego_cc)))

# ── 5. KEGG Enrichment ───────────────────────────────────────────────────
message("=== [5/6] KEGG Pathway Analysis ===")

kegg <- enrichKEGG(gene          = gene_entrez$ENTREZID,
                   universe      = all_entrez$ENTREZID,  # corrected background
                   organism      = "hsa",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)

kegg <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
write.csv(as.data.frame(kegg), "results/tables/KEGG_Enrichment.csv")
message(paste("KEGG pathways (raw):", nrow(kegg)))

# ── 5b. Filter non-cancer-relevant KEGG pathways ─────────────────────────
# Some pathways appear due to shared gene families (e.g. histones in
# Alcoholism/Lupus, collagens in Protein digestion) rather than true
# biological relevance to breast cancer. We remove these and save both
# the raw and filtered versions for full transparency.
non_cancer_terms <- paste(
  "Neuroactive", "Alcoholism", "Salivary", "Olfactory",
  "lupus", "Cornified", "digestion", "absorption",
  sep = "|"
)

kegg_filtered <- kegg
kegg_filtered@result <- kegg_filtered@result %>%
  filter(!grepl(non_cancer_terms, Description, ignore.case = TRUE))

write.csv(as.data.frame(kegg_filtered),
          "results/tables/KEGG_Enrichment_filtered.csv")
message(paste("KEGG pathways (cancer-relevant):", nrow(kegg_filtered)))

# ── 6. Visualization ─────────────────────────────────────────────────────
message("=== [6/6] Generating plots ===")

save_dotplot <- function(obj, title, filename, w = 10, h = 8) {
  p <- dotplot(obj, showCategory = 15) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  ggsave(file.path("results/figures", filename), p,
         width = w, height = h, dpi = 300)
  message(paste("Saved:", filename))
}

# GO plots
save_dotplot(ego_bp, "GO Biological Process - Top 15",  "GO_BP_Dotplot.png")
save_dotplot(ego_mf, "GO Molecular Function - Top 15",  "GO_MF_Dotplot.png")
save_dotplot(ego_cc, "GO Cellular Component - Top 15",  "GO_CC_Dotplot.png")

# KEGG — raw (all significant pathways)
save_dotplot(kegg,
             "KEGG Pathways - Top 15 (all)",
             "KEGG_Dotplot_raw.png")

# KEGG — filtered (cancer-relevant pathways only)
save_dotplot(kegg_filtered,
             "KEGG Pathways - Top 15 (Cancer-relevant)",
             "KEGG_Dotplot_filtered.png")

# Gene-concept network — top 5 BP terms
cnet <- cnetplot(ego_bp, showCategory = 5)
ggsave("results/figures/GO_BP_Cnetplot.png", cnet,
       width = 12, height = 10, dpi = 300)
message("Saved: GO_BP_Cnetplot.png")

message("=================================================")
message("Enrichment analysis complete.")
message("Tables  -> results/tables/")
message("Figures -> results/figures/")
message("=================================================")

