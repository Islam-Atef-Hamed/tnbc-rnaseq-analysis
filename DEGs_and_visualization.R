# =======================================================================
# Script 02: Differential Expression Analysis & Visualization
# Project : TNBC vs Normal RNA-seq Analysis
# Author  : Islam Atef Hamed
# Date    : 2025
# Input   : TCGA_TNBC_vs_Normal_RawCounts.rds (from script 01)
# Output  : DESeq2 results tables + PCA, Volcano, Heatmap plots
# =======================================================================

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(dplyr)

# ── 0. Create output directories ─────────────────────────────────────────
dir.create("results/tables",  recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# ── 1. Load Data ─────────────────────────────────────────────────────────
message("=== [1/7] Loading RNA-seq data ===")
data_se <- readRDS("TCGA_TNBC_vs_Normal_RawCounts.rds")

# Verify group labels
print(table(colData(data_se)$Group))

# ── 2. Build DESeq2 Object ───────────────────────────────────────────────
message("=== [2/7] Building DESeq2 object ===")

# Extract raw counts from the unstranded assay
countData <- assay(data_se, "unstranded")

# Build the experimental design table
colData_df        <- as.data.frame(colData(data_se))
colData_df$Group  <- factor(colData_df$Group, levels = c("Normal", "TNBC"))

dds <- DESeqDataSetFromMatrix(
  countData = countData,
  colData   = colData_df,
  design    = ~ Group
)

# ── 3. Pre-filtering ─────────────────────────────────────────────────────
message("=== [3/7] Pre-filtering low-count genes ===")

# Remove genes with fewer than 10 total counts across all samples
keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]
message(paste("Genes retained after filtering:", nrow(dds)))

# ── 4. Run DESeq2 ────────────────────────────────────────────────────────
message("=== [4/7] Running DESeq2 (may take 2-5 minutes) ===")
dds <- DESeq(dds)

# ── 5. Extract Results ───────────────────────────────────────────────────
message("=== [5/7] Extracting results ===")
res <- results(dds, contrast = c("Group", "TNBC", "Normal"))

# Map Ensembl IDs to Gene Symbols using rowData from SummarizedExperiment
gene_info <- as.data.frame(rowData(data_se))
if ("gene_name" %in% colnames(gene_info)) {
  res$gene_name <- gene_info$gene_name[match(rownames(res), rownames(gene_info))]
} else {
  # Fall back to Ensembl IDs if gene symbols are unavailable
  res$gene_name <- rownames(res)
}

# ── 6. Filter DEGs ───────────────────────────────────────────────────────
message("=== [6/7] Filtering significant DEGs ===")

# Thresholds: padj < 0.01, |log2FC| > 2
sig_up <- subset(res,
                 padj < 0.01 & log2FoldChange > 2 &
                 !is.na(padj) & !is.na(gene_name))
sig_up <- sig_up[order(sig_up$padj), ]

sig_down <- subset(res,
                   padj < 0.01 & log2FoldChange < -2 &
                   !is.na(padj) & !is.na(gene_name))
sig_down <- sig_down[order(sig_down$padj), ]

message(paste("Upregulated DEGs  :", nrow(sig_up)))
message(paste("Downregulated DEGs:", nrow(sig_down)))

# ── 7. Save Results ──────────────────────────────────────────────────────
message("=== [7/7] Saving results ===")

# Full DESeq2 results table
write.csv(as.data.frame(res),
          "results/tables/DESeq2_All_Results.csv")

# Upregulated DEGs — candidate drug targets for downstream analysis
write.csv(as.data.frame(sig_up),
          "results/tables/DEGs_Upregulated_TNBC.csv")

# Downregulated DEGs
write.csv(as.data.frame(sig_down),
          "results/tables/DEGs_Downregulated_TNBC.csv")

# Full DESeq2 object — required as input for script 03
saveRDS(dds, "DESeq2_Object.rds")

message("Results saved to results/tables/")

# =======================================================================
# Visualization
# =======================================================================

# ── A. Volcano Plot ──────────────────────────────────────────────────────
message("Generating Volcano Plot...")

# Prepare data frame
volcano_df           <- as.data.frame(res)
volcano_df$gene_name <- ifelse(!is.na(volcano_df$gene_name),
                               volcano_df$gene_name,
                               rownames(volcano_df))

# Classify each gene
volcano_df$significance <- "Not Significant"
volcano_df$significance[
  volcano_df$padj < 0.01 & volcano_df$log2FoldChange >  2] <- "Upregulated"
volcano_df$significance[
  volcano_df$padj < 0.01 & volcano_df$log2FoldChange < -2] <- "Downregulated"

volcano_df$significance <- factor(
  volcano_df$significance,
  levels = c("Downregulated", "Upregulated", "Not Significant")
)

# Select top 30 genes by adjusted p-value for labelling
top_labels <- volcano_df %>%
  filter(significance %in% c("Upregulated", "Downregulated")) %>%
  arrange(padj) %>%
  head(30) %>%
  pull(gene_name)

volcano_df$label <- ifelse(volcano_df$gene_name %in% top_labels,
                           volcano_df$gene_name, NA)

# Plot
ggplot(volcano_df,
       aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Downregulated"  = "blue",
                                "Upregulated"    = "red",
                                "Not Significant" = "grey50")) +
  geom_vline(xintercept = c(-2, 2),
             linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.01),
             linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_text(aes(label = label), size = 2.5, vjust = -0.5, hjust = 0.5,
            color = "black", na.rm = TRUE, check_overlap = TRUE) +
  labs(title = "TNBC vs Normal - Volcano Plot",
       x     = "Log2 Fold Change",
       y     = "-Log10 Adjusted P-value",
       color = "") +
  theme_bw() +
  theme(legend.position  = "right",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

ggsave("results/figures/Volcano_Plot.png", width = 12, height = 10, dpi = 300)
message("Saved: Volcano_Plot.png")

# ── B. Heatmap — Top 50 Upregulated Genes ───────────────────────────────
message("Generating Heatmap...")
top50_up <- head(sig_up, 50)

if (nrow(top50_up) > 0) {
  # Log-normalise counts for visualisation (not used in DE testing)
  ntd <- normTransform(dds)
  mat <- assay(ntd)[rownames(top50_up), ]

  # Order columns: Normal first, TNBC second
  col_order <- order(colData(dds)$Group)
  mat       <- mat[, col_order]

  # Use gene symbols as row labels
  row_labels              <- top50_up$gene_name
  names(row_labels)       <- rownames(top50_up)

  annotation_col          <- data.frame(Group = colData(dds)$Group[col_order])
  rownames(annotation_col) <- colnames(mat)

  pheatmap(mat,
           annotation_col = annotation_col,
           labels_row     = row_labels,
           show_colnames  = FALSE,
           cluster_cols   = FALSE,
           cluster_rows   = TRUE,
           scale          = "row",
           color          = colorRampPalette(c("blue", "white", "red"))(50),
           main           = "Top 50 Upregulated Genes in TNBC",
           filename       = "results/figures/Heatmap_Top50_Upregulated.png",
           width = 14, height = 12)

  message("Saved: Heatmap_Top50_Upregulated.png")
}

# ── C. PCA Plot ──────────────────────────────────────────────────────────
message("Generating PCA Plot...")

# Variance-stabilising transformation for PCA
vsd        <- vst(dds, blind = FALSE)
pcaData    <- plotPCA(vsd, intgroup = "Group", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color = Group)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  theme_bw() +
  ggtitle("PCA: TNBC vs Normal")

ggsave("results/figures/PCA_Plot.png", width = 8, height = 6, dpi = 300)
message("Saved: PCA_Plot.png")

message("=================================================")
message("DEG analysis and visualization complete.")
message("Tables  -> results/tables/")
message("Figures -> results/figures/")
message("=================================================")
