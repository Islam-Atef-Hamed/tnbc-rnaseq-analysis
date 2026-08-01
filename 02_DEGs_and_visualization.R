# =======================================================================
# DESeq2: Differential Expression Analysis - TNBC vs Normal
# =======================================================================

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(dplyr)

# 1. تحميل البيانات
# ----------------------------------------------------------------------
message("=== [1/7] تحميل بيانات RNA-Seq ===")
data_se <- readRDS("TCGA_TNBC_vs_Normal_RawCounts.rds")

# التأكد من المجموعات
print(table(colData(data_se)$Group))

# 2. إعداد كائن DESeq2
# ----------------------------------------------------------------------
message("=== [2/7] إعداد كائن DESeq2 ===")

# نستخدم الـ counts الخام (unstranded assay)
countData <- assay(data_se, "unstranded")

# جدول التصميم
colData_df <- as.data.frame(colData(data_se))
colData_df$Group <- factor(colData_df$Group, levels = c("Normal", "TNBC"))

dds <- DESeqDataSetFromMatrix(
  countData = countData,
  colData = colData_df,
  design = ~ Group
)

# 3. الفلترة المسبقة (إزالة الجينات ذات العدد المنخفض)
# ----------------------------------------------------------------------
message("=== [3/7] الفلترة المسبقة ===")
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
message(paste("عدد الجينات بعد الفلترة:", nrow(dds)))

# 4. تشغيل DESeq2
# ----------------------------------------------------------------------
message("=== [4/7] تشغيل DESeq2 (ممكن ياخد 2-5 دقايق) ===")
dds <- DESeq(dds)

# 5. استخراج النتائج
# ----------------------------------------------------------------------
message("=== [5/7] استخراج النتائج ===")
res <- results(dds, contrast = c("Group", "TNBC", "Normal"))

# إضافة أسماء الجينات (من Ensembl ID لـ Gene Symbol)
# نستخدم البيانات المرفقة في SummarizedExperiment
gene_info <- as.data.frame(rowData(data_se))
if("gene_name" %in% colnames(gene_info)) {
  res$gene_name <- gene_info$gene_name[match(rownames(res), rownames(gene_info))]
} else {
  res$gene_name <- rownames(res)  # لو مفيش gene_name نستخدم الـ Ensembl ID
}

# 6. تصفية الجينات المرتفعة المعنوية (Upregulated DEGs)
# ----------------------------------------------------------------------
message("=== [6/7] تصفية DEGs ===")
# الشروط: padj < 0.01, log2FoldChange > 2 (Upregulated في TNBC)
sig_up <- subset(res, padj < 0.01 & log2FoldChange > 2 & !is.na(padj) & !is.na(gene_name))
sig_up <- sig_up[order(sig_up$padj), ]

# الجينات المنخفضة (Downregulated) — ممكن نحتاجها برضه
sig_down <- subset(res, padj < 0.01 & log2FoldChange < -2 & !is.na(padj) & !is.na(gene_name))
sig_down <- sig_down[order(sig_down$padj), ]

message(paste("✓ عدد الجينات المرتفعة (Upregulated):", nrow(sig_up)))
message(paste("✓ عدد الجينات المنخفضة (Downregulated):", nrow(sig_down)))

# 7. حفظ النتائج
# ----------------------------------------------------------------------
message("=== [7/7] حفظ الملفات ===")

# كل النتائج
write.csv(as.data.frame(res), "DESeq2_All_Results.csv")

# الجينات المرتفعة فقط (دي اللي هنستهدفها في الـ Drug Repurposing)
write.csv(as.data.frame(sig_up), "DEGs_Upregulated_TNBC.csv")

# الجينات المنخفضة
write.csv(as.data.frame(sig_down), "DEGs_Downregulated_TNBC.csv")

# كائن DESeq2 كامل للتحليلات اللاحقة
saveRDS(dds, "DESeq2_Object.rds")






# =======================================================================
# Visualization
# =======================================================================
# =======================================================================
# Visualization - ALTERNATIVE (بدون EnhancedVolcano)
# =======================================================================

library(ggplot2)
library(pheatmap)
library(dplyr)

# A. Volcano Plot يدوي
# ----------------------------------------------------------------------
message("رسم Volcano Plot...")

# تحضير البيانات
volcano_df <- as.data.frame(res)
volcano_df$gene_name <- ifelse(!is.na(volcano_df$gene_name), 
                               volcano_df$gene_name, 
                               rownames(volcano_df))

# تصنيف النقاط
volcano_df$significance <- "Not Significant"
volcano_df$significance[volcano_df$padj < 0.01 & abs(volcano_df$log2FoldChange) > 2] <- "Significant"
volcano_df$significance[volcano_df$padj < 0.01 & volcano_df$log2FoldChange > 2] <- "Upregulated"
volcano_df$significance[volcano_df$padj < 0.01 & volcano_df$log2FoldChange < -2] <- "Downregulated"

volcano_df$significance <- factor(volcano_df$significance, 
                                  levels = c("Downregulated", "Upregulated", 
                                             "Significant", "Not Significant"))

# اختيار أسماء الجينات للتسمية (أعلى 20 upregulated + أعلى 20 downregulated)
top_labels <- volcano_df %>%
  filter(significance %in% c("Upregulated", "Downregulated")) %>%
  arrange(padj) %>%
  head(30) %>%
  pull(gene_name)

volcano_df$label <- ifelse(volcano_df$gene_name %in% top_labels, 
                           volcano_df$gene_name, NA)

# الرسم
ggplot(volcano_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Downregulated" = "blue", 
                                "Upregulated" = "red", 
                                "Significant" = "green", 
                                "Not Significant" = "grey50")) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_text(aes(label = label), size = 2.5, vjust = -0.5, hjust = 0.5, 
            color = "black", na.rm = TRUE, check_overlap = TRUE) +
  labs(title = "TNBC vs Normal - Volcano Plot",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value",
       color = "") +
  theme_bw() +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

ggsave("Volcano_Plot.png", width = 12, height = 10, dpi = 300)
message("✓ تم حفظ Volcano_Plot.png")

# B. Heatmap لأعلى 50 جين مرتفع
# ----------------------------------------------------------------------
message("رسم Heatmap...")
top50_up <- head(sig_up, 50)

if(nrow(top50_up) > 0) {
  ntd <- normTransform(dds)
  mat <- assay(ntd)[rownames(top50_up), ]
  
  col_order <- order(colData(dds)$Group)
  mat <- mat[, col_order]
  
  row_labels <- top50_up$gene_name
  names(row_labels) <- rownames(top50_up)
  
  annotation_col <- data.frame(
    Group = colData(dds)$Group[col_order]
  )
  rownames(annotation_col) <- colnames(mat)
  
  pheatmap(mat,
           annotation_col = annotation_col,
           labels_row = row_labels,
           show_colnames = FALSE,
           cluster_cols = FALSE,
           cluster_rows = TRUE,
           scale = "row",
           color = colorRampPalette(c("blue", "white", "red"))(50),
           main = "Top 50 Upregulated Genes in TNBC",
           filename = "Heatmap_Top50_Upregulated.png",
           width = 14, height = 12)
  
  message("✓ تم حفظ Heatmap_Top50_Upregulated.png")
}

# C. PCA Plot
# ----------------------------------------------------------------------
message("رسم PCA...")
vsd <- vst(dds, blind = FALSE)

pcaData <- plotPCA(vsd, intgroup = "Group", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color = Group)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  theme_bw() +
  ggtitle("PCA: TNBC vs Normal")

ggsave("PCA_Plot.png", width = 8, height = 6, dpi = 300)
message("✓ تم حفظ PCA_Plot.png")

message("========================================")
message("✓✓✓ تم الانتهاء من التحليل والرسومات!")
message("========================================")

