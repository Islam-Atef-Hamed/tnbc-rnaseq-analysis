# RNA-seq Differential Expression Analysis: TNBC vs Normal Breast Tissue

A complete end-to-end transcriptomic analysis of **Triple-Negative Breast Cancer (TNBC)** using RNA-seq data from TCGA, identifying key dysregulated genes and enriched biological pathways driving TNBC progression.

---

## Biological Question

Which genes and pathways are significantly dysregulated in Triple-Negative Breast Cancer compared to normal breast tissue, and what do they reveal about the molecular mechanisms driving this aggressive cancer subtype?

---

## Dataset

- **Source:** TCGA (The Cancer Genome Atlas) via GEO
- **Comparison:** TNBC tumor samples vs Normal breast tissue
- **Data type:** RNA-seq raw counts
- **Tool:** DESeq2 (Bioconductor)

---

## Analysis Workflow

```
TCGA Data Download
       ↓
Quality Control & Filtering (counts ≥ 10)
       ↓
DESeq2 Differential Expression
(padj < 0.01, |log2FC| > 2)
       ↓
Visualization: PCA · Volcano · Heatmap
       ↓
GO Enrichment (BP · MF · CC)
       ↓
KEGG Pathway Analysis
       ↓
STRING Network Preparation
```

---

## Key Results

### Top Upregulated Genes in TNBC

| Gene | Role in Cancer |
|------|---------------|
| **AURKA** | Mitotic kinase — overexpressed in aggressive breast cancers |
| **TPX2** | Spindle assembly — drives chromosomal instability |
| **PLK1** | Cell cycle regulator — established therapeutic target |
| **FOXM1** | Transcription factor — master regulator of TNBC proliferation |
| **EZH2** | Epigenetic silencer — drives TNBC metastasis |
| **CDK1** | Cell cycle kinase — essential for mitotic entry |
| **CCNB1/CCNB2** | Cyclin B — drives G2/M transition |

### Enrichment Summary

**GO Biological Process:** Dominated by cell division — organelle fission, nuclear division, chromosome segregation (padj < 1e-12)

**GO Cellular Component:** Chromosomal region, spindle apparatus, kinetochore

**KEGG Pathways:** Cell cycle, transcriptional misregulation in cancer, cytokine signaling

---

## Figures

| Figure | Description |
|--------|-------------|
| `PCA_Plot.png` | Clear separation of TNBC vs Normal (PC1: 43% variance) |
| `Volcano_Plot.png` | DEGs with top genes labeled |
| `Heatmap_Top50_Upregulated.png` | Expression patterns of top 50 upregulated genes |
| `GO_BP_Dotplot.png` | GO Biological Process enrichment |
| `GO_MF_Dotplot.png` | GO Molecular Function enrichment |
| `GO_CC_Dotplot.png` | GO Cellular Component enrichment |
| `GO_BP_Cnetplot.png` | Gene-concept network for BP |
| `KEGG_Dotplot.png` | KEGG pathway enrichment |

---

## Project Structure

```
tnbc-rnaseq-analysis/
│
├── scripts/
│   ├── 01_download_and_prepare_data.R     # TCGA data retrieval
│   ├── 02_DEGs_and_visualization.R        # DESeq2 + PCA/Volcano/Heatmap
│   ├── 03_enrichment_analysis.R           # GO + KEGG enrichment
│   └── 04_prepare_genes_list_for_string.R # STRING network prep
│
├── results/
│   ├── figures/                           # All plots (PNG, 300 DPI)
│   └── tables/
│       ├── DESeq2_All_Results.csv
│       ├── DEGs_Upregulated_TNBC.csv
│       ├── DEGs_Downregulated_TNBC.csv
│       └── Top50_Upregulated_for_STRING.txt
│
├── data/                                  # Raw data (not tracked by git)
└── README.md
```

---

## Requirements

```r
# Bioconductor
BiocManager::install(c("DESeq2", "clusterProfiler", "enrichplot",
                       "org.Hs.eg.db", "TCGAbiolinks"))

# CRAN
install.packages(c("ggplot2", "pheatmap", "dplyr", "EnhancedVolcano"))
```

R version ≥ 4.2.0 recommended.

---

## How to Reproduce

```r
# Run scripts in order:
source("scripts/01_download_and_prepare_data.R")
source("scripts/02_DEGs_and_visualization.R")
source("scripts/03_enrichment_analysis.R")
source("scripts/04_prepare_genes_list_for_string.R")
```

---

## Biological Interpretation

The transcriptomic landscape of TNBC is dominated by **cell cycle dysregulation**. The top upregulated genes (AURKA, PLK1, TPX2, FOXM1, CDK1) form a tightly interconnected network controlling mitotic progression. This signature is consistent with TNBC's hallmark of rapid, uncontrolled proliferation.

Notably, **EZH2** upregulation links proliferative signaling to epigenetic reprogramming — a known driver of TNBC metastasis and therapy resistance. The presence of **cytokine signaling** (KEGG) suggests active tumor-immune crosstalk, relevant to immunotherapy response prediction.

**Potential therapeutic targets identified:**
- AURKA → Alisertib (phase II trials in TNBC)
- PLK1 → Volasertib
- EZH2 → Tazemetostat (FDA approved)
- FOXM1 → Indirect targeting via CDK inhibitors

---

## Next Steps

- [ ] Survival analysis correlating DEG expression with TCGA clinical data
- [ ] Molecular docking of top drug candidates against AURKA/PLK1 structures
- [ ] Single-cell RNA-seq to resolve tumor heterogeneity
- [ ] Integrate with ATAC-seq for chromatin accessibility analysis

---

## Author

**Islam Atef Hamed**
Bioinformatics Diploma Student | Cancer Transcriptomics & Multi-omics
[GitHub](https://github.com/Islam-Atef-Hamed)
