# =======================================================================
# Script 01: TCGA-BRCA Data Download and TNBC Sample Preparation
# Project : TNBC vs Normal RNA-seq Analysis
# Author  : Islam Atef Hamed
# Date    : 2025
# Output  : TCGA_TNBC_vs_Normal_RawCounts.rds
# Notes   : Requires GDC access and stable internet connection.
#           Download may take 20-60 minutes depending on connection speed.
# =======================================================================

library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)

# ── 1. Download Clinical Data (BCR XML) ──────────────────────────────────
message("=== [1/5] Querying clinical XML files ===")

query_clin <- GDCquery(
  project       = "TCGA-BRCA",
  data.category = "Clinical",
  data.format   = "BCR XML"
)

GDCdownload(query_clin, method = "api", files.per.chunk = 10)

message("=== [2/5] Extracting clinical patient data ===")
clinical_patient <- GDCprepare_clinic(query_clin, clinical.info = "patient")

# ── 2. Identify TNBC Patients ────────────────────────────────────────────
message("=== [3/5] Filtering TNBC patients ===")

# TNBC definition: ER-negative, PR-negative, HER2-negative
# Column names are taken directly from the TCGA-BRCA BCR XML schema
tnbc_patients <- clinical_patient %>%
  filter(
    breast_carcinoma_estrogen_receptor_status              == "Negative" &
    breast_carcinoma_progesterone_receptor_status          == "Negative" &
    lab_proc_her2_neu_immunohistochemistry_receptor_status == "Negative"
  ) %>%
  pull(bcr_patient_barcode) %>%
  unique()

message(paste("TNBC patients identified:", length(tnbc_patients)))

# Safety check — stop early if no TNBC patients are found
if (length(tnbc_patients) == 0) {
  stop("No TNBC patients found. Please verify column names in the clinical data.")
}

# ── 3. Query RNA-seq Metadata ────────────────────────────────────────────
message("=== [4/5] Querying RNA-seq metadata ===")

query_rnaseq <- GDCquery(
  project       = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

rnaseq_metadata <- getResults(query_rnaseq)

# Keep all Normal samples + only Tumor samples from TNBC patients
target_cases <- rnaseq_metadata %>%
  filter(
    sample_type == "Solid Tissue Normal" |
    (sample_type == "Primary Tumor" &
     substr(cases, 1, 12) %in% tnbc_patients)
  ) %>%
  pull(cases)

# Summary of selected samples
normal_count <- sum(
  rnaseq_metadata$sample_type[rnaseq_metadata$cases %in% target_cases] ==
  "Solid Tissue Normal"
)
tumor_count <- sum(
  rnaseq_metadata$sample_type[rnaseq_metadata$cases %in% target_cases] ==
  "Primary Tumor"
)

message(paste("Total samples selected:", length(target_cases)))
message(paste("  Normal :", normal_count))
message(paste("  TNBC   :", tumor_count))

# Safety check — stop if no TNBC tumor samples were matched
if (tumor_count == 0) {
  stop("No TNBC tumor samples selected. Check patient barcode matching.")
}

# ── 4. Download and Prepare Final Dataset ────────────────────────────────
message("=== [5/5] Downloading and preparing final dataset ===")

query_final <- GDCquery(
  project       = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode       = target_cases
)

GDCdownload(query_final, method = "api", files.per.chunk = 10)
data_se <- GDCprepare(query_final)

# Add a Group column: "Normal" or "TNBC"
colData(data_se)$Group <- ifelse(
  colData(data_se)$sample_type == "Solid Tissue Normal",
  "Normal",
  "TNBC"
)

# Final sample summary
message("Final sample summary:")
print(table(colData(data_se)$Group))

# Save as RDS — input for script 02
saveRDS(data_se, "TCGA_TNBC_vs_Normal_RawCounts.rds")
message("Saved: TCGA_TNBC_vs_Normal_RawCounts.rds")
message("Ready for script 02: DEG analysis")
