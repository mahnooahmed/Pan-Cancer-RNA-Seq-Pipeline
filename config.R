# =============================================================================
#  config.R  —  Set all paths and parameters here before running any script
# =============================================================================
#  Edit the paths in this file to match your local directory structure.
#  All scripts source this file at the top; you should never need to edit
#  paths inside individual scripts.
# =============================================================================

# ── DIRECTORIES ──────────────────────────────────────────────────────────────
RNA_SEQ_DIR    <- "path/to/RNA-Seq"                   # root RNA-Seq folder
BY_CANCER_DIR  <- file.path(RNA_SEQ_DIR, "ByCancer")  # cancer-type subfolders
OUTPUT_DIR     <- "path/to/output"                    # DEG results and plots

# ── INPUT FILES ──────────────────────────────────────────────────────────────
MANIFEST_PATH  <- "path/to/gdc_manifest.txt"          # GDC download manifest
GENE_LIST_FILE <- "path/to/gene_list.xlsx"            # Excel with E1/E2/E3/DUB columns
MAF_FILE       <- "path/to/combined_pan_cancer.maf"   # somatic mutation MAF file

# ── DERIVED FILE PATHS (auto-built from RNA_SEQ_DIR; no need to edit) ────────
COUNTS_FILE    <- file.path(RNA_SEQ_DIR, "pan_cancer_counts_T.tsv")
TPM_FILE       <- file.path(RNA_SEQ_DIR, "pan_cancer_tpm_T.tsv")
GENE_ANNOT     <- file.path(RNA_SEQ_DIR, "gene_annotation.tsv")
METADATA_FILE  <- file.path(RNA_SEQ_DIR, "sample_metadata_clean.tsv")
UUID_MAP_FILE  <- file.path(RNA_SEQ_DIR, "uuid_to_barcode_map.tsv")
SAMPLE_TYPES   <- file.path(RNA_SEQ_DIR, "gdc_sample_types.tsv")
TPM_SUBSET_DIR <- file.path(RNA_SEQ_DIR, "TPM_subsets")

# ── COUNT COLUMN (GDC STAR-Counts TSV column to use) ─────────────────────────
#   "unstranded"       – strand-unaware raw counts  (recommended for most uses)
#   "stranded_first"   – first-strand protocol
#   "stranded_second"  – second-strand protocol
#   "tpm_unstranded"   – TPM values
COUNT_COLUMN <- "unstranded"

# ── MATRIX BUILDER TUNING ────────────────────────────────────────────────────
CHUNK_SIZE  <- 300L   # samples per write-chunk (lower if RAM is tight)
REPORT_N    <- 500L   # print progress every N files during build

# ── DESeq2 / ANALYSIS PARAMETERS ─────────────────────────────────────────────
PADJ_CUTOFF       <- 0.05
LOG2FC_CUTOFF     <- 1
MIN_COUNT_SUM     <- 10   # minimum total count across samples to keep a gene
TOP_HEATMAP_N     <- 50   # top N DEGs shown in heatmap
MIN_SAMPLES       <- 3    # minimum samples required per group to run DESeq2
MAX_SAMPLES_HEAT  <- 150  # max samples per group subsampled for heatmap

# ── UBL GENE CATEGORIES ──────────────────────────────────────────────────────
GENE_CATEGORIES <- c("E1", "E2", "E3", "DUB")
GENE_SHEET_NAME <- "COMBINED GENES"   # sheet name inside GENE_LIST_FILE

# ── SAMPLE TYPE LABELS (used in tumor/normal classification) ─────────────────
LOW_LIB_UUIDS <- character(0)
# Add sample UUIDs to exclude (low library size outliers), e.g.:
# LOW_LIB_UUIDS <- c("uuid-1", "uuid-2")
