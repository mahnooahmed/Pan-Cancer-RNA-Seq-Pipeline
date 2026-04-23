# =============================================================================
#  02_build_pan_cancer_matrices.R
#
#  Builds transposed pan-cancer counts and TPM matrices from GDC STAR-Counts
#  TSV files organised by cancer type. Performs three memory-efficient passes:
#    Pass 1 — identify non-zero genes across all samples
#    Pass 2 — write counts matrix   (samples × genes, integer)
#    Pass 3 — write TPM matrix      (samples × genes, float64)
#
#  Also handles sample-type annotation (Tumor / Normal / Metastatic etc.)
#  and removes flagged low-library-size samples defined in config.R.
#
#  Input  : ByCancer/{cancer}/*.augmented_star_gene_counts.tsv
#            gdc_sample_types.tsv   (from script 01)
#  Output : pan_cancer_counts_T.tsv
#            pan_cancer_tpm_T.tsv
#            gene_annotation.tsv
#            sample_metadata_clean.tsv
#
#  Run    : source("scripts/02_build_pan_cancer_matrices.R")
# =============================================================================

source("config.R")

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
  library(data.table)
})

STAR_ROWS <- c("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous")

# ── LOAD & PREPARE SAMPLE METADATA ───────────────────────────────────────────
cat("Loading sample metadata...\n")
sample_types <- fread(SAMPLE_TYPES, sep = "\t", header = TRUE)

sample_types[, sample_group := fcase(
  sample_type == "Primary Tumor",                                   "Primary_Tumor",
  sample_type == "Solid Tissue Normal",                             "Normal",
  sample_type == "Metastatic",                                      "Metastatic",
  sample_type == "Recurrent Tumor",                                 "Recurrent_Tumor",
  sample_type == "Primary Blood Derived Cancer - Peripheral Blood", "Primary_Tumor",
  sample_type == "Additional - New Primary",                        "Primary_Tumor",
  sample_type == "Additional Metastatic",                           "Metastatic",
  default = NA_character_
)]

cat("Sample group counts:\n")
print(table(sample_types$sample_group, useNA = "ifany"))

# Remove flagged low-library-size samples (configured in config.R)
if (length(LOW_LIB_UUIDS) > 0) {
  n_before <- nrow(sample_types)
  sample_types <- sample_types[!sample_uuid %in% LOW_LIB_UUIDS]
  cat(sprintf("Removed %d low-library samples.\n", n_before - nrow(sample_types)))
}

# Match samples to files on disk
all_files <- list.files(BY_CANCER_DIR, recursive = TRUE,
                        pattern = "augmented_star_gene_counts\\.tsv$",
                        full.names = TRUE)
file_lookup <- data.table(
  sample_uuid = sub("\\..*", "", basename(all_files)),
  filepath    = all_files
)
sample_types <- merge(sample_types, file_lookup, by = "sample_uuid", all.x = TRUE)

missing <- sum(is.na(sample_types$filepath))
if (missing > 0)
  cat(sprintf("WARNING: %d samples have no matching file — skipping.\n", missing))
sample_types <- sample_types[!is.na(filepath)]

n_samples <- nrow(sample_types)
filepaths <- sample_types$filepath
cat(sprintf("Samples ready to process: %d\n\n", n_samples))

# ── GENE METADATA FROM FIRST FILE ─────────────────────────────────────────────
ref <- fread(filepaths[1], sep = "\t", header = TRUE)
ref <- ref[!gene_id %in% STAR_ROWS]
ref[, gene_id := sub("\\..*", "", gene_id)]   # strip Ensembl version suffix

if (anyDuplicated(ref$gene_id))
  warning(sprintf("%d duplicate gene IDs after version stripping.", sum(duplicated(ref$gene_id))))

gene_ids   <- ref$gene_id
gene_names <- ref$gene_name
gene_types <- ref$gene_type
n_genes    <- length(gene_ids)
rm(ref); gc()
cat(sprintf("Genes in reference file: %d\n", n_genes))

# ── PASS 1/3 — IDENTIFY NON-ZERO GENES ───────────────────────────────────────
# Use numeric() (float64) to avoid integer overflow when summing across >10k samples
cat(sprintf("\n── Pass 1/3: Scanning %d files for expressed genes...\n", n_samples))
row_sums <- numeric(n_genes)

for (j in seq_len(n_samples)) {
  f <- fread(filepaths[j], sep = "\t", header = TRUE,
             select = c("gene_id", "unstranded"), showProgress = FALSE)
  f <- f[!gene_id %in% STAR_ROWS]
  row_sums <- row_sums + as.numeric(f$unstranded)

  if (j %% REPORT_N == 0L || j == n_samples)
    cat(sprintf("  [%d / %d]\n", j, n_samples))
}

keep   <- !is.na(row_sums) & row_sums > 0
n_keep <- sum(keep)
cat(sprintf("Keeping %d / %d genes (removed %d all-zero)\n\n",
            n_keep, n_genes, n_genes - n_keep))

gene_ids_k   <- gene_ids[keep]
gene_names_k <- gene_names[keep]
gene_types_k <- gene_types[keep]
rm(row_sums, gene_ids, gene_names, gene_types); gc()

# Save gene annotation
fwrite(
  data.table(gene_id = gene_ids_k, gene_name = gene_names_k, gene_type = gene_types_k),
  GENE_ANNOT, sep = "\t"
)
cat(sprintf("Gene annotation saved: %s\n", GENE_ANNOT))

# Shared header and metadata columns (reused in passes 2 and 3)
col_header <- paste(c("sample_id", "cancer_type", "sample_group", gene_ids_k),
                    collapse = "\t")
meta_cols <- sample_types[, .(
  sample_id    = sample_uuid,
  cancer_type  = sub("TCGA-", "", project_id),
  sample_group = sample_group
)]

make_chunks <- function(n, k) split(seq_len(n), ceiling(seq_len(n) / k))
chunks <- make_chunks(n_samples, CHUNK_SIZE)

# ── PASS 2/3 — WRITE COUNTS MATRIX ───────────────────────────────────────────
cat(sprintf("── Pass 2/3: Writing counts [chunk=%d, %d chunks]...\n",
            CHUNK_SIZE, length(chunks)))
writeLines(col_header, COUNTS_FILE)

for (ci in seq_along(chunks)) {
  idx       <- chunks[[ci]]
  chunk_mat <- matrix(0L, nrow = length(idx), ncol = n_keep)

  for (ji in seq_along(idx)) {
    f <- fread(filepaths[idx[ji]], sep = "\t", header = TRUE,
               select = c("gene_id", "unstranded"), showProgress = FALSE)
    f <- f[!gene_id %in% STAR_ROWS]
    chunk_mat[ji, ] <- f$unstranded[keep]
  }

  fwrite(cbind(meta_cols[idx], as.data.table(chunk_mat)),
         COUNTS_FILE, sep = "\t", append = TRUE, col.names = FALSE)
  rm(chunk_mat); gc()
  cat(sprintf("  Chunk %d / %d  (%d samples)\n", ci, length(chunks), length(idx)))
}
cat(sprintf("Counts matrix saved: %s\n\n", COUNTS_FILE))

# ── PASS 3/3 — WRITE TPM MATRIX ──────────────────────────────────────────────
cat(sprintf("── Pass 3/3: Writing TPM [chunk=%d, %d chunks]...\n",
            CHUNK_SIZE, length(chunks)))
writeLines(col_header, TPM_FILE)

for (ci in seq_along(chunks)) {
  idx       <- chunks[[ci]]
  chunk_mat <- matrix(0.0, nrow = length(idx), ncol = n_keep)

  for (ji in seq_along(idx)) {
    f <- fread(filepaths[idx[ji]], sep = "\t", header = TRUE,
               select = c("gene_id", "tpm_unstranded"), showProgress = FALSE)
    f <- f[!gene_id %in% STAR_ROWS]
    chunk_mat[ji, ] <- f$tpm_unstranded[keep]
  }

  fwrite(cbind(meta_cols[idx], as.data.table(chunk_mat)),
         TPM_FILE, sep = "\t", append = TRUE, col.names = FALSE)
  rm(chunk_mat); gc()
  cat(sprintf("  Chunk %d / %d  (%d samples)\n", ci, length(chunks), length(idx)))
}
cat(sprintf("TPM matrix saved: %s\n\n", TPM_FILE))

# ── SAVE CLEAN METADATA ───────────────────────────────────────────────────────
fwrite(meta_cols, METADATA_FILE, sep = "\t")
cat(sprintf("Metadata saved: %s\n\n", METADATA_FILE))

cat("═══ DONE ═══\n")
cat(sprintf("Dimensions : %d samples × %d genes\n", n_samples, n_keep))
cat(sprintf("Counts     : %s\n", COUNTS_FILE))
cat(sprintf("TPM        : %s\n", TPM_FILE))
cat(sprintf("Annotation : %s\n", GENE_ANNOT))
cat(sprintf("Metadata   : %s\n\n", METADATA_FILE))
cat("To reload as gene × sample matrix:\n")
cat("  cnt <- fread(COUNTS_FILE)\n")
cat("  mat <- t(as.matrix(cnt[, -(1:3)]))\n")
cat("  colnames(mat) <- cnt$sample_id\n")
