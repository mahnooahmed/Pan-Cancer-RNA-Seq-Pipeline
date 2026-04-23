# =============================================================================
#  03_verify_counts_matrix.R
#
#  Quality-checks the pan-cancer counts matrix produced by script 02.
#  Reads the file in memory-efficient chunks to avoid loading the full
#  matrix into RAM.
#
#  Checks performed:
#    1. Matrix dimensions vs metadata
#    2. Duplicate gene IDs
#    3. Missing values (NAs)
#    4. Negative values
#    5. All-zero sample columns
#    6. All-zero gene rows
#    7. Library size distribution + low-library outliers
#    8. Housekeeping gene spot-check (ACTB)
#
#  Input  : pan_cancer_counts_T.tsv + sample_metadata_clean.tsv
#  Output : console report (PASS / FAIL per check)
#
#  Run    : source("scripts/03_verify_counts_matrix.R")
# =============================================================================

source("config.R")

suppressPackageStartupMessages({
  library(data.table)
})

READ_CHUNK <- 2000L   # genes per pass chunk; lower if RAM is tight

# ── LOAD METADATA ─────────────────────────────────────────────────────────────
meta <- fread(METADATA_FILE, sep = "\t", header = TRUE)

# ── READ HEADER ONLY ──────────────────────────────────────────────────────────
cat("Reading header...\n")
header     <- fread(COUNTS_FILE, sep = "\t", nrows = 0L)
all_cols   <- colnames(header)
# Identify count columns: everything after the three metadata columns
meta_cols  <- c("sample_id", "cancer_type", "sample_group")
count_cols <- setdiff(all_cols, meta_cols)
n_samples  <- length(count_cols)
rm(header); gc()

# ── COUNT TOTAL ROWS ──────────────────────────────────────────────────────────
cat("Counting gene rows...\n")
con     <- file(COUNTS_FILE, "r")
n_lines <- 0L
repeat {
  batch <- readLines(con, n = 10000L, warn = FALSE)
  if (length(batch) == 0) break
  n_lines <- n_lines + length(batch)
}
close(con)
n_genes <- n_lines - 1L   # subtract header row
cat(sprintf("Matrix: %d genes × %d samples\n\n", n_genes, n_samples))

# ── INITIALISE ACCUMULATORS ───────────────────────────────────────────────────
col_sums    <- numeric(n_samples)
names(col_sums) <- count_cols
row_zero    <- 0L
na_total    <- 0L
neg_total   <- 0L
gene_ids    <- character(n_genes)
genes_seen  <- 0L
actb_counts <- NULL
skip_rows   <- 1L   # start after header

# ── CHUNKED SCAN ──────────────────────────────────────────────────────────────
cat("Scanning matrix in chunks...\n")
repeat {
  chunk <- fread(COUNTS_FILE, sep = "\t", header = FALSE,
                 skip       = skip_rows,
                 nrows      = READ_CHUNK,
                 col.names  = all_cols,
                 showProgress = FALSE)
  if (nrow(chunk) == 0L) break

  idx_from <- genes_seen + 1L
  idx_to   <- genes_seen + nrow(chunk)
  gene_ids[idx_from:idx_to] <- chunk$sample_id   # first col = gene_id in gene-row files

  # Extract numeric block only
  cmat <- as.matrix(chunk[, ..count_cols])

  na_total  <- na_total  + sum(is.na(cmat))
  neg_total <- neg_total + sum(cmat < 0L, na.rm = TRUE)
  col_sums  <- col_sums  + colSums(cmat, na.rm = TRUE)
  row_zero  <- row_zero  + sum(rowSums(cmat, na.rm = TRUE) == 0L)

  # ACTB spot-check (ENSG00000075624)
  if (is.null(actb_counts)) {
    actb_row <- grep("ENSG00000075624", chunk[[1]])
    if (length(actb_row) > 0)
      actb_counts <- as.numeric(cmat[actb_row[1], ])
  }

  skip_rows  <- skip_rows  + nrow(chunk)
  genes_seen <- genes_seen + nrow(chunk)
  cat(sprintf("  [%d / %d] genes scanned\n", genes_seen, n_genes))
  rm(chunk, cmat); gc()
}

# ── REPORT ────────────────────────────────────────────────────────────────────
cat("\n══════════════════════════════════════\n")
cat("  COUNTS MATRIX QUALITY REPORT\n")
cat("══════════════════════════════════════\n")

cat(sprintf("\n[1] Shape\n"))
cat(sprintf("    Matrix   : %d genes × %d samples\n", n_genes, n_samples))
cat(sprintf("    Metadata : %d samples\n", nrow(meta)))
cols_match <- all(count_cols == meta$sample_id)
cat(sprintf("    Column order matches metadata: %s\n", cols_match))

cat(sprintf("\n[2] Duplicate gene IDs\n"))
dup_genes <- gene_ids[duplicated(gene_ids)]
if (length(dup_genes) == 0) cat("    OK — none\n") else {
  cat(sprintf("    WARNING: %d duplicates\n", length(dup_genes)))
  print(head(dup_genes, 20))
}

cat(sprintf("\n[3] Missing values (NAs)\n"))
cat(sprintf("    %d %s\n", na_total, ifelse(na_total == 0, "(OK)", "<- WARNING")))

cat(sprintf("\n[4] Negative values\n"))
cat(sprintf("    %d %s\n", neg_total, ifelse(neg_total == 0, "(OK)", "<- WARNING")))

cat(sprintf("\n[5] All-zero sample columns\n"))
zero_cols <- names(col_sums[col_sums == 0])
if (length(zero_cols) == 0) cat("    OK — none\n") else {
  cat(sprintf("    WARNING: %d all-zero samples\n", length(zero_cols)))
  print(head(zero_cols, 20))
}

cat(sprintf("\n[6] All-zero gene rows\n"))
cat(sprintf("    %d / %d (%.1f%%)\n", row_zero, n_genes, 100 * row_zero / n_genes))
if (row_zero / n_genes > 0.3)
  cat("    WARNING: >30% all-zero genes is unusually high\n")

cat(sprintf("\n[7] Library size distribution\n"))
cat(sprintf("    Min    : %s\n", format(min(col_sums),    big.mark = ",")))
cat(sprintf("    Median : %s\n", format(median(col_sums), big.mark = ",")))
cat(sprintf("    Max    : %s\n", format(max(col_sums),    big.mark = ",")))
low_lib <- names(col_sums[col_sums < 0.1 * median(col_sums)])
if (length(low_lib) == 0) {
  cat("    OK — no low-library-size outliers\n")
} else {
  cat(sprintf("    WARNING: %d low-library samples (< 10%% of median):\n", length(low_lib)))
  low_df <- data.frame(
    sample_id    = low_lib,
    library_size = col_sums[low_lib],
    cancer_type  = meta$cancer_type[match(low_lib, meta$sample_id)]
  )
  print(low_df[order(low_df$library_size), ], row.names = FALSE)
}

cat(sprintf("\n[8] Housekeeping gene spot-check (ACTB — ENSG00000075624)\n"))
if (is.null(actb_counts)) {
  cat("    ACTB not found — verify gene ID format in matrix\n")
} else {
  cat(sprintf("    Median count  : %s\n", format(median(actb_counts), big.mark = ",")))
  cat(sprintf("    Zero samples  : %d\n", sum(actb_counts == 0)))
  if (median(actb_counts) < 100)
    cat("    WARNING: ACTB median very low — check COUNT_COLUMN in config.R\n")
}

cat("\n══════ SUMMARY ══════\n")
checks <- c(
  "Shape OK"           = cols_match,
  "No duplicate genes" = length(dup_genes) == 0,
  "No NAs"             = na_total  == 0,
  "No negatives"       = neg_total == 0,
  "No zero samples"    = length(zero_cols) == 0,
  "Library size OK"    = length(low_lib) == 0
)
for (nm in names(checks))
  cat(sprintf("  [%s] %s\n", ifelse(checks[nm], "PASS", "FAIL"), nm))
cat("\n")
