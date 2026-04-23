# =============================================================================
#  04_build_tpm_subsets.R
#
#  Extracts TPM values for gene categories defined in the Excel gene list
#  (E1, E2, E3, DUB) and saves one TSV per category.  The full TPM matrix
#  is read in a single pass so each gene file is loaded only once regardless
#  of how many categories are requested.
#
#  Ensembl IDs in the matrix are matched to gene symbols via gene_annotation.tsv.
#  Output columns use HGNC gene symbols; a lookup table is also saved.
#
#  Input  : pan_cancer_tpm_T.tsv
#            gene_annotation.tsv
#            gene_list.xlsx  (sheet defined by GENE_SHEET_NAME in config.R)
#  Output : TPM_subsets/tpm_{E1,E2,E3,DUB}.tsv
#            TPM_subsets/gene_category_lookup.tsv
#
#  Run    : source("scripts/04_build_tpm_subsets.R")
# =============================================================================

source("config.R")

suppressPackageStartupMessages({
  for (pkg in c("data.table", "readxl", "stringr")) {
    if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  }
  library(data.table)
  library(readxl)
  library(stringr)
})

dir.create(TPM_SUBSET_DIR, showWarnings = FALSE, recursive = TRUE)

META_COLS <- c("sample_id", "cancer_type", "sample_group")

# ── LOAD GENE ANNOTATION ──────────────────────────────────────────────────────
cat("Loading gene annotation...\n")
ann <- fread(GENE_ANNOT, sep = "\t", header = TRUE)
ensembl_to_symbol <- setNames(ann$gene_name, ann$gene_id)

# ── READ TPM COLUMN NAMES (no data loaded yet) ────────────────────────────────
cat("Reading TPM matrix header...\n")
hdr          <- fread(TPM_FILE, nrows = 0L)
avail_genes  <- setdiff(names(hdr), META_COLS)
rm(hdr); gc()
cat(sprintf("  %d gene columns available\n\n", length(avail_genes)))

# ── LOAD GENE LISTS FROM EXCEL ────────────────────────────────────────────────
cat("Loading gene lists from Excel...\n")
xl <- read_excel(GENE_LIST_FILE, sheet = GENE_SHEET_NAME)

gene_lists <- lapply(GENE_CATEGORIES, function(col) {
  toupper(str_trim(na.omit(as.character(xl[[col]]))))
})
names(gene_lists) <- GENE_CATEGORIES

# ── RESOLVE SYMBOLS -> ENSEMBL IDs ───────────────────────────────────────────
cat("Resolving gene symbols to Ensembl IDs:\n")
ensembl_lists <- lapply(GENE_CATEGORIES, function(cat_name) {
  syms    <- gene_lists[[cat_name]]
  ensembl <- ann[toupper(gene_name) %in% syms, gene_id]
  found   <- intersect(ensembl, avail_genes)
  missing <- setdiff(syms, toupper(ann$gene_name))

  cat(sprintf("  %-6s : %d symbols -> %d Ensembl IDs -> %d in TPM matrix",
              cat_name, length(syms), length(ensembl), length(found)))
  if (length(missing) > 0)
    cat(sprintf("  [not in annotation: %s]", paste(missing, collapse = ", ")))
  cat("\n")
  found
})
names(ensembl_lists) <- GENE_CATEGORIES

# ── LOAD NEEDED GENES IN A SINGLE PASS ───────────────────────────────────────
all_ensembl <- unique(unlist(ensembl_lists))
cat(sprintf("\nUnique genes to load: %d\n", length(all_ensembl)))
cat("Loading TPM subset (single pass)...\n")
tpm_all <- fread(TPM_FILE, sep = "\t", header = TRUE,
                 select = c(META_COLS, all_ensembl))
cat(sprintf("  Loaded: %d samples × %d genes\n\n", nrow(tpm_all), length(all_ensembl)))

# ── SAVE ONE FILE PER CATEGORY ────────────────────────────────────────────────
cat("Saving subsets:\n")
for (cat_name in GENE_CATEGORIES) {
  cat_ens <- ensembl_lists[[cat_name]]
  if (length(cat_ens) == 0) {
    cat(sprintf("  %-6s : SKIPPED (no genes found)\n", cat_name))
    next
  }

  sub_dt    <- tpm_all[, c(META_COLS, cat_ens), with = FALSE]
  sym_names <- ensembl_to_symbol[cat_ens]
  sym_names[is.na(sym_names)] <- cat_ens[is.na(sym_names)]  # fallback to Ensembl ID
  setnames(sub_dt, old = cat_ens, new = sym_names)

  out_file <- file.path(TPM_SUBSET_DIR, sprintf("tpm_%s.tsv", cat_name))
  fwrite(sub_dt, out_file, sep = "\t")
  cat(sprintf("  %-6s : %d genes -> %s\n", cat_name, length(cat_ens), basename(out_file)))
}

# ── SAVE GENE LOOKUP TABLE ────────────────────────────────────────────────────
lookup_dt <- rbindlist(lapply(GENE_CATEGORIES, function(cat_name) {
  ens <- ensembl_lists[[cat_name]]
  data.table(
    category   = cat_name,
    ensembl_id = ens,
    gene_name  = ensembl_to_symbol[ens]
  )
}))
lookup_file <- file.path(TPM_SUBSET_DIR, "gene_category_lookup.tsv")
fwrite(lookup_dt, lookup_file, sep = "\t")
cat(sprintf("\nGene lookup table saved: %s\n", lookup_file))

cat("\n═══ DONE ═══\n")
cat(sprintf("Subsets saved to: %s\n\n", TPM_SUBSET_DIR))
for (f in list.files(TPM_SUBSET_DIR, pattern = "\\.tsv$"))
  cat(sprintf("  %-35s  %.1f MB\n", f,
              file.size(file.path(TPM_SUBSET_DIR, f)) / 1e6))
