# =============================================================================
#  05_deseq2_analysis.R
#
#  Runs DESeq2 differential expression analysis for two comparison categories:
#    Category 1 — Tumor vs Normal
#    Category 2 — Mutated vs Wildtype  (requires MAF + UUID-barcode map)
#
#  For each gene group (E1, E2, E3, DUB, Combined), produces:
#    *_all_genes.csv          — full DESeq2 results table
#    *_significant_DEGs.csv   — filtered to padj < PADJ_CUTOFF, |LFC| > LOG2FC_CUTOFF
#    *_volcano.pdf            — volcano plot with top-15 gene labels
#    *_PCA.pdf                — PCA of VST-normalised counts
#    *_heatmap.tiff           — ComplexHeatmap of top N DEGs (Z-score scaled)
#    session_info.txt         — R session info for reproducibility
#
#  All thresholds and paths are set in config.R.
#
#  Input  : pan_cancer_counts_T.tsv
#            gene_annotation.tsv
#            sample_metadata_clean.tsv
#            uuid_to_barcode_map.tsv
#            gene_list.xlsx
#            combined_pan_cancer.maf
#  Output : OUTPUT_DIR/{Tumor_vs_Normal,Mutated_vs_Wildtype}/{E1,E2,E3,DUB,Combined}/
#
#  Run    : source("scripts/05_deseq2_analysis.R")
# =============================================================================

source("config.R")

# ── 0. INSTALL / LOAD DEPENDENCIES ───────────────────────────────────────────
required_cran <- c("ggplot2", "RColorBrewer", "dplyr", "readr",
                   "readxl", "stringr", "tibble", "ggrepel",
                   "data.table", "circlize")
required_bioc <- c("DESeq2", "ComplexHeatmap")

missing_cran <- required_cran[!required_cran %in% installed.packages()[, "Package"]]
if (length(missing_cran) > 0) install.packages(missing_cran, quiet = TRUE)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
missing_bioc <- required_bioc[!required_bioc %in% installed.packages()[, "Package"]]
if (length(missing_bioc) > 0) BiocManager::install(missing_bioc, ask = FALSE)

suppressPackageStartupMessages({
  library(DESeq2); library(ComplexHeatmap); library(circlize)
  library(ggplot2); library(ggrepel); library(RColorBrewer)
  library(dplyr); library(readr); library(readxl)
  library(stringr); library(tibble); library(data.table)
})
ht_opt$message <- FALSE

# ── 1. CREATE OUTPUT FOLDER STRUCTURE ────────────────────────────────────────
cat("Creating output folders...\n")
for (cat_folder in c("Tumor_vs_Normal", "Mutated_vs_Wildtype")) {
  for (grp in c(GENE_CATEGORIES, "Combined")) {
    dir.create(file.path(OUTPUT_DIR, cat_folder, grp),
               recursive = TRUE, showWarnings = FALSE)
  }
}

# ── 2. LOAD GENE LISTS ───────────────────────────────────────────────────────
cat("Loading gene lists from Excel...\n")
gene_sheet <- read_excel(GENE_LIST_FILE, sheet = GENE_SHEET_NAME)

missing_cols <- setdiff(GENE_CATEGORIES, colnames(gene_sheet))
if (length(missing_cols) > 0)
  stop("Missing columns in gene list: ", paste(missing_cols, collapse = ", "))

gene_groups <- lapply(setNames(GENE_CATEGORIES, GENE_CATEGORIES), function(col) {
  g <- gene_sheet[[col]]
  unique(trimws(as.character(g[!is.na(g) & g != ""])))
})
for (g in names(gene_groups))
  cat(sprintf("  %-4s : %d genes\n", g, length(gene_groups[[g]])))
all_ubl_genes <- unique(unlist(gene_groups))

# ── 3. LOAD GENE ANNOTATION ──────────────────────────────────────────────────
cat("Loading gene annotation...\n")
gene_annot        <- read_tsv(GENE_ANNOT, show_col_types = FALSE)
ensembl_to_symbol <- setNames(gene_annot$gene_name, gene_annot$gene_id)
cat(sprintf("  %d genes in annotation\n", nrow(gene_annot)))

# ── 4. LOAD SAMPLE METADATA ──────────────────────────────────────────────────
cat("Loading sample metadata...\n")
metadata <- read_tsv(METADATA_FILE, show_col_types = FALSE) %>%
  mutate(Category1 = case_when(
    sample_group == "Normal"                                         ~ "Normal",
    grepl("Tumor|tumor|Cancer|Metastatic|Recurrent", sample_group)  ~ "Tumor",
    TRUE                                                             ~ NA_character_
  ))
cat(sprintf("  Total samples: %d\n", nrow(metadata)))
print(table(metadata$Category1, useNA = "ifany"))

# ── 5. LOAD UUID -> TCGA BARCODE MAP ────────────────────────────────────────
cat("Loading UUID to TCGA barcode map...\n")
uuid_map <- read_tsv(UUID_MAP_FILE, show_col_types = FALSE)
cat(sprintf("  %d entries\n", nrow(uuid_map)))

# ── 6. LOAD MAF & BUILD MUTATION LABELS ──────────────────────────────────────
cat("Loading MAF file...\n")
maf_raw         <- readLines(MAF_FILE, n = 200, warn = FALSE)
n_comment_lines <- sum(startsWith(maf_raw, "#"))
maf_data <- fread(MAF_FILE,
                  select       = c("Hugo_Symbol", "Tumor_Sample_Barcode"),
                  skip         = n_comment_lines,
                  showProgress = FALSE)
setDF(maf_data)
cat(sprintf("  MAF rows: %d\n", nrow(maf_data)))

maf_ubl <- maf_data %>% filter(Hugo_Symbol %in% all_ubl_genes)
maf_ubl$barcode_short  <- substr(maf_ubl$Tumor_Sample_Barcode, 1, 15)
uuid_map$barcode_short <- substr(uuid_map$tcga_barcode,        1, 15)
mutated_barcodes_short <- unique(maf_ubl$barcode_short)

# Auto-detect which ID column in uuid_map aligns with metadata$sample_id
n_match_file    <- sum(metadata$sample_id %in% uuid_map$file_id[uuid_map$barcode_short %in% mutated_barcodes_short])
n_match_uuid    <- sum(metadata$sample_id %in% uuid_map$sample_uuid[uuid_map$barcode_short %in% mutated_barcodes_short])
n_match_barcode <- sum(metadata$sample_id %in% uuid_map$tcga_barcode[uuid_map$barcode_short %in% mutated_barcodes_short])

if (n_match_file >= n_match_uuid && n_match_file >= n_match_barcode) {
  per_group_fn    <- function(bs) uuid_map$file_id[uuid_map$barcode_short %in% bs]
  mutated_ids_all <- per_group_fn(mutated_barcodes_short)
  cat("  Using file_id to identify mutated samples\n")
} else if (n_match_uuid >= n_match_barcode) {
  per_group_fn    <- function(bs) uuid_map$sample_uuid[uuid_map$barcode_short %in% bs]
  mutated_ids_all <- per_group_fn(mutated_barcodes_short)
  cat("  Using sample_uuid to identify mutated samples\n")
} else {
  per_group_fn    <- function(bs) uuid_map$tcga_barcode[uuid_map$barcode_short %in% bs]
  mutated_ids_all <- per_group_fn(mutated_barcodes_short)
  cat("  Using tcga_barcode to identify mutated samples\n")
}

mutated_file_ids_by_group <- lapply(gene_groups, function(genes) {
  bs <- unique(maf_ubl$barcode_short[maf_ubl$Hugo_Symbol %in% genes])
  per_group_fn(bs)
})

metadata <- metadata %>%
  mutate(Category2 = case_when(
    Category1 != "Tumor"           ~ NA_character_,
    sample_id %in% mutated_ids_all ~ "Mutated",
    Category1 == "Tumor"           ~ "Wildtype",
    TRUE                           ~ NA_character_
  ))
cat("  Category2 (tumors only):\n")
print(table(metadata$Category2, useNA = "ifany"))

# ── 7. BUILD GENE ROW INDEX (read first column only) ─────────────────────────
cat("Scanning count matrix gene index...\n")
gene_id_index       <- fread(COUNTS_FILE, select = 1L, showProgress = FALSE)[[1]]
gene_id_index_clean <- sub("\\..*$", "", gene_id_index)

ubl_ensembl_all <- names(ensembl_to_symbol)[ensembl_to_symbol %in% all_ubl_genes]
ubl_rows        <- which(gene_id_index_clean %in% ubl_ensembl_all)
if (length(ubl_rows) == 0)
  stop("No UBL genes matched. Check gene symbols in Excel vs gene_annotation.tsv.")

gene_rows_by_group <- lapply(gene_groups, function(genes) {
  ens <- names(ensembl_to_symbol)[ensembl_to_symbol %in% genes]
  which(gene_id_index_clean %in% ens)
})
for (g in names(gene_rows_by_group))
  cat(sprintf("  %-4s : %d rows\n", g, length(gene_rows_by_group[[g]])))

# ── 8. RAM-EFFICIENT COUNT LOADER ─────────────────────────────────────────────
load_count_subset <- function(row_indices, sample_ids) {
  all_cols     <- colnames(fread(COUNTS_FILE, nrows = 0))
  cols_to_read <- c("gene_id", intersect(sample_ids, all_cols))
  if (length(cols_to_read) < 2) {
    cat("    *** No matching sample columns in count matrix\n"); return(NULL)
  }
  dt  <- fread(COUNTS_FILE, select = cols_to_read, showProgress = FALSE)
  dt  <- dt[row_indices, ]
  mat <- as.matrix(dt[, -1, with = FALSE])
  rownames(mat) <- sub("\\..*$", "", dt[[1]])
  mat[is.na(mat)] <- 0L
  storage.mode(mat) <- "integer"
  mat
}

# ── 9. CORE ANALYSIS FUNCTIONS ────────────────────────────────────────────────

run_deseq2 <- function(count_mat, ref_level, test_level,
                       sample_ids_ref, sample_ids_test) {
  valid_ref  <- intersect(sample_ids_ref,  colnames(count_mat))
  valid_test <- intersect(sample_ids_test, colnames(count_mat))
  cat(sprintf("    %s: %d | %s: %d\n",
              ref_level, length(valid_ref), test_level, length(valid_test)))

  if (length(valid_ref) < MIN_SAMPLES || length(valid_test) < MIN_SAMPLES) {
    cat(sprintf("    *** < %d samples in a group — skipping\n", MIN_SAMPLES))
    return(NULL)
  }

  sub_mat <- count_mat[, c(valid_ref, valid_test), drop = FALSE]
  coldata <- data.frame(
    Condition = factor(c(rep(ref_level, length(valid_ref)),
                         rep(test_level, length(valid_test))),
                       levels = c(ref_level, test_level)),
    row.names = c(valid_ref, valid_test)
  )

  sub_mat <- sub_mat[rowSums(sub_mat) >= MIN_COUNT_SUM, , drop = FALSE]
  cat(sprintf("    Genes after low-count filter: %d\n", nrow(sub_mat)))
  if (nrow(sub_mat) < 5) { cat("    *** Too few genes\n"); return(NULL) }

  dds <- tryCatch(DESeqDataSetFromMatrix(sub_mat, coldata, ~ Condition),
                  error = function(e) { cat("    *** DESeqDataSetFromMatrix:", conditionMessage(e), "\n"); NULL })
  if (is.null(dds)) return(NULL)

  dds <- tryCatch(DESeq(dds, quiet = TRUE),
                  error = function(e) { cat("    *** DESeq():", conditionMessage(e), "\n"); NULL })
  if (is.null(dds)) return(NULL)

  res    <- results(dds, contrast = c("Condition", test_level, ref_level))
  res_df <- as.data.frame(res)
  res_df$Ensembl  <- rownames(res_df)
  res_df$GeneName <- ensembl_to_symbol[res_df$Ensembl]
  res_df$GeneName[is.na(res_df$GeneName)] <- res_df$Ensembl[is.na(res_df$GeneName)]
  res_df$Significant <- case_when(
    !is.na(res_df$padj) & res_df$padj < PADJ_CUTOFF & res_df$log2FoldChange >  LOG2FC_CUTOFF ~ "Upregulated",
    !is.na(res_df$padj) & res_df$padj < PADJ_CUTOFF & res_df$log2FoldChange < -LOG2FC_CUTOFF ~ "Downregulated",
    TRUE ~ "Not Significant"
  )
  cat(sprintf("    Up: %d | Down: %d\n",
              sum(res_df$Significant == "Upregulated"),
              sum(res_df$Significant == "Downregulated")))
  list(dds = dds, res_df = res_df)
}

make_volcano <- function(res_df, title, out_path) {
  pd      <- res_df[!is.na(res_df$padj), ]
  pd$y    <- pmin(-log10(pd$padj), 100)
  top_lab <- pd %>% filter(Significant != "Not Significant") %>% arrange(padj) %>% head(15)

  p <- ggplot(pd, aes(log2FoldChange, y, color = Significant)) +
    geom_point(alpha = 0.45, size = 1.1) +
    scale_color_manual(values = c(Upregulated       = "#D73027",
                                  Downregulated     = "#4575B4",
                                  `Not Significant` = "grey72")) +
    geom_vline(xintercept = c(-LOG2FC_CUTOFF, LOG2FC_CUTOFF),
               linetype = "dashed", linewidth = 0.4) +
    geom_hline(yintercept = -log10(PADJ_CUTOFF),
               linetype = "dashed", linewidth = 0.4) +
    geom_text_repel(data = top_lab, aes(label = GeneName),
                    size = 2.6, max.overlaps = 20, box.padding = 0.3) +
    labs(title = title, x = "Log2 Fold Change",
         y = "-Log10 Adjusted P-value", color = NULL) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  pdf(out_path, width = 10, height = 8); print(p); dev.off()
  cat(sprintf("    Saved: %s\n", basename(out_path)))
}

make_pca <- function(dds, title, out_path) {
  tryCatch({
    vsd <- tryCatch(vst(dds, blind = FALSE),
                    error = function(e) varianceStabilizingTransformation(dds, blind = FALSE))
    p <- plotPCA(vsd, intgroup = "Condition") +
      theme_bw(base_size = 11) + ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    pdf(out_path, width = 8, height = 6); print(p); dev.off()
    cat(sprintf("    Saved: %s\n", basename(out_path)))
  }, error = function(e) cat("    *** PCA failed:", conditionMessage(e), "\n"))
}

make_heatmap <- function(dds, res_df, title, out_path,
                         n = TOP_HEATMAP_N,
                         max_per_group = MAX_SAMPLES_HEAT) {
  tryCatch({
    top <- res_df %>% filter(Significant != "Not Significant") %>%
      arrange(padj) %>% head(n)
    if (nrow(top) < 2) { cat("    *** < 2 significant genes — skipping heatmap\n"); return(invisible(NULL)) }

    norm <- counts(dds, normalized = TRUE)
    ids  <- intersect(top$Ensembl, rownames(norm))
    if (length(ids) < 2) { cat("    *** < 2 genes in normalised counts — skipping\n"); return(invisible(NULL)) }

    mat      <- norm[ids, , drop = FALSE]
    row_vars <- apply(mat, 1, var, na.rm = TRUE)
    mat      <- mat[row_vars > 0, , drop = FALSE]
    if (nrow(mat) < 2) { cat("    *** All genes zero-variance — skipping\n"); return(invisible(NULL)) }

    sample_condition <- setNames(as.character(dds$Condition), colnames(dds))
    groups           <- levels(dds$Condition)

    keep_cols <- unlist(lapply(groups, function(grp) {
      grp_cols <- intersect(names(sample_condition)[sample_condition == grp], colnames(mat))
      if (length(grp_cols) > max_per_group) { set.seed(42); sort(sample(grp_cols, max_per_group)) }
      else sort(grp_cols)
    }))
    mat <- mat[, keep_cols, drop = FALSE]

    scaled <- t(scale(t(mat)))
    scaled[!is.finite(scaled)] <- 0
    scaled <- pmax(pmin(scaled, 3), -3)

    rlabs <- top$GeneName[match(rownames(scaled), top$Ensembl)]
    rlabs[is.na(rlabs)] <- rownames(scaled)[is.na(rlabs)]
    rownames(scaled) <- make.unique(rlabs)

    cond_vec   <- factor(sample_condition[keep_cols], levels = groups)
    all_colors <- c(Normal = "#2196F3", Tumor = "#F44336",
                    Wildtype = "#4CAF50", Mutated = "#FF9800")
    cond_colors <- all_colors[names(all_colors) %in% groups]

    col_fun <- colorRamp2(c(-3, -1.5, 0, 1.5, 3),
                          c("#4575B4", "#91BFDB", "white", "#FC8D59", "#D73027"))
    top_ann <- HeatmapAnnotation(
      Condition = cond_vec,
      col       = list(Condition = cond_colors),
      annotation_name_gp   = gpar(fontsize = 11, fontface = "bold"),
      annotation_name_side = "left",
      simple_anno_size     = unit(6, "mm"),
      show_legend          = TRUE
    )

    fig_height <- max(8,  (nrow(scaled) * 8 + 40) / 25.4)
    fig_width  <- max(16, (ncol(scaled) * 1.5 + 60) / 25.4)

    ht <- Heatmap(
      scaled,
      name                 = "Z-score",
      col                  = col_fun,
      top_annotation       = top_ann,
      cluster_rows         = TRUE,
      cluster_columns      = FALSE,
      show_column_names    = FALSE,
      show_row_names       = TRUE,
      row_names_side       = "right",
      row_names_gp         = gpar(fontsize = 11, fontface = "italic"),
      row_dend_side        = "left",
      row_dend_width       = unit(15, "mm"),
      column_split         = cond_vec,
      column_gap           = unit(4, "mm"),
      column_title         = NULL,
      rect_gp              = gpar(col = "white", lwd = 0.5),
      heatmap_legend_param = list(
        title         = " ",
        labels_gp     = gpar(fontsize = 10),
        legend_height = unit(50, "mm"),
        at            = c(-3, -2, -1, 0, 1, 2, 3)
      ),
      use_raster     = TRUE,
      raster_quality = 8
    )

    tiff(out_path,
         width = as.integer(fig_width * 200), height = as.integer(fig_height * 200),
         res = 200, compression = "lzw")
    draw(ht,
         column_title           = title,
         column_title_gp        = gpar(fontsize = 14, fontface = "bold"),
         heatmap_legend_side    = "right",
         annotation_legend_side = "right",
         padding                = unit(c(12, 20, 12, 12), "mm"))
    dev.off()
    cat(sprintf("    Saved: %s\n", basename(out_path)))
  }, error = function(e) cat("    *** Heatmap failed:", conditionMessage(e), "\n"))
}

save_csvs <- function(res_df, out_dir, prefix) {
  sorted <- res_df[order(res_df$padj, na.last = TRUE), ]
  write.csv(sorted, file.path(out_dir, paste0(prefix, "_all_genes.csv")),        row.names = FALSE)
  degs <- sorted[!is.na(sorted$Significant) & sorted$Significant != "Not Significant", ]
  write.csv(degs,   file.path(out_dir, paste0(prefix, "_significant_DEGs.csv")), row.names = FALSE)
  cat(sprintf("    CSVs saved — %d significant DEGs\n", nrow(degs)))
  invisible(degs)
}

run_group_analysis <- function(group_name, row_indices,
                               ref_ids, test_ids,
                               ref_level, test_level,
                               category_folder) {
  cat(sprintf("\n-- %s | %s vs %s --\n", group_name, test_level, ref_level))
  count_mat <- load_count_subset(row_indices, unique(c(ref_ids, test_ids)))
  if (is.null(count_mat)) return(invisible(NULL))

  out <- run_deseq2(count_mat, ref_level, test_level, ref_ids, test_ids)
  if (is.null(out)) return(invisible(NULL))

  out_dir <- file.path(OUTPUT_DIR, category_folder, group_name)
  pfx     <- sprintf("%s_%s_vs_%s", group_name, test_level, ref_level)

  save_csvs(out$res_df, out_dir, pfx)
  make_volcano(out$res_df, sprintf("Volcano: %s | %s vs %s", group_name, test_level, ref_level),
               file.path(out_dir, paste0(pfx, "_volcano.pdf")))
  make_pca(out$dds, sprintf("PCA: %s | %s vs %s", group_name, test_level, ref_level),
           file.path(out_dir, paste0(pfx, "_PCA.pdf")))
  make_heatmap(out$dds, out$res_df,
               sprintf("Top%d DEGs: %s | %s vs %s", TOP_HEATMAP_N, group_name, test_level, ref_level),
               file.path(out_dir, paste0(pfx, "_heatmap.tiff")))
  invisible(out$res_df)
}

# ── 10. CATEGORY 1 — TUMOR vs NORMAL ─────────────────────────────────────────
cat("\n\n╔══════════════════════════════════════╗\n")
cat("║   CATEGORY 1 : TUMOR vs NORMAL      ║\n")
cat("╚══════════════════════════════════════╝\n")

normal_ids <- metadata$sample_id[!is.na(metadata$Category1) & metadata$Category1 == "Normal"]
tumor_ids  <- metadata$sample_id[!is.na(metadata$Category1) & metadata$Category1 == "Tumor"]
cat(sprintf("  Normal: %d | Tumor: %d\n", length(normal_ids), length(tumor_ids)))

cat1_res_list <- list()
for (grp in names(gene_groups)) {
  cat1_res_list[[grp]] <- run_group_analysis(
    group_name = grp, row_indices = gene_rows_by_group[[grp]],
    ref_ids = normal_ids, test_ids = tumor_ids,
    ref_level = "Normal", test_level = "Tumor",
    category_folder = "Tumor_vs_Normal"
  )
}

cat("\n-- Combined Top-50 (Tumor vs Normal) --\n")
all_degs_cat1 <- bind_rows(cat1_res_list) %>%
  filter(Significant != "Not Significant") %>%
  arrange(padj) %>% distinct(Ensembl, .keep_all = TRUE) %>%
  head(TOP_HEATMAP_N)
if (nrow(all_degs_cat1) >= 2) {
  run_group_analysis("Combined", ubl_rows, normal_ids, tumor_ids,
                     "Normal", "Tumor", "Tumor_vs_Normal")
}

# ── 11. CATEGORY 2 — MUTATED vs WILDTYPE ─────────────────────────────────────
cat("\n\n╔══════════════════════════════════════╗\n")
cat("║   CATEGORY 2 : MUTATED vs WILDTYPE  ║\n")
cat("╚══════════════════════════════════════╝\n")

wildtype_ids <- metadata$sample_id[!is.na(metadata$Category2) & metadata$Category2 == "Wildtype"]
mutated_ids  <- metadata$sample_id[!is.na(metadata$Category2) & metadata$Category2 == "Mutated"]
cat(sprintf("  Wildtype: %d | Mutated: %d\n", length(wildtype_ids), length(mutated_ids)))

if (length(wildtype_ids) < MIN_SAMPLES || length(mutated_ids) < MIN_SAMPLES) {
  cat("\n  *** Not enough samples for Category 2 (need >= MIN_SAMPLES per group)\n")
  cat("  *** Verify MAF barcodes match tcga_barcode in uuid_to_barcode_map.tsv\n")
} else {
  cat2_res_list <- list()
  for (grp in names(gene_groups)) {
    grp_mut_ids <- intersect(mutated_file_ids_by_group[[grp]], tumor_ids)
    cat2_res_list[[grp]] <- run_group_analysis(
      group_name = grp, row_indices = gene_rows_by_group[[grp]],
      ref_ids = wildtype_ids, test_ids = grp_mut_ids,
      ref_level = "Wildtype", test_level = "Mutated",
      category_folder = "Mutated_vs_Wildtype"
    )
  }

  cat("\n-- Combined Top-50 (Mutated vs Wildtype) --\n")
  all_degs_cat2 <- bind_rows(cat2_res_list) %>%
    filter(Significant != "Not Significant") %>%
    arrange(padj) %>% distinct(Ensembl, .keep_all = TRUE) %>%
    head(TOP_HEATMAP_N)
  if (nrow(all_degs_cat2) >= 2) {
    run_group_analysis("Combined", ubl_rows, wildtype_ids, mutated_ids,
                       "Wildtype", "Mutated", "Mutated_vs_Wildtype")
  }
}

# ── 12. SAVE SESSION INFO ─────────────────────────────────────────────────────
sink(file.path(OUTPUT_DIR, "session_info.txt"))
sessionInfo()
sink()

cat("\n\n╔══════════════════════════════════════╗\n")
cat("║          ANALYSIS COMPLETE           ║\n")
cat("╚══════════════════════════════════════╝\n")
cat("Output saved to:", OUTPUT_DIR, "\n")
cat("Each subfolder contains:\n")
cat("  *_all_genes.csv  |  *_significant_DEGs.csv\n")
cat("  *_volcano.pdf    |  *_PCA.pdf  |  *_heatmap.tiff\n")
cat("  session_info.txt\n\n")
