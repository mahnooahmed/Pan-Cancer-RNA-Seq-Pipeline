# =============================================================================
#  01_fetch_gdc_metadata.R
#
#  Queries the GDC API to retrieve sample type and project ID for every
#  file UUID in the GDC manifest, then saves the result as
#  gdc_sample_types.tsv for use by downstream scripts.
#
#  Input  : GDC manifest (.txt)
#  Output : gdc_sample_types.tsv  (file_id | project_id | sample_type | sample_uuid)
#
#  Run    : source("scripts/01_fetch_gdc_metadata.R")
# =============================================================================

source("config.R")

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
  if (!requireNamespace("httr",       quietly = TRUE)) install.packages("httr")
  if (!requireNamespace("jsonlite",   quietly = TRUE)) install.packages("jsonlite")
  library(data.table)
  library(httr)
  library(jsonlite)
})

# ── LOAD MANIFEST ─────────────────────────────────────────────────────────────
cat("Reading manifest...\n")
manifest  <- fread(MANIFEST_PATH, sep = "\t", header = TRUE)
cat(sprintf("  Columns : %s\n", paste(colnames(manifest), collapse = ", ")))
cat(sprintf("  Files   : %d\n\n", nrow(manifest)))

all_ids   <- manifest$id   # 'id' column = GDC file UUID
n_batches <- ceiling(length(all_ids) / 100)
results   <- vector("list", n_batches)

# ── BATCH QUERY GDC API ───────────────────────────────────────────────────────
cat(sprintf("Fetching sample types for %d files in %d batches...\n",
            length(all_ids), n_batches))

for (i in seq(1, length(all_ids), by = 100)) {
  b     <- ceiling(i / 100)
  batch <- all_ids[i:min(i + 99, length(all_ids))]

  body <- list(
    filters = list(
      op      = "in",
      content = list(field = "file_id", value = as.list(batch))
    ),
    fields = "file_id,cases.samples.sample_type,cases.project.project_id",
    size   = length(batch)
  )

  resp <- POST(
    "https://api.gdc.cancer.gov/files",
    body        = toJSON(body, auto_unbox = TRUE),
    add_headers("Content-Type" = "application/json", "Accept" = "application/json"),
    timeout(60)
  )

  if (status_code(resp) != 200) {
    cat(sprintf("  Batch %d — API error %d, skipping\n", b, status_code(resp)))
    next
  }

  hits <- fromJSON(
    content(resp, as = "text", encoding = "UTF-8"),
    simplifyVector = TRUE
  )$data$hits

  results[[b]] <- rbindlist(lapply(seq_len(nrow(hits)), function(j) {
    data.table(
      file_id     = hits$file_id[j],
      project_id  = tryCatch(hits$cases[[j]]$project$project_id[1],
                             error = function(e) NA_character_),
      sample_type = tryCatch(hits$cases[[j]]$samples[[1]]$sample_type[1],
                             error = function(e) NA_character_)
    )
  }), fill = TRUE)

  if (b %% 10 == 0 || b == n_batches)
    cat(sprintf("  Batch %d / %d done\n", b, n_batches))
  Sys.sleep(0.3)   # be polite to the API
}

# ── COMBINE & ANNOTATE ────────────────────────────────────────────────────────
sample_types <- rbindlist(results, fill = TRUE)

# Map file UUID -> sample UUID (the prefix of the filename in the manifest)
sample_types[, sample_uuid := sub("\\..*", "",
                                  manifest$filename[match(file_id, manifest$id)])]

cat(sprintf("\nTotal records parsed: %d\n", nrow(sample_types)))
cat("\nSample type breakdown:\n")
print(table(sample_types$sample_type))
cat("\nSamples per cancer type:\n")
print(sample_types[, .N, by = .(project_id, sample_type)][order(project_id)])

# ── SAVE ──────────────────────────────────────────────────────────────────────
fwrite(sample_types, SAMPLE_TYPES, sep = "\t")
cat(sprintf("\nSaved -> %s\n", SAMPLE_TYPES))
