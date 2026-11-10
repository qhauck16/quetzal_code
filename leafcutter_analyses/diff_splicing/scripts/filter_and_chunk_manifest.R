#!/usr/bin/env Rscript
# Filter a raw GDC manifest to RNA-Seq tumor BAMs (STAR 2-Pass Genome) using
# the companion JSON metadata, then split into chunks for parallel
# gdc-client downloads.
#
# Usage:
#   Rscript scripts/filter_and_chunk_manifest.R \
#       --manifest download/gdc_manifest.<date>.txt \
#       --json     download/metadata.cohort.<date>.json \
#       [--chunk_size 25]
#
# Outputs (in download/):
#   gdc_manifest.filtered.tsv         -- RNA-Seq tumor BAMs only
#   manifest_chunks/chunk_NN.tsv      -- chunk_size rows each (gdc-client format)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
get  <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  args[i + 1]
}
MANIFEST_PATH <- "download/gdc_manifest.2026-04-20.180256.txt"
JSON_PATH     <- "download/metadata.cohort.2026-04-20.json"
CHUNK         <- 50

if (is.null(MANIFEST_PATH) || is.null(JSON_PATH)) {
  stop("Both --manifest and --json are required.")
}

manifest <- read_tsv(MANIFEST_PATH, show_col_types = FALSE)
json     <- fromJSON(JSON_PATH, flatten = TRUE)

cat(sprintf("Loaded manifest: %d files\n", nrow(manifest)))
cat(sprintf("Loaded JSON metadata: %d files\n", nrow(json)))

# -- Filter JSON to RNA-Seq tumor BAMs ---------------------------------------
keep <- json %>%
  filter(data_format == "BAM",
         experimental_strategy == "RNA-Seq",
         analysis.workflow_type == "STAR 2-Pass Genome")

ent <- bind_rows(keep$associated_entities)
ent$file_id <- keep$file_id
ent <- ent %>%
  mutate(sample_type_code = substr(entity_submitter_id, 14, 15)) %>%
  filter(sample_type_code == "01")  # Primary Tumor

keep_ids <- unique(ent$file_id)
cat(sprintf("After RNA-Seq + STAR + tumor filter: %d files\n", length(keep_ids)))

# -- Filter manifest ---------------------------------------------------------
manifest_filt <- manifest %>% filter(id %in% keep_ids)
cat(sprintf("Manifest after intersect: %d files\n", nrow(manifest_filt)))

write_tsv(manifest_filt, "download/gdc_manifest.filtered.tsv")
cat("Wrote download/gdc_manifest.filtered.tsv\n")

# -- Chunk -------------------------------------------------------------------
out_dir <- "download/manifest_chunks"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

old <- list.files(out_dir, pattern = "^chunk_\\d+\\.tsv$", full.names = TRUE)
if (length(old)) file.remove(old)

n_chunks <- ceiling(nrow(manifest_filt) / CHUNK)
manifest_filt$chunk_idx <- rep(seq_len(n_chunks), each = CHUNK,
                               length.out = nrow(manifest_filt))

for (i in seq_len(n_chunks)) {
  rows <- manifest_filt %>% filter(chunk_idx == i) %>% dplyr::select(-chunk_idx)
  out  <- file.path(out_dir, sprintf("chunk_%02d.tsv", i))
  write_tsv(rows, out)
}

cat(sprintf("Split %d files into %d chunks (size %d) under %s\n",
            nrow(manifest_filt), n_chunks, CHUNK, out_dir))
cat("Run e.g.:\n")
cat(sprintf("  for f in %s/chunk_*.tsv; do gdc-client download -m $f -d download/raw/$(basename $f .tsv) & done; wait\n",
            out_dir))
