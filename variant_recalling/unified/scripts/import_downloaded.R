
# After running gdc-client download into <gene_dir>/download/<file_id>/,
# this script moves+renames each BAM/BAI into <gene_dir>/bams/<target_name>
# according to to_download.tsv.

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(fs)
})

args <- commandArgs(trailingOnly = TRUE)
gene_dir <- args[1]

dl_dir   <- file.path(gene_dir, "download")
bams_dir <- file.path(gene_dir, "bams")
dir.create(bams_dir, showWarnings = FALSE, recursive = TRUE)

to_dl <- fread(file.path(gene_dir, "to_download.tsv"),
               data.table = FALSE) %>% as_tibble()

n_moved <- 0L
n_missing <- 0L
for (i in seq_len(nrow(to_dl))) {
  fid    <- to_dl$file_id[i]
  target <- to_dl$target_name[i]
  src_dir <- file.path(dl_dir, fid)
  if (!dir.exists(src_dir)) {
    n_missing <- n_missing + 1L
    next
  }
  bams <- dir_ls(src_dir, regexp = "\\.bam$", type = "file")
  if (length(bams) != 1) {
    warning("Expected 1 BAM in ", src_dir, "; found ", length(bams))
    next
  }
  bais <- dir_ls(src_dir, regexp = "\\.bai$", type = "file")
  if (length(bais) != 1) {
    warning("Expected 1 BAI in ", src_dir, "; found ", length(bais))
  }
  tgt_bam <- file.path(bams_dir, target)
  tgt_bai <- sub("\\.bam$", ".bai", tgt_bam)
  if (!file.exists(tgt_bam)) file_move(bams, tgt_bam)
  if (length(bais) == 1 && !file.exists(tgt_bai)) file_move(bais, tgt_bai)
  n_moved <- n_moved + 1L
}

cat("Moved", n_moved, "BAMs into", bams_dir, "\n")
cat("Still missing (download dir not present):", n_missing, "\n")
