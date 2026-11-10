
# Concatenate VEP-annotated VCFs from one (gene, source) pair into a flat TSV.

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
gene     <- args[1]                          # nkap / ylpm1
source   <- args[2]                          # wxs / rna
in_dir   <- args[3]                          # vep dir
out_path <- args[4]                          # output tsv

vcfs <- list.files(in_dir, pattern = "\\.vep\\.vcf$", full.names = TRUE)

parse_csq_header <- function(lines) {
  csq_line <- grep('^##INFO=<ID=CSQ,', lines, value = TRUE)
  if (length(csq_line) == 0) return(character(0))
  fmt <- sub('.*Format: ([^"]+)".*', '\\1', csq_line[1])
  strsplit(fmt, "\\|")[[1]]
}

parse_one <- function(p) {
  lines <- readLines(p)
  csq_fields <- parse_csq_header(lines)
  body <- lines[!grepl("^#", lines)]
  if (length(body) == 0) return(NULL)
  cols <- do.call(rbind, strsplit(body, "\t", fixed = TRUE))
  if (ncol(cols) < 8) return(NULL)
  df <- tibble(
    CHROM  = cols[, 1],
    POS    = as.integer(cols[, 2]),
    ID     = cols[, 3],
    REF    = cols[, 4],
    ALT    = cols[, 5],
    QUAL   = cols[, 6],
    FILTER = cols[, 7],
    INFO   = cols[, 8]
  )
  df$rail_id <- sub("\\.vep\\.vcf$", "", basename(p))
  csq <- str_match(df$INFO, "CSQ=([^;]+)")[, 2]
  csq[is.na(csq)] <- ""
  # Take only the first transcript annotation for a flat table
  first_csq <- sapply(strsplit(csq, ","), `[`, 1)
  if (length(csq_fields) > 0) {
    csq_split <- strsplit(first_csq, "\\|", fixed = FALSE)
    n <- length(csq_fields)
    csq_mat <- do.call(rbind, lapply(csq_split, function(x) {
      length(x) <- n; x
    }))
    colnames(csq_mat) <- csq_fields
    df <- bind_cols(df, as_tibble(csq_mat))
  }
  df
}

dat <- bind_rows(lapply(vcfs, parse_one))
if (nrow(dat) == 0) {
  message("No variants found in ", in_dir)
  fwrite(tibble(), out_path, sep = "\t")
  quit(status = 0)
}

dat <- dat %>%
  mutate(gene = gene, source = source) %>%
  dplyr::select(gene, source, rail_id, CHROM, POS, REF, ALT,
                FILTER, everything(), -INFO)

fwrite(dat, out_path, sep = "\t")
cat("Wrote", out_path, "(", nrow(dat), "variant rows from",
    length(vcfs), "VCFs)\n")
