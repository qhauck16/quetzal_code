#!/usr/bin/env Rscript
# One-shot: load the flashier fit from RDS, extract the sample x factor
# loading matrix, and write a small TSV that select_and_match_bam.R can read
# cheaply on every snakemake invocation.
#
# Rerun whenever the RDS changes (not per gene / per comparison).

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

RDS_PATH <- "../../quetzal_runs/genome_wide/softimpute_flash_300_qc_filtered.RDS"
OUT_PATH <- "flashier_loadings.tsv"

cat(sprintf("Loading %s ...\n", RDS_PATH))
out <- readRDS(RDS_PATH)

L <- out$L_pm
if (is.null(rownames(L))) stop("L_pm has no rownames — cannot map to rail_id")

loadings <- as.data.frame(L)
colnames(loadings) <- paste0("k", seq_len(ncol(L)))
loadings <- tibble::rownames_to_column(loadings, var = "rail_id") %>%
  mutate(rail_id = as.integer(rail_id))

write_tsv(loadings, OUT_PATH)
cat(sprintf("Wrote %s: %d samples x %d factors\n",
            OUT_PATH, nrow(loadings), ncol(L)))
