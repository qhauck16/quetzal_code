
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
# args: [psi_tsv, rpm_tsv, output_tsv]
psi_path <- args[1]
rpm_path <- args[2]
out_path <- args[3]

psi <- fread(psi_path, data.table = FALSE)
rpm <- fread(rpm_path, data.table = FALSE)

pool <- inner_join(psi, rpm, by = c("gene_name", "start", "end")) %>%
  group_by(gene_name, start, end) %>%
  slice_max(median_psi, n = 1, with_ties = FALSE) %>%
  ungroup()

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
write_tsv(pool, out_path)
