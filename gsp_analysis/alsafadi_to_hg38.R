# Alsafadi et al. 2016 (Nat. Commun. 10.1038/ncomms10615) Supplementary Data 1: junctions
# differentially spliced in SF3B1-mutant uveal melanoma, lifted hg19 -> hg38.
# Extracted from spec_flash_fact.Rmd (chunk "alsafadi 2016 SF3B1 supplement -> hg38
# dataframe"), code unchanged apart from paths.
# Writes Alsafadi2016_supp1_hg38.tsv into OUT (default: this directory).
#   Rscript alsafadi_to_hg38.R [OUT]     (run from gsp_analysis/)
suppressPackageStartupMessages({
  library(tidyverse); library(readxl); library(rtracklayer); library(GenomicRanges)
})
args <- commandArgs(TRUE); OUT <- if (length(args)) args[1] else "."

ALSAFADI_XLSX <- "../inputs/Alsafadi2016_supp1.xlsx"
CHAIN_GZ      <- "../inputs/alphamissense/hg19ToHg38.over.chain.gz"
CHAIN_FILE    <- file.path(tempdir(), "hg19ToHg38.over.chain")
writeLines(readLines(gzfile(CHAIN_GZ)), CHAIN_FILE)
chain <- import.chain(CHAIN_FILE)

# Header is on row 2; row 1 is the supplementary-data title.
alsafadi_hg19 <- read_excel(ALSAFADI_XLSX, sheet = 1, skip = 1)
names(alsafadi_hg19) <- c("chromosome", "start_hg19", "end_hg19", "strand",
                         "gene_name", "norm_rc_mean", "log2fc",
                         "padj", "type_alsafadi")
alsafadi_hg19$start_hg19 <- as.integer(alsafadi_hg19$start_hg19)
alsafadi_hg19$end_hg19   <- as.integer(alsafadi_hg19$end_hg19)
alsafadi_hg19$row_id     <- seq_len(nrow(alsafadi_hg19))

gr_hg19 <- GRanges(
  seqnames = alsafadi_hg19$chromosome,
  ranges   = IRanges(start = alsafadi_hg19$start_hg19, end = alsafadi_hg19$end_hg19),
  strand   = alsafadi_hg19$strand,
  row_id   = alsafadi_hg19$row_id
)
lifted_grl <- liftOver(gr_hg19, chain)
gr_hg38 <- unlist(lifted_grl)
lifted_df <- data.frame(
  row_id          = gr_hg38$row_id,
  chromosome_hg38 = as.character(seqnames(gr_hg38)),
  start_hg38      = start(gr_hg38),
  end_hg38        = end(gr_hg38),
  strand_hg38     = as.character(strand(gr_hg38)),
  stringsAsFactors = FALSE
)
alsafadi_hg38 <- alsafadi_hg19 %>% left_join(lifted_df, by = "row_id")
write_tsv(alsafadi_hg38, file.path(OUT, "Alsafadi2016_supp1_hg38.tsv"))
