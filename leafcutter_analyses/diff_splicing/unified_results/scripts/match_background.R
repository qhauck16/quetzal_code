
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(MatchIt)
  library(GenomicRanges)
})

args <- commandArgs(trailingOnly = TRUE)
# args: [df_final, psi, rpm, pool, gencode_rds, output, splice_types]
# splice_types: comma-separated list of `type` values to keep from
# df_final BEFORE matching (e.g. "alt 3'ss" or "alt 3'ss,alt 5'ss").
# Empty / missing -> keep everything.
df_path       <- args[1]
psi_path      <- args[2]
rpm_path      <- args[3]
pool_path     <- args[4]
gencode_path  <- args[5]
out_path      <- args[6]
splice_types  <- if (length(args) >= 7) args[7] else ""

df_final <- read_tsv(df_path, show_col_types = FALSE)
pool     <- fread(pool_path, data.table = FALSE)
psi      <- fread(psi_path,  data.table = FALSE)
rpm      <- fread(rpm_path,  data.table = FALSE)
gencode  <- readRDS(gencode_path)

if (nrow(df_final) == 0) {
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  write_tsv(tibble(), out_path)
  quit(save = "no", status = 0)
}

# Filter test set to splice types of interest (per-gene).
if (nzchar(splice_types)) {
  keep_types <- trimws(strsplit(splice_types, ",")[[1]])
  df_final <- df_final %>% filter(type %in% keep_types)
  cat("match_background: keeping splice types ",
      paste(keep_types, collapse = ", "), " -> ",
      nrow(df_final), " test rows\n", sep = "")
  if (nrow(df_final) == 0) {
    dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
    write_tsv(tibble(), out_path)
    quit(save = "no", status = 0)
  }
}

critical <- df_final %>%
  left_join(psi, by = c("gene_name",
                        "mane_intron_start" = "start",
                        "mane_intron_end"   = "end")) %>%
  left_join(rpm, by = c("gene_name",
                        "mane_intron_start" = "start",
                        "mane_intron_end"   = "end")) %>%
  dplyr::rename(annot_median_rpm = median_rpm,
                annot_median_psi = median_psi) %>%
  group_by(gene_name, start, end) %>%
  slice_max(annot_median_psi, with_ties = FALSE, na_rm = TRUE) %>%
  ungroup()

pool$condition <- "background"

test_rows <- critical %>%
  filter(!is.na(annot_median_psi)) %>%
  dplyr::select(gene_name, mane_intron_start, mane_intron_end,
                annot_median_psi, annot_median_rpm) %>%
  dplyr::rename(median_rpm = annot_median_rpm,
                median_psi = annot_median_psi,
                start      = mane_intron_start,
                end        = mane_intron_end) %>%
  mutate(condition = "test")

if (nrow(test_rows) == 0) {
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  write_tsv(tibble(), out_path)
  quit(save = "no", status = 0)
}

both <- test_rows %>%
  bind_rows(pool) %>%
  group_by(across(-condition)) %>%
  filter(!(condition == "background" & any(condition == "test"))) %>%
  ungroup() %>%
  filter(!is.na(median_rpm))

both$condition <- factor(both$condition, levels = c("background", "test"))

m.out <- matchit(condition ~ median_psi + median_rpm,
                 data     = both,
                 method   = "nearest",
                 ratio    = 1,
                 distance = "mahalanobis")

matched <- match.data(m.out)

# Annotate each row with strand + chromosome from gencode
g_exons <- gencode[gencode$type == "exon"]
gene_tbl <- data.frame(
    gene_name  = g_exons$gene_name,
    strand     = as.character(strand(g_exons)),
    chromosome = as.character(seqnames(g_exons))
  ) %>%
  distinct(gene_name, .keep_all = TRUE)

matched <- matched %>%
  left_join(gene_tbl, by = "gene_name") %>%
  filter(!is.na(strand), !is.na(chromosome))

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
write_tsv(matched, out_path)
