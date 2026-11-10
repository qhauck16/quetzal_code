
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
# args: [full_leaf_sig, introns, mane, dpsi, fc, output]
full_leaf_sig_path <- args[1]
introns_path       <- args[2]
mane_path          <- args[3]
dpsi_thresh        <- as.numeric(args[4])
fc_thresh          <- as.numeric(args[5])
out_path           <- args[6]

full_leaf_sig <- read_tsv(full_leaf_sig_path, show_col_types = FALSE)

introns <- fread(introns_path, data.table = FALSE)
mane_tx <- fread(mane_path,   data.table = FALSE)
mane_introns <- introns %>%
  filter(transcript_id %in% mane_tx$Ensembl_nuc) %>%
  dplyr::select(start, end, strand, gene_name, seqnames) %>%
  unique()

hits <- full_leaf_sig %>% filter(deltapsi > dpsi_thresh)
if (fc_thresh > 0) {
  hits <- hits %>% filter(test / control > fc_thresh)
}

make_id <- function(df) paste0(df$gene_name, ":", df$start, "-", df$end)

top <- hits %>%
  mutate(start = as.numeric(start),
         end   = as.numeric(end)) %>%
  filter(gene_name %in% mane_introns$gene_name)

canonical <- top %>%
  inner_join(mane_introns %>% dplyr::select(gene_name, strand, start, end) %>%
               rename(mane_start = start, mane_end = end),
             by = c("gene_name", "strand", "start" = "mane_start", "end" = "mane_end")) %>%
  mutate(type = "canonical",
         mane_intron_start = NA_integer_,
         mane_intron_end   = NA_integer_)
canonical_ids <- make_id(canonical)

skipped <- top %>%
  inner_join(mane_introns %>% dplyr::select(gene_name, strand, start, end) %>%
               rename(up_start = start, up_end = end),
             by = c("gene_name", "strand", "start" = "up_start")) %>%
  inner_join(mane_introns %>% dplyr::select(gene_name, strand, start, end) %>%
               rename(down_start = start, down_end = end),
             by = c("gene_name", "strand", "end" = "down_end")) %>%
  filter(up_end < down_start) %>%
  mutate(type = "skipped exon")
skipped_ids <- make_id(skipped)

assigned_ids <- c(canonical_ids, skipped_ids)

start_matches <- top %>%
  filter(!make_id(.) %in% assigned_ids) %>%
  inner_join(mane_introns %>% dplyr::select(gene_name, strand, start, end) %>%
               rename(mane_start = start, mane_end = end),
             by = c("gene_name", "strand", "start" = "mane_start")) %>%
  mutate(type = if_else(strand == "+", "alt 3'ss", "alt 5'ss"),
         mane_intron_start = start,
         mane_intron_end   = mane_end)

end_matches <- top %>%
  filter(!make_id(.) %in% assigned_ids) %>%
  inner_join(mane_introns %>% dplyr::select(gene_name, strand, start, end) %>%
               rename(mane_start = start, mane_end = end),
             by = c("gene_name", "strand", "end" = "mane_end")) %>%
  mutate(type = if_else(strand == "+", "alt 5'ss", "alt 3'ss"),
         mane_intron_start = mane_start,
         mane_intron_end   = end)

annotated <- bind_rows(canonical, skipped, start_matches, end_matches) %>%
  mutate(junc_id = make_id(.)) %>%
  dplyr::select(junc_id, type, mane_intron_start, mane_intron_end)

df_final <- top %>%
  mutate(junc_id = make_id(.)) %>%
  left_join(annotated, by = "junc_id") %>%
  mutate(type = if_else(is.na(type), "other", type)) %>%
  group_by(gene_name, start, end) %>%
  filter(n() < 2) %>%
  ungroup()

crypt_dist <- df_final %>%
  filter(type == "alt 3'ss") %>%
  mutate(crypt_dist = if_else(strand == "+",
                              mane_intron_end - end,
                              start - mane_intron_start)) %>%
  pull(crypt_dist)

df_final$distance_to_cryptic <- NA_integer_
df_final$distance_to_cryptic[df_final$type == "alt 3'ss"] <- crypt_dist

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
write_tsv(df_final, out_path)
