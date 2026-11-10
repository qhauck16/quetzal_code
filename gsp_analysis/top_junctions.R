# Top unannotated junction per top-loaded feature for GSP #10 (SF3B1), #13 (NKAP) and
# #30 (YLPM1), with splice-type assignment against MANE introns.
# Extracted from spec_flash_fact.Rmd (chunks "pre-processing" and "top 250 GSPs 10 13 30
# unannotated junctions + splice types"), code unchanged apart from paths.
# Writes fact_<k>_top_junctions.tsv into OUT (default: this directory).
#   Rscript top_junctions.R [OUT]        (run from gsp_analysis/)
suppressPackageStartupMessages({
  library(data.table); library(tidyverse); library(flashier); library(fastTopics)
})
first <- dplyr::first
args <- commandArgs(TRUE); OUT <- if (length(args)) args[1] else "."

out <- readRDS('../quetzal_runs/genome_wide/softimpute_flash_300_qc_filtered.RDS')
obj_dir <- '../quetzal_runs/genome_wide/'
gencode <- readRDS('../inputs/gencode_v46_granges.RDS')

ldf <- flashier::ldf(out, type = 'i')
norm_L <- ldf$L
norm_F <- t(diag(ldf$D) %*% t(ldf$F))   # absorb D into F

gene_chr_map <- as.data.frame(gencode) %>%
  filter(type == "exon") %>%
  distinct(gene_name, seqnames, strand) %>%
  group_by(gene_name) %>%
  summarise(chr    = as.character(first(seqnames)),
             strand = as.character(first(strand)),
             .groups = "drop")

FACTS  <- c(10L, 13L, 30L)
N_TOP  <- 250L

intron_df        <- fread('../leafcutter_analyses/diff_splicing/gencode_v46_introns.tsv', data.table = F)
mane_transcripts <- fread('../inputs/MANE.GRCh38.v1.4.summary.txt', data.table = F)
mane_introns <- intron_df %>%
  filter(transcript_id %in% mane_transcripts$Ensembl_nuc) %>%
  dplyr::select(start, end, strand, gene_name, seqnames) %>%
  unique()

snaptron_cols <- c('gene_name',
                    colnames(fread('../inputs/snaptron/tcga/chr1/snaptron_output/ABCA4_snaptron.tsv',
                                    nrows = 0L)))

for (fact in FACTS) {
message(">>> fact = ", fact)

top_features <- rownames(norm_F)[order(norm_F[, fact], decreasing = TRUE)][seq_len(N_TOP)]
parts        <- str_split_fixed(top_features, "\\.", 2)
genes_top    <- parts[, 1]
lfacts_top   <- parts[, 2]

top_junc_list <- vector("list", length(top_features))
for (j in seq_along(top_features)) {
  gene_name_j <- genes_top[j]
  pnmf_fact   <- lfacts_top[j]

  chr_match <- gene_chr_map$chr[gene_chr_map$gene_name == gene_name_j]
  if (!length(chr_match)) next
  chr_j <- chr_match[1]

  res_path      <- paste0(obj_dir, '/', chr_j, '/FastTopics_output/', gene_name_j, '/res.RDS')
  snaptron_path <- paste0('../inputs/snaptron/tcga/', chr_j, '/snaptron_output/', gene_name_j, '_snaptron.tsv')
  if (!file.exists(res_path) || !file.exists(snaptron_path)) next

  res        <- readRDS(res_path)
  gene_table <- fread(snaptron_path, data.table = F)
  res        <- poisson2multinom(res)
  if (!(pnmf_fact %in% colnames(res$F))) next

  unn_juncs <- paste0(gene_table$start[gene_table$annotated == 0], '-',
                       gene_table$end  [gene_table$annotated == 0])
  cand_idx  <- rownames(res$F) %in% unn_juncs
  if (!any(cand_idx)) next

  loads_cand <- res$F[cand_idx, pnmf_fact]
  top_junc   <- rownames(res$F)[cand_idx][which.max(loads_cand)]

  jc       <- as.numeric(str_split_fixed(top_junc, '-', 2))
  matching <- gene_table[gene_table$start == jc[1] & gene_table$end == jc[2], , drop = FALSE]
  if (!nrow(matching)) next

  to_keep <- as.data.frame(c(list(gene_name = gene_name_j), as.list(matching[1, ])),
                           stringsAsFactors = FALSE)
  names(to_keep) <- snaptron_cols
  top_junc_list[[j]] <- to_keep
}

top_junc_df <- bind_rows(top_junc_list)
top_junc_df$start <- as.numeric(top_junc_df$start)
top_junc_df$end   <- as.numeric(top_junc_df$end)

make_id <- function(df) paste0(df$gene_name, ":", df$start, "-", df$end)

canonical <- top_junc_df %>%
  inner_join(mane_introns %>% dplyr::select(gene_name, strand, start, end) %>%
               dplyr::rename(mane_start = start, mane_end = end),
             by = c("gene_name", "strand", "start" = "mane_start", "end" = "mane_end")) %>%
  mutate(type = "canonical", mane_intron_start = NA_integer_, mane_intron_end = NA_integer_)
canonical_ids <- make_id(canonical)

skipped_exon <- top_junc_df %>%
  inner_join(mane_introns %>% dplyr::select(gene_name, strand, start, end) %>%
               dplyr::rename(up_start = start, up_end = end),
             by = c("gene_name", "strand", "start" = "up_start")) %>%
  inner_join(mane_introns %>% dplyr::select(gene_name, strand, start, end) %>%
               dplyr::rename(down_start = start, down_end = end),
             by = c("gene_name", "strand", "end" = "down_end")) %>%
  filter(up_end < down_start) %>%
  mutate(type = "skipped exon")
skipped_ids  <- make_id(skipped_exon)
assigned_ids <- c(canonical_ids, skipped_ids)

start_matches <- top_junc_df %>%
  filter(!make_id(.) %in% assigned_ids) %>%
  inner_join(mane_introns %>% dplyr::select(gene_name, strand, start, end) %>%
               dplyr::rename(mane_start = start, mane_end = end),
             by = c("gene_name", "strand", "start" = "mane_start")) %>%
  mutate(type = if_else(strand == "+", "alt 3'ss", "alt 5'ss"),
         mane_intron_start = start, mane_intron_end = mane_end)

end_matches <- top_junc_df %>%
  filter(!make_id(.) %in% assigned_ids) %>%
  inner_join(mane_introns %>% dplyr::select(gene_name, strand, start, end) %>%
               dplyr::rename(mane_start = start, mane_end = end),
             by = c("gene_name", "strand", "end" = "mane_end")) %>%
  mutate(type = if_else(strand == "+", "alt 5'ss", "alt 3'ss"),
         mane_intron_start = mane_start, mane_intron_end = end)

annotated <- bind_rows(canonical, skipped_exon, start_matches, end_matches) %>%
  mutate(junc_id = make_id(.)) %>%
  dplyr::select(junc_id, type, mane_intron_start, mane_intron_end)

df_final <- top_junc_df %>%
  mutate(junc_id = make_id(.)) %>%
  left_join(annotated, by = "junc_id") %>%
  mutate(type = if_else(is.na(type), "other", type)) %>%
  group_by(gene_name, start, end) %>%
  filter(n() < 2) %>%
  ungroup()

crypt_dist <- df_final %>%
  filter(type == "alt 3'ss") %>%
  mutate(crypt_dist = if_else(strand == "+", mane_intron_end - end, start - mane_intron_start)) %>%
  pull(crypt_dist)
df_final$distance_to_cryptic <- NA_integer_
df_final$distance_to_cryptic[df_final$type == "alt 3'ss"] <- crypt_dist

write_tsv(df_final, file.path(OUT, paste0('fact_', fact, '_top_junctions.tsv')))
}
