.libPaths('/scratch/midway2/qhauck/conda_env/rstudio-server/lib/R/library/')
# Top-250 driving unannotated junctions per GSP, with splice types. A verbatim
# port of spec_flash_fact.Rmd chunk "top 250 GSPs 10 13 30 unannotated
# junctions + splice types" (build_gsp_top + assign_types), parameterised by
# GSP. Writes gsp_top/fact_<N>_top_junctions.tsv.
# Usage: Rscript build_gsp_top.R 1 12 10      (10 = check vs the published set)
if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
suppressPackageStartupMessages({
  library(data.table); library(tidyverse); library(fastTopics); library(flashier)
})
FACTS <- as.integer(commandArgs(TRUE)); N_TOP <- 250L
OUT <- file.path(QROOT, "mutation_ds/gsp_top")
first <- dplyr::first; rename <- dplyr::rename

out <- readRDS(file.path(QROOT, "quetzal_runs/genome_wide/softimpute_flash_300_qc_filtered.RDS"))
ldf <- flashier::ldf(out, type = "i")
norm_F <- t(diag(ldf$D) %*% t(ldf$F)); rm(out); gc()
gencode <- readRDS(file.path(QROOT, "inputs/gencode_v46_granges.RDS"))
gene_chr_map <- as.data.frame(gencode) %>% filter(type == "exon") %>%
  distinct(gene_name, seqnames, strand) %>% group_by(gene_name) %>%
  summarise(chr = as.character(first(seqnames)), strand = as.character(first(strand)), .groups = "drop")
obj_dir <- file.path(QROOT, "quetzal_runs/genome_wide")
AG <- file.path(QROOT, "inputs/snaptron/tcga")
intron_df        <- fread(file.path(QROOT, "leafcutter_analyses/diff_splicing/gencode_v46_introns.tsv"), data.table = F)
mane_transcripts <- fread(file.path(QROOT, "inputs/MANE.GRCh38.v1.4.summary.txt"), data.table = F)
mane_introns <- intron_df %>% filter(transcript_id %in% mane_transcripts$Ensembl_nuc) %>%
  dplyr::select(start, end, strand, gene_name, seqnames) %>% unique()
snaptron_cols <- c("gene_name", colnames(fread(file.path(AG, "chr1/snaptron_output/ABCA4_snaptron.tsv"), nrows = 0L)))

for (fact in FACTS) {
  message(">>> fact = ", fact)
  top_features <- rownames(norm_F)[order(norm_F[, fact], decreasing = TRUE)][seq_len(N_TOP)]
  parts <- str_split_fixed(top_features, "\\.", 2)
  top_junc_list <- vector("list", length(top_features))
  for (j in seq_along(top_features)) {
    g <- parts[j, 1]; pf <- parts[j, 2]
    chr_match <- gene_chr_map$chr[gene_chr_map$gene_name == g]
    if (!length(chr_match)) next
    res_path <- file.path(obj_dir, chr_match[1], "FastTopics_output", g, "res.RDS")
    snap     <- file.path(AG, chr_match[1], "snaptron_output", paste0(g, "_snaptron.tsv"))
    if (!file.exists(res_path) || !file.exists(snap)) next
    res <- poisson2multinom(readRDS(res_path)); gt <- fread(snap, data.table = F)
    if (!(pf %in% colnames(res$F))) next
    unn  <- paste0(gt$start[gt$annotated == 0], "-", gt$end[gt$annotated == 0])
    cand <- rownames(res$F) %in% unn
    if (!any(cand)) next
    top_junc <- rownames(res$F)[cand][which.max(res$F[cand, pf])]
    jc <- as.numeric(str_split_fixed(top_junc, "-", 2))
    m  <- gt[gt$start == jc[1] & gt$end == jc[2], , drop = FALSE]
    if (!nrow(m)) next
    tk <- as.data.frame(c(list(gene_name = g), as.list(m[1, ])), stringsAsFactors = FALSE)
    names(tk) <- snaptron_cols; top_junc_list[[j]] <- tk
  }
  top_junc_df <- bind_rows(top_junc_list) %>% mutate(start = as.numeric(start), end = as.numeric(end))
  cat(sprintf("fact %d: top unannotated junctions resolved %d / %d\n", fact, nrow(top_junc_df), N_TOP))

  make_id <- function(df) paste0(df$gene_name, ":", df$start, "-", df$end)
  mi <- mane_introns %>% dplyr::select(gene_name, strand, start, end)
  canonical <- top_junc_df %>%
    inner_join(mi %>% rename(mane_start = start, mane_end = end),
               by = c("gene_name", "strand", "start" = "mane_start", "end" = "mane_end")) %>%
    mutate(type = "canonical", mane_intron_start = NA_integer_, mane_intron_end = NA_integer_)
  skipped_exon <- top_junc_df %>%
    inner_join(mi %>% rename(up_start = start, up_end = end), by = c("gene_name", "strand", "start" = "up_start")) %>%
    inner_join(mi %>% rename(down_start = start, down_end = end), by = c("gene_name", "strand", "end" = "down_end")) %>%
    filter(up_end < down_start) %>% mutate(type = "skipped exon")
  assigned_ids <- c(make_id(canonical), make_id(skipped_exon))
  start_matches <- top_junc_df %>% filter(!make_id(.) %in% assigned_ids) %>%
    inner_join(mi %>% rename(mane_start = start, mane_end = end), by = c("gene_name", "strand", "start" = "mane_start")) %>%
    mutate(type = if_else(strand == "+", "alt 3'ss", "alt 5'ss"), mane_intron_start = start, mane_intron_end = mane_end)
  end_matches <- top_junc_df %>% filter(!make_id(.) %in% assigned_ids) %>%
    inner_join(mi %>% rename(mane_start = start, mane_end = end), by = c("gene_name", "strand", "end" = "mane_end")) %>%
    mutate(type = if_else(strand == "+", "alt 5'ss", "alt 3'ss"), mane_intron_start = mane_start, mane_intron_end = end)
  annotated <- bind_rows(canonical, skipped_exon, start_matches, end_matches) %>%
    mutate(junc_id = make_id(.)) %>% dplyr::select(junc_id, type, mane_intron_start, mane_intron_end)
  df_final <- top_junc_df %>% mutate(junc_id = make_id(.)) %>% left_join(annotated, by = "junc_id") %>%
    mutate(type = if_else(is.na(type), "other", type)) %>%
    group_by(gene_name, start, end) %>% filter(n() < 2) %>% ungroup()
  a3 <- df_final$type == "alt 3'ss"
  df_final$distance_to_cryptic <- NA_integer_
  df_final$distance_to_cryptic[a3] <- with(df_final[a3, ], if_else(strand == "+", mane_intron_end - end, start - mane_intron_start))
  print(table(df_final$type))
  write_tsv(df_final, file.path(OUT, sprintf("fact_%d_top_junctions.tsv", fact)))
}
