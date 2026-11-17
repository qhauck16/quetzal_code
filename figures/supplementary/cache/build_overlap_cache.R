# Cache for supplementary figures S6 (Alsafadi 2016) and S9 (GSP vs
# LeafCutter). The top-250 driving-junction sets are read from
# fact_{10,13,30}_top_junctions.tsv rather than rebuilt from the 3.7 GB
# flashier RDS and 250 per-gene fits.
if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
suppressPackageStartupMessages({ library(tidyverse); library(data.table); library(MatchIt) })
OUT <- file.path(QROOT, "figures/supplementary/cache")
set.seed(1)

gencode <- readRDS(file.path(QROOT, "inputs/gencode_v46_granges.RDS"))
gene_chr_map <- as.data.frame(gencode) %>% filter(type == "exon") %>%
  distinct(gene_name, seqnames, strand) %>% group_by(gene_name) %>%
  summarise(chr = as.character(dplyr::first(seqnames)), strand = as.character(dplyr::first(strand)),
            .groups = "drop")
rpm_psi <- fread(file.path(QROOT, "leafcutter_analyses/diff_splicing/unified_results/shared/all_junction_rpm_and_psi.tsv"),
                 data.table = FALSE)

build_sets <- function(fact, types) {
  df <- fread(file.path(QROOT, "gsp_analysis", sprintf("fact_%d_top_junctions.tsv", fact)), data.table = FALSE)
  pr <- df %>% filter(!is.na(mane_intron_start), !is.na(mane_intron_end)) %>%
    left_join(rpm_psi, by = c("gene_name", "mane_intron_start" = "start",
                              "mane_intron_end" = "end")) %>%
    filter(!is.na(median_psi), !is.na(median_rpm))
  keys <- paste(pr$gene_name, pr$mane_intron_start, pr$mane_intron_end, sep = "|")
  pool <- rpm_psi %>% mutate(key = paste(gene_name, start, end, sep = "|")) %>%
    filter(!key %in% keys, !is.na(median_psi), !is.na(median_rpm))
  both <- bind_rows(
    pr %>% transmute(gene_name, start = mane_intron_start, end = mane_intron_end,
                     median_psi, median_rpm, type, treat = 1L),
    pool %>% transmute(gene_name, start, end, median_psi, median_rpm,
                       type = NA_character_, treat = 0L))
  m <- as.data.frame(match.data(matchit(treat ~ median_psi + median_rpm,
        data = as.data.frame(both), method = "nearest", distance = "mahalanobis",
        ratio = 1L, replace = FALSE)))
  mt <- m %>% filter(treat == 1L) %>% left_join(gene_chr_map, by = "gene_name") %>% filter(!is.na(chr))
  mb <- m %>% filter(treat == 0L) %>% left_join(gene_chr_map, by = "gene_name") %>% filter(!is.na(chr))
  subc <- mt %>% filter(type %in% types) %>% distinct(subclass) %>% pull(subclass)
  list(gsp = df %>% filter(type %in% types) %>%
         transmute(chr = chromosome, strand, start = as.integer(start), end = as.integer(end)) %>%
         distinct(),
       mb  = mb %>% filter(subclass %in% subc) %>%
         transmute(chr, strand, start = as.integer(start), end = as.integer(end)) %>% distinct())
}
count_exact <- function(q, t) nrow(inner_join(distinct(q), distinct(t),
                                              by = c("chr", "strand", "start", "end")))
cells_for <- function(sets, up, down, label) {
  ng <- nrow(sets$gsp); nm <- nrow(sets$mb)
  c4 <- c(count_exact(sets$gsp, up), count_exact(sets$gsp, down),
          count_exact(sets$mb, up),  count_exact(sets$mb, down))
  ft <- fisher.test(matrix(c(c4[1], ng - c4[1], c4[3], nm - c4[3]), nrow = 2))
  data.frame(panel = label, set = rep(c("gsp", "matched"), each = 2),
             direction = rep(c("upregulated", "downregulated"), 2),
             n_overlap = c4, n_set = rep(c(ng, nm), each = 2), fisher_p = ft$p.value)
}

## ---- S6: Alsafadi 2016 ------------------------------------------------
als <- fread(file.path(QROOT, "gsp_analysis/Alsafadi2016_supp1_hg38.tsv"), data.table = FALSE) %>%
  filter(!is.na(start_hg38), !is.na(end_hg38)) %>% mutate(start_hg38 = start_hg38 + 1L) %>%
  distinct(chromosome_hg38, strand_hg38, start_hg38, end_hg38, .keep_all = TRUE) %>%
  mutate(chr = chromosome_hg38, strand = strand_hg38, start = start_hg38, end = end_hg38)
s10 <- build_sets(10, "alt 3'ss")
s6 <- cells_for(s10, als %>% filter(log2fc > 0) %>% select(chr, strand, start, end),
                als %>% filter(log2fc < 0) %>% select(chr, strand, start, end), "alsafadi")
# Control row for S6a: the near-zero-weight junction groups of GSP #10
# (alt 3'ss, none in the top-250 set), as in S9.
source(file.path(QROOT, "figures/supplementary/cache/null_junction_set.R"))
nul10 <- null_junction_set(10, "alt 3'ss",
                           fread(file.path(QROOT, "gsp_analysis/fact_10_top_junctions.tsv"), data.table = FALSE))
als_up <- als %>% filter(log2fc > 0) %>% select(chr, strand, start, end)
als_dn <- als %>% filter(log2fc < 0) %>% select(chr, strand, start, end)
s6 <- bind_rows(s6, data.frame(panel = "alsafadi", set = "null",
                               direction = c("upregulated", "downregulated"),
                               n_overlap = c(count_exact(nul10, als_up), count_exact(nul10, als_dn)),
                               n_set = nrow(nul10), fisher_p = NA_real_))
fwrite(s6, file.path(OUT, "s6_alsafadi_cells.tsv"), sep = "\t")
up_key  <- als %>% filter(log2fc > 0, padj < 0.05) %>% transmute(k = paste(chr, strand, start, end)) %>% distinct() %>% pull(k)
gsp_key <- s10$gsp %>% transmute(k = paste(chr, strand, start, end)) %>% pull(k)
fwrite(data.frame(only_alsafadi = length(setdiff(up_key, gsp_key)),
                  shared = length(intersect(up_key, gsp_key)),
                  only_gsp = length(setdiff(gsp_key, up_key))),
       file.path(OUT, "s6_venn.tsv"), sep = "\t")
print(s6); print(fread(file.path(OUT, "s6_venn.tsv")))

## ---- S7: GSP vs LeafCutter ---------------------------------------------
load_leaf <- function(gene, comp) {
  base <- file.path(QROOT, "leafcutter_analyses/diff_splicing/results", paste0(gene, "_", comp))
  ds  <- fread(file.path(base, "leafcutter_ds_effect_sizes.txt"), data.table = FALSE)
  sig <- fread(file.path(base, "leafcutter_ds_cluster_significance.txt"), data.table = FALSE) %>%
    separate(cluster, into = c("chr_sg", "cluster_id"), sep = ":", remove = FALSE)
  ds <- ds %>% separate(intron, into = c("chr", "start", "end", "cluster_id"), sep = ":",
                        convert = TRUE, remove = FALSE) %>%
    mutate(strand = if_else(grepl("_\\+$", cluster_id), "+", "-"),
           start = as.integer(start) + 1L, end = as.integer(end) - 1L) %>%
    left_join(sig %>% dplyr::select(cluster_id, p.adjust), by = "cluster_id") %>%
    filter(!is.na(p.adjust), p.adjust < 0.05)
  list(up = ds %>% filter(deltapsi > 0) %>% transmute(chr, strand, start, end) %>% distinct(),
       down = ds %>% filter(deltapsi < 0) %>% transmute(chr, strand, start, end) %>% distinct())
}
specs <- list(list(10, "sf3b1", "alt 3'ss"), list(13, "nkap", "alt 3'ss"),
              list(30, "ylpm1", c("alt 3'ss", "alt 5'ss")))
# S9: RPM/PSI-matched MANE introns are not like-for-like with unannotated
# driving junctions, so the baseline is the near-zero-weight junction groups
# of the same GSP (below). Each panel reports the up- and down-regulated LeafCutter
# hits in each set and tests direction among the hits:
# binom.test(up, up + down, p = 0.5), separately per set.
source(file.path(QROOT, "figures/supplementary/cache/null_junction_set.R"))
top_tbl <- function(fact)
  fread(file.path(QROOT, "gsp_analysis", sprintf("fact_%d_top_junctions.tsv", fact)), data.table = FALSE)
gsp_set <- function(fact, types) {
  top_tbl(fact) %>%
    filter(type %in% types) %>%
    transmute(chr = chromosome, strand, start = as.integer(start), end = as.integer(end)) %>%
    distinct()
}
# Baseline: top unannotated junctions of the 250 junction groups with weight
# closest to zero on the GSP, same splice types, none in the GSP's top set
# (null_junction_set.R).
rnd_sets <- lapply(specs, function(sp) null_junction_set(sp[[1]], sp[[3]], top_tbl(sp[[1]])))
s7 <- bind_rows(lapply(seq_along(specs), function(i) {
  sp <- specs[[i]]
  g <- gsp_set(sp[[1]], sp[[3]]); r <- rnd_sets[[i]]
  bind_rows(lapply(c("10v30", "20v60", "30v90"), function(cmp) {
    lf <- load_leaf(sp[[2]], cmp)
    bind_rows(lapply(list(gsp = g, null = r), function(q) {
      up <- count_exact(q, lf$up); dn <- count_exact(q, lf$down)
      data.frame(direction = c("upregulated", "downregulated"), n_overlap = c(up, dn),
                 n_set = nrow(q), binom_p = if (up + dn > 0) binom.test(up, up + dn, p = 0.5)$p.value else NA_real_)
    }), .id = "set") %>%
      mutate(panel = paste(toupper(sp[[2]]), cmp), fact = sp[[1]], .before = 1)
  }))
}))
fwrite(s7, file.path(OUT, "s7_leafcutter_cells.tsv"), sep = "\t")
print(s7 %>% select(panel, set, direction, n_overlap, n_set, binom_p))
