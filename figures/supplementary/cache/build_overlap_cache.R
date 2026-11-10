# Cache for supplementary figures S6 (Alsafadi 2016) and S7 (GSP vs
# LeafCutter). Logic is spec_flash_fact.Rmd's ("alsafadi GSP10 overlap",
# "matched background junctions by RPM and PSI", "4-cell overlap
# heatmap ...", "venn diagram ...", "flashier-vs-leafcutter GSP
# comparison panels"), except that the top-250 driving-junction sets are
# read from fact_{10,13,30}_top_junctions.tsv -- the Rmd's own saved
# output of build_gsp_top() + assign_types() -- instead of being rebuilt
# from the 3.7 GB flashier RDS and 250 per-gene fits. Counts check out
# against the published panels (GSP #10 alt 3'ss n = 168, GSP #13 203,
# GSP #30 alt 3'/5'ss 188).
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
# S9 (2026-11-10, Quinn): no matched baseline -- the RPM/PSI-matched MANE introns
# were not like-for-like with unannotated driving junctions. Each panel reports the
# driving junctions that are LeafCutter up- and down-regulated hits, and tests
# direction among the hits: binom.test(up, up + down, p = 0.5).
gsp_set <- function(fact, types) {
  fread(file.path(QROOT, "gsp_analysis", sprintf("fact_%d_top_junctions.tsv", fact)), data.table = FALSE) %>%
    filter(type %in% types) %>%
    transmute(chr = chromosome, strand, start = as.integer(start), end = as.integer(end)) %>%
    distinct()
}
s7 <- bind_rows(lapply(specs, function(sp) {
  g <- gsp_set(sp[[1]], sp[[3]])
  bind_rows(lapply(c("10v30", "20v60", "30v90"), function(cmp) {
    lf <- load_leaf(sp[[2]], cmp)
    up <- count_exact(g, lf$up); dn <- count_exact(g, lf$down)
    data.frame(panel = paste(toupper(sp[[2]]), cmp), set = "gsp",
               direction = c("upregulated", "downregulated"), n_overlap = c(up, dn),
               n_set = nrow(g), binom_p = binom.test(up, up + dn, p = 0.5)$p.value,
               fact = sp[[1]])
  }))
}))
fwrite(s7, file.path(OUT, "s7_leafcutter_cells.tsv"), sep = "\t")
print(s7 %>% select(panel, direction, n_overlap, n_set, binom_p))
