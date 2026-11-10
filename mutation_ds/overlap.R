.libPaths('/scratch/midway2/qhauck/conda_env/rstudio-server/lib/R/library/')
# Overlap of each mutation comparison's LeafCutter DS junctions with its GSP's
# top-250 driving junctions -- the S7 test (figures/supplementary/cache/
# build_overlap_cache.R), unchanged: GSP alt-ss junctions vs a background of
# junctions matched 1:1 on the MANE intron's median PSI + RPM (Mahalanobis),
# Fisher's exact test on the overlap with LeafCutter-upregulated junctions
# (p.adjust < 0.05, dPSI > 0). Matching needs a MANE partner intron, so only
# alt 3'ss / alt 5'ss junctions can be tested (as in S7).
# Out: out/overlap_cells.tsv (+ printed summary).
if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
suppressPackageStartupMessages({ library(tidyverse); library(data.table); library(MatchIt) })
MD  <- file.path(QROOT, "mutation_ds")
dir.create(file.path(MD, "out"), showWarnings = FALSE)
set.seed(1)

gencode <- readRDS(file.path(QROOT, "inputs/gencode_v46_granges.RDS"))
gene_chr_map <- as.data.frame(gencode) %>% filter(type == "exon") %>%
  distinct(gene_name, seqnames, strand) %>% group_by(gene_name) %>%
  summarise(chr = as.character(dplyr::first(seqnames)), strand = as.character(dplyr::first(strand)),
            .groups = "drop")
rpm_psi <- fread(file.path(QROOT, "leafcutter_analyses/diff_splicing/unified_results/shared/all_junction_rpm_and_psi.tsv"),
                 data.table = FALSE)

# ---- build_sets / cells_for: copied from build_overlap_cache.R ----------
build_sets <- function(fact, types) {
  f  <- file.path(MD, "gsp_top", sprintf("fact_%d_top_junctions.tsv", fact))
  if (!file.exists(f)) f <- file.path(QROOT, "gsp_analysis", sprintf("fact_%d_top_junctions.tsv", fact))  # S7's published sets
  df <- fread(f, data.table = FALSE)
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
             n_overlap = c4, n_set = rep(c(ng, nm), each = 2), fisher_p = ft$p.value,
             odds_ratio = unname(ft$estimate))
}
# min_dpsi: optional effect-size floor on |deltapsi| (S7 itself uses 0).
# min_fc: optional PSI fold-change floor in the junction's direction
# (up: test / control; down: control / test), as in the diff-splicing
# df_final tables (dPSI > 0.05, FC > 2.5). The effect-size file names its PSI
# columns by group, so test / control is right whatever the group order.
load_leaf <- function(base, min_dpsi = 0, min_fc = 1) {
  ds  <- fread(file.path(base, "leafcutter_ds_effect_sizes.txt"), data.table = FALSE)
  # leafcutter_ds.R encodes groups by order of first appearance in the
  # phenotype file (its log: "Encoding as <first> =0, <second> =1"), and
  # deltapsi = group 1 - group 0. Orient it as mutant - control either way.
  first_grp <- fread(file.path(base, "leaf_pheno_data.txt"), header = FALSE, nrows = 1)$V2
  if (first_grp == "test") ds$deltapsi <- -ds$deltapsi
  sig <- fread(file.path(base, "leafcutter_ds_cluster_significance.txt"), data.table = FALSE) %>%
    separate(cluster, into = c("chr_sg", "cluster_id"), sep = ":", remove = FALSE)
  ds <- ds %>% separate(intron, into = c("chr", "start", "end", "cluster_id"), sep = ":",
                        convert = TRUE, remove = FALSE) %>%
    mutate(strand = if_else(grepl("_\\+$", cluster_id), "+", "-"),
           start = as.integer(start) + 1L, end = as.integer(end) - 1L) %>%
    left_join(sig %>% dplyr::select(cluster_id, p.adjust), by = "cluster_id") %>%
    filter(!is.na(p.adjust), p.adjust < 0.05, abs(deltapsi) > min_dpsi) %>%
    filter(ifelse(deltapsi > 0, test / control, control / test) >= min_fc)
  list(up = ds %>% filter(deltapsi > 0) %>% transmute(chr, strand, start, end) %>% distinct(),
       down = ds %>% filter(deltapsi < 0) %>% transmute(chr, strand, start, end) %>% distinct())
}

A35 <- c("alt 3'ss", "alt 5'ss")
TYPES <- list(`alt 3'ss` = "alt 3'ss", `alt 5'ss` = "alt 5'ss", `alt 3'/5'ss` = A35)
# Mutation comparisons (all three type sets) + the nine S7 GSP comparisons with
# S7's own type sets, so the min_dpsi = 0 rows re-derive the published S7.
COMPS <- bind_rows(
  tibble(comp = c("IDH1_R132", rep(c("SPOP_F133", "SPOP_W131", "SPOP_F102", "SPOP_Y87"), 2)),
         gsp = c(1L, rep(12L, 4), rep(40L, 4)), base = file.path(MD, "results", comp),
         types = list(names(TYPES))),
  expand_grid(g = c("sf3b1", "nkap", "ylpm1"), cmp = c("10v30", "20v60", "30v90")) %>%
    transmute(comp = paste(toupper(g), cmp), gsp = c(sf3b1 = 10L, nkap = 13L, ylpm1 = 30L)[g],
              base = file.path(QROOT, "leafcutter_analyses/diff_splicing/results", paste0(g, "_", cmp)),
              types = ifelse(g == "ylpm1", list("alt 3'/5'ss"), list("alt 3'ss"))))
FILTERS <- tribble(~min_dpsi, ~min_fc,  0, 1,  0.05, 1,  0, 2.5,  0.05, 2.5)

sets <- list()
res <- bind_rows(lapply(seq_len(nrow(COMPS)), function(i) bind_rows(lapply(seq_len(nrow(FILTERS)), function(f) {
  md <- FILTERS$min_dpsi[f]; mf <- FILTERS$min_fc[f]
  lf <- load_leaf(COMPS$base[i], md, mf)
  bind_rows(lapply(COMPS$types[[i]], function(tn) {
    key <- paste(COMPS$gsp[i], tn)
    if (is.null(sets[[key]])) sets[[key]] <<- build_sets(COMPS$gsp[i], TYPES[[tn]])
    cells_for(sets[[key]], lf$up, lf$down, COMPS$comp[i]) %>%
      mutate(gsp = COMPS$gsp[i], types = tn, min_dpsi = md, min_fc = mf,
             n_ds_up = nrow(lf$up), n_ds_down = nrow(lf$down))
  }))
}))))
fwrite(res, file.path(MD, "out", "overlap_cells.tsv"), sep = "\t")
res %>% filter(set == "gsp", direction == "upregulated") %>%
  left_join(res %>% filter(set == "matched", direction == "upregulated") %>%
              select(panel, gsp, types, min_dpsi, min_fc, bg_overlap = n_overlap, bg_n = n_set),
            by = c("panel", "gsp", "types", "min_dpsi", "min_fc")) %>%
  transmute(comparison = panel, gsp, types, min_dpsi, min_fc, n_ds_up, gsp_hits = sprintf("%d/%d", n_overlap, n_set),
            bg_hits = sprintf("%d/%d", bg_overlap, bg_n), odds_ratio = signif(odds_ratio, 3),
            fisher_p = signif(fisher_p, 3)) %>% as.data.frame() %>% print()

# ---- per-junction dPSI for the figure ------------------------------------
# Every GSP / matched-background junction that LeafCutter tested (any p),
# with its mutant - control dPSI and cluster p.adjust. Type set per row as
# above (alt 3'/5'ss for the mutation comparisons).
load_dpsi <- function(base) {
  ds  <- fread(file.path(base, "leafcutter_ds_effect_sizes.txt"), data.table = FALSE)
  if (fread(file.path(base, "leaf_pheno_data.txt"), header = FALSE, nrows = 1)$V2 == "test")
    ds$deltapsi <- -ds$deltapsi
  sig <- fread(file.path(base, "leafcutter_ds_cluster_significance.txt"), data.table = FALSE) %>%
    separate(cluster, into = c("chr_sg", "cluster_id"), sep = ":", remove = FALSE)
  ds %>% separate(intron, into = c("chr", "start", "end", "cluster_id"), sep = ":",
                  convert = TRUE, remove = FALSE) %>%
    mutate(strand = if_else(grepl("_\\+$", cluster_id), "+", "-"),
           start = as.integer(start) + 1L, end = as.integer(end) - 1L) %>%
    left_join(sig %>% dplyr::select(cluster_id, p.adjust), by = "cluster_id") %>%
    filter(!is.na(p.adjust)) %>% dplyr::select(chr, strand, start, end, deltapsi, psi_test = test,
                                               psi_control = control, p.adjust)
}
jd <- bind_rows(lapply(seq_len(nrow(COMPS)), function(i) {
  tn <- if ("alt 3'/5'ss" %in% COMPS$types[[i]]) "alt 3'/5'ss" else COMPS$types[[i]][1]
  s  <- sets[[paste(COMPS$gsp[i], tn)]]; d <- load_dpsi(COMPS$base[i])
  bind_rows(inner_join(s$gsp, d, by = c("chr", "strand", "start", "end")) %>% mutate(set = "gsp"),
            inner_join(s$mb,  d, by = c("chr", "strand", "start", "end")) %>% mutate(set = "matched")) %>%
    mutate(comparison = COMPS$comp[i], gsp = COMPS$gsp[i], types = tn)
}))
fwrite(jd, file.path(MD, "out", "junction_dpsi.tsv"), sep = "\t")
