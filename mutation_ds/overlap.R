# Overlap of each mutation comparison's LeafCutter DS junctions (cluster
# p.adjust < 0.05) with its GSP's top-250 driving junctions, and with a control
# set: the top unannotated junctions of the 250 junction groups whose weight on
# the GSP is closest to zero (null_junction_set.R). Per set and direction, the
# number of junctions LeafCutter calls up / down, and a two-sided binomial test of
# up vs down (p0 = 0.5). Used by Supplementary Note 1 Fig. 3b.
# Out: out/overlap_cells.tsv, out/junction_dpsi.tsv (+ printed summary).
if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
suppressPackageStartupMessages({ library(tidyverse); library(data.table) })
MD  <- file.path(QROOT, "mutation_ds")
dir.create(file.path(MD, "out"), showWarnings = FALSE)
set.seed(1)

# Near-zero-weight baseline -- shared with Fig. S9.
source(file.path(QROOT, "figures/supplementary/cache/null_junction_set.R"))
# ---- driving + control sets, and the overlap cells ------------------------
build_sets <- function(fact, types) {
  f  <- file.path(MD, "gsp_top", sprintf("fact_%d_top_junctions.tsv", fact))
  if (!file.exists(f)) f <- file.path(QROOT, "gsp_analysis", sprintf("fact_%d_top_junctions.tsv", fact))  # S7's published sets
  df <- fread(f, data.table = FALSE)
  list(gsp  = df %>% filter(type %in% types) %>%
         transmute(chr = chromosome, strand, start = as.integer(start), end = as.integer(end)) %>%
         distinct(),
       null = null_junction_set(fact, types, exclude_df = df))
}
count_exact <- function(q, t) nrow(inner_join(distinct(q), distinct(t),
                                              by = c("chr", "strand", "start", "end")))
cells_for <- function(sets, up, down, label) {
  bind_rows(lapply(c("gsp", "null"), function(st) {
    q <- sets[[st]]; u <- count_exact(q, up); d <- count_exact(q, down)
    data.frame(panel = label, set = st, direction = c("upregulated", "downregulated"),
               n_overlap = c(u, d), n_set = nrow(q),
               binom_p = if (u + d > 0) binom.test(u, u + d, p = 0.5)$p.value else NA_real_)
  }))
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
# The five mutation comparisons of Supplementary Note 1 (all three type sets).
COMPS <- tibble(comp = c("IDH1_R132", "SPOP_F133", "SPOP_W131", "SPOP_F102", "SPOP_Y87"),
                gsp = c(1L, rep(12L, 4)), base = file.path(MD, "results", comp),
                types = list(names(TYPES)))
FILTERS <- tribble(~min_dpsi, ~min_fc,  0, 1)   # the setting SN1 Fig. 3b uses

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
print(res %>% filter(types == "alt 3'/5'ss") %>%
        select(panel, set, direction, n_overlap, n_set, binom_p) %>% as.data.frame())

# ---- per-junction dPSI for the figure ------------------------------------
# Every GSP-driving / control junction that LeafCutter tested (any p),
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
            inner_join(s$null, d, by = c("chr", "strand", "start", "end")) %>% mutate(set = "null")) %>%
    mutate(comparison = COMPS$comp[i], gsp = COMPS$gsp[i], types = tn)
}))
fwrite(jd, file.path(MD, "out", "junction_dpsi.tsv"), sep = "\t")
