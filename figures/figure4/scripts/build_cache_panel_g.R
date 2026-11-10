# ------------------------------------------------------------------
# Figure 4, panel g -- cache builder.
#
# Everything expensive for the YLPM1 expression ECDF lives here:
#   1. the Leaf2-unproductive target gene set,
#   2. a ROW SLICE of the 7.7 G scaled-RPM matrix (all 11,173 columns,
#      only the ~600 target rows) pulled out with awk,
#   3. the 100-iteration stratified bootstrap,
#   4. the Wilcoxon p-value and the highlight-gene positions.
#
# Writes small artifacts to figure4/cache/ that panel_g/build.R reads.
# Re-run only when the upstream data changes:
#   Rscript --vanilla \
#       scripts/build_cache_panel_g.R
# ------------------------------------------------------------------

if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tibble)
})
source(file.path(QROOT, "figures/scripts/vt_paths.R"))  # vt_path()
HERE     <- normalizePath(file.path(dirname(sub(
  "--file=", "", commandArgs(FALSE)[grep("--file=", commandArgs(FALSE))])), ".."),
  mustWork = FALSE)
if (!nzchar(HERE) || !dir.exists(HERE)) HERE <- getwd()
CACHE <- file.path(HERE, "cache")
dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)

# ---- constants ---------------------------------------------------
N_TOP       <- 20L
N_BOOTSTRAP <- 100L
N_GRID      <- 200L
PSI         <- 0.01
FC          <- 2.5
SEED_BASE   <- 1000L * N_TOP + 1L          # bootstrap_random_cdfs() call site
HIGHLIGHTS  <- c("ATM", "WRN", "FH", "STAG2", "SPOP",
                 "LATS1", "MSH6", "TP53BP1", "XPC", "PBRM1")

EXPR_TSV <- file.path(QROOT, "expression/combined_scaled_rpm.tsv")
SLICE    <- file.path(CACHE, "panel_g_expr.tsv")

# ------------------------------------------------------------------
# 1. Candidate target genes -- Leaf2 unproductive hits.
# ------------------------------------------------------------------
leaf2_class <- fread(file.path(QROOT, "leafcutter_analyses/leaf2_classification/clustering/Leaf2_junction_classifications.txt"), data.table = FALSE)
stopifnot(is.logical(leaf2_class$Coding), is.logical(leaf2_class$UTR))

df <- fread(file.path(QROOT, "leafcutter_analyses/diff_splicing/unified_results/ylpm1/20v60/dpsi0.01_fc2.5/df_final.tsv"),
            data.table = FALSE) |>
  # NB the start-1 offset: df_final is 1-based, Leaf2's Intron_coord is not.
  mutate(start = as.numeric(start) - 1,
         Intron_coord = paste0(chrom, ":", start, "-", end)) |>
  left_join(leaf2_class |>
              select(Intron_coord, Annot, Coding, UTR, GencodePC),
            by = "Intron_coord")

# Alt 3'ss + alt 5'ss only, matching Fig. 4f and the rest of the paper.
hits <- df |> filter(deltapsi > PSI, test / control > FC, type %in% c("alt 3'ss", "alt 5'ss"))
cand <- hits |> filter(!Coding, !UTR) |> pull(gene_name) |> unique()
cat("[cache g] candidate unproductive genes:", length(cand), "\n")

writeLines(cand, file.path(CACHE, "panel_g_cand_genes.txt"))

# ------------------------------------------------------------------
# 2. Row-slice the 7.7 G matrix with awk.
#
# A candidate gene that has no expression row simply never matches, so
# the sliced file's rows ARE the final target set. All ~11k columns are
# kept because the bootstrap samples from the whole column pool.
# ------------------------------------------------------------------
if (!file.exists(SLICE)) {
  cat("[cache g] streaming row slice out of", EXPR_TSV, "...\n")
  cmd <- sprintf(
    "awk -F'\\t' 'NR==FNR{g[$1]=1; next} FNR==1{print; next} ($1 in g) && !seen[$1]++' %s %s > %s",
    shQuote(file.path(CACHE, "panel_g_cand_genes.txt")),
    shQuote(EXPR_TSV), shQuote(SLICE))
  st <- system(cmd)
  stopifnot(st == 0L)
}
cat("[cache g] slice:", round(file.size(SLICE) / 1e6), "MB\n")

expr_mat <- fread(SLICE, header = TRUE, data.table = FALSE) |>
  column_to_rownames("gene_name") |>
  as.matrix()
target <- rownames(expr_mat)
cat("[cache g] target genes in expr_mat:", length(target),
    "| columns:", ncol(expr_mat), "\n")

# ------------------------------------------------------------------
# 3. Top-20 samples on GSP #30.
# ------------------------------------------------------------------
full_info <- fread(vt_path("YLPM1_full_plus_variants"),
                   data.table = FALSE)
ordered <- full_info |> arrange(desc(k30))
top_ids <- as.character(ordered$rail_id[seq_len(N_TOP)])
top_ids <- top_ids[top_ids %in% colnames(expr_mat)]
cat("[cache g] top ids:", length(top_ids), "\n")
writeLines(top_ids, file.path(CACHE, "panel_g_top_ids.txt"))

# ------------------------------------------------------------------
# 4. The red curve: per-gene mean over the top-20 columns.
# ------------------------------------------------------------------
m <- apply(expr_mat[target, top_ids, drop = FALSE], 1, mean, na.rm = TRUE)
red <- tibble(gene_name = names(m), x = as.numeric(m))
fwrite(red, file.path(CACHE, "panel_g_red.tsv"), sep = "\t")

# ------------------------------------------------------------------
# 5. Bootstrap -- 100 stratified random sample groups.
# ------------------------------------------------------------------
stratified_random <- function(top_ids, full_info, seed) {
  proj <- full_info |>
    filter(rail_id %in% as.numeric(top_ids)) |>
    dplyr::count(project_id, name = "n")
  pool <- full_info |>
    mutate(rid = as.character(rail_id)) |>
    filter(rid %in% colnames(expr_mat)) |>
    dplyr::select(rid, project_id)
  set.seed(seed)
  ids <- character()
  for (j in seq_len(nrow(proj))) {
    pp <- pool$rid[pool$project_id == proj$project_id[j]]
    k  <- min(proj$n[j], length(pp))
    if (k > 0) ids <- c(ids, sample(pp, k))
  }
  ids
}

rand_means <- lapply(seq_len(N_BOOTSTRAP), function(it) {
  rid <- stratified_random(top_ids, full_info, seed = SEED_BASE + it)
  if (length(rid) == 0) return(NULL)
  apply(expr_mat[target, rid, drop = FALSE], 1, mean, na.rm = TRUE)
})
rand_means <- Filter(Negate(is.null), rand_means)
all_vals   <- unlist(rand_means)
all_vals   <- all_vals[is.finite(all_vals)]
x_grid     <- seq(min(all_vals), max(all_vals), length.out = N_GRID)
mat        <- vapply(rand_means, function(v) ecdf(v)(x_grid), numeric(N_GRID))

ribbon <- data.frame(
  x        = x_grid,
  mean_cdf = rowMeans(mat),
  lower    = apply(mat, 1, quantile, probs = 0.025, na.rm = TRUE),
  upper    = apply(mat, 1, quantile, probs = 0.975, na.rm = TRUE))
fwrite(ribbon, file.path(CACHE, "panel_g_ribbon.tsv"), sep = "\t")

# ------------------------------------------------------------------
# 6. One-sided Wilcoxon: target-in-top-20 < target-in-random-samples.
# ------------------------------------------------------------------
tx <- red$x[is.finite(red$x)]
pv <- wilcox.test(tx, all_vals, alternative = "less")$p.value
p_lab <- sprintf("p = %.2g", pv)

# ------------------------------------------------------------------
# 7. Highlight positions -- percentile ON THE RED CURVE only.
# ------------------------------------------------------------------
red_rank <- red |>
  filter(is.finite(x)) |>
  mutate(ecdf_y = rank(x, ties.method = "first") / n())
hl <- red_rank |> filter(gene_name %in% HIGHLIGHTS) |>
  select(gene_name, x, ecdf_y) |> arrange(ecdf_y)
fwrite(hl, file.path(CACHE, "panel_g_highlights.tsv"), sep = "\t")

meta <- list(n_target = length(target), n_top = length(top_ids),
             p_value = pv, p_label = p_lab,
             n_highlight = nrow(hl),
             missing_highlight = setdiff(HIGHLIGHTS, hl$gene_name),
             n_boot = length(rand_means))
saveRDS(meta, file.path(CACHE, "panel_g_meta.rds"))

cat("[cache g] n target =", length(target), " ", p_lab,
    " labelled", nrow(hl), "of", length(HIGHLIGHTS),
    if (length(meta$missing_highlight))
      paste("| missing:", paste(meta$missing_highlight, collapse = ", ")) else "",
    "\n")
