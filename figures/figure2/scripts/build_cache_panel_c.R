.libPaths('/scratch/midway2/qhauck/conda_env/rstudio-server/lib/R/library/')

# Cache builder for Figure 2 panel c (GSP x covariate correlation heatmap).
#
# Panel c's plotting code is cheap; getting to its inputs is not -- the
# flashier fit is a 3.5 G RDS and the AlphaMissense mutation table is
# 536 M. This script runs that expensive chain ONCE and persists only
# what the plot needs, so compose.R never touches the big files.
#
# Source of truth: tcga/genome_wide/flashier/flash_fact_char.Rmd, chunks
#   "loading re-used files", "correlating flashier factors with mutation
#   status", "setting up covariates", "pearson correlation heatmap - all
#   covariates", "correlation heatmap", and "correlation heatmap top 50
#   tissue-ordered".
# The upstream chain is reproduced VERBATIM, quirks included (see the
# max_expl note below). Do not "clean it up" -- the goal is output
# identical to the Rmd.
#
#   Rscript --vanilla build_cache_panel_c.R
#
# Runs for several minutes and needs tens of GB of RAM on a cold run.

if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(flashier)
})

# ---------------------------------------------------------------- paths
FLASH_RDS   <- file.path(QROOT, "quetzal_runs/genome_wide/softimpute_flash_300_qc_filtered.RDS")
META_TSV    <- file.path(QROOT, "inputs/tcga_v2_samples.tsv")
AM_MUT_TSV  <- file.path(QROOT, "mutation_processing/merged_grch38_with_am.tsv")
COVAR_RDS   <- file.path(QROOT, "support_tables/standardized_covariates.RDS")
CACHE_DIR   <- file.path(QROOT, "figures/figure2/cache")

OUT_MATRIX <- file.path(CACHE_DIR, "panel_c_matrix.tsv")
OUT_COLS   <- file.path(CACHE_DIR, "panel_c_cols.tsv")
OUT_ROWS   <- file.path(CACHE_DIR, "panel_c_rows.tsv")

N_TOP      <- 50
TISSUE_THR <- 0.4
TECH_THR   <- 0.4

stage <- function(...) cat(sprintf("[panel_c] %s\n", paste0(...)))
mem   <- function(tag) {
  g <- gc(verbose = FALSE)
  stage(sprintf("  mem after %s: used %.1f MB, max used %.1f MB",
                tag, sum(g[, "used"   ] * c(8, 8)) / 1024,
                     sum(g[, "max used"] * c(8, 8)) / 1024))
}

# ------------------------------------------------------- freshness guard
outs <- c(OUT_MATRIX, OUT_COLS, OUT_ROWS)
if (all(file.exists(outs)) && file.exists(FLASH_RDS) &&
    min(file.info(outs)$mtime) > file.info(FLASH_RDS)$mtime) {
  stage("cache is present and newer than the flashier RDS -- nothing to do.")
  stage("delete the cache/*.tsv files to force a rebuild.")
  quit(save = "no", status = 0)
}

dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)
for (f in c(FLASH_RDS, META_TSV, AM_MUT_TSV, COVAR_RDS)) {
  if (!file.exists(f)) stop("missing input: ", f)
}

# ------------------------------------------- 1. flashier fit -> norm_L
stage("loading flashier RDS (3.5 G, this is the slow part) ...")
out <- readRDS(FLASH_RDS)
stage("computing ldf(type = 'i') ...")
ldf <- flashier::ldf(out, type = 'i')
norm_L <- ldf$L
rm(out, ldf); mem("dropping flashier fit")
stage(sprintf("norm_L: %d samples x %d GSPs", nrow(norm_L), ncol(norm_L)))

# --------------------------------------------------------- 2. metadata
stage("reading TCGA metadata ...")
tcga_metadata <- fread(META_TSV, data.table = FALSE)
tcga_metadata <- tcga_metadata %>%
  mutate(sampleId = sapply(strsplit(tcga_barcode, "-"), function(x) {
    x[4] <- sub("[A-Za-z]+$", "", x[4])   # drop trailing letters from 4th piece
    paste(x[1:4], collapse = "-")
  }))

# ------------------------------------------------ 3. mutation -> TMB
# AlphaMissense-filtered calls: keep am_pathogenicity > 0.9, and keep NA
# (non-missense events -- indels, frameshifts, splice -- are unscored).
stage("reading AlphaMissense mutation table (536 M) ...")
am_mut_data <- fread(AM_MUT_TSV, data.table = FALSE)
mut_data_am09 <- am_mut_data %>%
  filter(is.na(am_pathogenicity) | am_pathogenicity > 0.9)
rm(am_mut_data); mem("dropping raw mutation table")

stage("pivoting sample x gene mutation status ...")
mut_status <- mut_data_am09 %>%
  dplyr::select(hugoGeneSymbol, sampleId) %>%
  mutate(value = 1) %>%
  unique() %>%
  pivot_wider(names_from = hugoGeneSymbol, values_from = value, values_fill = 0) %>%
  right_join(dplyr::select(tcga_metadata,
                           c(rail_id, sampleId, gdc_cases.project.project_id)),
             by = "sampleId") %>%
  dplyr::select(-sampleId) %>%
  column_to_rownames('rail_id') %>%
  na.omit()
rm(mut_data_am09); mem("dropping filtered mutation table")

# Drop gene columns that became all-zero after the join + na.omit. Uses
# name-based selection deliberately: the positional form silently
# recycles values across columns when any gene is zero-only.
.meta_col   <- "gdc_cases.project.project_id"
.gene_cols  <- setdiff(colnames(mut_status), .meta_col)
.keep_genes <- .gene_cols[colSums(mut_status[, .gene_cols, drop = FALSE]) > 0]
mut_status  <- mut_status[, c(.keep_genes, .meta_col), drop = FALSE]
stage(sprintf("mut_status: %d samples x %d genes",
              nrow(mut_status), length(.keep_genes)))

# ------------------------------------------------ 4. standardized covars
stage("reading standardized covariates ...")
pre_model <- readRDS(COVAR_RDS) %>% mutate(rail_id = as.numeric(rail_id))

# -------------------------------------------- 5. covariate matrix build
stage("building covariate matrix ...")
covar_df <- tcga_metadata %>%
  filter(rail_id %in% as.numeric(rownames(norm_L))) %>%
  mutate(rail_id_chr = as.character(rail_id))

tissue_dummies <- model.matrix(~ 0 + gdc_cases.project.project_id, data = covar_df)
colnames(tissue_dummies) <- gsub("gdc_cases.project.project_id", "",
                                 colnames(tissue_dummies))
rownames(tissue_dummies) <- covar_df$rail_id_chr

tech_df <- pre_model %>%
  transmute(
    rail_id_chr      = as.character(rail_id),
    is_ffpe          = as.numeric(gdc_cases.samples.is_ffpe == "TRUE" |
                                  gdc_cases.samples.is_ffpe == TRUE),
    is_male          = as.numeric(gdc_cases.demographic.gender == "male"),
    avg_len          = avg_len,
    junction_count   = junction_count,
    junction_avg_cov = junction_avg_coverage,
    pct_C            = pct_C,
    pct_G            = pct_G,
    RIN              = analytes.rna_integrity_number
  ) %>%
  column_to_rownames("rail_id_chr")

common <- Reduce(intersect, list(rownames(norm_L), rownames(tissue_dummies),
                                 rownames(tech_df)))
L_aligned      <- norm_L[common, ]
tissue_aligned <- tissue_dummies[common, ]
tech_aligned   <- as.matrix(tech_df[common, ])
stage(sprintf("aligned on %d common samples", length(common)))

# TMB. Note tmb[common] yields NA for samples absent from mut_status;
# cor(use = "complete.obs") drops those pairs. Matches the Rmd.
tmb    <- rowSums(mut_status[, setdiff(colnames(mut_status), .meta_col)],
                  na.rm = TRUE)
tmb_df <- data.frame(TMB = tmb[common], row.names = common)

all_covars <- cbind(tech_aligned, tmb_df, tissue_aligned)
rm(norm_L, tissue_dummies, tech_df, mut_status, tcga_metadata, pre_model, covar_df)
mem("dropping intermediates")

# ------------------------------------------------ 6. correlation matrix
n_gsps <- ncol(L_aligned)
stage(sprintf("correlating %d GSPs x %d covariates ...",
              n_gsps, ncol(all_covars)))
corr_mat <- matrix(NA, nrow = n_gsps, ncol = ncol(all_covars))
for (k in seq_len(n_gsps)) {
  for (j in seq_len(ncol(all_covars))) {
    corr_mat[k, j] <- cor(L_aligned[, k], all_covars[, j], use = "complete.obs")
  }
}
rownames(corr_mat) <- seq_len(n_gsps)
colnames(corr_mat) <- colnames(all_covars)

tissue_cols <- colnames(tissue_aligned)
tech_cols   <- c(colnames(tech_aligned), "TMB")

# Top1/top2 loading ratio -> "single sample GSP" flag
top_ratios <- data.frame(gsp = seq_len(n_gsps))
top_ratios$ratio_1v2 <- sapply(seq_len(n_gsps), function(k) {
  sorted <- sort(L_aligned[, k], decreasing = TRUE)
  sorted[1] / sorted[2]
})
rm(L_aligned, all_covars, tissue_aligned, tech_aligned); mem("post-correlation")

# ------------------------------------------------- 7. summary columns
corr_df <- as.data.frame(corr_mat)
corr_df$gsp           <- seq_len(n_gsps)
corr_df$max_tissue    <- apply(corr_mat[, tissue_cols, drop = FALSE], 1,
                               function(x) max(abs(x)))
corr_df$max_technical <- apply(corr_mat[, tech_cols, drop = FALSE], 1,
                               function(x) max(abs(x)))
corr_df$`top1/2>5`    <- as.numeric(top_ratios$ratio_1v2 > 5)

summary_cols <- c("max_tissue", "max_technical", "top1/2>5", "max_expl")
# Verbatim from the Rmd: max_expl is taken over the tissue and technical
# columns PLUS the other three summary columns, so it partly summarises
# summaries. Preserved deliberately -- changing it would shift both the
# "Max explained" column and the [Unknown] row labels.
corr_df$max_expl <- apply(corr_df[, c(tissue_cols, tech_cols,
                                      setdiff(summary_cols, "max_expl"))],
                          1, function(x) max(abs(x)))

stopifnot(all(tissue_cols %in% names(corr_df)),
          all(tech_cols   %in% names(corr_df)))

display_cols <- c(tissue_cols, tech_cols, summary_cols)
heatmap_mat  <- as.matrix(corr_df[, display_cols])
rownames(heatmap_mat) <- seq_len(n_gsps)
heatmap_mat_trunc <- pmin(pmax(heatmap_mat, -0.7), 0.7)   # truncate at +/-0.7

# ---------------------------------------------------- 8. row labels
stage("deriving row labels ...")
gsps       <- as.integer(rownames(heatmap_mat_trunc))
row_labels <- rep("", length(gsps))

tissue_cor      <- corr_df[, tissue_cols, drop = FALSE]
tissue_max_abs  <- apply(abs(tissue_cor), 1, max)
tissue_top_idx  <- apply(abs(tissue_cor), 1, which.max)
tissue_top_name <- gsub("TCGA-", "", tissue_cols[tissue_top_idx])

tech_cor      <- corr_df[, tech_cols, drop = FALSE]
tech_max_abs  <- apply(abs(tech_cor), 1, max)
tech_top_idx  <- apply(abs(tech_cor), 1, which.max)
tech_top_name <- tech_cols[tech_top_idx]

unexplained <- corr_df$max_expl < TISSUE_THR

# Whichever of tissue/technical has the higher |r| names the GSP.
for (i in seq_along(gsps)) {
  if (tissue_max_abs[i] >= tech_max_abs[i] && tissue_max_abs[i] > TISSUE_THR) {
    row_labels[i] <- sprintf("GSP %d [%s]", gsps[i], tissue_top_name[i])
  } else if (tech_max_abs[i] > TECH_THR) {
    row_labels[i] <- sprintf("GSP %d [%s]", gsps[i], tech_top_name[i])
  }
}
# Unexplained GSPs, without overwriting an assigned label.
for (i in seq_along(gsps)) {
  if (row_labels[i] == "" && unexplained[i]) {
    row_labels[i] <- sprintf("GSP %d [Unknown]", gsps[i])
  }
}

# ------------------------------------- 9. top 50, tissue-ordered columns
heatmap_mat_top <- heatmap_mat_trunc[seq_len(N_TOP), , drop = FALSE]
row_labels_top  <- row_labels[seq_len(N_TOP)]

# Order tissue columns by which GSP explains them best; ties go to the
# lower GSP index (which.max takes the first).
tissue_mat          <- heatmap_mat_top[, tissue_cols, drop = FALSE]
best_gsp_idx        <- apply(abs(tissue_mat), 2, which.max)
tissue_order        <- order(best_gsp_idx)
tissue_cols_ordered <- tissue_cols[tissue_order]

display_cols_o    <- c(tissue_cols_ordered, tech_cols, summary_cols)
heatmap_mat_top_o <- heatmap_mat_top[, display_cols_o, drop = FALSE]
col_group_o <- ifelse(display_cols_o %in% tissue_cols, "Tissue",
               ifelse(display_cols_o %in% tech_cols,   "Technical", "Summary"))

pretty_col_label <- function(x) {
  lookup <- c(
    is_ffpe          = "FFPE status",
    is_male          = "Sex",
    avg_len          = "Avg read length",
    junction_count   = "Junction count",
    junction_avg_cov = "Avg junction cov",
    pct_C            = "%C",
    pct_G            = "%G",
    RIN              = "RIN",
    TMB              = "TMB",
    max_tissue       = "Max tissue |r|",
    max_technical    = "Max technical |r|",
    `top1/2>5`       = "Single sample GSP",
    max_expl         = "Max explained"
  )
  ifelse(x %in% names(lookup), lookup[x], x)
}
# Tissue columns show the bare TCGA code; everything else gets the
# human-readable name.
display_labels_o <- vapply(display_cols_o, function(x) {
  if (x %in% tissue_cols) sub("^TCGA-", "", x) else pretty_col_label(x)
}, character(1), USE.NAMES = FALSE)

# ------------------------------------------------------- 10. write cache
stage("writing cache ...")

mat_out <- as.data.frame(heatmap_mat_top_o, check.names = FALSE)
mat_out <- cbind(gsp = as.integer(rownames(heatmap_mat_top_o)), mat_out)
readr::write_tsv(mat_out, OUT_MATRIX)

readr::write_tsv(tibble(
  column        = display_cols_o,
  group         = col_group_o,
  label         = display_labels_o,
  display_order = seq_along(display_cols_o)
), OUT_COLS)

readr::write_tsv(tibble(
  gsp           = as.integer(rownames(heatmap_mat_top_o)),
  row_label     = row_labels_top,
  display_order = seq_len(N_TOP)
), OUT_ROWS)

stage(sprintf("wrote %s  (%d x %d)", basename(OUT_MATRIX),
              nrow(mat_out), ncol(mat_out) - 1L))
stage(sprintf("wrote %s  (%d columns: %d tissue / %d technical / %d summary)",
              basename(OUT_COLS), length(display_cols_o),
              sum(col_group_o == "Tissue"), sum(col_group_o == "Technical"),
              sum(col_group_o == "Summary")))
stage(sprintf("wrote %s  (%d GSPs, %d labelled)", basename(OUT_ROWS),
              N_TOP, sum(row_labels_top != "")))
mem("final")
stage("done.")
