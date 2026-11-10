.libPaths('/scratch/midway2/qhauck/conda_env/rstudio-server/lib/R/library/')

# Cache builder for Figure 2 panel g ("Mutation enrichment discovers
# additional mutation-GSP relationships").
#
# Panel g's CMH half is already cached on disk
# (cmh_mutation_enrichment_am0_9.tsv), but its LINEAR half is not: it
# comes from a MatrixEQTL run that needs a 3.7 GB flashier RDS and a
# 536 MB AlphaMissense mutation table. This script reproduces the
# **AM > 0.9 pass** of flash_fact_char.Rmd (the pass whose output,
# paper_figures/mutation_test_am0_9.pdf, matches the published panel --
# the 0.34 and 0.564 passes do not) and persists the final joined table
# so compose.R can draw the scatter in seconds.
#
# This is a faithful port of flash_fact_char.Rmd, NOT an improvement:
#   - chunk "loading re-used files"                    (lines ~219-230)
#   - chunk "correlating flashier factors ..."         (lines ~234-272)
#   - build_mut_status_am()                            (lines ~289-317)
#   - chunk "finding GSPs we don't care about"         (lines ~356-370)
#   - chunk "third pass: AM threshold 0.9"             (lines ~1533-1600)
#   - both_tests_09 assembly                           (lines ~1686-1713)
# Quirks below that look like bugs are deliberate matches to the Rmd;
# they are flagged with "Rmd parity:".
#
# Runtime: ~20-40 min cold (dominated by the RDS load + MatrixEQTL).
# No plotting happens here.
#
#   Rscript --vanilla build_cache_panel_g.R

if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(readr)
  library(MatrixEQTL)
  library(flashier)
})

# ---- absolute paths (runnable from any cwd) ------------------------
FIG2_DIR     <- file.path(QROOT, "figures/figure2")

RDS_FLASH  <- file.path(QROOT, "quetzal_runs/genome_wide/softimpute_flash_300_qc_filtered.RDS")
AM_MUT_TSV <- file.path(QROOT, "mutation_processing/merged_grch38_with_am.tsv")
CMH_TSV    <- file.path(QROOT, "gsp_analysis/cmh_mutation_enrichment_am0_9.tsv")
META_TSV   <- file.path(QROOT, "inputs/tcga_v2_samples.tsv")
COVAR_RDS  <- file.path(QROOT, "support_tables/standardized_covariates.RDS")

CACHE_DIR    <- file.path(FIG2_DIR, "cache")
OUT_TESTS    <- file.path(CACHE_DIR, "panel_g_both_tests.tsv")
OUT_THRESH   <- file.path(CACHE_DIR, "panel_g_thresholds.tsv")

dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- guard: skip the expensive load if the cache is already fresh ---
if (file.exists(OUT_TESTS) && file.exists(OUT_THRESH) &&
    file.info(OUT_TESTS)$mtime  > file.info(RDS_FLASH)$mtime &&
    file.info(OUT_THRESH)$mtime > file.info(RDS_FLASH)$mtime) {
  cat("[panel_g] Cache is newer than the flashier RDS -- nothing to do.\n")
  cat("[panel_g]   ", OUT_TESTS,  "\n")
  cat("[panel_g]   ", OUT_THRESH, "\n")
  cat("[panel_g] Delete those files to force a rebuild.\n")
  quit(save = "no", status = 0)
}

for (f in c(RDS_FLASH, AM_MUT_TSV, CMH_TSV, META_TSV, COVAR_RDS)) {
  if (!file.exists(f)) stop("Missing required input: ", f)
}

t0 <- Sys.time()
step <- function(msg) cat(sprintf("[panel_g][%5.1f min] %s\n",
  as.numeric(difftime(Sys.time(), t0, units = "mins")), msg))

# ---- 1. flashier loadings -----------------------------------------
step("Loading flashier RDS (3.7 GB, slow) ...")
out <- readRDS(RDS_FLASH)
step("Computing ldf() ...")
ldf <- flashier::ldf(out, type = 'i')
norm_L <- ldf$L
rm(out, ldf); gc(verbose = FALSE)
step(sprintf("norm_L: %d samples x %d GSPs. RDS freed.",
             nrow(norm_L), ncol(norm_L)))

# ---- 2. metadata ---------------------------------------------------
step("Reading TCGA metadata ...")
tcga_metadata <- fread(META_TSV, data.table = FALSE)
tcga_metadata <- tcga_metadata %>%
  mutate(sampleId = sapply(strsplit(tcga_barcode, "-"), function(x) {
    x[4] <- sub("[A-Za-z]+$", "", x[4])   # drop trailing letters from 4th piece
    paste(x[1:4], collapse = "-")
  }))

# ---- 3. AlphaMissense-filtered sample x gene matrix at AM > 0.9 -----
step("Reading AlphaMissense mutation table (536 MB) ...")
am_mut_data <- fread(AM_MUT_TSV, data.table = FALSE)

# Verbatim build_mut_status_am() from the Rmd. A mutation counts if
# am_pathogenicity is NA (indels/frameshifts/splice -- unscored by
# AlphaMissense) or exceeds the threshold.
build_mut_status_am <- function(thresh) {
  m <- am_mut_data %>%
    filter(is.na(am_pathogenicity) | am_pathogenicity > thresh) %>%
    dplyr::select(hugoGeneSymbol, sampleId) %>%
    mutate(value = 1) %>%
    unique() %>%
    pivot_wider(names_from = hugoGeneSymbol,
                values_from = value, values_fill = 0) %>%
    right_join(dplyr::select(tcga_metadata,
                             c(rail_id, sampleId, gdc_cases.project.project_id)),
               by = "sampleId") %>%
    dplyr::select(-sampleId) %>%
    column_to_rownames('rail_id') %>%
    na.omit()
  # Drop all-zero gene columns BY NAME. The Rmd notes that the earlier
  # positional reassignment silently recycled column data whenever any
  # gene was zero-only; keep the name-based form.
  meta_col <- "gdc_cases.project.project_id"
  gene_cols <- setdiff(colnames(m), meta_col)
  keep_genes <- gene_cols[colSums(m[, gene_cols, drop = FALSE]) > 0]
  m[, c(keep_genes, meta_col), drop = FALSE]
}

step("Building mut_status_am_09 (AM > 0.9) ...")
mut_status_am <- build_mut_status_am(0.9)
rm(am_mut_data); gc(verbose = FALSE)
step(sprintf("mut_status_am_09: %d samples x %d cols. Mutation table freed.",
             nrow(mut_status_am), ncol(mut_status_am)))

# ---- 4. single-sample GSPs to drop ---------------------------------
# A GSP whose top loading is > 5x its third-highest is driven by one or
# two samples; those columns are dropped from the phenotype matrix.
selected_cols <- apply(norm_L, 2, function(col) {
  sorted_vals <- sort(col, decreasing = TRUE)
  max_val <- sorted_vals[1]
  third_val <- if (length(sorted_vals) >= 3) sorted_vals[3] else sorted_vals[length(sorted_vals)]
  max_val > 5 * third_val
})
single_sample_factors <- seq(1:ncol(norm_L))[selected_cols]
step(sprintf("Dropping %d single-sample GSP(s).", length(single_sample_factors)))

# ---- 5. AM>0.9 pass: hard TMB cap, then per-gene count filters ------
# Rmd parity: the 0.9 pass uses a hard TMB cap of 500 rather than the
# mean+3*SD rule the 0.34/0.564 passes use, and applies the per-gene
# count filters AFTER that cap.
tmb_09 <- rowSums(mut_status_am[, setdiff(colnames(mut_status_am),
                                          "gdc_cases.project.project_id"),
                                drop = FALSE], na.rm = TRUE)
keep_samples_09 <- names(tmb_09)[tmb_09 <= 500]
step(sprintf("TMB cap = 500; dropping %d sample(s); keeping %d",
             sum(tmb_09 > 500), length(keep_samples_09)))
mut_status_am <- mut_status_am[keep_samples_09, , drop = FALSE]

# Rmd parity: these two lines index with a logical vector over the gene
# columns plus a trailing TRUE that keeps the project_id meta column.
mut_status_am <- mut_status_am[, c(colSums(mut_status_am[, 1:(ncol(mut_status_am) - 1)]) > 10, T)]
mut_status_am <- mut_status_am[, c(colSums(mut_status_am[, 1:(ncol(mut_status_am) - 1)]) < 1000, T)]
step(sprintf("After 10 < colSums < 1000: %d samples x %d cols",
             nrow(mut_status_am), ncol(mut_status_am)))

# ---- 6. align samples, build MatrixEQTL inputs ----------------------
common_samples <- intersect(rownames(norm_L), rownames(mut_status_am))
L_pm_filtered <- norm_L[common_samples, , drop = FALSE]
mut_status_ordered <- mut_status_am[common_samples, , drop = FALSE]
L_pm_filtered <- L_pm_filtered[, -single_sample_factors]

geno <- as.matrix(mut_status_ordered[, setdiff(colnames(mut_status_ordered),
                                               "gdc_cases.project.project_id")])

step("Reading standardized covariates ...")
pre_model <- readRDS(COVAR_RDS) %>% mutate(rail_id = as.numeric(rail_id))

full_covar <- pre_model %>%
  filter(rail_id %in% common_samples) %>%
  filter(complete.cases(.)) %>%
  column_to_rownames('rail_id')
full_covar <- full_covar[intersect(common_samples, rownames(full_covar)),
                         , drop = FALSE]

phenos <- as.matrix(L_pm_filtered)
colnames(phenos) <- setdiff(seq_len(ncol(norm_L)), single_sample_factors)

covar_matrix <- model.matrix(~ ., data = full_covar)
qr_decomposition <- qr(covar_matrix)
matrix_rank <- qr_decomposition$rank
num_cols <- ncol(covar_matrix)
if (matrix_rank < num_cols) {
  dependent_indices <- qr_decomposition$pivot[(matrix_rank + 1):num_cols]
  cat("[panel_g] Removing colinear covariate columns: ",
      paste(colnames(covar_matrix)[dependent_indices], collapse = ", "), "\n")
  covar_matrix_final <- covar_matrix[, -dependent_indices]
} else {
  cat("[panel_g] No colinear covariate columns detected.\n")
  covar_matrix_final <- covar_matrix
}
covar_matrix_final <- covar_matrix_final[, 2:ncol(covar_matrix_final)]  # drop intercept

covar_sliced <- SlicedData$new()
covar_sliced$CreateFromMatrix(t(covar_matrix_final))

filt_geno   <- geno[rownames(geno) %in% rownames(covar_matrix_final), ]
filt_phenos <- phenos[rownames(phenos) %in% rownames(covar_matrix_final), ]
snps <- SlicedData$new(); snps$CreateFromMatrix(t(filt_geno))
gene <- SlicedData$new(); gene$CreateFromMatrix(t(filt_phenos))
step(sprintf("MatrixEQTL inputs: %d genes(GSPs) x %d snps(mutated genes) x %d samples",
             ncol(filt_phenos), ncol(filt_geno), nrow(filt_geno)))

# ---- 7. linear model ------------------------------------------------
step("Running Matrix_eQTL_main (slow) ...")
meqtl_res_09 <- Matrix_eQTL_main(
  snps = snps, gene = gene, cvrt = covar_sliced,
  output_file_name = tempfile(),
  pvOutputThreshold = 1, useModel = modelLINEAR,
  errorCovariance = numeric(),
  verbose = TRUE, pvalue.hist = "qqplot"
)
rm(snps, gene, covar_sliced, filt_geno, filt_phenos, geno, phenos); gc(verbose = FALSE)
step("MatrixEQTL done.")

# ---- 8. join to the cached CMH table --------------------------------
# Rmd parity: ncol(mut_status_am) here INCLUDES the trailing
# project_id meta column, so the Bonferroni denominator is
# 50 * (n_genes + 1). Kept as-is to match the published panel.
bonf_thresh_09 <- 0.05 / (50 * ncol(mut_status_am))
neg_log10_bonf_09 <- -log10(bonf_thresh_09)
step(sprintf("Bonferroni threshold = %.3g  (-log10 = %.3f)",
             bonf_thresh_09, neg_log10_bonf_09))

cmh_enrich_09 <- fread(CMH_TSV, data.table = FALSE)

both_tests_09 <- meqtl_res_09$all$eqtls %>%
  filter(statistic > 0) %>%
  filter(as.numeric(gene) < 51) %>%
  filter(pvalue < bonf_thresh_09) %>%
  mutate(linear_pval = pvalue, id = paste0(gene, '-', snps)) %>%
  dplyr::select(id, linear_pval) %>%
  left_join(cmh_enrich_09 %>%
              mutate(cmh_pval = pvalue, id = paste0(gsp, '-', gene)) %>%
              dplyr::select(id, cmh_pval), by = 'id')
both_tests_09$cmh_pval[is.na(both_tests_09$cmh_pval)] <- 1
# Rmd parity: the 0.9 pass has NO `filter(!grepl(44, id))` -- that MSI
# GSP drop appears only in the 0.34 pass.

step(sprintf("both_tests_09: %d GSP-gene pairs past Bonferroni.",
             nrow(both_tests_09)))

# ---- 9. n_mut_top25 --------------------------------------------------
step("Computing n_mut_top25 ...")
both_tests_09 <- both_tests_09 %>%
  mutate(gsp_num = as.numeric(str_split_i(id, '-', 1)),
         gene_name = str_split_i(id, '-', 2)) %>%
  rowwise() %>%
  mutate(
    n_mut_top25 = {
      top25_ids <- rownames(norm_L)[order(norm_L[, gsp_num], decreasing = TRUE)[1:25]]
      top25_ids <- top25_ids[top25_ids %in% rownames(mut_status_am)]
      if (gene_name %in% colnames(mut_status_am))
        sum(mut_status_am[top25_ids, gene_name])
      else 0L
    }
  ) %>%
  ungroup()

# ---- 10. persist -----------------------------------------------------
# The plot draws `both_tests_09 %>% filter(n_mut_top25 > 2)`; that
# filter is left to compose.R so the cache keeps the full table.
write_tsv(both_tests_09, OUT_TESTS)
write_tsv(data.frame(
  quantity  = c("bonf_thresh_09", "neg_log10_bonf_09",
                "n_gsps", "n_mut_cols_incl_meta"),
  value     = c(bonf_thresh_09, neg_log10_bonf_09,
                50, ncol(mut_status_am))
), OUT_THRESH)

step(sprintf("Wrote %s (%d rows)", OUT_TESTS, nrow(both_tests_09)))
step(sprintf("Wrote %s", OUT_THRESH))
step("Done.")
