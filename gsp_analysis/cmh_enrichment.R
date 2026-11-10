# GSP x mutated-gene enrichment among the top 25 samples of each GSP, Cochran-Mantel-Haenszel
# stratified by TCGA cohort, AlphaMissense > 0.9 (or unscored) mutations (Fig. 2d).
# Extracted from flash_fact_char.Rmd (chunks "loading re-used files", "correlating flashier
# factors with mutation status", "AlphaMissense-filtered mutation status ...", "setting up
# covariates", the GSP-alignment / top1-top2 ratio part of "pearson correlation heatmap", and
# the CMH half of "third pass: AM threshold 0.9"), code unchanged apart from paths.
# Writes cmh_mutation_enrichment_am0_9.tsv into OUT (default: this directory).
#   Rscript cmh_enrichment.R [OUT]       (run from gsp_analysis/)
suppressPackageStartupMessages({ library(data.table); library(tidyverse); library(flashier) })
args <- commandArgs(TRUE); OUT <- if (length(args)) args[1] else "."

out <- readRDS('../quetzal_runs/genome_wide/softimpute_flash_300_qc_filtered.RDS')
tcga_metadata <- fread('../inputs/tcga_v2_samples.tsv', data.table = F)
ldf <- flashier::ldf(out, type = 'i')
norm_L <- ldf$L

am_mut_data <- fread('../mutation_processing/merged_grch38_with_am.tsv', data.table = FALSE)
tcga_metadata <- tcga_metadata %>%
  mutate(sampleId = sapply(strsplit(tcga_barcode, "-"), function(x) {
    x[4] <- sub("[A-Za-z]+$", "", x[4])   # drop trailing letters from 4th piece
    paste(x[1:4], collapse = "-")
  }))

build_mut_status_am <- function(thresh) {
  m <- am_mut_data %>%
    filter(is.na(am_pathogenicity) | am_pathogenicity > thresh) %>%
    dplyr::select(hugoGeneSymbol, sampleId) %>%
    mutate(value = 1) %>%
    unique() %>%
    pivot_wider(names_from = hugoGeneSymbol, values_from = value, values_fill = 0) %>%
    right_join(dplyr::select(tcga_metadata, c(rail_id, sampleId, gdc_cases.project.project_id)),
               by = "sampleId") %>%
    dplyr::select(-sampleId) %>%
    column_to_rownames('rail_id') %>%
    na.omit()
  meta_col <- "gdc_cases.project.project_id"
  gene_cols <- setdiff(colnames(m), meta_col)
  keep_genes <- gene_cols[colSums(m[, gene_cols, drop = FALSE]) > 0]
  m[, c(keep_genes, meta_col), drop = FALSE]
}
mut_status_am_09 <- build_mut_status_am(0.9)

pre_model <- readRDS('../support_tables/standardized_covariates.RDS') %>%
  mutate(rail_id = as.numeric(rail_id))

# GSPs not driven by a single sample: top1/top2 loading ratio <= 5 over the samples that
# have metadata and covariates (L_aligned in flash_fact_char.Rmd).
covar_df <- tcga_metadata %>%
  filter(rail_id %in% as.numeric(rownames(norm_L))) %>%
  mutate(rail_id_chr = as.character(rail_id))
common <- Reduce(intersect, list(rownames(norm_L), covar_df$rail_id_chr, as.character(pre_model$rail_id)))
L_aligned <- norm_L[common, ]
ratio_1v2 <- sapply(seq_len(ncol(L_aligned)), function(k) {
  sorted <- sort(L_aligned[, k], decreasing = TRUE); sorted[1] / sorted[2]
})
non_singleton_gsps <- seq_len(ncol(L_aligned))[as.numeric(ratio_1v2 > 5) == 0]

# ---- third pass, AM threshold 0.9 (CMH half) ----
mut_status_am <- mut_status_am_09
tmb_09 <- rowSums(mut_status_am[, setdiff(colnames(mut_status_am), "gdc_cases.project.project_id"),
                                drop = FALSE], na.rm = TRUE)
keep_samples_09 <- names(tmb_09)[tmb_09 <= 500]
mut_status_am <- mut_status_am[keep_samples_09, , drop = FALSE]
mut_status_am <- mut_status_am[, c(colSums(mut_status_am[, 1:(ncol(mut_status_am) - 1)]) > 10, T)]
mut_status_am <- mut_status_am[, c(colSums(mut_status_am[, 1:(ncol(mut_status_am) - 1)]) < 1000, T)]

all_common      <- intersect(rownames(norm_L), rownames(mut_status_am))
mut_genes_only  <- as.matrix(mut_status_am[all_common,
                                           setdiff(colnames(mut_status_am), "gdc_cases.project.project_id")])
gene_counts <- colSums(mut_genes_only)
keep_genes  <- names(gene_counts[gene_counts > 10 & gene_counts < 1000])
mut_mat     <- mut_genes_only[, keep_genes]
n_genes     <- ncol(mut_mat)
cancer_type <- mut_status_am[all_common, "gdc_cases.project.project_id"]

strata    <- split(seq_len(length(all_common)), cancer_type)
strata    <- strata[sapply(strata, length) >= 2]
strata_n  <- sapply(strata, length)
strata_m1 <- do.call(rbind, lapply(strata, function(idx) colSums(mut_mat[idx, , drop = FALSE])))

top_n <- 25
cmh_results <- vector("list", length(non_singleton_gsps))
for (ki in seq_along(non_singleton_gsps)) {
  k       <- non_singleton_gsps[ki]
  l_col   <- norm_L[all_common, k]
  is_top  <- integer(length(all_common))
  is_top[order(l_col, decreasing = TRUE)[1:top_n]] <- 1L

  sum_a <- sum_E <- sum_V <- sum_ad_over_n <- sum_bc_over_n <- numeric(n_genes)
  n_valid_per_gene <- integer(n_genes)
  for (s_name in names(strata)) {
    idx <- strata[[s_name]]
    ns  <- strata_n[s_name]
    n1s <- sum(is_top[idx])
    n0s <- ns - n1s
    if (n1s == 0 || n0s == 0) next
    m1 <- strata_m1[s_name, ]
    m0 <- ns - m1
    gene_valid <- m1 > 0 & m0 > 0
    a <- as.numeric(crossprod(is_top[idx], mut_mat[idx, , drop = FALSE]))
    b <- n1s - a; c <- m1 - a; d <- n0s - c
    sum_a         <- sum_a + a
    sum_E         <- sum_E + ifelse(gene_valid, n1s * m1 / ns, 0)
    sum_V         <- sum_V + ifelse(gene_valid, n1s * n0s * m1 * m0 / (ns^2 * (ns - 1)), 0)
    sum_ad_over_n <- sum_ad_over_n + ifelse(gene_valid, a * d / ns, 0)
    sum_bc_over_n <- sum_bc_over_n + ifelse(gene_valid, b * c / ns, 0)
    n_valid_per_gene <- n_valid_per_gene + as.integer(gene_valid)
  }
  testable <- n_valid_per_gene >= 1 & sum_a >= 5 & sum_V > 0
  hits     <- which(testable)
  if (length(hits) == 0) next
  z     <- (sum_a[hits] - sum_E[hits]) / sqrt(sum_V[hits])
  pval  <- pnorm(z, lower.tail = FALSE)
  or_mh <- sum_ad_over_n[hits] / pmax(sum_bc_over_n[hits], 1e-10)
  cmh_results[[ki]] <- data.frame(
    gsp = k, gene = keep_genes[hits], n_top_mutated = as.integer(sum_a[hits]),
    n_strata = n_valid_per_gene[hits], odds_ratio = round(or_mh, 2), pvalue = pval)
}
n_total_tests <- length(non_singleton_gsps) * n_genes
cmh_enrich_09 <- bind_rows(cmh_results) %>%
  mutate(padj = pmin(pvalue * n_total_tests, 1)) %>%
  arrange(pvalue)
fwrite(cmh_enrich_09, file.path(OUT, "cmh_mutation_enrichment_am0_9.tsv"), sep = "\t")
