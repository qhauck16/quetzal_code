#!/usr/bin/env Rscript
# Figure 1b: percent-recovery vs |deltaPSI| for 3 pairwise TCGA
# comparisons (KICH vs LUAD, KICH vs TGCT, TGCT vs LUAD), FC >= 1.5
# threshold only.
#
# Extracted verbatim from tcga/model_eval/quetzal_benchmarking.Rmd
# (setup + `run_comparison` + FC helper + the three "FC 1.5" chunks).
# Reads inputs through ../../../paths.yaml so it is portable.
#
# Output: out/{KICH_vs_LUAD,KICH_vs_TGCT,TGCT_vs_LUAD}_recovery_fc1.5.pdf

if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
suppressPackageStartupMessages({
  library(tidyverse)
  library(rtracklayer)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(data.table)
  library(MatchIt)
})
select <- dplyr::select
rename <- dplyr::rename

# ------------------------------------------------------------------
# Path resolution -- read paths.yaml from the repo root.
# ------------------------------------------------------------------
this_dir <- dirname(sys.frame(1)$ofile %||% getwd())
`%||%` <- function(a, b) if (is.null(a) || !nzchar(a)) b else a
repo_root <- normalizePath(file.path(this_dir, "..", "..", ".."))

out_dir <- file.path(this_dir, "out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cache_dir <- file.path(this_dir, "recovery_cache")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------
# Intron reference (needed to name junctions -> genes).
# ------------------------------------------------------------------
gencode <- readRDS(file.path(QROOT, "inputs/gencode_v46_granges.RDS"))

txdb <- makeTxDbFromGRanges(gencode)
introns_by_tx <- intronsByTranscript(txdb, use.names = TRUE)
introns <- unlist(introns_by_tx, use.names = FALSE)

tx2gene <- unique(
  mcols(gencode)[gencode$type == "transcript",
             c("transcript_id", "gene_id", "gene_name")]
)

introns$transcript_id <- names(introns_by_tx)[
  rep(seq_along(introns_by_tx), elementNROWS(introns_by_tx))
]

mcols(introns) <- cbind(
  mcols(introns),
  tx2gene[match(introns$transcript_id, tx2gene$transcript_id), ]
)

intron_df <- as.data.frame(introns)

intron_to_gene <- intron_df %>%
  dplyr::select(start, end, strand, gene_name) %>%
  unique()

tcga_metadata <- fread(file.path(QROOT, "inputs/tcga_v2_samples.tsv"), data.table = FALSE)

# ------------------------------------------------------------------
# Load LeafCutter DS output for one comparison, filter to significant
# junctions with |deltapsi| > threshold.
# ------------------------------------------------------------------
load_leafcutter <- function(comparison, deltapsi_thresh = 0.1) {
  sig_file <- file.path(file.path(QROOT, "benchmarking/leafcutter"),
                        paste0(comparison, "_ds_cluster_significance.txt"))
  eff_file <- file.path(file.path(QROOT, "benchmarking/leafcutter"),
                        paste0(comparison, "_ds_effect_sizes.txt"))

  leaf_sig <- fread(sig_file, data.table = FALSE) %>%
    filter(status == "Success") %>%
    separate(cluster, into = c("chrom", "cluster_id"), sep = ":") %>%
    select(-chrom)

  leaf_effect_size <- fread(eff_file, data.table = FALSE)

  leaf_to_gene_name <- leaf_effect_size %>%
    separate(intron, into = c("chrom", "start", "end", "cluster_id"), sep = ":") %>%
    mutate(start = as.numeric(start) + 1, end = as.numeric(end) - 1)
  leaf_to_gene_name$strand <- str_split_i(leaf_to_gene_name$cluster_id, "_", 3)

  leaf_to_gene_name <- leaf_to_gene_name %>%
    left_join(intron_to_gene, by = c("start", "end", "strand"))

  cluster_to_gene_name <- leaf_to_gene_name %>%
    group_by(cluster_id) %>%
    mutate(
      gene_name = {
        x <- gene_name[!is.na(gene_name)]
        if (length(x) == 0) NA_character_
        else names(which.max(table(x)))
      }
    ) %>%
    ungroup()

  cluster_to_gene_name %>%
    left_join(leaf_sig, by = "cluster_id") %>%
    filter(p.adjust <= 0.05) %>%
    filter(abs(deltapsi) > deltapsi_thresh) %>%
    unique()
}

# ------------------------------------------------------------------
# factor_loading_summary + run_comparison + run_recovery{,_fc}
# -- lifted verbatim from quetzal_benchmarking.Rmd.
# ------------------------------------------------------------------
factor_loading_summary <- function(leaf_sig_df, chr, cohorts) {
  chr_hits <- leaf_sig_df %>%
    filter(chrom %in% chr) %>%
    filter(!is.na(gene_name)) %>%
    arrange(desc(abs(deltapsi))) %>%
    mutate(junc_id = paste0(start, "-", end)) %>%
    distinct(junc_id, .keep_all = TRUE)

  relevant_rail_ids <- tcga_metadata %>%
    filter(gdc_cases.project.project_id %in% cohorts) %>%
    pull(rail_id)

  gene_dir <- file.path(QROOT, "quetzal_runs/gene_level")
  res_cache <- list()
  results   <- list()
  skipped   <- list(no_gene_dir = c(), junc_not_in_res = c(),
                    too_few_samples = c(), no_factor_match = c(),
                    matched = c())

  for (i in seq_len(nrow(chr_hits))) {
    if (i %% 100 == 0) print(i)
    gene <- chr_hits$gene_name[i]
    junc <- chr_hits$junc_id[i]
    junc_chr <- chr_hits$chrom[i]

    if (is.null(res_cache[[gene]])) {
      res_path <- file.path(gene_dir, junc_chr, gene, "res.RDS")
      if (!file.exists(res_path) || file.size(res_path) < 100) {
        skipped$no_gene_dir <- c(skipped$no_gene_dir, paste0(i, ":", gene))
        next
      }
      res_cache[[gene]] <- readRDS(res_path)
    }
    res <- res_cache[[gene]]

    if (!junc %in% rownames(res$F)) {
      skipped$junc_not_in_res <- c(skipped$junc_not_in_res, paste0(i, ":", gene))
      next
    }

    L_ids <- as.numeric(rownames(res$L)[rownames(res$L) %in% relevant_rail_ids])
    L_groups <- tcga_metadata$gdc_cases.project.project_id[
      match(L_ids, tcga_metadata$rail_id)
    ]

    if (any(table(L_groups)[cohorts] < 10 |
            any(!cohorts %in% names(table(L_groups)[cohorts])),
            na.rm = TRUE)) {
      skipped$too_few_samples <- c(skipped$too_few_samples, paste0(i, ":", gene))
      next
    }

    found_factor <- FALSE
    for (k in seq_len(ncol(res$F))) {
      junc_val <- res$F[junc, k]
      max_val <- max(res$F[, k])
      if (max_val == 0 || junc_val < max_val / 2) next

      found_factor <- TRUE
      L_sub <- res$L[rownames(res$L) %in% relevant_rail_ids, k, drop = FALSE]
      groups <- tcga_metadata$gdc_cases.project.project_id[
        match(as.numeric(rownames(L_sub)), tcga_metadata$rail_id)
      ]

      col_name <- paste0(i, ":", gene, "-k", k)
      for (cohort in cohorts) {
        results[[col_name]][[cohort]] <- mean(L_sub[groups == cohort, 1])
      }
    }
    if (found_factor) {
      skipped$matched <- c(skipped$matched, paste0(i, ":", gene))
    } else {
      skipped$no_factor_match <- c(skipped$no_factor_match, paste0(i, ":", gene))
    }
  }

  if (length(results) == 0) {
    return(list(summary = data.frame(row.names = cohorts), skipped = skipped))
  }
  summary_df <- do.call(cbind, lapply(results, function(x) {
    sapply(cohorts, function(c) x[[c]])
  })) %>% as.data.frame()

  list(summary = summary_df, skipped = skipped)
}

# Background PSI/RPM tables for the matched-sampling step.
aggregated_psi_medians <- fread(file.path(QROOT, "support_tables/aggregated_psi_medians.tsv"), data.table = FALSE)
aggregated_rpm_medians <- fread(file.path(QROOT, "support_tables/aggregated_junction_rpms.tsv"), data.table = FALSE)
all_junction_and_psi   <- fread(file.path(QROOT, "leafcutter_analyses/diff_splicing/unified_results/shared/all_junction_rpm_and_psi.tsv"))

gene_to_chrom <- as.data.frame(mcols(gencode)) %>%
  filter(type == "gene") %>%
  dplyr::select(gene_name) %>%
  mutate(chrom = as.character(seqnames(gencode[gencode$type == "gene"]))) %>%
  distinct(gene_name, .keep_all = TRUE)

chroms <- c(paste0("chr", 1:22), "chrX")

run_recovery <- function(sp, hits_ref, sk, psi_vals) {
  storage <- matrix(0, nrow = length(psi_vals), ncol = 2)
  storage[, 1] <- psi_vals
  for (i in seq_along(psi_vals)) {
    psi <- psi_vals[i]
    cutoff_number <- sum(abs(hits_ref$deltapsi) > psi) + 1
    col_ranks <- as.numeric(str_split_i(colnames(sp), ":", 1))
    in_cutoff <- col_ranks < cutoff_number
    candidates <- sp[, in_cutoff, drop = FALSE]
    cand_signs <- sign(hits_ref$deltapsi[col_ranks[in_cutoff]])
    successes <- sum(
      (cand_signs > 0 & candidates[2, ] > candidates[1, ]) |
      (cand_signs < 0 & candidates[1, ] > candidates[2, ])
    )
    skip_ranks <- as.numeric(str_split_i(sk$no_factor_match, ":", 1))
    failures <- ncol(candidates) - successes + sum(skip_ranks < cutoff_number)
    storage[i, 2] <- successes / (successes + failures)
  }
  storage
}

run_comparison <- function(leaf_sig, comparison, cohorts) {
  match_test_file  <- file.path(cache_dir, paste0(comparison, "_matched_test.tsv"))
  match_bg_file    <- file.path(cache_dir, paste0(comparison, "_matched_background.tsv"))
  fls_test_file    <- file.path(cache_dir, paste0(comparison, "_fls_test.RDS"))
  fls_matched_file <- file.path(cache_dir, paste0(comparison, "_fls_matched.RDS"))

  if (file.exists(match_test_file) && file.exists(match_bg_file)) {
    matched_test       <- fread(match_test_file, data.table = FALSE)
    matched_background <- fread(match_bg_file, data.table = FALSE)
  } else {
    chr_hits_ref_full <- leaf_sig %>%
      filter(!is.na(gene_name)) %>%
      arrange(desc(abs(deltapsi))) %>%
      mutate(junc_id = paste0(start, "-", end)) %>%
      distinct(junc_id, .keep_all = TRUE)

    test_juncs <- chr_hits_ref_full %>%
      mutate(junc_coord = junc_id) %>%
      left_join(all_junction_and_psi, by = c("gene_name", "start", "end")) %>%
      filter(!is.na(median_psi), !is.na(median_rpm)) %>%
      mutate(condition = "test") %>%
      dplyr::select(gene_name, start, end, median_psi, median_rpm,
                    condition, deltapsi, junc_coord)

    background_juncs <- all_junction_and_psi %>%
      anti_join(test_juncs, by = c("gene_name", "start", "end")) %>%
      mutate(condition = "background", deltapsi = NA,
             junc_coord = paste0(start, "-", end))

    match_input <- bind_rows(test_juncs, background_juncs) %>%
      mutate(condition = factor(condition,
                                levels = c("background", "test")))

    m.out <- matchit(
      condition ~ median_psi + median_rpm,
      data = match_input,
      method = "nearest",
      distance = "mahalanobis",
      ratio = 1
    )
    matched_data <- match.data(m.out)

    matched_background <- matched_data %>%
      filter(condition == "background") %>%
      arrange(subclass)
    matched_test <- matched_data %>%
      filter(condition == "test") %>%
      arrange(subclass)

    fwrite(matched_test,       match_test_file, sep = "\t")
    fwrite(matched_background, match_bg_file,    sep = "\t")
  }

  matched_test_for_fls <- matched_test %>%
    left_join(gene_to_chrom, by = "gene_name") %>%
    filter(!is.na(chrom))

  if (file.exists(fls_test_file)) {
    all_chr_result <- readRDS(fls_test_file)
  } else {
    all_chr_result <- factor_loading_summary(matched_test_for_fls, chroms, cohorts)
    saveRDS(all_chr_result, fls_test_file)
  }

  results <- all_chr_result$summary
  junc_ids <- str_split_i(colnames(results), "-", 1)
  summed_props <- matrix(0, nrow = 2, ncol = length(unique(junc_ids)))
  colnames(summed_props) <- unique(junc_ids)
  for (junc_id in unique(junc_ids)) {
    duplicates <- which(junc_ids == junc_id)
    summed_props[, junc_id] <- rowSums(results[duplicates])
  }
  skipped <- all_chr_result$skipped

  chr_hits_ref <- matched_test_for_fls %>%
    arrange(desc(abs(deltapsi))) %>%
    mutate(junc_id = paste0(start, "-", end)) %>%
    distinct(junc_id, .keep_all = TRUE)

  matched_bg_for_fls <- matched_background %>%
    dplyr::select(-deltapsi) %>%
    left_join(matched_test %>% dplyr::select(subclass, deltapsi),
              by = "subclass") %>%
    left_join(gene_to_chrom, by = "gene_name") %>%
    filter(!is.na(chrom)) %>%
    arrange(desc(abs(deltapsi)))

  if (file.exists(fls_matched_file)) {
    matched_chr_result <- readRDS(fls_matched_file)
  } else {
    matched_chr_result <- factor_loading_summary(matched_bg_for_fls, chroms, cohorts)
    saveRDS(matched_chr_result, fls_matched_file)
  }

  matched_results <- matched_chr_result$summary
  matched_junc_ids <- str_split_i(colnames(matched_results), "-", 1)
  matched_summed_props <- matrix(0, nrow = 2,
                                  ncol = length(unique(matched_junc_ids)))
  colnames(matched_summed_props) <- unique(matched_junc_ids)
  for (junc_id in unique(matched_junc_ids)) {
    duplicates <- which(matched_junc_ids == junc_id)
    matched_summed_props[, junc_id] <- rowSums(matched_results[duplicates])
  }
  matched_skipped <- matched_chr_result$skipped

  matched_chr_hits_ref <- matched_bg_for_fls %>%
    mutate(junc_id = paste0(start, "-", end)) %>%
    distinct(junc_id, .keep_all = TRUE)

  psi_vals <- c(50, 47.5, 45, 42.5, 40, 39:10) / 100
  storage         <- run_recovery(summed_props,         chr_hits_ref,         skipped,         psi_vals)
  matched_storage <- run_recovery(matched_summed_props, matched_chr_hits_ref, matched_skipped, psi_vals)

  recovery_df <- data.frame(
    psi_cutoff = psi_vals,
    leafcutter = storage[, 2],
    matched    = matched_storage[, 2]
  ) %>%
    pivot_longer(c(leafcutter, matched), names_to = "set",
                 values_to = "recovery")

  list(
    recovery_df          = recovery_df,
    storage              = storage,
    matched_storage      = matched_storage,
    test_summed_props    = summed_props,
    test_hits_ref        = chr_hits_ref,
    test_skipped         = skipped,
    matched_summed_props = matched_summed_props,
    matched_hits_ref     = matched_chr_hits_ref,
    matched_skipped      = matched_skipped
  )
}

# FC-thresholded recovery (a candidate only counts if its cohort-
# proportion FC exceeds fc_thresh AND its sign matches deltaPSI).
run_recovery_fc <- function(sp, hits_ref, sk, psi_vals, fc_thresh) {
  storage <- matrix(0, nrow = length(psi_vals), ncol = 3)
  storage[, 1] <- psi_vals
  storage[, 3] <- fc_thresh
  for (i in seq_along(psi_vals)) {
    psi <- psi_vals[i]
    cutoff_number <- sum(abs(hits_ref$deltapsi) > psi) + 1
    col_ranks <- as.numeric(str_split_i(colnames(sp), ":", 1))
    in_cutoff <- col_ranks < cutoff_number
    candidates <- sp[, in_cutoff, drop = FALSE]
    cand_signs <- sign(hits_ref$deltapsi[col_ranks[in_cutoff]])

    c1 <- candidates[1, ]; c2 <- candidates[2, ]
    hi <- pmax(c1, c2); lo <- pmin(c1, c2)
    ratio <- ifelse(lo > 0, hi / lo, ifelse(hi > 0, Inf, 1))
    above_fc <- ratio >= fc_thresh

    direction_correct <- (cand_signs > 0 & c2 > c1) |
                          (cand_signs < 0 & c1 > c2)
    successes <- sum(above_fc & direction_correct)

    skip_ranks <- as.numeric(str_split_i(sk$no_factor_match, ":", 1))
    failures <- ncol(candidates) - successes + sum(skip_ranks < cutoff_number)
    storage[i, 2] <- successes / (successes + failures)
  }
  storage
}

recovery_df_fc1.5 <- function(pkg) {
  psi_vals <- pkg$storage[, 1]
  test_storage <- run_recovery_fc(pkg$test_summed_props,
                                    pkg$test_hits_ref,
                                    pkg$test_skipped,
                                    psi_vals, 1.5)
  matched_storage <- run_recovery_fc(pkg$matched_summed_props,
                                      pkg$matched_hits_ref,
                                      pkg$matched_skipped,
                                      psi_vals, 1.5)
  data.frame(
    psi_cutoff = psi_vals,
    leafcutter = test_storage[, 2],
    matched    = matched_storage[, 2]
  ) %>%
    pivot_longer(c(leafcutter, matched), names_to = "set",
                 values_to = "recovery")
}

# ------------------------------------------------------------------
# Run all three comparisons.
# ------------------------------------------------------------------
COMPARISONS <- list(
  KICH_vs_LUAD = c("TCGA-KICH", "TCGA-LUAD"),
  KICH_vs_TGCT = c("TCGA-KICH", "TCGA-TGCT"),
  TGCT_vs_LUAD = c("TCGA-TGCT", "TCGA-LUAD")
)

for (comp_name in names(COMPARISONS)) {
  cat(sprintf("\n== %s ==\n", comp_name))
  leaf_sig <- load_leafcutter(comp_name)
  pkg      <- run_comparison(leaf_sig, comp_name, COMPARISONS[[comp_name]])
  fc_df    <- recovery_df_fc1.5(pkg)

  ref_x <- range(pkg$recovery_df$psi_cutoff)

  p <- ggplot(fc_df, aes(x = psi_cutoff, y = recovery, color = set)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_x_reverse(limits = c(ref_x[2], ref_x[1])) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual(values = c(leafcutter = "#d95f02", matched = "#7570b3")) +
    labs(x = "deltaPSI cutoff", y = "Percent recovery", color = NULL) +
    theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank())

  out_pdf <- file.path(out_dir, sprintf("%s_recovery_fc1.5.pdf", comp_name))
  ggsave(out_pdf, p, width = 6, height = 4)
  cat("Wrote", out_pdf, "\n")
}
