# Top-N alt-splice-site PSI matrices and the SVM "phenotype axis" score per
# TCGA sample, for SF3B1 (GSP #10), NKAP (GSP #13) and YLPM1 (GSP #30).
#
# Ported 2026-10-06 from tcga/genome_wide/flashier/diff_splicing/
# {sf3b1,nkap,ylpm1}/top500_alt3_fc_vs_loading.Rmd (chunks "shared inputs",
# "load <comp> comparison", "assign splice types", "pick top ...", "build
# PSI matrix over selected junctions", "phenotype axis svm setup / fit",
# "phenotype axis plots" -> phen_df), so quetzal_code builds every table
# figures 2f, S13, S17 and S18 read. Logic is the Rmds'; the three copies
# differed only in gene, comparison, splice types and loading table.
#
# One deliberate fix: the YLPM1 Rmd wrote ylpm1_top500_alt3_5_fc5_psi.tsv
# but trained the SVM on ylpm1_top500_alt3_fc5_psi.tsv, a stale April file
# (same 467 junctions then, so no difference in the published numbers). Here
# the SVM always uses the matrix built in the same run.
#
# Drivers (one per figure, each writes its own cache/):
#   figure2/scripts/build_cache_sf3b1_phenotype.R  (30v90, alt 3'ss, k10)
#   figure3/scripts/build_cache_nkap_phenotype.R   (20v60, alt 3'ss, k13)
#   figure4/scripts/build_cache_ylpm1_phenotype.R  (10v30, alt 3'ss + 5'ss, k30)
.libPaths('/scratch/midway2/qhauck/conda_env/rstudio-server/lib/R/library/')
if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
suppressPackageStartupMessages({
  library(yaml); library(data.table); library(tidyverse); library(LiblineaR)
})
select <- dplyr::select; filter <- dplyr::filter; rename <- dplyr::rename

PA_FIG_ROOT <- file.path(QROOT, "figures")

# ---- "shared inputs" -----------------------------------------------------
pa_annotation <- function() {
  intron_df <- fread(file.path(QROOT, "leafcutter_analyses/diff_splicing/gencode_v46_introns.tsv"), data.table = FALSE)
  mane_transcripts <- fread(file.path(QROOT, "inputs/MANE.GRCh38.v1.4.summary.txt"), data.table = FALSE)
  list(intron_to_gene = intron_df %>% select(start, end, strand, gene_name) %>% unique(),
       mane_introns = intron_df %>%
         filter(transcript_id %in% mane_transcripts$Ensembl_nuc) %>%
         select(start, end, strand, gene_name, seqnames) %>% unique())
}

# ---- "load <comp> comparison" --------------------------------------------
pa_load_comparison <- function(results_dir, A) {
  leaf_sig <- read_tsv(file.path(results_dir, "leafcutter_ds_cluster_significance.txt"),
                       show_col_types = FALSE)
  leaf_effect_size <- read_tsv(file.path(results_dir, "leafcutter_ds_effect_sizes.txt"),
                               show_col_types = FALSE)
  leaf_to_gene_name <- leaf_effect_size %>%
    separate(intron, into = c("chrom", "start", "end", "cluster_id"), sep = ":") %>%
    mutate(start = as.numeric(start) + 1, end = as.numeric(end) - 1)
  leaf_to_gene_name$strand <- str_split_i(leaf_to_gene_name$cluster_id, "_", 3)
  leaf_to_gene_name <- leaf_to_gene_name %>%
    left_join(A$intron_to_gene, by = c("start", "end", "strand"))
  cluster_to_gene_name <- leaf_to_gene_name %>%
    group_by(cluster_id) %>%
    mutate(gene_name = {
      x <- gene_name[!is.na(gene_name)]
      if (length(x) == 0) NA_character_ else names(which.max(table(x)))
    }) %>% ungroup()
  leaf_sig <- leaf_sig %>% separate(cluster, into = c("chrom", "cluster_id"), sep = ":")
  cluster_to_gene_name %>%
    left_join(leaf_sig, by = c("chrom", "cluster_id")) %>%
    filter(p.adjust <= 0.05) %>% unique()
}

# ---- "assign splice types" (dpsi > 0, FC > 5, MANE-intron classes) --------
pa_assign_splice_types <- function(full_leaf_sig, A, dpsi_thresh = 0, fc_thresh = 5) {
  mane_introns <- A$mane_introns
  high_psi_hits <- full_leaf_sig %>% filter(deltapsi > dpsi_thresh, test / control > fc_thresh)
  make_id <- function(df) paste0(df$gene_name, ":", df$start, "-", df$end)
  top_junc_df <- high_psi_hits %>%
    mutate(start = as.numeric(start), end = as.numeric(end)) %>%
    filter(gene_name %in% mane_introns$gene_name)
  canonical <- top_junc_df %>%
    inner_join(mane_introns %>% select(gene_name, strand, start, end) %>%
                 rename(mane_start = start, mane_end = end),
               by = c("gene_name", "strand", "start" = "mane_start", "end" = "mane_end")) %>%
    mutate(type = "canonical", mane_intron_start = NA_integer_, mane_intron_end = NA_integer_)
  skipped_exon <- top_junc_df %>%
    inner_join(mane_introns %>% select(gene_name, strand, start, end) %>%
                 rename(up_start = start, up_end = end),
               by = c("gene_name", "strand", "start" = "up_start")) %>%
    inner_join(mane_introns %>% select(gene_name, strand, start, end) %>%
                 rename(down_start = start, down_end = end),
               by = c("gene_name", "strand", "end" = "down_end")) %>%
    filter(up_end < down_start) %>% mutate(type = "skipped exon")
  assigned_ids <- c(make_id(canonical), make_id(skipped_exon))
  start_matches <- top_junc_df %>% filter(!make_id(.) %in% assigned_ids) %>%
    inner_join(mane_introns %>% select(gene_name, strand, start, end) %>%
                 rename(mane_start = start, mane_end = end),
               by = c("gene_name", "strand", "start" = "mane_start")) %>%
    mutate(type = if_else(strand == "+", "alt 3'ss", "alt 5'ss"),
           mane_intron_start = start, mane_intron_end = mane_end)
  end_matches <- top_junc_df %>% filter(!make_id(.) %in% assigned_ids) %>%
    inner_join(mane_introns %>% select(gene_name, strand, start, end) %>%
                 rename(mane_start = start, mane_end = end),
               by = c("gene_name", "strand", "end" = "mane_end")) %>%
    mutate(type = if_else(strand == "+", "alt 5'ss", "alt 3'ss"),
           mane_intron_start = mane_start, mane_intron_end = end)
  annotated <- bind_rows(canonical, skipped_exon, start_matches, end_matches) %>%
    mutate(junc_id = make_id(.)) %>%
    select(junc_id, type, mane_intron_start, mane_intron_end)
  top_junc_df %>%
    mutate(junc_id = make_id(.)) %>%
    left_join(annotated, by = "junc_id") %>%
    mutate(type = if_else(is.na(type), "other", type)) %>%
    group_by(gene_name, start, end) %>% filter(n() < 2) %>% ungroup()
}

# ---- "pick top ..." + "build PSI matrix over selected junctions" ----------
pa_top_junctions <- function(df_final, types, n = 500) {
  df_final %>% filter(type %in% types) %>% slice_max(deltapsi, n = n) %>%
    mutate(start = as.integer(start), end = as.integer(end))
}
pa_psi_matrix <- function(top) {
  psi_base <- file.path(file.path(QROOT, "leafcutter_analyses/psi_files"), "tcga")
  result_list <- list()
  suppressWarnings(for (gene in unique(top$gene_name)) {
    if (grepl("ENSG", gene)) next
    chrs <- unique(top$chrom[top$gene_name == gene])
    chr <- chrs[grepl("chr", chrs)][1]
    psi_path <- file.path(psi_base, chr, gene, "psi_table.tsv")
    if (!file.exists(psi_path) || file.size(psi_path) < 100) next
    header <- names(fread(psi_path, nrows = 0))
    col_coords <- sub("^[^:]+:", "", header)
    wanted_se <- paste0(top$start[top$gene_name == gene], "-", top$end[top$gene_name == gene])
    col_idx <- which(col_coords %in% wanted_se)
    if (length(col_idx) == 0) next
    matched <- col_coords[col_idx]
    col_idx <- col_idx[!duplicated(matched)]; matched <- matched[!duplicated(matched)]
    dt <- fread(psi_path, select = c(1L, col_idx))
    setnames(dt, c("rail_id", paste0(gene, ":", matched)))
    result_list[[gene]] <- dt
  })
  Reduce(function(x, y) merge(x, y, by = "rail_id", all = TRUE), result_list)
}

# ---- "phenotype axis svm setup / fit / plots" -> phen_df ------------------
pa_extract_scores <- function(fit, X) {
  dv <- predict(fit, X, decisionValues = TRUE)$decisionValues
  useful_col <- which.max(colSums(abs(dv)))
  s <- dv[, useful_col]
  if (colnames(dv)[useful_col] != "P") s <- -s
  s
}
pa_phenotype_axis <- function(psi_mat, info, kcol, P_thresh = 0.1, N_thresh = 1e-4) {
  psi_mat <- as.data.frame(psi_mat)
  shared <- psi_mat %>% inner_join(info %>% select(rail_id, all_of(kcol)), by = "rail_id")
  X_raw <- as.matrix(shared[, setdiff(colnames(shared), c("rail_id", kcol))])
  X_raw[is.na(X_raw)] <- 0
  k_vec <- shared[[kcol]]
  eps <- 1e-4
  psi_clip  <- pmin(pmax(X_raw, 0), 1)
  psi_logit <- log((psi_clip + eps) / (1 - psi_clip + eps))
  psi_scaled <- scale(psi_logit)
  keep <- apply(psi_scaled, 2, function(v) all(is.finite(v)))
  psi_scaled <- psi_scaled[, keep, drop = FALSE]
  group <- case_when(k_vec >= P_thresh ~ "P", k_vec <= N_thresh ~ "N", TRUE ~ "intermediate")
  train_idx <- which(group %in% c("P", "N"))
  y_train <- factor(group[train_idx], levels = c("N", "P"))
  n_p <- sum(y_train == "P"); n_n <- sum(y_train == "N")
  wi <- c("N" = (n_p + n_n) / (2 * n_n), "P" = (n_p + n_n) / (2 * n_p))
  fp_candidates <- 0:n_p
  target_fp <- fp_candidates[which.min(abs(fp_candidates / (fp_candidates + n_p) - 0.10))]
  score_at_cost <- function(co) {
    fit_co <- LiblineaR(data = psi_scaled[train_idx, , drop = FALSE], target = y_train,
                        type = 2, cost = co, bias = 1, wi = wi)
    tibble(group = factor(group[train_idx], levels = c("N", "P")),
           score = as.numeric(pa_extract_scores(fit_co, psi_scaled[train_idx, , drop = FALSE])),
           cost = co)
  }
  fp_from_scores <- function(sc) sc %>% group_by(cost) %>%
    summarize(FP = sum(group == "N" & score > 0), .groups = "drop") %>% arrange(cost)
  cost_grid <- c(1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3); max_iter <- 8
  all_fp_df <- fp_from_scores(bind_rows(lapply(cost_grid, score_at_cost)))
  exact_hit <- any(all_fp_df$FP == target_fp, na.rm = TRUE); iter <- 1
  while (!exact_hit && iter < max_iter) {
    below_idx <- which(all_fp_df$FP <= target_fp)
    if (length(below_idx) == 0 || below_idx[1] == 1) break
    i <- below_idx[1]; lo <- all_fp_df$cost[i - 1]; hi <- all_fp_df$cost[i]
    if (hi / lo < 1.01) break
    iter <- iter + 1
    fine <- 10^seq(log10(lo), log10(hi), length.out = 11)[2:10]
    fine_fp <- fp_from_scores(bind_rows(lapply(fine, score_at_cost)))
    all_fp_df <- bind_rows(all_fp_df, fine_fp) %>% distinct(cost, .keep_all = TRUE) %>% arrange(cost)
    exact_hit <- any(fine_fp$FP == target_fp, na.rm = TRUE)
  }
  optimal_c <- if (exact_hit) min(all_fp_df$cost[all_fp_df$FP == target_fp]) else {
    el <- all_fp_df %>% filter(FP <= target_fp)
    if (nrow(el) == 0) all_fp_df$cost[which.min(all_fp_df$FP)] else min(el$cost)
  }
  fit <- LiblineaR(data = psi_scaled[train_idx, , drop = FALSE], target = y_train,
                   type = 2, cost = optimal_c, bias = 1, wi = wi)
  cat(sprintf("[phenotype_axis] %s: %d samples x %d junctions; P=%d N=%d; target FP %d; optimal_c %g\n",
              kcol, nrow(psi_scaled), ncol(psi_scaled), n_p, n_n, target_fp, optimal_c))
  tibble(rail_id = shared$rail_id, !!kcol := k_vec,
         group = factor(group, levels = c("N", "intermediate", "P")),
         score = as.numeric(pa_extract_scores(fit, psi_scaled)))
}

# One gene end to end. Writes <cache_dir>/<psi_name>.tsv and
# <cache_dir>/<gene>_phenotype_axis_scores.tsv.
pa_build <- function(gene, comparison, types, info, kcol, psi_name, cache_dir,
                     results_root = file.path(dirname(file.path(QROOT, "leafcutter_analyses/diff_splicing/gencode_v46_introns.tsv")), "results")) {
  A <- pa_annotation()
  fls <- pa_load_comparison(file.path(results_root, paste0(gene, "_", comparison)), A)
  df_final <- pa_assign_splice_types(fls, A)
  top <- pa_top_junctions(df_final, types)
  cat(sprintf("[phenotype_axis] %s %s: %d significant junctions -> %d selected (%s)\n",
              gene, comparison, nrow(fls), nrow(top), paste(types, collapse = " + ")))
  psi <- pa_psi_matrix(top)
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  fwrite(psi, file.path(cache_dir, paste0(psi_name, ".tsv")), sep = "\t")
  phen <- pa_phenotype_axis(psi, info, kcol)
  fwrite(phen, file.path(cache_dir, paste0(gene, "_phenotype_axis_scores.tsv")), sep = "\t")
  cat(sprintf("[phenotype_axis] wrote %s/%s.tsv (%d x %d) and %s_phenotype_axis_scores.tsv\n",
              cache_dir, psi_name, nrow(psi), ncol(psi) - 1, gene))
  invisible(list(psi = psi, phen = phen))
}
