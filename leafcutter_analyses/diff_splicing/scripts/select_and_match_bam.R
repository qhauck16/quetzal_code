#!/usr/bin/env Rscript
# BAM-flow version of select_and_match.R.
#
# Same logic as the snaptron-flow version, but the "available" pool comes
# from download/bam_inventory.tsv (rail_id <-> bam_path) rather than from
# .junc files in a snaptron-derived directory.
#
# Usage:
#   Rscript select_and_match_bam.R \
#     --loadings flashier_loadings.tsv \
#     --factor_num 20 \
#     --control_threshold 1e-4 \
#     --n_test 10 \
#     --control_ratio 3 \
#     --metadata ../../../tcga_v2_samples.tsv \
#     --inventory download/bam_inventory.tsv \
#     --gene nkap \
#     --out_dir results/nkap_10v30 \
#     [--control_mode]
#
# Ranking column is built from --factor_num (so factor_num=20 -> column "k20"),
# which is read from flashier_loadings.tsv (produced by
# scripts/extract_loadings_from_rds.R from the canonical qc_filtered RDS).

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(MatchIt)
})

args <- commandArgs(trailingOnly = TRUE)
parse_arg <- function(flag) {
  idx <- which(args == flag)
  if (length(idx) == 0) stop(paste("Missing argument:", flag))
  args[idx + 1]
}

loadings_file     <- parse_arg("--loadings")
factor_num        <- as.integer(parse_arg("--factor_num"))
factor_col        <- paste0("k", factor_num)
control_threshold <- as.numeric(parse_arg("--control_threshold"))
n_test            <- as.integer(parse_arg("--n_test"))
control_ratio     <- as.integer(parse_arg("--control_ratio"))
metadata_file     <- parse_arg("--metadata")
inventory_file    <- parse_arg("--inventory")
gene              <- parse_arg("--gene")
out_dir           <- parse_arg("--out_dir")
control_mode      <- "--control_mode" %in% args

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tcga_metadata <- fread(metadata_file, data.table = FALSE)
loadings      <- read_tsv(loadings_file, show_col_types = FALSE)
inventory     <- read_tsv(inventory_file, show_col_types = FALSE)

if (!(factor_col %in% colnames(loadings))) {
  stop(sprintf("Column %s not found in %s (available: %s)",
               factor_col, loadings_file,
               paste(colnames(loadings), collapse = ", ")))
}

gene_info <- loadings %>%
  dplyr::select(rail_id, all_of(factor_col))

available_rail_ids <- unique(inventory$rail_id)

covariates_to_match <- c("rail_id", "gdc_cases.project.project_id", "avg_len",
                         "gdc_cases.samples.is_ffpe")

full_info <- tcga_metadata[, covariates_to_match] %>%
  filter(rail_id %in% available_rail_ids) %>%
  filter(!is.na(gdc_cases.project.project_id),
         !is.na(avg_len)) %>%
  mutate(gdc_cases.samples.is_ffpe = ifelse(
    is.na(gdc_cases.samples.is_ffpe), FALSE, gdc_cases.samples.is_ffpe
  )) %>%
  inner_join(gene_info, by = "rail_id")

# Per-project cap: for the control-mode run to be feasible, each project P
# needs (control_ratio + 1) * n_from_P samples drawn from the control-eligible
# pool (factor_col < control_threshold) — n_from_P as pseudo-test plus
# control_ratio * n_from_P as pseudo-control. Cap n_from_P at
# floor(n_P_low_loaders / (control_ratio + 1)) and backfill the test set with
# the next-best loaders from uncapped projects.
proj_col <- "gdc_cases.project.project_id"
project_caps <- full_info %>%
  filter(.data[[factor_col]] < control_threshold) %>%
  count(.data[[proj_col]], name = "n_low_loaders") %>%
  mutate(cap = n_low_loaders %/% (control_ratio + 1)) %>%
  { setNames(.$cap, .[[proj_col]]) }

ranked <- full_info %>% arrange(desc(!!sym(factor_col)))
proj_used <- setNames(integer(length(project_caps)), names(project_caps))
top_ids <- integer(0)
capped_out <- 0L
for (i in seq_len(nrow(ranked))) {
  if (length(top_ids) >= n_test) break
  p <- ranked[[proj_col]][i]
  if (is.na(p) || is.na(project_caps[p])) next
  if (proj_used[p] >= project_caps[p]) {
    capped_out <- capped_out + 1L
    next
  }
  top_ids <- c(top_ids, ranked$rail_id[i])
  proj_used[p] <- proj_used[p] + 1L
}
if (capped_out > 0) {
  cat(sprintf(
    "Per-project cap active: skipped %d higher-loading samples whose cohort hit cap.\n",
    capped_out))
  cat("Per-project caps used (project: used/cap):\n")
  for (p in names(project_caps)) {
    if (proj_used[p] > 0) {
      cat(sprintf("  %s: %d/%d\n", p, proj_used[p], project_caps[p]))
    }
  }
}

full_info$group <- NA_character_
full_info$group[full_info$rail_id %in% top_ids] <- "test"
full_info$group[is.na(full_info$group) &
                  full_info[[factor_col]] < control_threshold] <- "control"
full_info <- full_info %>% filter(!is.na(group))

n_controls_available <- sum(full_info$group == "control")
actual_test <- sum(full_info$group == "test")
if (n_controls_available < actual_test * control_ratio) {
  max_test <- floor(n_controls_available / control_ratio)
  cat(sprintf("WARNING: Only %d controls available, capping test to %d (from %d)\n",
              n_controls_available, max_test, actual_test))
  test_ids <- full_info %>%
    filter(group == "test") %>%
    arrange(desc(!!sym(factor_col))) %>%
    slice_head(n = max_test) %>%
    pull(rail_id)
  full_info$group[full_info$group == "test" & !(full_info$rail_id %in% test_ids)] <- NA
  full_info <- full_info %>% filter(!is.na(group))
}

full_info$group <- factor(full_info$group, levels = c("control", "test"))

m.out <- matchit(
  group ~ avg_len + gdc_cases.samples.is_ffpe,
  data = full_info,
  exact = ~ gdc_cases.project.project_id,
  method = "nearest",
  ratio = control_ratio
)
matched_df <- match.data(m.out)
n_test_final <- sum(matched_df$group == "test")
n_ctrl_final <- sum(matched_df$group == "control")
cat(sprintf("Test: %d, Control: %d\n", n_test_final, n_ctrl_final))

if (control_mode) {
  true_test_ids <- matched_df %>% filter(group == "test") %>% pull(rail_id)

  # Per-project counts in the TRUE test set. The fake test set mirrors the
  # cohort (project) proportions exactly; matchit below then picks fake
  # controls on project + avg_len + FFPE the same way the TRUE pipeline does.
  test_counts <- matched_df %>%
    filter(group == "test") %>%
    count(.data[[proj_col]], name = "n_needed")

  # Same NA filtering / FFPE defaulting as the TRUE pipeline, and exclude the
  # actual TRUE test rail_ids so the control doesn't recycle the mutants.
  # Restricted to samples in the flashier loadings (same universe as the
  # TRUE pipeline's full_info). Until 2026-09-25 this drew from every BAM
  # in the inventory, so rail_ids outside the flashier object -- incl.
  # normal-labelled ones from the patient-level sort bug -- could enter.
  candidate_pool <- tcga_metadata[, covariates_to_match] %>%
    filter(rail_id %in% available_rail_ids,
           rail_id %in% gene_info$rail_id,
           !is.na(gdc_cases.project.project_id),
           !is.na(avg_len)) %>%
    mutate(gdc_cases.samples.is_ffpe = ifelse(
      is.na(gdc_cases.samples.is_ffpe), FALSE, gdc_cases.samples.is_ffpe
    )) %>%
    filter(!(rail_id %in% true_test_ids))

  set.seed(as.integer(paste0(n_test_final, n_ctrl_final)))

  fake_test <- bind_rows(lapply(seq_len(nrow(test_counts)), function(i) {
    p <- test_counts[[proj_col]][i]
    n_needed <- test_counts$n_needed[i]
    avail <- candidate_pool %>% filter(.data[[proj_col]] == p)
    if (nrow(avail) < n_needed) {
      stop(sprintf(
        "Control mode: project %s has %d candidates but needs %d",
        p, nrow(avail), n_needed))
    }
    avail %>% slice_sample(n = n_needed)
  }))

  remaining_pool <- candidate_pool %>%
    filter(!(rail_id %in% fake_test$rail_id))

  matchit_input <- bind_rows(
    fake_test %>% mutate(group = "test"),
    remaining_pool %>% mutate(group = "control")
  )
  matchit_input$group <- factor(matchit_input$group, levels = c("control", "test"))

  m.ctrl <- matchit(
    group ~ avg_len + gdc_cases.samples.is_ffpe,
    data = matchit_input,
    exact = ~ gdc_cases.project.project_id,
    method = "nearest", ratio = control_ratio
  )
  matched_df <- match.data(m.ctrl)

  cat(sprintf("Control mode: matched draw (Test: %d, Control: %d)\n",
              sum(matched_df$group == "test"),
              sum(matched_df$group == "control")))
}

# -- Outputs ----------------------------------------------------------------
write_tsv(data.frame(sample_id = matched_df$rail_id),
          file.path(out_dir, "leaf_final_rail_ids.tsv"))

# pheno labels each sample by its rail_id (suffix .junc to match leafcutter
# convention; the .junc files are produced by the regtools rule)
write.table(
  data.frame(
    sample_id = paste0(matched_df$rail_id, ".junc"),
    condition = matched_df$group
  ),
  file.path(out_dir, "leaf_pheno_data.txt"),
  sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE
)

# juncfiles.txt: absolute paths to the per-sample .junc files.
# Ensure junc_files/ exists so normalizePath returns a real absolute path
# (otherwise it silently returns the relative input when the dir is missing).
dir.create("junc_files", recursive = TRUE, showWarnings = FALSE)
JUNC_DIR <- normalizePath("junc_files")
junc_paths <- file.path(JUNC_DIR, paste0(matched_df$rail_id, ".junc"))
writeLines(junc_paths, file.path(out_dir, "juncfiles.txt"))

# Also write the BAM list this run depends on (Snakefile turns this into
# regtools rule inputs)
bam_lookup <- inventory %>% filter(rail_id %in% matched_df$rail_id) %>%
  dplyr::select(rail_id, home_bam)
write_tsv(bam_lookup, file.path(out_dir, "bams_used.tsv"))

cat(sprintf("Wrote outputs to %s\n", out_dir))
