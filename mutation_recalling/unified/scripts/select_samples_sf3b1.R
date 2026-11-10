
# SF3B1 case+control selection. Same framework as NKAP/YLPM1
# (loading-based) with SF3B1-specific exclusions on top:
#
# Cases: norm_L[, factor_key] > loading_thresh
#        AND no SUGP1 mutation
#        AND no SF3B1 hotspot mutation
#        AND no SF3B1 mutation near a hotspot (positions 620-800)
#        AND no SF3B1 p.R549C
#
# Controls: per-project proportion-matched from samples with
#          norm_L[, factor_key] < control_cutoff and the SAME
#          (SUGP1 / SF3B1-hotspot / near / R549C) exclusions.
#
# Output schema (samples.tsv) matches select_samples.R so
# build_manifests.R / import_downloaded.R consume it unchanged.

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(flashier)
})

args <- commandArgs(trailingOnly = TRUE)
gene            <- args[1]                   # "sf3b1"
factor_key      <- args[2]                   # "k10"
rds_path        <- args[3]
tcga_metadata   <- args[4]
loading_thresh  <- as.numeric(args[5])
control_cutoff  <- as.numeric(args[6])
seed_val        <- as.integer(args[7])
out_dir         <- args[8]
mut_data_paths  <- args[-(1:8)]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- SF3B1 exclusion criteria ----
SF3B1_HOTSPOTS <- c(
  "R625H", "R625C", "K700E",
  "K666T", "K666E", "K666N",
  "G740E", "G740V",
  "K741E", "K741N", "K741Q"
)
NEAR_HOTSPOT_LO <- 620
NEAR_HOTSPOT_HI <- 800
EXTRA_EXCLUDE   <- c("R549C")

# ---- Load RDS, get norm_L[, factor_key] -------------------------
cat("Loading RDS:", rds_path, "\n")
out <- readRDS(rds_path)

ldf <- flashier::ldf(out, type = "i")
L   <- ldf$L
colnames(L) <- paste0("k", seq_len(ncol(L)))
stopifnot(factor_key %in% colnames(L))

sample_ids <- as.integer(rownames(L))
loadings   <- L[, factor_key]

cat("Total samples in RDS:", length(sample_ids), "\n")
cat("Cases (loading >", loading_thresh, "):",
    sum(loadings > loading_thresh, na.rm = TRUE), "\n")
cat("Eligible controls (loading <", control_cutoff, "):",
    sum(loadings < control_cutoff, na.rm = TRUE), "\n")

# ---- Metadata + sampleId construction (matches flash_fact_char.Rmd) ----
meta <- fread(tcga_metadata, data.table = FALSE) %>%
  dplyr::select(rail_id, tcga_barcode,
                project_id   = gdc_cases.project.project_id,
                submitter_id = gdc_cases.submitter_id) %>%
  mutate(sampleId = sapply(strsplit(tcga_barcode, "-"), function(x) {
    if (length(x) >= 4) x[4] <- sub("[A-Za-z]+$", "", x[4])
    paste(x[seq_len(min(4, length(x)))], collapse = "-")
  }))

# ---- Combined mutation table, build per-rail_id exclusion set ---
stopifnot(length(mut_data_paths) > 0)
mut <- bind_rows(lapply(mut_data_paths, function(p) {
  fread(p, data.table = FALSE) %>%
    dplyr::select(any_of(c("hugoGeneSymbol", "sampleId", "proteinChange")))
})) %>%
  distinct()

sf3b1_mut <- mut %>%
  filter(hugoGeneSymbol == "SF3B1") %>%
  mutate(
    pos = suppressWarnings(as.numeric(str_extract(proteinChange, "\\d+"))),
    is_hotspot      = proteinChange %in% SF3B1_HOTSPOTS,
    is_near_hotspot = !is.na(pos) & pos > NEAR_HOTSPOT_LO &
                      pos < NEAR_HOTSPOT_HI,
    is_extra_excl   = proteinChange %in% EXTRA_EXCLUDE,
    is_excluded     = is_hotspot | is_near_hotspot | is_extra_excl
  )

excluded_sampleids <- unique(c(
  mut$sampleId[mut$hugoGeneSymbol == "SUGP1"],
  sf3b1_mut$sampleId[sf3b1_mut$is_excluded]
))

excluded_rail_ids <- meta %>%
  filter(sampleId %in% excluded_sampleids) %>%
  pull(rail_id) %>%
  unique() %>%
  as.integer()

cat("Excluded rail_ids (SUGP1 / SF3B1 hotspot / near / R549C):",
    length(excluded_rail_ids), "\n")

# ---- Cases: high loading minus exclusions -----------------------
case_df <- tibble(rail_id = sample_ids, loading = loadings) %>%
  filter(loading > loading_thresh,
         !rail_id %in% excluded_rail_ids) %>%
  inner_join(meta, by = "rail_id")

case_counts <- case_df %>%
  count(project_id, name = "n_cases")

cat("\nCases per project_id:\n")
print(case_counts)

# ---- Controls: low loading minus exclusions, per-project matched ----
control_pool <- tibble(rail_id = sample_ids, loading = loadings) %>%
  filter(loading < control_cutoff,
         !rail_id %in% excluded_rail_ids) %>%
  inner_join(meta, by = "rail_id") %>%
  semi_join(case_counts, by = "project_id") %>%
  anti_join(case_df, by = "rail_id")

set.seed(seed_val)
controls_used <- control_pool %>%
  inner_join(case_counts, by = "project_id") %>%
  group_by(project_id) %>%
  group_modify(~ dplyr::slice_sample(.x, n = min(.x$n_cases[1], nrow(.x)))) %>%
  ungroup() %>%
  dplyr::select(rail_id, loading, tcga_barcode, project_id, submitter_id)

# ---- Final samples.tsv (same schema as select_samples.R) --------
samples <- bind_rows(
  case_df %>%
    dplyr::select(rail_id, loading, tcga_barcode,
                  project_id, submitter_id) %>%
    mutate(role = "case"),
  controls_used %>% mutate(role = "control")
) %>%
  distinct(rail_id, .keep_all = TRUE) %>%
  arrange(role, project_id, desc(loading))

samples_path  <- file.path(out_dir, "samples.tsv")
barcodes_path <- file.path(out_dir, "barcodes_to_request.txt")

fwrite(samples, samples_path, sep = "\t")
writeLines(samples$tcga_barcode, barcodes_path)

cat("\nWrote", samples_path, "(", nrow(samples), "samples:",
    sum(samples$role == "case"),    "cases /",
    sum(samples$role == "control"), "controls)\n")
cat("Wrote", barcodes_path, "\n")
