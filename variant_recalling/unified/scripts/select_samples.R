
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(flashier)
})

args <- commandArgs(trailingOnly = TRUE)
gene              <- args[1]                   # "nkap" / "ylpm1"
factor_key        <- args[2]                   # "k13" / "k30"
rds_path          <- args[3]
tcga_metadata     <- args[4]
loading_thresh    <- as.numeric(args[5])
control_cutoff    <- as.numeric(args[6])
seed_val          <- as.integer(args[7])
out_dir           <- args[8]
exclude_mut_genes <- args[9]                   # comma-separated HUGO symbols (or "")
mut_data_paths    <- args[-(1:9)]              # remaining args = compiled mutation tsvs

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("Loading RDS:", rds_path, "\n")
out <- readRDS(rds_path)

ldf  <- flashier::ldf(out, type = 'i')
L    <- ldf$L
colnames(L) <- paste0('k', seq_len(ncol(L)))

if (!factor_key %in% colnames(L)) {
  stop("factor_key ", factor_key, " not in norm_L columns (have: ",
       paste(head(colnames(L), 5), collapse = ", "), " ... )")
}

sample_ids <- as.integer(rownames(L))
loadings   <- L[, factor_key]

cat("Total samples in RDS:", length(sample_ids), "\n")
cat("Cases (loading >", loading_thresh, "):",
    sum(loadings > loading_thresh, na.rm = TRUE), "\n")
cat("Eligible controls (loading <", control_cutoff, "):",
    sum(loadings < control_cutoff, na.rm = TRUE), "\n")

# Bring in tcga metadata
meta <- fread(tcga_metadata, data.table = FALSE) %>%
  dplyr::select(rail_id, tcga_barcode,
                project_id   = gdc_cases.project.project_id,
                submitter_id = gdc_cases.submitter_id) %>%
  mutate(sampleId = sapply(strsplit(tcga_barcode, "-"), function(x) {
    if (length(x) >= 4) x[4] <- sub("[A-Za-z]+$", "", x[4])
    paste(x[seq_len(min(4, length(x)))], collapse = "-")
  }))

# ---- Optional: exclude rail_ids with a pre-existing variant in any of
#      `exclude_mut_genes` (HUGO symbols, comma-separated). Builds the
#      same sample-x-gene matrix logic as flash_fact_char.Rmd.
exclude_genes <- if (nzchar(exclude_mut_genes)) {
  trimws(strsplit(exclude_mut_genes, ",")[[1]])
} else {
  character(0)
}

mutated_rail_ids <- integer(0)
if (length(exclude_genes) > 0 && length(mut_data_paths) > 0) {
  cat("Excluding rail_ids with a variant in:",
      paste(exclude_genes, collapse = ", "), "\n")
  mut_data <- bind_rows(lapply(mut_data_paths, function(p) {
    fread(p, data.table = FALSE) %>%
      dplyr::select(any_of(c("hugoGeneSymbol", "sampleId")))
  })) %>%
    distinct() %>%
    filter(hugoGeneSymbol %in% exclude_genes)

  mutated_rail_ids <- meta %>%
    inner_join(mut_data, by = "sampleId") %>%
    pull(rail_id) %>%
    unique() %>%
    as.integer()

  cat("  -> excluding", length(mutated_rail_ids),
      "rail_ids carrying a known variant\n")
}

# Cases: every sample with loading above the threshold and no excluded variant
case_df <- tibble(rail_id = sample_ids, loading = loadings) %>%
  filter(loading > loading_thresh,
         !rail_id %in% mutated_rail_ids) %>%
  inner_join(meta, by = "rail_id")

case_counts <- case_df %>%
  count(project_id, name = "n_cases")

cat("Cases per project_id:\n")
print(case_counts)

# Controls: matched per project_id, also excluding mutated rail_ids
control_pool <- tibble(rail_id = sample_ids, loading = loadings) %>%
  filter(loading < control_cutoff,
         !rail_id %in% mutated_rail_ids) %>%
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

samples <- bind_rows(
  case_df       %>% dplyr::select(rail_id, loading, tcga_barcode,
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
    sum(samples$role == "case"), "cases /",
    sum(samples$role == "control"), "controls)\n")
cat("Wrote", barcodes_path,
    "-- paste into GDC portal, filter for WXS + RNA-Seq STAR 2-Pass Genome BAMs,\n",
    "  drop the manifest + metadata at:", out_dir,
    "as gdc_manifest.txt and metadata.cohort.json.\n")
