# Standardized per-sample TCGA covariate table for downstream pipelines
# (cancer_specific_factors.R, flash_fact_char.Rmd).
#
# Numeric covariates are z-scored across the full TCGA population using
# means and SDs computed on observed values only; NAs are left as NA --
# downstream analyses are expected to drop incomplete samples via
# complete.cases() (or to exclude columns they can't carry, as
# cancer_specific_factors.R does with RIN).
#
# Run interactively in RStudio with the project's working directory set
# to .../spliceNMF/tcga/. Edit the paths below if anything moves.


suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(bit64)
})

tcga_metadata_path <- "../inputs/tcga_v2_samples.tsv"
rin_score_file     <- "../inputs/analyte.tsv"
out_path           <- "standardized_covariates.RDS"

tcga_metadata <- fread(tcga_metadata_path, data.table = FALSE)

rin_scores <- read_tsv(rin_score_file, show_col_types = FALSE) %>%
  filter(analytes.rna_integrity_number != "'--") %>%
  dplyr::select(analytes.submitter_id, analytes.rna_integrity_number) %>%
  mutate(analytes.rna_integrity_number =
           as.numeric(analytes.rna_integrity_number))

columns_to_keep <- c(
  "gdc_cases.project.project_id", "rail_id",
  "gdc_cases.demographic.gender", "gdc_platform",
  "analytes.rna_integrity_number", "gdc_cases.samples.is_ffpe",
  "%C", "%G", "junction_count", "junction_avg_coverage",
  "avg_len", "plate"
)

covar_tbl <- tcga_metadata %>%
  mutate(plate = str_split_i(tcga_barcode, "-", 6)) %>%
  left_join(rin_scores,
            by = c("gdc_cases.samples.portions.analytes.submitter_id" =
                     "analytes.submitter_id")) %>%
  dplyr::select(all_of(columns_to_keep)) %>%
  # `%C` / `%G` get mangled by base data-frame coercion; rename now so
  # consumers don't need backticks.
  dplyr::rename(pct_C = `%C`, pct_G = `%G`)

# integer64 -> numeric for downstream compatibility
covar_tbl[] <- lapply(covar_tbl, function(x) {
  if (inherits(x, "integer64")) as.numeric(x) else x
})

# rail_id stays an identifier; pull it out before z-scoring so it isn't
# itself scaled.
rail_ids <- covar_tbl$rail_id
covar_tbl$rail_id <- NULL

# NA-safe z-score: mean/sd ignore NAs, NAs in input stay NA in output.
scale_keep_na <- function(x) {
  if (!is.numeric(x)) return(x)
  m <- mean(x, na.rm = TRUE)
  s <- sd(x,   na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - m) / s
}
covar_tbl <- as.data.frame(lapply(covar_tbl, scale_keep_na))

covar_tbl$rail_id <- rail_ids
covar_tbl <- covar_tbl[, c("rail_id",
                           setdiff(colnames(covar_tbl), "rail_id"))]

saveRDS(covar_tbl, out_path)
tsv_path <- sub("\\.RDS$", ".tsv", out_path, ignore.case = TRUE)
write.table(covar_tbl, tsv_path,
            sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("Wrote %d rows x %d cols to %s\n",
            nrow(covar_tbl), ncol(covar_tbl), out_path))
