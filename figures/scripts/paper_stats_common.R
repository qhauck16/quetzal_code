# Shared helpers for figures/figure*/paper_stats.Rmd: paths, and a claim()
# recorder that sets each recomputed number next to the value printed in the
# manuscript text, so every Rmd ends with an "old vs new" table.
#
# Each paper_stats.Rmd takes params$inputs:
#   "current"    -- today's pipeline outputs (default)
#   "pre_rerun"  -- the backed-up inputs from before the 2026-10 re-runs
#                   (diff-splicing df_final, Rmd-era phenotype scores), used
#                   to check that the code reproduces the text's numbers.
if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
suppressPackageStartupMessages({
  library(yaml); library(data.table); library(tidyverse)
})
select <- dplyr::select; filter <- dplyr::filter; rename <- dplyr::rename

PS_FIG  <- file.path(QROOT, "figures")
source(file.path(PS_FIG, "scripts/vt_paths.R"))

ps_inputs <- function() if (exists("params") && !is.null(params$inputs)) params$inputs else "current"

# df_final for a diff-splicing comparison, current or pre-rerun backup
uni_path <- function(gene, comp, thr, file = "df_final.tsv") {
  root <- file.path(QROOT, "leafcutter_analyses/diff_splicing/unified_results")
  if (ps_inputs() == "pre_rerun") root <- file.path(root, ".pre_fix_2026-10-05")
  file.path(root, gene, comp, thr, file)
}
# phenotype-axis scores: quetzal_code cache (current) or the Rmd-era file
phen_path <- function(gene) {
  fig <- c(sf3b1 = "figure2", nkap = "figure3", ylpm1 = "figure4")[[gene]]
  if (ps_inputs() == "pre_rerun")
    return(file.path(QROOT, "leafcutter_analyses/diff_splicing", gene, paste0(gene, "_phenotype_axis_scores.tsv")))
  file.path(PS_FIG, fig, "cache", paste0(gene, "_phenotype_axis_scores.tsv"))
}
# phenotype-positive = score > 0 among P and intermediates (P samples the SVM
# scores <= 0 are misclassified and excluded), as in the original Rmds.
phen_pos <- function(phen) phen %>% filter(group %in% c("intermediate", "P"), score > 0)

tcga_meta <- function() fread(file.path(QROOT, "inputs/tcga_v2_samples.tsv"), data.table = FALSE) %>%
  transmute(rail_id, project_id = gdc_cases.project.project_id, tcga_barcode)

is_blank <- function(x) is.na(x) | x == "" | x == "NA"

# ---- claim recorder -------------------------------------------------
CLAIMS <- tibble(section = character(), claim = character(),
                 text = character(), now = character())
claim <- function(section, claim, text, now) {
  CLAIMS <<- add_row(CLAIMS, section = section, claim = claim,
                     text = as.character(text), now = as.character(now))
  cat(sprintf("[%s] %s: text = %s | now = %s\n", section, claim, text, now))
}
fmt_p <- function(p) formatC(p, format = "e", digits = 1)
claims_table <- function() {
  CLAIMS %>% mutate(changed = ifelse(text == now, "", "**CHANGED**"))
}
