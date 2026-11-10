
# Shared helpers for the supplementary_plots pipeline.

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(yaml)
  library(fastTopics)
  library(ggsci)
})

# Null-coalescing operator (avoid rlang dependency).
`%||%` <- function(x, y) if (is.null(x)) y else x

# Resolve which group (read_length / tcga_cohort) a gene belongs to.
lookup_entry <- function(cfg, gene) {
  if (gene %in% names(cfg$read_length)) {
    return(list(group = "read_length", entry = cfg$read_length[[gene]]))
  }
  if (gene %in% names(cfg$tcga_cohort)) {
    return(list(group = "tcga_cohort", entry = cfg$tcga_cohort[[gene]]))
  }
  stop("Gene not present in either config group: ", gene)
}

# Build the per-sample group label used to facet the structure plot.
#   read_length group: paste0(<cohort sans TCGA->, "-", <read_len>, "bp")
#   tcga_cohort group: paste0(<cohort sans TCGA->)
make_group_labels <- function(meta_df, group) {
  short <- sub("^TCGA-", "", meta_df$project_id)
  if (identical(group, "read_length")) {
    rl <- round(suppressWarnings(as.numeric(meta_df$avg_len)))
    paste0(short, "-", rl, "bp")
  } else {
    short
  }
}

# Read TCGA metadata, restrict to (rail_id, project_id, avg_len).
read_tcga_meta <- function(path) {
  m <- fread(path, data.table = FALSE)
  if (!"avg_len" %in% names(m)) m$avg_len <- NA_real_
  m %>%
    dplyr::select(rail_id,
                   project_id = gdc_cases.project.project_id,
                   avg_len) %>%
    dplyr::mutate(rail_id = as.integer(rail_id))
}

# Read GTEx metadata, restrict to (rail_id, study).
read_gtex_meta <- function(path) {
  m <- fread(path, data.table = FALSE)
  m %>%
    dplyr::select(rail_id, study) %>%
    dplyr::mutate(rail_id = as.integer(rail_id))
}

# Compute the DE junctions for a gene under (factor, thresholds).
de_junction_names <- function(res, de_res, fact, thresholds) {
  fact_idx <- as.integer(sub("^k", "", fact))
  rownames(res$F)[
    de_res$svalue[, fact_idx]   < thresholds$sval_thresh &
    de_res$postmean[, fact_idx] > thresholds$lfc_thresh &
    res$F[, fact_idx]           > thresholds$abs_thresh
  ]
}

# Reuse the parent showcase_genes palette so structure-plot colours stay
# in lock-step with showcase_genes outputs. We sys.source the parent lib
# into a private env and lift out just `factor_palette` -- the parent
# helpers (parse_args, load_gene_data, ...) are intentionally not
# re-exported.
.supp_this_dir   <- dirname(sub("--file=", "",
  commandArgs(trailingOnly = FALSE)[
    grep("--file=", commandArgs(trailingOnly = FALSE))]))
.supp_parent_lib <- normalizePath(
  file.path(.supp_this_dir, "..", "..", "scripts", "lib.R"),
  mustWork = TRUE)
.supp_parent_env <- new.env(parent = baseenv())
sys.source(.supp_parent_lib, envir = .supp_parent_env)
factor_palette <- .supp_parent_env$factor_palette
rm(.supp_this_dir, .supp_parent_lib, .supp_parent_env)
