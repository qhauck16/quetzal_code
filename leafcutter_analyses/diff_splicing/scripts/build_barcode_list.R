#!/usr/bin/env Rscript
# Build the TCGA barcode list for GDC manifest download.
#
# Output:
#   download/barcodes_to_request.tsv   -- one TCGA barcode per line; submit to GDC
#                                        portal -> get back a manifest file
#   download/sample_plan.tsv           -- rail_id, role(s), gene(s); the
#                                        bookkeeping used post-download to route
#                                        BAMs into the correct folder
#   download/already_have.tsv          -- rail_ids skipped because their BAM is
#                                        already on disk somewhere
#
# Selection rules (per gene):
#   * top 50 by factor_col       -> case set (gene-specific)
#   * shared control candidates  -> samples with all three gene factors
#                                   below their respective thresholds AND
#                                   not in any case set
#   * per-project control budget -> max across genes of
#                                       (n_test_30 in that project) * 3 (ratio)
#                                       * 1.5 (control_run draw)
#                                       * 1.1 (10% buffer)
#                                   = 1.65x sized per project
#
# Existing BAMs (currently in {gene}/bams/) are excluded from the download
# request -- they're already on disk.

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(yaml)
})

set.seed(123)

# Anchor all paths to diff_splicing/ regardless of cwd.
find_diff_dir <- function() {
  # 1. Try to locate via Rscript --file=
  args0 <- commandArgs(trailingOnly = FALSE)
  fa <- args0[grep("--file=", args0)]
  if (length(fa)) {
    sd <- dirname(normalizePath(sub("--file=", "", fa)))
    cand <- normalizePath(file.path(sd, ".."), mustWork = FALSE)
    if (file.exists(file.path(cand, "config.yaml"))) return(cand)
  }
  # 2. Try cwd, parent, grandparent
  for (d in c(getwd(), file.path(getwd(), ".."),
              file.path(getwd(), "../.."))) {
    d <- normalizePath(d, mustWork = FALSE)
    if (file.exists(file.path(d, "config.yaml")) &&
        dir.exists(file.path(d, "scripts"))) return(d)
  }
  stop("Could not locate diff_splicing/ (no config.yaml found near cwd or script).")
}
DIFF_DIR <- find_diff_dir()
setwd(DIFF_DIR)
cat("Working from:", DIFF_DIR, "\n")

CONFIG_PATH       <- file.path(DIFF_DIR, "config.yaml")
METADATA_PATH     <- normalizePath(file.path(DIFF_DIR, "../../../tcga_v2_samples.tsv"))
EXISTING_BAM_DIRS <- c("nkap/bams", "sf3b1/bams", "ylpm1/bams",
                       "spop/bams", "idh1/bams",
                       "nkap_bams", "sf3b1_bams", "ylpm1_bams", "control_bams")
OUT_DIR           <- "download"
N_CASES_PER_GENE  <- 50
OVERSAMPLE        <- 1.65   # 1.5x for control-mode draw * 1.1 buffer

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cfg   <- read_yaml(CONFIG_PATH)
genes <- cfg$genes
comps <- cfg$comparisons
max_n_test <- max(sapply(comps, `[[`, "n_test"))   # 30
max_ratio  <- max(sapply(comps, `[[`, "control_ratio"))   # 3

# -- 1. Load metadata + per-gene info ----------------------------------------
meta <- fread(METADATA_PATH, data.table = FALSE) %>%
  dplyr::select(rail_id, tcga_barcode, gdc_cases.submitter_id,
                gdc_cases.project.project_id,
                avg_len, gdc_cases.samples.is_ffpe) %>%
  filter(!is.na(rail_id), !is.na(tcga_barcode),
         !is.na(gdc_cases.project.project_id))

gene_info <- lapply(names(genes), function(g) {
  fc <- genes[[g]]$factor_col
  read_tsv(genes[[g]]$info_file, show_col_types = FALSE) %>%
    dplyr::select(rail_id, all_of(fc)) %>%
    rename(!!paste0("loading_", g) := !!fc)
})
names(gene_info) <- names(genes)

# Wide table: rail_id + loading_<gene> for each gene, joined with metadata
loadings <- purrr::reduce(gene_info, dplyr::full_join, by = "rail_id")
samples  <- inner_join(meta, loadings, by = "rail_id")

# -- 2. Pick top N_CASES_PER_GENE for each gene ------------------------------
case_sets <- lapply(names(genes), function(g) {
  samples %>%
    arrange(desc(.data[[paste0("loading_", g)]])) %>%
    slice_head(n = N_CASES_PER_GENE) %>%
    pull(rail_id)
})
names(case_sets) <- names(genes)

all_cases <- unique(unlist(case_sets))

# -- 3. Greedy shared control pool across genes ------------------------------
# Per gene g, target per-project budget = (top_n_test_in_p) * max_ratio *
# OVERSAMPLE. We accumulate a single shared pool: for each (gene, project),
# first count members already in the pool that ALSO satisfy gene g's
# threshold and aren't a case for g — those count toward g's budget for
# free. Only draw additional samples to fill the residual gap.
#
# A control may serve multiple genes if its loadings pass each gene's
# threshold; the controls_by_gene.tsv diagnostic records which gene(s)
# each control is eligible for (== which gene's budget it counted toward).

thresholds <- sapply(genes, `[[`, "control_threshold")

gene_proj_budget <- lapply(names(genes), function(g) {
  loading_col <- paste0("loading_", g)
  proj_dist <- samples %>%
    arrange(desc(.data[[loading_col]])) %>%
    slice_head(n = max_n_test) %>%
    count(gdc_cases.project.project_id, name = "n_test_proj") %>%
    mutate(budget = ceiling(n_test_proj * max_ratio * OVERSAMPLE))
  cat(sprintf("\n[%s] top-%d project distribution -> budget per project:\n",
              g, max_n_test))
  print(proj_dist %>% arrange(desc(budget)))
  proj_dist
})
names(gene_proj_budget) <- names(genes)

# Accumulator: rail_ids in the shared pool, plus per-gene assignment ledger
pool_ids <- integer(0)
assignments <- list()  # per gene: tibble(rail_id, gene)

for (g in names(genes)) {
  loading_col <- paste0("loading_", g)
  thr         <- thresholds[[g]]
  proj_dist   <- gene_proj_budget[[g]]

  # Eligibility for gene g: low loading AND not a case for g
  eligible_g <- samples %>%
    filter(!rail_id %in% case_sets[[g]]) %>%
    filter(.data[[loading_col]] < thr)

  picked_for_g <- list()

  for (p in proj_dist$gdc_cases.project.project_id) {
    n_want <- proj_dist$budget[proj_dist$gdc_cases.project.project_id == p]
    proj_eligible <- eligible_g %>%
      filter(gdc_cases.project.project_id == p)

    already <- proj_eligible %>% filter(rail_id %in% pool_ids)
    n_already <- nrow(already)

    if (n_already >= n_want) {
      # Existing pool already covers this gene's project budget; reuse all
      kept <- already %>% slice_head(n = n_want)
      picked_for_g[[p]] <- kept$rail_id
      next
    }

    n_need <- n_want - n_already
    new_pool <- proj_eligible %>% filter(!rail_id %in% pool_ids)
    n_have_new <- nrow(new_pool)
    if (n_have_new < n_need) {
      cat(sprintf("  [%s] %s: wanted %d, have %d in pool + %d new available (capped)\n",
                  g, p, n_want, n_already, n_have_new))
    }
    drawn <- new_pool %>% slice_sample(n = min(n_need, n_have_new))
    pool_ids <- c(pool_ids, drawn$rail_id)
    picked_for_g[[p]] <- c(already$rail_id, drawn$rail_id)
  }

  assignments[[g]] <- tibble(
    rail_id = unlist(picked_for_g, use.names = FALSE),
    gene    = g
  )
  cat(sprintf("[%s] assigned %d controls (target ~%d, shared pool now %d)\n",
              g, nrow(assignments[[g]]), sum(proj_dist$budget),
              length(unique(pool_ids))))
}

controls_long <- bind_rows(assignments)
controls      <- unique(pool_ids)

cat(sprintf("\n=== Total unique controls (union across genes): %d ===\n",
            length(controls)))

# -- 6. Assemble plan --------------------------------------------------------
# Cases are gene-specific (route to {gene}_bams/), controls go to control_bams/.
# Controls list per gene preserved in controls_long for diagnostic write-out.
plan_rows <- bind_rows(
  lapply(names(case_sets), function(g) {
    tibble(rail_id = case_sets[[g]], role = "case", gene = g)
  }),
  tibble(rail_id = controls, role = "control", gene = NA_character_)
)

# A rail_id may appear more than once (case for one gene, control elsewhere).
# Routing rule: case > control. If case for >1 gene, list both (and the
# router will pick the first alphabetically by gene as physical home).
plan <- plan_rows %>%
  distinct(rail_id, role, gene) %>%
  left_join(meta %>% dplyr::select(rail_id, tcga_barcode), by = "rail_id")

# -- 7. Identify already-downloaded BAMs --------------------------------------
existing_rail_ids <- unique(unlist(lapply(EXISTING_BAM_DIRS, function(d) {
  if (!dir.exists(d)) return(integer(0))
  bams <- list.files(d, pattern = "\\.bam$", full.names = FALSE)
  as.integer(sub("_.*", "", bams))
})))
existing_rail_ids <- existing_rail_ids[!is.na(existing_rail_ids)]

cat(sprintf("\nAlready-downloaded rail_ids on disk: %d\n", length(existing_rail_ids)))

needed_plan <- plan %>% filter(!rail_id %in% existing_rail_ids)

barcodes_to_request <- needed_plan %>%
  distinct(rail_id, tcga_barcode) %>%
  filter(!is.na(tcga_barcode))

# -- 8. Write outputs ---------------------------------------------------------
write_tsv(plan, file.path(OUT_DIR, "sample_plan.tsv"))

# Diagnostic: which gene each control was sampled for (a control may have
# been sampled by >1 gene's pool draw, in which case it appears once per gene).
write_tsv(controls_long, file.path(OUT_DIR, "controls_by_gene.tsv"))

write_tsv(tibble(rail_id = existing_rail_ids), file.path(OUT_DIR, "already_have.tsv"))

write_tsv(barcodes_to_request,
          file.path(OUT_DIR, "barcodes_to_request.tsv"))

writeLines(barcodes_to_request$tcga_barcode,
           file.path(OUT_DIR, "barcodes_only.txt"))

cat(sprintf("\n--- Summary ---\nUnique rail_ids in plan:        %d\n",
            n_distinct(plan$rail_id)))
cat(sprintf("Already on disk (skip):         %d\n", length(existing_rail_ids)))
cat(sprintf("New rail_ids to download:       %d\n",
            n_distinct(needed_plan$rail_id)))
cat(sprintf("Wrote: %s/{sample_plan,barcodes_to_request,barcodes_only,already_have}.tsv\n",
            OUT_DIR))
