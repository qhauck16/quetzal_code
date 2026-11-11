# Cache for supplementary figure S8 (LeafCutter summary). Runs from
# leafcutter_analyses/diff_splicing/unified_results so relative paths
# resolve. Writes cache/s8_counts.tsv.
if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
suppressPackageStartupMessages({ library(tidyverse); library(data.table) })
setwd(file.path(QROOT, "leafcutter_analyses/diff_splicing/unified_results"))
GENES        <- c("sf3b1", "nkap", "ylpm1")
COMPARISONS  <- c("10v30", "20v60", "30v90")
DPSI_VALUES  <- c(0.01, 0.05)
FC_VALUES    <- c(1, 2.5, 5)

# Paths (relative to unified_results).
REAL_DIR     <- "../results"
CTRL_DIR     <- "../control_results"
OUT_PDF      <- "leafcutter_meta_analysis.pdf"

SPLICE_TYPES <- c("alt 3'ss", "alt 5'ss", "canonical", "skipped exon", "other")
GENE_LABELS  <- c(sf3b1 = "SF3B1", nkap = "NKAP", ylpm1 = "YLPM1")

# Order condition factor so each gene sits next to its control on the
# x-axis: SF3B1, SF3B1 ctrl, NKAP, NKAP ctrl, YLPM1, YLPM1 ctrl.
COND_LEVELS  <- c("SF3B1",  "SF3B1 ctrl",
                  "NKAP",   "NKAP ctrl",
                  "YLPM1",  "YLPM1 ctrl")

source(file.path(QROOT, "figures/supplementary/cache/leafcutter_counts_lib.R"))

combos <- expand_grid(
  gene       = GENES,
  condition  = c("real", "control"),
  comparison = COMPARISONS
) %>%
  mutate(dir_path = if_else(
    condition == "real",
    file.path(REAL_DIR, paste0(gene, "_", comparison)),
    file.path(CTRL_DIR, paste0(gene, "_", comparison))
  ))

combos$leaf <- map(combos$dir_path, function(p) {
  cat(sprintf("  loading %s\n", p))
  load_full_leaf_sig(p)
})

combos %>% mutate(n_sig = map_int(leaf, ~ if (is.null(.x)) NA_integer_ else nrow(.x))) %>%
  dplyr::select(gene, condition, comparison, n_sig) %>%
  arrange(gene, comparison, condition)

threshold_grid <- expand_grid(dpsi = DPSI_VALUES, fc = FC_VALUES)

counts_list <- list()
for (i in seq_len(nrow(combos))) {
  cm   <- combos[i, ]
  leaf <- cm$leaf[[1]]
  for (j in seq_len(nrow(threshold_grid))) {
    th <- threshold_grid[j, ]
    df_local <- assign_splice_types(leaf, th$dpsi, th$fc)
    if (nrow(df_local) > 0L) {
      cnt <- df_local %>% dplyr::count(type)
    } else {
      cnt <- tibble(type = character(0), n = integer(0))
    }
    cnt <- tibble(type = SPLICE_TYPES) %>%
      left_join(cnt, by = "type") %>%
      mutate(n = replace_na(n, 0L))
    cnt$gene       <- cm$gene
    cnt$condition  <- cm$condition
    cnt$comparison <- cm$comparison
    cnt$dpsi       <- th$dpsi
    cnt$fc         <- th$fc
    counts_list[[length(counts_list) + 1L]] <- cnt
  }
}
counts <- bind_rows(counts_list)
counts

fwrite(counts, file.path(QROOT, "figures/supplementary/cache/s8_counts.tsv"), sep = "\t")
