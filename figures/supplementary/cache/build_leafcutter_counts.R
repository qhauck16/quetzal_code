# Cache for supplementary figure S8 (LeafCutter summary). The code below
# is leafcutter_meta_analysis.Rmd's config / mane / helpers /
# load_leaf_sigs / counts chunks, VERBATIM (extracted by
# figures/supplementary build notes), run from the Rmd's own directory so
# its relative paths resolve. Writes cache/s8_counts.tsv.
.libPaths('/scratch/midway2/qhauck/conda_env/rstudio-server/lib/R/library/')
if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
suppressPackageStartupMessages({ library(tidyverse); library(data.table) })
setwd(file.path(QROOT, "leafcutter_analyses/diff_splicing/unified_results"))
GENES        <- c("sf3b1", "nkap", "ylpm1")
COMPARISONS  <- c("10v30", "20v60", "30v90")
DPSI_VALUES  <- c(0.01, 0.05)
FC_VALUES    <- c(1, 2.5, 5)

# Paths (relative to this Rmd inside diff_splicing/unified_results).
REAL_DIR     <- "../results"
CTRL_DIR     <- "../control_results"
INTRONS_PATH <- "../gencode_v46_introns.tsv"
MANE_PATH    <- "../MANE.GRCh38.v1.4.summary.txt"
OUT_PDF      <- "leafcutter_meta_analysis.pdf"

SPLICE_TYPES <- c("alt 3'ss", "alt 5'ss", "canonical", "skipped exon", "other")
GENE_LABELS  <- c(sf3b1 = "SF3B1", nkap = "NKAP", ylpm1 = "YLPM1")

# Order condition factor so each gene sits next to its control on the
# x-axis: SF3B1, SF3B1 ctrl, NKAP, NKAP ctrl, YLPM1, YLPM1 ctrl.
COND_LEVELS  <- c("SF3B1",  "SF3B1 ctrl",
                  "NKAP",   "NKAP ctrl",
                  "YLPM1",  "YLPM1 ctrl")

introns <- fread(INTRONS_PATH, data.table = FALSE)
mane_tx <- fread(MANE_PATH,    data.table = FALSE)

mane_introns <- introns %>%
  filter(transcript_id %in% mane_tx$Ensembl_nuc) %>%
  dplyr::select(start, end, strand, gene_name, seqnames) %>%
  unique()

intron_to_gene <- introns %>%
  dplyr::select(start, end, strand, gene_name) %>%
  unique()

# Mirrors unified_results/scripts/load_comparison.R: join effect_sizes +
# cluster_significance, attach gene_name via the introns table, keep
# clusters with p.adjust <= 0.05.
load_full_leaf_sig <- function(dir_path) {
  es_path <- file.path(dir_path, "leafcutter_ds_effect_sizes.txt")
  sg_path <- file.path(dir_path, "leafcutter_ds_cluster_significance.txt")
  if (!file.exists(es_path) || !file.exists(sg_path)) return(NULL)

  leaf_eff <- read_tsv(es_path, show_col_types = FALSE)
  leaf_sig <- read_tsv(sg_path, show_col_types = FALSE)

  leaf_to_gene <- leaf_eff %>%
    separate(intron,
              into = c("chrom", "start", "end", "cluster_id"),
              sep  = ":") %>%
    mutate(start  = as.numeric(start) + 1,
            end    = as.numeric(end)   - 1,
            strand = str_split_i(cluster_id, "_", 3),
            # Normalize sign: leafcutter writes columns test/control in
            # whatever order leaf_pheno_data.txt has them; the test cohort
            # is swapped to col-1 for control comparisons, which flips
            # the sign of the supplied deltapsi. Always recompute as
            # `test - control` so positive == "more in test" everywhere.
            deltapsi = test - control) %>%
    left_join(intron_to_gene, by = c("start", "end", "strand"))

  cluster_to_gene <- leaf_to_gene %>%
    group_by(cluster_id) %>%
    mutate(gene_name = {
      x <- gene_name[!is.na(gene_name)]
      if (length(x) == 0L) NA_character_
      else names(which.max(table(x)))
    }) %>%
    ungroup()

  leaf_sig <- leaf_sig %>%
    separate(cluster, into = c("chrom", "cluster_id"), sep = ":")

  cluster_to_gene %>%
    left_join(leaf_sig, by = c("chrom", "cluster_id")) %>%
    filter(!is.na(p.adjust), p.adjust <= 0.05) %>%
    unique()
}

# Mirrors unified_results/scripts/assign_splice_types.R: filter by
# deltapsi + fold-change thresholds and classify each junction.
assign_splice_types <- function(full_leaf_sig, dpsi_thresh, fc_thresh) {
  if (is.null(full_leaf_sig) || nrow(full_leaf_sig) == 0L) {
    return(tibble(type = character(0)))
  }
  hits <- full_leaf_sig %>% filter(deltapsi > dpsi_thresh)
  if (fc_thresh > 0) hits <- hits %>% filter(test / control > fc_thresh)
  if (nrow(hits) == 0L) return(tibble(type = character(0)))

  top <- hits %>%
    mutate(start = as.numeric(start), end = as.numeric(end)) %>%
    filter(!is.na(gene_name), gene_name %in% mane_introns$gene_name)
  if (nrow(top) == 0L) return(tibble(type = character(0)))

  mid <- function(d) paste0(d$gene_name, ":", d$start, "-", d$end)

  canonical <- top %>%
    inner_join(mane_introns %>%
                  dplyr::select(gene_name, strand, start, end) %>%
                  rename(mane_start = start, mane_end = end),
                by = c("gene_name", "strand",
                       "start" = "mane_start", "end" = "mane_end")) %>%
    mutate(type = "canonical")
  canonical_ids <- mid(canonical)

  skipped <- top %>%
    inner_join(mane_introns %>%
                  dplyr::select(gene_name, strand, start, end) %>%
                  rename(up_start = start, up_end = end),
                by = c("gene_name", "strand", "start" = "up_start")) %>%
    inner_join(mane_introns %>%
                  dplyr::select(gene_name, strand, start, end) %>%
                  rename(down_start = start, down_end = end),
                by = c("gene_name", "strand", "end" = "down_end")) %>%
    filter(up_end < down_start) %>%
    mutate(type = "skipped exon")
  skipped_ids <- mid(skipped)
  assigned    <- c(canonical_ids, skipped_ids)

  start_m <- top %>% filter(!mid(.) %in% assigned) %>%
    inner_join(mane_introns %>%
                  dplyr::select(gene_name, strand, start, end) %>%
                  rename(mane_start = start, mane_end = end),
                by = c("gene_name", "strand", "start" = "mane_start")) %>%
    mutate(type = if_else(strand == "+", "alt 3'ss", "alt 5'ss"))

  end_m <- top %>% filter(!mid(.) %in% assigned) %>%
    inner_join(mane_introns %>%
                  dplyr::select(gene_name, strand, start, end) %>%
                  rename(mane_start = start, mane_end = end),
                by = c("gene_name", "strand", "end" = "mane_end")) %>%
    mutate(type = if_else(strand == "+", "alt 5'ss", "alt 3'ss"))

  ann <- bind_rows(canonical, skipped, start_m, end_m) %>%
    mutate(junc_id = mid(.)) %>%
    dplyr::select(junc_id, type)

  top %>%
    mutate(junc_id = mid(.)) %>%
    left_join(ann, by = "junc_id") %>%
    mutate(type = if_else(is.na(type), "other", type)) %>%
    group_by(gene_name, start, end) %>%
    filter(n() < 2) %>%
    ungroup()
}

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
