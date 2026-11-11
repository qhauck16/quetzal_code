# Shared LeafCutter junction counting for Fig. S8 (cache/build_leafcutter_counts.R)
# and Supplementary Note 1 Fig. 4 (note1/fig4.R), so both count the same junctions:
# significant clusters (adj. P <= 0.05), dPSI = test - control, junctions assigned to a
# gene with a MANE Select transcript, splice type relative to MANE introns.
if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
suppressPackageStartupMessages({ library(tidyverse); library(data.table) })
LC_DIR <- file.path(QROOT, "leafcutter_analyses/diff_splicing")
introns <- fread(file.path(LC_DIR, "gencode_v46_introns.tsv"), data.table = FALSE)
mane_tx <- fread(file.path(QROOT, "inputs/MANE.GRCh38.v1.4.summary.txt"),data.table = FALSE)

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
