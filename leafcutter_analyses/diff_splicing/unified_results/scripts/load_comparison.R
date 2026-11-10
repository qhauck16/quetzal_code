
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
# args: [cluster_sig, effect_size, introns, output]
cluster_sig_path <- args[1]
effect_size_path <- args[2]
introns_path     <- args[3]
out_path         <- args[4]

introns <- fread(introns_path, data.table = FALSE)
intron_to_gene <- introns %>%
  dplyr::select(start, end, strand, gene_name) %>%
  unique()

leaf_sig <- read_tsv(cluster_sig_path, show_col_types = FALSE)
leaf_eff <- read_tsv(effect_size_path, show_col_types = FALSE)

leaf_to_gene <- leaf_eff %>%
  separate(intron, into = c("chrom", "start", "end", "cluster_id"), sep = ":") %>%
  mutate(start = as.numeric(start) + 1,
         end   = as.numeric(end)   - 1)
leaf_to_gene$strand <- str_split_i(leaf_to_gene$cluster_id, "_", 3)

leaf_to_gene <- leaf_to_gene %>%
  left_join(intron_to_gene, by = c("start", "end", "strand"))

cluster_to_gene <- leaf_to_gene %>%
  group_by(cluster_id) %>%
  mutate(
    gene_name = {
      x <- gene_name[!is.na(gene_name)]
      if (length(x) == 0) NA_character_ else names(which.max(table(x)))
    }
  ) %>%
  ungroup()

leaf_sig <- leaf_sig %>%
  separate(cluster, into = c("chrom", "cluster_id"), sep = ":")

full_leaf_sig <- cluster_to_gene %>%
  left_join(leaf_sig, by = c("chrom", "cluster_id")) %>%
  filter(p.adjust <= 0.05) %>%
  unique()

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
write_tsv(full_leaf_sig, out_path)
