
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
# args: [features, scored_3ss, scored_5ss, branchpoints, output]
features_path     <- args[1]
scored_3ss_path   <- args[2]
scored_5ss_path   <- args[3]
branchpoints_path <- args[4]
out_path          <- args[5]

feats <- read_tsv(features_path, show_col_types = FALSE)

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

if (nrow(feats) == 0) {
  write_tsv(tibble(), out_path)
  quit(save = "no", status = 0)
}

read_scores <- function(path) {
  if (file.info(path)$size == 0) return(rep(NA_real_, nrow(feats)))
  s <- read_tsv(path, col_names = c("sequence", "score"), show_col_types = FALSE)
  if (nrow(s) != nrow(feats)) {
    stop(sprintf("Row-count mismatch between MaxEnt scores (%d) and features (%d): %s",
                 nrow(s), nrow(feats), path))
  }
  suppressWarnings(as.numeric(s$score))
}

feats$score_3ss <- read_scores(scored_3ss_path)
feats$score_5ss <- read_scores(scored_5ss_path)

# ---- branchpoint: nearest branchpoint UPSTREAM of the 3'ss --------------
# Branchpoints sit in the intron just upstream of the 3'ss. Upstream
# direction depends on strand:
#   + strand: branchpoint_position < three_ss; signed distance = three_ss - bp
#   - strand: branchpoint_position > three_ss; signed distance = bp - three_ss
# A simple abs(bp - three_ss) used to pick "the nearest branchpoint in
# either direction", which on minus-strand introns picked the nearest
# downstream branchpoint (i.e., one belonging to a *different* intron).
bp <- fread(branchpoints_path, data.table = FALSE)
colnames(bp)[colnames(bp) == "Chromosome"] <- "chromosome"
colnames(bp)[colnames(bp) == "Strand"]     <- "strand"

joined <- feats %>%
  inner_join(bp, by = c("chromosome", "strand")) %>%
  mutate(bp_distance_to_3ss = case_when(
    strand == "+" ~ three_ss - Branchpoint_position,
    strand == "-" ~ Branchpoint_position - three_ss,
    TRUE          ~ NA_real_
  )) %>%
  filter(bp_distance_to_3ss > 0) %>%
  group_by(gene_name, chromosome, start, end, strand) %>%
  slice_min(bp_distance_to_3ss, n = 1, with_ties = FALSE) %>%
  ungroup()

bp_cols <- joined %>%
  dplyr::select(gene_name, chromosome, start, end, strand,
                Branchpoint_score, bp_distance_to_3ss)

feats <- feats %>%
  left_join(bp_cols, by = c("gene_name", "chromosome", "start", "end", "strand"))

write_tsv(feats, out_path)
