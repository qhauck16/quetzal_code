
# Per (gsp, plot_gene, region) build:
#   sample_bigwigs.tsv   (rail_id, abs bigwig path, group, strand)
#   groups_file.txt      (Group_label, color, Bedgz, Supergroup)
#   junctions.bed        (mean PSI per junction across the union of the
#                         condition + matched-controls sample set)

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
samples_tsv     <- args[1]
plot_gene       <- args[2]
chr             <- args[3]
strand          <- args[4]
condition_label <- args[5]   # e.g. "High_k13"
condition_color <- args[6]
tcga_bw_dir     <- args[7]   # already absolute or rooted in pygenome/
tcga_psi_dir    <- args[8]
out_bigwigs_tsv <- args[9]
out_groups_txt  <- args[10]
out_junc_bed    <- args[11]

dir.create(dirname(out_bigwigs_tsv), showWarnings = FALSE, recursive = TRUE)

samples <- fread(samples_tsv, data.table = FALSE) %>% as_tibble()
stopifnot(all(c("rail_id", "role", "project_short") %in% colnames(samples)))

bw_root <- normalizePath(tcga_bw_dir, mustWork = FALSE)

sample_bigwigs <- samples %>%
  mutate(
    bigwigPath  = file.path(bw_root, project_short, paste0(rail_id, ".bw")),
    Group_label = if_else(role == "condition",
                          condition_label, "Matched_Controls"),
    strand      = "."
  ) %>%
  dplyr::select(sample = rail_id, bigwigPath, Group_label, strand)

write.table(sample_bigwigs, out_bigwigs_tsv,
            sep = "\t", quote = FALSE, row.names = FALSE)

groups <- tibble(
  Group_label = c(condition_label, "Matched_Controls"),
  color       = c(condition_color, "#969696"),
  Bedgz       = "",
  Supergroup  = ""
)
write.table(groups, out_groups_txt,
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- Junctions BED -------------------------------------------------
psi_path <- file.path(tcga_psi_dir, chr, plot_gene, "psi_table.tsv")
if (!file.exists(psi_path)) {
  stop("PSI table not found: ", psi_path)
}

suppressWarnings({
  psi <- fread(psi_path, data.table = FALSE) %>%
    filter(V1 %in% sample_bigwigs$sample) %>%
    column_to_rownames("V1")
})
psi[is.na(psi)] <- 0

make_bed_df <- function(mat, chr, strand) {
  df_vals <- as.data.frame(t(mat * 100))
  df_vals$junc_id <- rownames(df_vals)

  df_meta <- df_vals %>%
    separate(junc_id,
             into = c("clu", "range"),
             sep = ":", remove = FALSE) %>%
    separate(range,
             into = c("start", "end"),
             sep = "-", convert = TRUE)

  df_meta$`#Chrom` <- chr
  df_meta$strand   <- strand
  df_meta$start    <- df_meta$start - 1
  df_meta$junc <- paste0(df_meta$`#Chrom`, ":", df_meta$start, ":",
                         df_meta$end, ":clu_", df_meta$clu, "_",
                         df_meta$strand)
  df_meta$gid  <- paste0(df_meta$`#Chrom`, "_clu_", df_meta$clu, "_",
                         df_meta$strand)

  df_meta %>%
    relocate(`#Chrom`, start, end, junc, gid, strand) %>%
    dplyr::select(-clu, -junc_id, -junc)
}

bed <- make_bed_df(psi, chr, strand) %>%
  group_by(start, end) %>%
  slice_head(n = 1) %>%
  ungroup()

num_cols <- sapply(bed, is.numeric)
bed[, num_cols] <- lapply(bed[, num_cols, drop = FALSE],
                         function(x) round(as.double(x), 4))

write.table(bed, out_junc_bed,
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Wrote:\n  ", out_bigwigs_tsv,
    "\n  ", out_groups_txt,
    "\n  ", out_junc_bed,
    "\n")
