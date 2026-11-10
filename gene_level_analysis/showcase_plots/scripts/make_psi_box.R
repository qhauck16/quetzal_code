#!/usr/bin/env Rscript
# Per-sample PSI box-plot of the TOP DE junction of a read-length gene,
# split by sample read-length (<=50 bp vs >=75 bp). Read-length scan
# already showed these factors load differently across 48bp / 76bp
# samples within a cohort -- this plot verifies the per-junction PSI
# difference for the top driver junction.

script_dir <- dirname(sub("--file=", "",
  commandArgs(trailingOnly = FALSE)[
    grep("--file=", commandArgs(trailingOnly = FALSE))]))
source(file.path(script_dir, "lib_supp.R"))

suppressPackageStartupMessages({
  library(ggplot2)
})

args        <- commandArgs(trailingOnly = TRUE)
gene        <- args[1]
config_path <- args[2]
out_pdf     <- args[3]

dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)

cfg <- yaml::read_yaml(config_path)
le  <- lookup_entry(cfg, gene)
if (!identical(le$group, "read_length")) {
  stop("psi_box is only defined for read_length genes; ", gene,
       " is in the ", le$group, " group.")
}
ent     <- le$entry
chr     <- ent$chr
fact    <- if (is.null(ent$factor)) "" else as.character(ent$factor)
cohorts <- as.character(unlist(ent$cohorts %||% list()))
if (!nzchar(fact)) {
  stop("config[read_length][", gene, "]: factor must be set.")
}
if (length(cohorts) == 0L) {
  stop("config[read_length][", gene, "]: cohorts must be set.")
}

cat(sprintf("[psi_box] %s | %s | %s\n", gene, chr, fact))

# ---- DE junctions and the top one by postmean -----------------------------

res_path    <- file.path(cfg$paths$tcga_chr_dir, chr, gene, "res.RDS")
de_res_path <- file.path(cfg$paths$tcga_chr_dir, chr, gene, "de_res.RDS")
res    <- readRDS(res_path)
de_res <- readRDS(de_res_path)

fact_idx <- as.integer(sub("^k", "", fact))
pass <- (de_res$svalue[, fact_idx]   < cfg$thresholds$sval_thresh
         & de_res$postmean[, fact_idx] > cfg$thresholds$lfc_thresh
         & res$F[, fact_idx]           > cfg$thresholds$abs_thresh)
if (!any(pass)) {
  stop("No DE junctions for ", gene, " ", fact, " under current thresholds.")
}
ranked    <- order(de_res$postmean[, fact_idx] * pass, decreasing = TRUE)
top_junc  <- rownames(res$F)[ranked[1L]]
top_post  <- de_res$postmean[ranked[1L], fact_idx]
cat(sprintf("[psi_box] top DE junction: %s (postmean = %.3f)\n",
            top_junc, top_post))

# ---- PSI table + match the top junction -----------------------------------

psi_path <- file.path(cfg$paths$tcga_psi_dir, chr, gene, "psi_table.tsv")
if (!file.exists(psi_path)) stop("PSI table missing: ", psi_path)
psi <- fread(psi_path, sep = "\t", data.table = FALSE)
rownames(psi) <- psi$V1
psi$V1 <- NULL

cols    <- names(psi)
hits    <- cols[endsWith(cols, paste0(":", top_junc))]
if (length(hits) == 0L) {
  stop("Top junction ", top_junc, " not found as a column in the PSI table.")
}
col_key <- hits[[1L]]

psi_df <- data.frame(
  rail_id = as.integer(rownames(psi)),
  psi     = suppressWarnings(as.numeric(psi[[col_key]]))
)
psi_df <- psi_df[is.finite(psi_df$psi), ]

# ---- Join metadata, bucket read length, restrict to plot cohorts ----------

md <- read_tcga_meta(cfg$paths$tcga_metadata)
md <- md[md$project_id %in% cohorts, ]
md$cohort   <- sub("^TCGA-", "", md$project_id)
md$read_len <- round(suppressWarnings(as.numeric(md$avg_len)))
md$rl_group <- dplyr::case_when(
  md$read_len <= 50L ~ "<= 50 bp",
  md$read_len >= 75L ~ ">= 75 bp",
  TRUE               ~ NA_character_
)
md <- md[!is.na(md$rl_group), c("rail_id", "cohort", "rl_group")]

plot_df <- merge(psi_df, md, by = "rail_id")
plot_df$rl_group <- factor(plot_df$rl_group, levels = c("<= 50 bp", ">= 75 bp"))
plot_df$cohort   <- factor(plot_df$cohort,   levels = sub("^TCGA-", "", cohorts))

n_per <- aggregate(rail_id ~ rl_group, plot_df, length)
cat("[psi_box] per-group sample counts:\n")
print(n_per)
if (nrow(plot_df) == 0L) {
  stop("No samples left after read-length bucket + cohort filter.")
}

# ---- Box plot --------------------------------------------------------------

p <- ggplot(plot_df, aes(x = rl_group, y = psi)) +
  geom_boxplot(width = 0.5, fill = "white", colour = "black",
                outlier.colour = "black", outlier.size = 1) +
  labs(title    = sprintf("%s (%s) top DE junction: %s",
                          gene, fact, top_junc),
       subtitle = sprintf("read-length stratified PSI (%d samples total)",
                          nrow(plot_df)),
       x        = "Sample read length",
       y        = "PSI") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

ggplot2::ggsave(out_pdf, p, width = 7, height = 4.5, device = cairo_pdf)
cat("[psi_box] wrote ", out_pdf, "\n", sep = "")
