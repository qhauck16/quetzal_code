#!/usr/bin/env Rscript
# Per-region prep for pyGenomeTracks: load the cached gd RDS and write
# junctions.bed, sample_bigwigs.tsv, groups_file.txt, region.txt into out_dir.
#
# An empty --region "" causes the default (auto-computed from gd) region to
# be used.

source(file.path(dirname(sub("--file=", "",
       commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))])),
       "lib.R"))

args <- parse_args()
gd_rds  <- require_arg(args, "gd_rds")
out_dir <- require_arg(args, "out_dir")
region  <- args[["region"]]            # may be NULL or ""

cfg   <- yaml::read_yaml(require_arg(args, "config"))
paths <- cfg$paths

gd <- readRDS(gd_rds)
# Backfill paths in case gd was generated under a different config invocation.
gd$paths <- paths

cat(sprintf("[prep_pygenome] %s | region='%s'\n", gd$gene_name,
            ifelse(is.null(region), "<auto>", region)))

res <- prep_pygenome_inputs(gd, out_dir = out_dir, region = region)

# MANE gene body PDF, with the resolved region highlighted.
mane_p <- plot_mane_gene_model_static(gd, region = res$region)
ggplot2::ggsave(file.path(out_dir, "mane_gene_model.pdf"),
                mane_p, width = 8, height = 1, device = cairo_pdf)

cat("[prep_pygenome] done -> ", out_dir, "\n")
