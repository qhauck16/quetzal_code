#!/usr/bin/env Rscript
# Build the per-(sample, gene) outputs:
#   junctions_classified.tsv       DE junctions + Leaf2 classification + glm stats
#   junction_arcs_mane.pdf         factor arcs over the MANE gene model
#   psi_cdfs_tcga_vs_gtex.pdf      4-cohort ECDF per DE junction
#   gene_data.RDS                  cached gd for the per-region rules

source(file.path(dirname(sub("--file=", "",
       commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))])),
       "lib.R"))

args <- parse_args()
gene        <- require_arg(args, "gene")
chr         <- require_arg(args, "chr")
factor_id   <- require_arg(args, "factor")
cancer_type <- require_arg(args, "cancer_type")
config_path <- require_arg(args, "config")
out_dir     <- require_arg(args, "out_dir")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cfg        <- yaml::read_yaml(config_path)
paths      <- cfg$paths
thresholds <- cfg$thresholds

cat(sprintf("[make_gene_outputs] %s | %s | %s | %s\n",
            gene, chr, factor_id, cancer_type))

gd <- load_gene_data(
  gene_name   = gene,
  chr         = chr,
  fact        = factor_id,
  cancer_type = cancer_type,
  paths       = paths,
  thresholds  = thresholds
)

cat(sprintf("[make_gene_outputs] %d DE junctions\n",
            length(gd$junctions_to_compare)))

# 1. junctions tsv with Leafcutter2 classifications
jc <- junction_classifications(gd)
write.table(jc, file.path(out_dir, "junctions_classified.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# 2. MANE-overlay arc plot (PDF)
arcs_p <- plot_factor_junctions_static(gd)
ggplot2::ggsave(file.path(out_dir, "junction_arcs_mane.pdf"),
       arcs_p, width = 10, height = 4, device = cairo_pdf)

# 3. PSI ECDFs (4-cohort)
n_juncs <- length(gd$junctions_to_compare)
ncol_cdf <- min(2, max(1, n_juncs))
cdfs_p <- plot_psi_cdfs_full_static(gd, ncol = ncol_cdf)
n_rows <- ceiling(max(1, n_juncs) / ncol_cdf)
ggplot2::ggsave(file.path(out_dir, "psi_cdfs_tcga_vs_gtex.pdf"),
       cdfs_p, width = 5 * ncol_cdf, height = 4 * n_rows,
       device = cairo_pdf, limitsize = FALSE)

# 4. Factor-loading violin (highlighted cohort vs Other)
violin_p <- plot_loading_violins(gd)
ggplot2::ggsave(file.path(out_dir, "loading_violin.pdf"),
       violin_p, width = 4, height = 4, device = cairo_pdf)

# 5. cached gd for the per-region rules
saveRDS(gd, file.path(out_dir, "gene_data.RDS"))

cat("[make_gene_outputs] done\n")
