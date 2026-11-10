#!/usr/bin/env Rscript
# Build the per-showcase outputs (no GTEx, multi-cohort):
#   junctions_classified.tsv  DE junctions + Leaf2 classification
#   junction_arcs_mane.pdf    factor arcs over the MANE gene model
#   structure_plot.pdf        fastTopics structure plot over the chosen cohorts
#   gene_data.RDS             cached gd for the per-region rules

source(file.path(dirname(sub("--file=", "",
       commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))])),
       "lib.R"))

args <- parse_args()
sample_id   <- require_arg(args, "sample")
config_path <- require_arg(args, "config")
out_dir     <- require_arg(args, "out_dir")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cfg        <- yaml::read_yaml(config_path)
paths      <- cfg$paths
thresholds <- cfg$thresholds
entry      <- cfg$showcase[[sample_id]]
if (is.null(entry)) stop(sprintf("No showcase entry for sample '%s'", sample_id))

gene        <- entry$gene
chr         <- entry$chr
factor_id   <- entry$factor
cohorts     <- entry$cohorts

cat(sprintf("[make_gene_outputs] %s | %s | %s | cohorts=%s\n",
            gene, chr, factor_id, paste(cohorts, collapse = ",")))

gd <- load_gene_data(
  gene_name  = gene,
  chr        = chr,
  fact       = factor_id,
  cohorts    = cohorts,
  paths      = paths,
  thresholds = thresholds
)

cat(sprintf("[make_gene_outputs] %d DE junctions, %d cohort samples\n",
            length(gd$junctions_to_compare), nrow(gd$cohort_lookup)))

# 1. junctions tsv with Leafcutter2 classifications
jc <- junction_classifications(gd)
write.table(jc, file.path(out_dir, "junctions_classified.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# 2. MANE-overlay arc plot (PDF)
arcs_p <- plot_factor_junctions_static(gd)
ggplot2::ggsave(file.path(out_dir, "junction_arcs_mane.pdf"),
       arcs_p, width = 10, height = 4, device = cairo_pdf)

# 3. Structure plot over the specified cohorts (PDF)
struc_p <- plot_structure_static(gd)
ggplot2::ggsave(file.path(out_dir, "structure_plot.pdf"),
       struc_p, width = 10, height = 3, device = cairo_pdf)

# 4. cached gd for the per-region rules
saveRDS(gd, file.path(out_dir, "gene_data.RDS"))

cat("[make_gene_outputs] done\n")
