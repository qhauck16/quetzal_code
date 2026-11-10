#!/usr/bin/env Rscript
# DE-junction arcs over the MANE gene body. Reuses the parent
# showcase_genes/scripts/lib.R::load_gene_data + ::plot_factor_junctions_static
# so the figure stays identical in style to the main showcase pipeline.

script_dir <- dirname(sub("--file=", "",
  commandArgs(trailingOnly = FALSE)[
    grep("--file=", commandArgs(trailingOnly = FALSE))]))
# Parent lib.R brings in load_gene_data, plot_factor_junctions_static,
# factor_palette, and all required libraries (tidyverse, ggbio, ...).
source(file.path(script_dir, "..", "..", "scripts", "lib.R"))
# lib_supp.R adds lookup_entry, %||%, and the supp-specific helpers.
source(file.path(script_dir, "lib_supp.R"))

args        <- commandArgs(trailingOnly = TRUE)
gene        <- args[1]
config_path <- args[2]
out_pdf     <- args[3]

dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)

cfg <- yaml::read_yaml(config_path)
le  <- lookup_entry(cfg, gene)
ent <- le$entry

chr     <- ent$chr
fact    <- if (is.null(ent$factor)) "" else as.character(ent$factor)
cohorts <- as.character(unlist(ent$cohorts %||% list()))

if (!nzchar(fact)) {
  stop("config[", le$group, "][", gene, "]: `factor` must be set.")
}
if (length(cohorts) == 0L) {
  stop("config[", le$group, "][", gene, "]: `cohorts` must be set.")
}

cat(sprintf("[arcs_mane] %s | %s | %s | cohorts=%s\n",
            gene, chr, fact, paste(cohorts, collapse = ",")))

gd <- load_gene_data(
  gene_name  = gene,
  chr        = chr,
  fact       = fact,
  cohorts    = cohorts,
  paths      = cfg$paths,
  thresholds = cfg$thresholds
)
cat(sprintf("[arcs_mane] %d DE junctions\n",
            length(gd$junctions_to_compare)))

p <- plot_factor_junctions_static(gd)
ggplot2::ggsave(out_pdf, p, width = 10, height = 4, device = cairo_pdf)
cat("[arcs_mane] wrote ", out_pdf, "\n", sep = "")
