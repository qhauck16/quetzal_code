#!/usr/bin/env Rscript
# Structure plot for a supplementary_plots gene. For read-length genes the
# grouping is cohort + read length (e.g. COAD-48bp); for tcga_cohort genes
# the grouping is the TCGA cohort short name.

script_dir <- dirname(sub("--file=", "",
  commandArgs(trailingOnly = FALSE)[
    grep("--file=", commandArgs(trailingOnly = FALSE))]))
source(file.path(script_dir, "lib_supp.R"))

args        <- commandArgs(trailingOnly = TRUE)
gene        <- args[1]
config_path <- args[2]
out_pdf     <- args[3]

dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)

cfg <- yaml::read_yaml(config_path)
le  <- lookup_entry(cfg, gene)
ent <- le$entry; group <- le$group

chr     <- ent$chr
fact    <- if (is.null(ent$factor)) "" else as.character(ent$factor)
cohorts <- as.character(unlist(ent$cohorts %||% list()))
if (length(cohorts) == 0L) {
  stop("config['",
       group, "']['", gene, "']['cohorts'] is empty; fill it in")
}

cat(sprintf("[structure_plot] %s | %s | group=%s | cohorts=%s\n",
            gene, chr, group, paste(cohorts, collapse = ",")))

res_path <- file.path(cfg$paths$tcga_chr_dir, chr, gene, "res.RDS")
res      <- readRDS(res_path)

tcga_md <- read_tcga_meta(cfg$paths$tcga_metadata)
md_sub <- tcga_md %>%
  dplyr::filter(rail_id %in% as.integer(rownames(res$L)),
                project_id %in% cohorts) %>%
  dplyr::mutate(group_label = make_group_labels(., group))

if (nrow(md_sub) == 0L) {
  stop("No samples matched cohorts ", paste(cohorts, collapse = ","))
}

# Subset every per-sample component of the fit consistently.
keep    <- rownames(res$L) %in% as.character(md_sub$rail_id)
res_sub <- res
res_sub$L  <- res_sub$L[keep, , drop = FALSE]
if (!is.null(res_sub$Ln)) res_sub$Ln <- res_sub$Ln[keep, , drop = FALSE]
if (!is.null(res_sub$s) && length(res_sub$s) == length(keep)) {
  res_sub$s <- res_sub$s[keep]
}

grp <- md_sub$group_label[
  match(rownames(res_sub$L), as.character(md_sub$rail_id))
]
# Stable ordering: keep cohort order from config, then numeric read length.
if (identical(group, "read_length")) {
  cohort_order <- sub("^TCGA-", "", cohorts)
  parse_cohort <- function(lbl) sub("-\\d+bp$", "", lbl)
  parse_rl     <- function(lbl) as.integer(sub(".*-(\\d+)bp$", "\\1", lbl))
  lvl <- unique(grp)
  lvl <- lvl[order(match(parse_cohort(lvl), cohort_order),
                    parse_rl(lvl))]
  grp <- factor(grp, levels = lvl)
} else {
  grp <- factor(grp, levels = sub("^TCGA-", "", cohorts))
}

num_factors <- ncol(res_sub$F)
pal         <- factor_palette(num_factors)

title <- if (nzchar(fact)) {
  paste0("TCGA factorisation of ", gene, " (", fact, ")")
} else {
  paste0("TCGA factorisation of ", gene)
}

struc_p <- fastTopics::structure_plot(
  res_sub,
  grouping = grp,
  gap      = 100,
  colors   = pal,
  n        = 5000
) +
  theme(aspect.ratio = 1 / 5) +
  labs(title = title)

ggplot2::ggsave(out_pdf, struc_p,
                 width = 12, height = 3, device = cairo_pdf)
cat("[structure_plot] wrote ", out_pdf, "\n", sep = "")
