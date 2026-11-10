# Supplementary figure: SVM method for estimating phenotype prevalence.
# SVM-derived phenotype-axis scores vs GSP loading for GSPs #10, #13 and
# #30. Labelled positives / negatives trained the decision boundary
# (dotted blue, score = 0); intermediates were scored by it.
# (Panel b, the C* mutation rate near the boundary, was dropped 2026-10-26;
# that test lives in figure3/paper_stats.Rmd only.)
# Source: diff_splicing/svm_plots.Rmd; data: quetzal_code phenotype-axis
# caches + NKAP variant table (vt_path).
if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
source("common.R")

DS <- file.path(QROOT, "leafcutter_analyses/diff_splicing")
COL <- c(P = "#D62728", intermediate = "#F0A400", N = "grey65")
LAB <- c(P = "positives", intermediate = "intermediates", N = "negatives")

panel <- function(gene, fact, ylab = TRUE) {
  d <- fread(vt_path(paste0(tolower(gene), "_phenotype_axis_scores")),
             data.table = FALSE)
  d$l <- log(d[[paste0("k", fact)]])
  d <- d %>% filter(is.finite(l)) %>%
    mutate(group = factor(group, levels = c("N", "intermediate", "P"))) %>%
    arrange(group)                          # positives drawn last, on top
  ggplot(d, aes(l, score, colour = group)) +
    geom_hline(yintercept = 0, colour = "steelblue", linetype = "dotted",
               linewidth = 0.5) +
    geom_point(size = 0.9, alpha = 0.75, stroke = 0) +
    scale_colour_manual(values = COL, labels = LAB, name = NULL,
                        breaks = c("P", "intermediate", "N")) +
    labs(x = paste0("log(GSP #", fact, ")"), y = if (ylab) "Phenotype-axis SVM score" else NULL,
         title = gene_md(paste0("GSP #", fact, " (", gene, ")"))) +
    theme_publication() +
    theme(plot.title = md_title()) +
    guides(colour = guide_legend(override.aes = list(size = 2.2, alpha = 1)))
}
# y title once (left panel); each panel keeps its own y scale.
p_a <- (panel("SF3B1", 10) | panel("NKAP", 13, FALSE) | panel("YLPM1", 30, FALSE)) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.box.spacing = unit(4, "pt"),
        legend.margin = margin(0, 0, 0, 0))

fig <- p_a
save_supp(fig, "figS20", 8.5, 3.4)
