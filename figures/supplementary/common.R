# Shared setup for every supplementary figure. Each figure lives in
# figS<N>/build.R, sources this file, and writes out/figS<N>.pdf.
#
#   cd quetzal_code/figures/supplementary
#   /scratch/midway2/qhauck/conda_env/rstudio-server/bin/Rscript --vanilla figS7/build.R
#   ./build_all.sh            # every figure
#
# Same conventions as the main figures (../figure2/compose.md): FONT_FAMILY
# and SIZE_* tokens from ../theme.R, gene symbols italic via gene_md() /
# gene_pm(), x tick labels angled 45 deg, page sizes in absolute inches.
.libPaths('/scratch/midway2/qhauck/conda_env/rstudio-server/lib/R/library/')

if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(patchwork)
  library(cowplot)
  library(ggrepel)
  library(scales)
  library(RColorBrewer)
})
select <- dplyr::select; rename <- dplyr::rename; filter <- dplyr::filter

FIG_ROOT <- file.path(QROOT, "figures")
source(file.path(FIG_ROOT, "theme.R"))
source(file.path(FIG_ROOT, "scripts", "coverage_plot.R"))
source(file.path(FIG_ROOT, "scripts", "pie_chart.R"))
OUT_DIR  <- file.path(FIG_ROOT, "supplementary", "out")
dir.create(OUT_DIR, showWarnings = FALSE)

# Write a supplementary figure via ggsave_pub() (../theme.R: svglite ->
# cairosvg, exact glyph spacing). Text measurement still uses cairo (the device every text
# measurement in this project is calibrated against).
save_supp <- function(p, name, width, height) {
  f <- file.path(OUT_DIR, paste0(name, ".pdf"))
  ggsave_pub(f, p, width = width, height = height, units = "in")
  cat("[supp] wrote", f, "\n")
  invisible(f)
}

# 2 x 2 exact-overlap heatmap used by S6 and S7: rows = driving GSP
# junctions / RPM-PSI-matched controls, columns = up / down in the
# reference set; each cell "count\n(percent of its row set)". `cells`
# has set (gsp | matched), direction, n_overlap, n_set, fisher_p.
overlap_heatmap <- function(cells, gsp_label, title = NULL, xlab = NULL,
                            pct_digits = 0) {
  rows <- c(gsp = gsp_label, matched = "RPM/PSI matched\njunctions")
  cells <- cells %>%
    mutate(frac = n_overlap / pmax(n_set, 1),
           row  = factor(rows[set], levels = rev(rows)),
           direction = factor(direction, levels = c("upregulated", "downregulated")),
           lab = sprintf(paste0("%d\n(%.", pct_digits, "f%%)"), n_overlap, 100 * frac))
  p <- cells$fisher_p[1]
  ggplot(cells, aes(direction, row, fill = frac)) +
    geom_tile(colour = "white", linewidth = 1.2) +
    geom_text(aes(label = lab, colour = frac > 0.55),
              size = SIZE_BASE / .pt, family = FONT_FAMILY, lineheight = 0.95) +
    scale_colour_manual(values = c(`FALSE` = "grey10", `TRUE` = "white"), guide = "none") +
    scale_fill_gradient(low = "#DEEBF7", high = "#08519C", limits = c(0, 1),
                        labels = scales::percent_format(accuracy = 1),
                        name = "Overlap") +
    scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
    labs(x = xlab, y = NULL, title = title,
         caption = if (p < 1e-3) sprintf("p = %.1e", p) else sprintf("p = %.2g", p)) +
    theme_publication() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
          plot.caption = element_text(size = SIZE_BASE * 0.9, hjust = 1,
                                      family = FONT_FAMILY),
          plot.title = element_text(size = SIZE_BASE, face = "bold", hjust = 0))
}
