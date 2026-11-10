# Figure 4 (YLPM1 / GSP #30) assembly.
#
# Panels:
#   a  YLPM1 LoF enrichment CDF + zoom   -- flash_fact_char.Rmd
#   b  deep-intronic splice sites, 2 loci -- coverage_tracks/fig4b_tracks (KDM4C, TRIM38)
#   c  intron-feature beta heatmap        -- diff_splicing/unified_results
#   d  LeafCutter event-type pie          -- unified_results/ylpm1/20v60/dpsi0.01_fc2.5
#   e  four mis-spliced tumour suppressors -- pygenome/YLPM1
#                                            (TP53BP1, PBRM1, MSH6, ASXL1)
#   f  GIST PSI heatmap + dendrogram      -- diff_splicing/ylpm1/leafcutter_prep_and_analysis.Rmd
#   g  expression ECDF + bootstrap ribbon -- combined_expression_analysis.Rmd
#
# The d/e/g panels reuse patterns already solved in figures 2 and 3
# (pie-with-leaders, sashimi quadrants, expression ECDF) -- start from
# those rather than from the Rmds.
#
#   Rscript compose.R
#
# Same conventions as ../figure2/compose.R -- read ../figure2/compose.md
# before changing layout, it documents the traps (absolute-inch ROW_H,
# scale limits vs coord, title alignment by measurement).


# ------------------------------------------------------------------
# 0. Setup
# ------------------------------------------------------------------
HERE <- normalizePath(dirname(sub("--file=", "",
    commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))])),
    mustWork = FALSE)
if (!nzchar(HERE) || !dir.exists(HERE)) HERE <- getwd()
setwd(HERE)

source("../theme.R")
source("../scripts/raster_utils.R")
source("../scripts/coverage_plot.R")
source("../scripts/pie_chart.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(patchwork)
  library(cowplot)
  library(ggrepel)
  library(scales)
  library(RColorBrewer)
})
select <- dplyr::select
rename <- dplyr::rename
filter <- dplyr::filter

dir.create("out",          showWarnings = FALSE)
dir.create("raster_cache", showWarnings = FALSE)
CACHE    <- "cache"

have_cache <- function(...) all(file.exists(file.path(CACHE, c(...))))

# ------------------------------------------------------------------
# Row heights in INCHES; PAGE_H is their sum. Relative heights let one
# row's change silently rescale every other row's derived text
# geometry -- see ../figure2/compose.md.
# Measured off the reference render at 100 dpi:
#   a|b 2.50 in, c|d 2.00 in, e 3.50 in, f|g 3.00 in
# ------------------------------------------------------------------
PAGE_W <- 8.5
# Row 1 2.50 -> 2.85 in: KDM4C has four arcs sharing one foot, and at 2.50
# the top two (4.87, 1.92) ran almost together (Quinn: make them distinct).
# 2026-11-09 (Quinn): citation order. Rows: a|b, panel_e/ (full width, tagged c), panel_d/ | panel_c/
# (flipped, tagged d | e), f|g. Rows and columns are NAMED by source panel so each keeps its size.
ROW_H  <- c(ab = 2.85, cd = 2.00, e = 3.50, fg = 3.00)
PAGE_H <- sum(ROW_H)

COL_W <- list(ab = c(1.15, 1), cd = c(c = 1.15, d = 1), e = 1, fg = c(1.09, 1))

PANEL_DIR <- "panels"
dir.create(PANEL_DIR, showWarnings = FALSE)
.args    <- commandArgs(trailingOnly = TRUE)
REBUILD  <- if (any(.args == "--all")) "ALL" else .args

panel_cached <- function(name, expr) {
  f <- file.path(PANEL_DIR, paste0(name, ".rds"))
  if (!identical(REBUILD, "ALL") && !(name %in% REBUILD) && file.exists(f)) {
    cat("[panel]", name, "(cached)\n")
    return(readRDS(f))
  }
  cat("[panel]", name, "building ...\n")
  p <- force(expr)
  saveRDS(p, f)
  p
}

missing_panel <- function(what, script) {
  ggplot() +
    annotate("text", x = 0, y = 0.6, label = what,
             size = SIZE_TITLE / .pt, fontface = "bold",
             family = FONT_FAMILY, color = "grey35") +
    annotate("text", x = 0, y = 0.2,
             label = paste0("not built yet - ", script),
             size = SIZE_BASE / .pt, family = FONT_FAMILY, color = "grey55") +
    coord_cartesian(xlim = c(-1, 1), ylim = c(0, 1)) +
    theme_void() +
    theme(panel.border = element_rect(color = "grey85", fill = NA,
                                       linetype = "dashed"))
}

# ------------------------------------------------------------------
# Panels. One file each; each must leave panel_<x> defined.
# ------------------------------------------------------------------
for (.p in c('a', 'b', 'c', 'd', 'e', 'f', 'g')) {
  src <- file.path(paste0('panel_', .p), 'build.R')
  if (file.exists(src)) source(src, local = FALSE)
  else assign(paste0('panel_', .p),
              missing_panel(paste(.p, " (no build.R yet)"), src))
}

# ------------------------------------------------------------------
# Compose.
# ------------------------------------------------------------------
figure4 <- ((panel_a | panel_b) + plot_layout(widths = COL_W$ab)) /
           panel_e /
           ((panel_d | panel_c) + plot_layout(widths = unname(COL_W$cd[c("d", "c")]))) /
           ((panel_f | panel_g) + plot_layout(widths = COL_W$fg)) +
  plot_layout(heights = unname(ROW_H[c("ab", "e", "cd", "fg")])) +
  plot_annotation_publication()

out_pdf <- file.path("out", "figure4.pdf")
cat("[compose] writing", out_pdf, "\n")
ggsave_pub(out_pdf, figure4, width = PAGE_W, height = PAGE_H,
                units = "in")
cat("[compose] done.\n")
