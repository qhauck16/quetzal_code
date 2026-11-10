# Figure 3 (NKAP / GSP #13) assembly.
#
# Panels:
#   a  somatic-mutation enrichment CDF + zoom
#   b  NKAP / PRPF8 domain lollipops
#   c  four mis-splicing loci, sashimi + zoom
#   d  cryptic 3'ss displacement + pie
#   e  expression ECDF by consequence
#   f  C* spliceosome structure renders (raster)
#
#   Rscript compose.R              # reuse cached panels/*.rds
#   Rscript compose.R a c          # rebuild only panels a and c
#   Rscript compose.R --all        # rebuild every panel


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
# Row heights in INCHES; PAGE_H is their sum. Absolute heights keep a
# change to one row from rescaling the text geometry of the others.
# ------------------------------------------------------------------
PAGE_W <- 8.5
# Rows 1 (a|b) and 3 (d|e) are sized so panel a's zoom labels and panel
# e's label ladders fit at SIZE_BASE: stretch the plot, not the type.
# Row order is a|b, structure (panel_f/), d|e, tracks (panel_c/), to match
# the text's citation order. Rows are NAMED so panels look their own
# height up by name (panel_c uses ROW_H[["c"]]); the structure row's
# height (F_ROW_H) is derived from the art's aspect below.
ROW_H  <- c(ab = 2.95, de = 2.75, c = 3.35)

# Column fractions per row. Panels d and e carry width-dependent geometry
# constants; re-tune them if COL_W$de changes.
COL_W <- list(ab = c(0.58, 0.42), c = 1, de = c(0.47, 0.53))

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
for (.p in c('a', 'b', 'c', 'd', 'e', 'f')) {
  src <- file.path(paste0('panel_', .p), 'build.R')
  if (file.exists(src)) source(src, local = FALSE)
  else assign(paste0('panel_', .p),
              missing_panel(paste(.p, " (no build.R yet)"), src))
}

# ------------------------------------------------------------------
# Compose.
# ------------------------------------------------------------------
# Structure-row height: f's art at the full page width (less its 6+6 pt
# side margins) at its own aspect, plus the title band (title pt x 1.2 +
# the 4 pt title margin + 4/2 pt plot margins).
F_DIM   <- dim(png::readPNG(file.path("raster_cache", "panel_f_improved.png")))
F_ROW_H <- (PAGE_W - 12 / 72) * F_DIM[1] / F_DIM[2] +
           (SIZE_TITLE * 1.2 + 4 + 6) / 72
ROW_H   <- c(ROW_H, f = F_ROW_H)
PAGE_H  <- sum(ROW_H)
cat(sprintf("[compose] row heights %s in; page %.2f in\n",
            paste(round(ROW_H, 2), collapse = " / "), PAGE_H))

figure3 <- ((panel_a | panel_b) + plot_layout(widths = COL_W$ab)) /
           panel_f /
           ((panel_d | panel_e) + plot_layout(widths = COL_W$de)) /
           panel_c +
  plot_layout(heights = unname(ROW_H[c("ab", "f", "de", "c")])) +
  plot_annotation_publication()

out_pdf <- file.path("out", "figure3.pdf")
cat("[compose] writing", out_pdf, "\n")
ggsave_pub(out_pdf, figure3, width = PAGE_W, height = PAGE_H,
                units = "in")
cat("[compose] done.\n")
