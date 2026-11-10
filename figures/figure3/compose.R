# Figure 3 (NKAP / GSP #13) assembly.
#
# Panels:
#   a  somatic-mutation enrichment CDF + zoom  -- flash_fact_char.Rmd
#   b  NKAP / PRPF8 domain lollipops          -- NKAP_full_plus_variants
#   c  four mis-splicing loci, sashimi + zoom -- pygenome/NKAP
#   d  cryptic 3'ss displacement + pie        -- diff_splicing/unified_results
#   e  expression ECDF by consequence         -- combined_expression_analysis.Rmd
#   f  structure screenshot (raster, as-is)   -- extracted from the reference PDF
#
#   /scratch/midway2/qhauck/conda_env/rstudio-server/bin/Rscript compose.R
#
# Same conventions as ../figure2/compose.R -- read ../figure2/compose.md
# before changing layout, it documents the traps (absolute-inch ROW_H,
# scale limits vs coord, title alignment by measurement).

.libPaths('/scratch/midway2/qhauck/conda_env/rstudio-server/lib/R/library/')

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
# Provisional, measured off the reference render at 100 dpi:
#   a|b 2.55 in, c 3.35 in, d|e 2.45 in, f 2.45 in
# ------------------------------------------------------------------
PAGE_W <- 8.5
# Row 1 is 0.40 in taller than the reference's ink band: panel a's zoom
# labels 47 mutations in two columns, and at 2.55 in the left column was
# full -- the PRPF8/FAM50A labels at its top had nowhere to spread and
# overprinted. Stretch the plot rather than shrink the type.
# Row 3 is 0.50 in taller than the reference's ink band: panel e labels 16 genes
# in three ladders, and at 2.45 in they only fitted at 7 pt against the figure's
# 8 pt SIZE_BASE. Quinn's rule is stretch the plot rather than shrink the type.
# ROW_H is in ABSOLUTE INCHES, so this leaves rows 1, 2 and 4 untouched.
# THREE rows, not four. Panel f used to sit alone on row 4 as a centred
# square with white space either side; Quinn: "arrange 3 so that f is not
# on its own at the bottom, make stuff fit together in a rectangle".
# f now joins d and e on row 3, so the page is a filled rectangle.
# ROW_H is in ABSOLUTE INCHES; the page is their sum.
# 2026-10-05 (Quinn): f is new full-width artwork on its own bottom row;
# d | e share row 3 as two wider columns, so row 3 drops 3.05 -> 2.75 in.
# Row 4 (f) is derived from the art's aspect after panel f is built
# (F_ROW_H below); rows 1-2 are unchanged -- a, b and c carry geometry
# measured against them.
# 2026-11-09 (Quinn): row order is a|b, structure (panel_f/, tagged c), d|e, tracks (panel_c/, tagged f),
# to match the text's citation order. Rows are NAMED so panels look
# their own row up by name (panel_c uses ROW_H[["c"]]), not by position.
ROW_H  <- c(ab = 2.95, de = 2.75, c = 3.35)

# Panel a widened (Quinn: "3a should be wider"); panel b takes the
# difference and restructures its legend to suit.
# Row 3 carries d | e | f: 2.98 / 3.23 / 2.30 in. f's artwork is 2.30 in
# wide so it sets the right column. BOTH d and e must be re-tuned for
# these widths -- their geometry constants were measured against the old
# two-column row (d 4.42 in, e 3.40 in) and go stale here. That
# staleness is the known fragility of this design; see build.md files.
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
# Row 4 height: f's art at the full page width (less its 6+6 pt side
# margins) at its own aspect, plus the title band (title pt x 1.2 + the
# 4 pt title margin + 4/2 pt plot margins).
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
