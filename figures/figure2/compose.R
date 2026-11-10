# Figure 2 (SF3B1 / genome-wide GSPs) assembly.
#
# Panels:
#   a  genome-wide workflow schematic  -- vector art cropped from the
#      Illustrator reference (title re-set in R so it matches b-g)
#   c  GSP #1 neuronal-junction ECDF   -- panel_b/ (cache/panel_b_running_psi.tsv); shown as "c" since 2026-11-09
#   b  GSP x covariate correlation heatmap -- panel_c/ (cache/panel_c_*.tsv); shown as "b" since 2026-11-09
#   d  GSP #10 vs SF3B1 hotspot mutations  -- SF3B1_full_info.tsv
#   e  cryptic 3'ss distance + event-type pie -- fact_10_top_junctions.tsv
#   f  PSI fold-change vs binned GSP #10   -- sf3b1_top500_alt3_fc5_psi.tsv
#   g  mutation enrichment scatter     -- cache/panel_g_*.tsv
#
# The three cache/ inputs come from scripts/build_cache_panel_{b,c,g}.R,
# which are expensive (each loads a 3.5 G flashier RDS). Run those once;
# this script then rebuilds the figure in seconds.
#
#   Rscript compose.R


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
source("../scripts/pdf_font_subst.R")
source("../scripts/coverage_plot.R")

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
# Per-panel build cache.
#
# Rebuilding all seven panels takes ~60 s, most of it in panels d/f
# (which read a 19 M TSV) and the heatmap. When iterating on ONE
# panel's layout that is wasted time, so each panel's built object is
# serialised to panels/<name>.rds and reused unless asked for.
#
#   Rscript compose.R            # reuse every cached panel
#   Rscript compose.R c a        # rebuild only c and a
#   Rscript compose.R --all      # rebuild everything
#
# A panel is also rebuilt automatically if its .rds is missing. Note
# the serialised object carries its data, so a panel will NOT pick up
# changed upstream data until you name it (or --all).
# ------------------------------------------------------------------
# Row heights, in INCHES (patchwork treats them as relative, but the
# page height below is their sum, so they are literal). Absolute units
# matter here because several panels derive text-collision spacing from
# their own drawn height -- with relative units, growing one row
# silently rescaled every other row's geometry constants.
#
# ROW_H[1] is tuned so panel a's CELL aspect matches its art's 5.83:1
# (filling a 4.90:1 cell stretched the schematic's text ~19%).
# Row 3 is 0.70 in taller than the reference's ink band: panel d's zoom
# labels 52 hotspot mutations, and every class tops out at the same
# 100% row, so that band needs vertical room the reference bought in
# Illustrator by shrinking the type.
PAGE_W <- 8.5
ROW_H  <- c(1.5714, 3.9370, 3.2346, 2.9570)
PAGE_H <- sum(ROW_H)

# Column split within each row, also used by panels that must convert
# label widths into data units.
# Row 2 is ordered heatmap | ECDF since 2026-11-09 (the heatmap is cited first, so it is tagged b). Named so
# each panel's width stays tied to the panel, not to its slot: panel_b/ and panel_c/ index by name.
# Rows 3-4 reordered 2026-11-09 to follow citation order: row 3 = mutation scatter (panel_g/, tagged d) |
# GSP #10 CDF (panel_d/, tagged e; width unchanged); row 4 = 3'ss histogram + pie (panel_e/, tagged f; widened,
# it was cramped) | serial boxplots (panel_f/, tagged g).
COL_W <- list(bc = c(heatmap = 1.50, ecdf = 1.00),
              de = c(scatter = 1, cdf = 1.75),
              fg = c(hist = 1.30, box = 1))

# Panels b and e are drawn at this fraction of their row's height, with the
# remainder split as whitespace above and below (they read as too tall when
# filling the row next to c / d). panel_b/build.R reads PANEL_H_SCALE["b"]
# for its inset geometry, so change it here only.
PANEL_H_SCALE <- c(b = 0.85, e = 0.85)
v_shrink <- function(p, s) {
  wrap_plots(plot_spacer(), p, plot_spacer(), ncol = 1,
             heights = c((1 - s) / 2, s, (1 - s) / 2))
}

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

# A placeholder so the figure still composes while an expensive cache
# is still building -- keeps the layout reviewable end-to-end.
missing_panel <- function(what, script) {
  ggplot() +
    annotate("text", x = 0, y = 0.6, label = what,
             size = SIZE_TITLE / .pt, fontface = "bold",
             family = FONT_FAMILY, color = "grey35") +
    annotate("text", x = 0, y = 0.2,
             label = paste0("cache missing - run scripts/", script),
             size = SIZE_BASE / .pt, family = FONT_FAMILY, color = "grey55") +
    coord_cartesian(xlim = c(-1, 1), ylim = c(0, 1)) +
    theme_void() +
    theme(panel.border = element_rect(color = "grey85", fill = NA,
                                       linetype = "dashed"))
}


# ------------------------------------------------------------------
# Panels. One file each so they can be worked on independently; each
# defines panel_<x>. Order matters only for d before f (shared data).
# ------------------------------------------------------------------
sf3b1_info <- read_tsv(vt_path("SF3B1_full_info"),
                        show_col_types = FALSE)
sf3b1_info$proteinChange <- str_split_i(sf3b1_info$gene_mutation, ':', 2)

for (.p in c('a', 'd', 'e', 'f', 'b', 'g', 'c')) {
  source(file.path(paste0('panel_', .p), 'build.R'), local = FALSE)
}

# ------------------------------------------------------------------
# Compose. Row heights are taken from the reference figure's own ink
# bands (a 235 px, b/c 547, d/e 296, f/g 412 at 150 dpi).
# ------------------------------------------------------------------
figure2 <- panel_a /
           # heatmap (built in panel_c/) first -> auto-tagged "b"; ECDF + WDR11 inset (panel_b/) -> "c"
           ((panel_c | v_shrink(panel_b, PANEL_H_SCALE["b"])) +
              plot_layout(widths = COL_W$bc)) /
           ((panel_g | panel_d) + plot_layout(widths = COL_W$de)) /
           ((panel_e | panel_f) + plot_layout(widths = COL_W$fg)) +
  plot_layout(heights = ROW_H) +
  plot_annotation_publication()

out_pdf <- file.path("out", "figure2.pdf")
cat("[compose] writing", out_pdf, "\n")
ggsave_pub(out_pdf, figure2, width = PAGE_W, height = PAGE_H,
                units = "in")
cat("[compose] done.\n")

