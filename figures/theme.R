# Shared design tokens for ALL paper figures (figure1/, figure2/, ...).
# Each figure's compose.R does `source("../theme.R")`. Edit here to
# change fonts, colours, sizes or panel-label styling once and have it
# propagate across every figure -- that single source is what keeps
# typography consistent between figures.
.libPaths('/scratch/midway2/qhauck/conda_env/rstudio-server/lib/R/library/')

if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

# ------------------------------------------------------------------
# Typography
# ------------------------------------------------------------------
# Nimbus Sans = URW metric clone of Helvetica. Cairo/fontconfig on midway
# does NOT alias "Helvetica" to it -- it fell back to DejaVuSans in every
# PDF until 2026-09-24 -- so name the face directly.
FONT_FAMILY <- "Nimbus Sans"

# Point sizes for on-page text (converts to ggplot2's ~mm via /.pt)
SIZE_BASE        <- 8              # default axis text
SIZE_AXIS_TITLE  <- 9
SIZE_STRIP       <- 9
SIZE_TITLE       <- 9
SIZE_LEGEND      <- 8
SIZE_PANEL_LABEL <- 10             # a, b, c, d, e (small on purpose)

# ------------------------------------------------------------------
# Color palette (change these once, propagates everywhere)
# ------------------------------------------------------------------
PAL <- list(
  # Panel 1b recovery-curve series
  leafcutter = "#d95f02",
  matched    = "#7570b3",
  # Cohort accents. The TCGA highlight color is deliberately shared
  # across violins and CDFs so panel 1d TCGA-LAML looks the same in
  # both. Same idea for GTEx.
  tcga_test  = "#e39f4f",    # TCGA-LAML / TCGA-LUSC highlight (orange)
  tcga_other = "grey70",     # pooled "Other TCGA" (gray)
  gtex       = "#3d8f4a",    # GTEx-BLOOD / GTEx-LUNG (green)
  # Factor-index -> arc fill. Keyed by fact_idx used in plot_factor_junctions_static.
  # BDH1 k4 = red (like the reference); EZH2 k6 = navy; NF1 k1 = pink.
  factor = c(
    "1"  = "#d33f6a",         # NF1 k1
    "2"  = "#fc8d62",
    "3"  = "#8da0cb",
    "4"  = "#d13030",         # BDH1 k4 (red)
    "5"  = "#a6d854",
    "6"  = "#132b73",         # EZH2 k6 (navy)
    "7"  = "#e5c494",
    "8"  = "#b3b3b3",
    "9"  = "#1b9e77",
    "10" = "#d95f02"
  ),

  # ---- Figure 2 series ------------------------------------------
  # Panel b: cohort ECDF curves. Brain cohorts get the saturated
  # hues, "other" the muted ones, so the brain/non-brain split reads
  # before the GTEx/TCGA split.
  cohort = c(
    "GTEx-brain" = "#1f77b4",
    "TCGA-brain" = "#ff7f0e",
    "GTEx-other" = "#8c564b",
    "TCGA-other" = "grey55"
  ),
  # Panel d: SF3B1 mutation classes (shared by the CDF and the dot
  # plot so a class looks identical in both).
  sf3b1_mut = c(
    "SF3B1 hotspot mutation"      = "#C0392B",
    # Was #A6CEE3: too pale to read as label TEXT on white (QA review).
    "Mutation Near SF3B1 Hotspot" = "#5B9BD5",
    "Other SF3B1 Mutation"        = "#2E8B57",
    "No SF3B1 Mutation"           = "#34495E"
  ),
  # Panel e: splicing event classes (ColorBrewer Set2, pinned by name
  # so the mapping can't shift if a class is absent).
  splice_type = c(
    "alt 3'ss"     = "#66C2A5",
    "alt 5'ss"     = "#FC8D62",
    "other"        = "#8DA0CB",
    "skipped exon" = "#E78AC3"
  ),
  # Panels e/f bar + box fill
  fill_neutral = "steelblue"
)

# Diverging ramp for the panel c correlation heatmap (RdBu, reversed
# so positive r is red).
heatmap_diverging <- function() rev(RColorBrewer::brewer.pal(11, "RdBu"))

# ------------------------------------------------------------------
# Gene-structure arrow density
# ------------------------------------------------------------------
# ggbio::geom_alignment's gap.geom controls the intron ornament:
#   "arrow"   = strand arrows every ~few pixels (busy but directional)
#   "chevron" = light chevrons (cleaner but no direction)
#   NULL      = no ornament
GAP_GEOM <- "arrow"

# ------------------------------------------------------------------
# Publication theme -- apply as `p + theme_publication()`.
# ------------------------------------------------------------------
theme_publication <- function(base_size = SIZE_BASE) {
  theme_bw(base_size = base_size, base_family = FONT_FAMILY) +
    theme(
      # kill visual noise -- no gridlines anywhere per publication style
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_blank(),
      panel.border      = element_blank(),
      axis.line         = element_line(color = "grey20", linewidth = 0.35),
      axis.ticks        = element_line(color = "grey20", linewidth = 0.30),
      # text sizing
      plot.title        = element_text(size = SIZE_TITLE,       face = "bold",
                                        family = FONT_FAMILY),
      plot.subtitle     = element_text(size = SIZE_BASE,        family = FONT_FAMILY),
      axis.title        = element_text(size = SIZE_AXIS_TITLE,  family = FONT_FAMILY),
      axis.text         = element_text(size = SIZE_BASE,        family = FONT_FAMILY),
      strip.text        = element_text(size = SIZE_STRIP,       face = "bold",
                                        family = FONT_FAMILY),
      strip.background  = element_rect(fill = "grey95", color = NA),
      legend.title      = element_text(size = SIZE_LEGEND,      family = FONT_FAMILY),
      legend.text       = element_text(size = SIZE_LEGEND,      family = FONT_FAMILY),
      legend.key.size   = unit(3.2, "mm"),
      legend.margin     = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      # tighter margins so panels can be packed
      plot.margin       = margin(4, 6, 4, 6, unit = "pt"),
      plot.tag          = element_text(size = SIZE_PANEL_LABEL, face = "bold",
                                        family = FONT_FAMILY)
    )
}

# Angled x tick labels. THE DEFAULT whenever x labels are too long or
# too dense to sit horizontally -- use this rather than rotating them
# to vertical. Vertical labels force the reader to tilt their head and
# eat a lot of height; 45 degrees stays readable at a glance.
# `size` defaults to NULL so the calling theme's axis.text size wins.
angled_x_labels <- function(angle = 45, size = NULL) {
  theme(axis.text.x = element_text(angle = angle, hjust = 1, vjust = 1,
                                    size = size))
}

# Empty theme for the gene-track / arcs / structure sub-panels where
# axes just chew whitespace.
theme_track <- function() {
  theme_publication() +
    theme(
      panel.border   = element_blank(),
      panel.grid     = element_blank(),
      axis.text.y    = element_blank(),
      axis.ticks.y   = element_blank(),
      axis.title.y   = element_blank(),
      # No axis lines at all -- the horizontal x-axis line was reading
      # visually as a "separator" between the arcs and the gene body.
      axis.line      = element_blank(),
      axis.ticks.x   = element_line(color = "grey40", linewidth = 0.3)
    )
}

# Convenience -- pull a factor color from PAL$factor with a fallback.
factor_color <- function(k) {
  key <- as.character(k)
  if (key %in% names(PAL$factor)) unname(PAL$factor[key]) else "grey60"
}

# Panel-label styling for patchwork (applied at compose time)
plot_annotation_publication <- function() {
  plot_annotation(
    tag_levels = "a",
    theme = theme(
      plot.tag          = element_text(size = SIZE_PANEL_LABEL,
                                        face = "bold",
                                        family = FONT_FAMILY),
      plot.tag.position = c(0, 1)
    )
  )
}

# patchwork styles each tag from ITS OWN panel's theme, so the
# plot_annotation(theme=) above never reaches wrap_elements() panels or
# nested patchworks -- every "a" rendered plain until 2026-09-24. Set the
# tag style (and the text family, so no stray DejaVuSans) on the global
# default theme and geom defaults too; only these elements are touched.
theme_update(text     = element_text(family = FONT_FAMILY),
             plot.tag = element_text(size = SIZE_PANEL_LABEL, face = "bold",
                                     family = FONT_FAMILY))
for (.g in c("text", "label")) update_geom_defaults(.g, list(family = FONT_FAMILY))
rm(.g)

# ------------------------------------------------------------------
# Gene symbols in italics (Quinn: "Everywhere we say YLPM1 or NKAP or
# any gene name should be italicized across all figures").
#
# A title like "YLPM1 LoF confers ..." is MIXED roman/italic, which
# element_text cannot do, so such strings go through gene_md() and are
# drawn with a ggtext markdown element instead. This is applied ONLY
# where a gene name occurs, not globally: element_markdown does not
# honour "\n" (markdown folds it to a space), so switching every title
# to markdown would silently unwrap every multi-line title in the paper.
# gene_md() converts "\n" to "<br>" for exactly that reason.
#
# Uses HTML <i>, not markdown *...*, because mutation strings carry
# literal asterisks ("Q437*") that markdown would read as emphasis.
# Those are escaped here too.
#
# AFTER switching a title to md_title(), RE-MEASURE the row's title
# alignment (scripts/title_align.R, strip cropped from the row's top
# edge): gridtext sets line height slightly differently from grid, so
# a markdown title and a plain one in the same row can drift.
# ------------------------------------------------------------------
# ggtext must be ATTACHED, not just referenced as ggtext::element_markdown.
# A panel built fresh loads the namespace as a side effect of that call,
# but a panel read back from panels/*.rds never calls it -- and without
# the namespace loaded, ggtext's element_grob method is not registered,
# so ggplot falls back to plain element_text and prints the raw markup
# ("<i>NKAP</i>, <br>...") in the title. That is exactly what figure 3
# panel a's title did after a cached rebuild.
suppressPackageStartupMessages(library(ggtext))

GENE_SYMBOLS <- c("NKAP", "PRPF8", "PRP8", "FAM50A", "SLU7", "YLPM1",
                  "SF3B1", "WDR11", "IDH1", "SPOP", "NF1", "BDH1",
                  "SENP2", "PRKDC", "EZH2", "BRD9", "MSH6", "TP53BP1",
                  "PBRM1", "ASXL1", "KDM4C", "TRIM38", "CD44", "EGFR",
                  "SUGP1", "KAT5", "SBNO2", "TCF12", "SETD2", "SRC", "MVK",
                  "NHEJ1", "MAP2K2", "KLHDC4")

gene_md <- function(x, genes = GENE_SYMBOLS) {
  x <- gsub("\\*", "\\\\*", x)
  x <- gsub("_",   "\\\\_", x)
  # Longest first so PRPF8 is not half-matched as PRP8. The lookarounds
  # stop a match inside a longer token (e.g. "NF1" inside "SNF1").
  for (g in genes[order(-nchar(genes))])
    x <- gsub(paste0("(?<![A-Za-z0-9>])", g, "(?![A-Za-z0-9<])"),
              paste0("<i>", g, "</i>"), x, perl = TRUE)
  gsub("\n", "<br>", x, fixed = TRUE)
}

# Markdown twins of theme_publication()'s text elements. Same size,
# face and family; halign carries the multi-line centring that hjust
# alone gives element_text.
md_title <- function(hjust = 0.5, size = SIZE_TITLE, ...)
  ggtext::element_markdown(size = size, face = "bold",
                           family = FONT_FAMILY, lineheight = 0.9,
                           hjust = hjust, halign = hjust, ...)
md_text <- function(size = SIZE_BASE, ...)
  ggtext::element_markdown(size = size, family = FONT_FAMILY, ...)

# Evaluate `expr` with a cairo_pdf device current, so any text measured
# while grobs are BUILT uses the same metrics the figure is drawn with.
# Needed wherever a legend or other grob is extracted up front (e.g.
# cowplot::get_plot_component) and contains markdown text: gridtext lays
# out words at build time, and under the default device's narrower AFM
# metrics the words were placed too close and printed over each other
# once cairo drew them ("frameshift/early stop" -> "early" and "stop"
# overlapping in figure 4 panel a's legend).
with_cairo_metrics <- function(expr) {
  old <- grDevices::dev.cur()
  tf  <- tempfile(fileext = ".pdf")
  grDevices::cairo_pdf(tf, width = 8.5, height = 11)
  on.exit({
    grDevices::dev.off()
    if (old > 1) try(grDevices::dev.set(old), silent = TRUE)
    unlink(tf)
  }, add = TRUE)
  force(expr)
}

# Plotmath twin of gene_md(), for LEGEND labels. Returns a list of
# expressions with each gene symbol wrapped in italic(). Legends use
# this rather than markdown: a gridtext legend label measured its words
# in one font and drew them in another, so in figure 4 panel a
# "early" and "stop" printed on top of each other. Plotmath is laid out
# and drawn by grid itself, so measurement and drawing always agree.
gene_pm <- function(x, genes = GENE_SYMBOLS) {
  pat <- paste0("(?<![A-Za-z0-9])(",
                paste(genes[order(-nchar(genes))], collapse = "|"),
                ")(?![A-Za-z0-9])")
  lapply(x, function(s) {
    parts <- regmatches(s, gregexpr(pat, s, perl = TRUE), invert = NA)[[1]]
    parts <- parts[nzchar(parts)]
    if (!length(parts)) return(s)
    isg <- grepl(paste0("^", substring(pat, 1)), parts, perl = TRUE) &
           parts %in% genes
    # A space right after an italic gene collapses in the PDF (the slant eats
    # it); carry it inside the italic token, as figure 4a's legend does.
    nxt <- c(parts[-1], "")
    sp  <- isg & c(!isg[-1], FALSE) & startsWith(nxt, " ")
    parts[isg & sp] <- paste0(parts[isg & sp], " ")
    after <- c(FALSE, sp[-length(sp)])
    parts[after] <- sub("^ ", "", parts[after])
    keep <- nzchar(parts); parts <- parts[keep]; isg <- isg[keep]
    tok <- ifelse(isg, sprintf('italic("%s")', parts),
                  sprintf('"%s"', gsub('"', '\\\\"', parts)))
    str2lang(if (length(tok) == 1) tok else
             paste0("paste(", paste(tok, collapse = ", "), ")"))
  })
}

# Final-output PDF device. Plain cairo_pdf sets plotmath operators and
# spaces (the "=" and " " in italic(p)[adj] == ...) in the device's
# *symbol* family, which fontconfig resolves to DejaVuSans; point that
# at FONT_FAMILY too so the whole page is one typeface.
cairo_pdf_pub <- function(...) {
  grDevices::cairo_pdf(...,
    symbolfamily = grDevices::cairoSymbolFont(FONT_FAMILY, usePUA = FALSE))
}

# Final-output writer for every figure (2026-10-05). cairo_pdf lays text
# out through pango 1.50, which rounds every glyph position to a whole
# point: at 6-10 pt that opened visible gaps inside words ("Mutati ons",
# "Enconuc lease"). svglite writes each text run with its exact position
# and width and leaves glyph layout to the renderer; cairosvg then draws
# it through cairo's own text API (no pango, no rounding) and embeds
# Nimbus Sans. Same page size, fonts, plotmath and Unicode as before.
# Text MEASUREMENT for layout (with_cairo_metrics) is unchanged.
CAIROSVG <- Sys.getenv("CAIROSVG", unname(Sys.which("cairosvg")))   # cairosvg binary (set CAIROSVG if not on PATH)
ggsave_pub <- function(filename, plot, width, height, units = "in", ...) {
  svg <- tempfile(fileext = ".svg")
  on.exit(unlink(svg), add = TRUE)
  ggplot2::ggsave(svg, plot, width = width, height = height, units = units,
                  device = svglite::svglite, bg = "white",
                  system_fonts = list(sans = FONT_FAMILY, serif = FONT_FAMILY,
                                      mono = FONT_FAMILY, symbol = FONT_FAMILY), ...)
  st <- system2(CAIROSVG, c(shQuote(svg), "-o", shQuote(filename)))
  if (st != 0 || !file.exists(filename)) stop("cairosvg failed for ", filename)
  invisible(filename)
}


# vt_path("<table>"): home of each per-gene variant table (see the file).
source(file.path(QROOT, "figures/scripts/vt_paths.R"))
