# Figure 3, panel d. Sourced by ../compose.R. Must leave `panel_d`.
#
# Rebuild just this panel:   Rscript compose.R d
#
# d. GSP #13 phenotype captures consistent displacement of canonical
#    3'ss downstream. Histogram of cryptic-3'ss position relative to
#    the annotated 3'ss, with an event-type pie inset.
#
# Source: the diff_splicing/unified_results Snakemake pipeline --
#   scripts/plot_3ss_barplots.R  (zoom branch) and
#   scripts/plot_proportions.R
# both reading ONE file, df_final.tsv, for the nkap 20v60
# dpsi0.01_fc2.5 instance. That instance is identified by its row count
# (2806) and its type breakdown, which reproduce the published
# percentages exactly.
#
# NB the published title says "GSP #20". That is a stale Illustrator
# leftover from an older factor key (spec_flash_fact.Rmd still carries
# 'NKAP'=20); NKAP is k13 everywhere in current code and panel a of the
# same figure says #13. Quinn confirmed #13 is correct.

if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
D_DIR <- file.path(QROOT, "leafcutter_analyses/diff_splicing/unified_results/nkap/20v60/dpsi0.01_fc2.5")
D_ZOOM <- c(-8, 20)   # unified_results/config.yaml: barplot_zoom.nkap

panel_d <- panel_cached("d", local({
  df <- fread(file.path(D_DIR, "df_final.tsv"), data.table = FALSE)
  if (nrow(df) != 2806) message("[panel d] df_final has ", nrow(df), " rows (published: 2806)")

  # ---- histogram --------------------------------------------------
  # `pos` is the NEGATED distance: distance_to_cryptic is measured
  # cryptic-relative-to-annotated, and the published axis runs the
  # other way ("position of cryptic 3'ss relative to annotated").
  d_hist <- df %>%
    filter(type == "alt 3'ss", !is.na(distance_to_cryptic)) %>%
    mutate(pos = -distance_to_cryptic) %>%
    filter(pos >= D_ZOOM[1], pos <= D_ZOOM[2])

  # ---- pie ---------------------------------------------------------
  # Over ALL rows, not just alt 3'ss -- that is what makes 88.06%.
  d_pie <- df %>%
    count(type) %>%
    mutate(prop  = n / sum(n),
           label = sprintf("%.2f%%", prop * 100))

  d_lvl <- sort(unique(d_pie$type))
  d_pal <- setNames(RColorBrewer::brewer.pal(5, "Set2"), d_lvl)

  # Shared helper: cartesian wedges, every percentage ON its own wedge.
  # Four of the five wedges here are under 6% and three of them are
  # adjacent, so all four get a boxed label. The helper solves those in
  # TWO DIMENSIONS (angle as well as radius) and forbids a chip from
  # touching the pie disc, so the fan spreads across the free space
  # instead of stacking, and no chip covers a wedge it does not name.
  # See ../scripts/pie_chart.R.
  #
  # Drawn radius, MEASURED off the render at 400 dpi. It cannot be
  # computed from the inset's width: coord_equal makes the pie
  # HEIGHT-limited, and the panel's drawn height is only settled by
  # patchwork at layout time. Deriving it from the width once gave
  # 0.73 in against a real 0.31 in, which collapsed the label fan.
  # Measuring is a FIXED POINT -- a larger constant makes the label
  # boxes smaller in pie radii, which shrinks the computed limits and
  # grows the drawn pie -- so re-render and re-measure until it stops
  # moving. History: 0.31 (row 3 at 2.45 in) -> 0.475 -> 0.485 (row 3
  # restructured to d|e|f) -> 0.432 (2-D chip placement; the fan now
  # spends width as well as height, so the composition the pie has to
  # share its box with is wider and the pie itself lands smaller).
  # Convergence is not smooth -- the helper's placement search picks a
  # different discrete arrangement as the constant moves -- so take the
  # value that REPRODUCES ITSELF, not the one nearest your last guess.
  # 0.432 -> 0.4333 measured; 0.44 -> 0.4367; 0.425 -> 0.4133.
  # Re-measure if D_PIE_BOX, ang_tol_deg or ROW_H changes. The measuring
  # script is a saturation scan for the widest coloured run:
  #   sat > 0.12 & grey < 0.95, longest TRUE run per row, max over rows.
  #
  # bottom 0.30 -> 0.18: the 2-D fan made the whole composition wider,
  # which shrank the (height-limited) pie from 0.85 to 0.75 in across.
  # Giving the box the extra 0.12 of panel height puts it back to
  # 0.87 in, larger than it ever was. The two bars at pos 3-4 sit under
  # the box but are 14 counts tall, so the pie still clears them.
  D_PIE_BOX  <- c(left = 0.04, bottom = 0.16, right = 0.46, top = 0.99)   # 2026-10-05: off the y axis in the wider 2-column row
  D_PIE_R_IN <- 0.6450
  # Angular span of every wedge that gets a boxed chip (all but alt 3'ss).
  D_SMALL_DEG <- 360 * sum(d_pie$prop[d_pie$type != "alt 3'ss"])
  pie <- pie_labels_on_wedges(d_pie, "type", d_pal,
                              levels_order = d_lvl,
                              pie_r_in = D_PIE_R_IN,
                              label_pt = SIZE_BASE * 0.75,
                              # Chips sit ON the pie, filled with their
                              # wedge's colour -- the published panel's
                              # look (see ../scripts/pie_chart.R). No
                              # rotation: wedges start at 12 o'clock as
                              # in the reference.
                              family = FONT_FAMILY)
  pie_leg <- pie_legend_grob(d_pie, "type", d_pal, levels_order = d_lvl,
                             label_pt = SIZE_BASE * 0.75, family = FONT_FAMILY,
                             key_mm = 2.6)

  hist <- ggplot(d_hist, aes(x = pos)) +
    geom_histogram(breaks = seq(D_ZOOM[1] - 0.5, D_ZOOM[2] + 0.5, by = 1),
                   fill = PAL$fill_neutral, color = "white",
                   linewidth = 0.2) +
    scale_x_continuous(breaks = seq(D_ZOOM[1], D_ZOOM[2], by = 2)) +
    coord_cartesian(xlim = D_ZOOM) +
    # Both titles are wrapped for a 2.98 in CELL (row 3 became d|e|f).
    # At 2.836 in of drawn width the old one-line axis title (2.71 in
    # under cairo's Helvetica, and it is centred on the PANEL, not the
    # plot, so only ~2.20 in is usable) and the old two-line title both
    # ran off the cell's right edge into panel e.
    labs(x = "Position of cryptic 3'ss\nrelative to annotated 3'ss",
         y = "Count",
         title = paste0("GSP #13 phenotype captures\n",
                        "consistent displacement of\n",
                        "canonical 3'ss downstream")) +
    theme_publication() +
    theme(plot.title  = element_text(hjust = 0.5),
          axis.title  = element_text(size = SIZE_BASE),
          # ZERO top margin. Measured at 400 dpi with a strip that
          # starts ABOVE the row's titles, d/e/f sit at 22/21/21 px.
          #
          # HOW TO MEASURE THIS, because I got it wrong once: the strip
          # handed to ../scripts/title_align.R must begin above EVERY
          # title in the row. If it starts inside them, every cell
          # reports "ink from row 1" and the row looks perfectly
          # aligned when it is not -- that false pass is exactly what
          # hid d 75 / e 21 / f 86 (a 65 px, 0.16 in spread that Quinn
          # spotted by eye). Crop from the row's top edge:
          #   row 3 starts at sum(ROW_H[1:2]) in; splits at
          #   PAGE_W * cumsum(COL_W$def)/sum(COL_W$def).
          plot.margin = margin(0, 6, 4, 4, "pt")) +
    angled_x_labels(size = SIZE_BASE * 0.8)

  # Pie and its legend are inset into the histogram's empty left half
  # (all the mass is at pos >= 3). The legend is placed separately so
  # it centres at its natural size instead of stealing pie width; it
  # sits in the headroom above the bars, which top out at ~56% of the
  # y range. Its type dropped to 6 pt and its key to 2.6 mm to fit the
  # 2.98 in cell -- at 6.4 pt / 3.0 mm "skipped exon" pushed the legend
  # block past the panel's right edge.
  wrap_elements(full =
    hist +
      annotate("text", x = D_ZOOM[1] + 0.5, y = 0, hjust = 0, vjust = -0.4,
               label = paste0("N = ", nrow(df)),
               size = SIZE_BASE * 0.8 / .pt, family = FONT_FAMILY,
               color = "grey25") +
      inset_element(pie, left = D_PIE_BOX["left"],
                    bottom = D_PIE_BOX["bottom"],
                    right = D_PIE_BOX["right"], top = D_PIE_BOX["top"],
                    align_to = "panel") +
      inset_element(wrap_elements(full = pie_leg),
                    left = 0.56, bottom = 0.56,
                    right = 0.90, top = 0.99, align_to = "panel"))
}))
