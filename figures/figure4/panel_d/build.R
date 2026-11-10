# Figure 4, panel d. Sourced by ../compose.R. Must leave `panel_d`.
#
# Rebuild just this panel:   Rscript compose.R d
#
# d. LeafCutter results show upregulation of both alt 3'ss and 5'ss.
#    A single pie of splice-event types over every differential
#    junction, percentages ON the wedges.
#
# Instance: ylpm1 / 20v60 / dpsi0.01_fc2.5. Confirmed by its row count
# (2460 = the published "N = 2460") and by its type breakdown
# reproducing the published percentages to the decimal:
#   alt 3'ss 30.3, alt 5'ss 50.1, canonical 1.1, other 12.8, skipped 5.7
# NB Quinn's recollection was dPSI > 0.05; the data says 0.01. Same
# threshold family as figure 3 panel d.

if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
D_DIR <- file.path(QROOT, "leafcutter_analyses/diff_splicing/unified_results/ylpm1/20v60/dpsi0.01_fc2.5")

panel_d <- panel_cached("d", local({
  df <- fread(file.path(D_DIR, "df_final.tsv"), data.table = FALSE)
  if (nrow(df) != 2460) message("[panel d] df_final has ", nrow(df), " rows (published: 2460)")

  d_pie <- df %>%
    count(type) %>%
    mutate(prop  = n / sum(n),
           label = sprintf("%.1f%%", prop * 100))

  d_lvl <- sort(unique(d_pie$type))
  d_pal <- setNames(RColorBrewer::brewer.pal(5, "Set2"), d_lvl)

  # Drawn radius, MEASURED off the render. coord_equal makes the pie
  # HEIGHT-limited, and the panel's drawn height is only settled by
  # patchwork at layout time, so this cannot be computed from the box
  # width -- doing that in figure 3 gave a 2x error and collapsed the
  # label fan. Re-measure if ROW_H or the box below changes.
  # Re-measured after the labels moved onto the wedges: the plot's
  # limits are now computed from where the labels land, so the drawn
  # radius moved with them (0.62 -> 0.645 in). Measuring is a fixed
  # point -- a bigger constant makes the label boxes smaller in pie
  # radii, which shrinks the limits and grows the drawn pie -- so
  # re-render and re-measure until the number stops moving. It took
  # two passes here (0.685 -> 0.646 -> 0.645).
  D_PIE_R_IN <- 0.580

  # Percentages sit ON the wedges; the two slivers (5.7 / 1.1%) get a
  # boxed label straddling the rim on their own mid-ray. See
  # ../scripts/pie_chart.R.
  pie <- pie_labels_on_wedges(d_pie, "type", d_pal,
                              levels_order = d_lvl,
                              pie_r_in = D_PIE_R_IN,
                              label_pt = SIZE_BASE * 0.85,
                              family = FONT_FAMILY)
  pie_leg <- pie_legend_grob(d_pie, "type", d_pal, levels_order = d_lvl,
                             label_pt = SIZE_BASE * 0.85,
                             family = FONT_FAMILY)

  # Title lives on a void plot behind the pie so the pie itself keeps
  # coord_equal (a title on the pie would be pushed around by it).
  # An EMPTY ggplot has no scale range, so an annotate() on it centres
  # itself instead of honouring its coordinates -- coord_cartesian
  # pins a 0..1 frame to place against. labs(x/y = NULL) kills the
  # stray "x"/"y" axis titles theme_publication would otherwise draw.
  frame <- ggplot() +
    labs(title = paste0("LeafCutter results show upregulation\n",
                        "of both alt 3'ss and 5'ss"),
         x = NULL, y = NULL) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_publication() +
    theme(plot.title   = element_text(hjust = 0.5),
          axis.line    = element_blank(),
          axis.ticks   = element_blank(),
          axis.text    = element_blank(),
          panel.border = element_blank(),
          plot.margin  = margin(4, 6, 4, 4, "pt"))

  wrap_elements(full =
    frame +
      # Centred under the pie (whose inset spans x 0.02-0.66), not in
      # the cell's bottom-left corner, where it read as belonging to
      # panel c's colourbar (QA review).
      annotate("text", x = 0.34, y = 0.0, hjust = 0.5, vjust = 0,
               label = paste0("N = ", nrow(df)),
               size = SIZE_BASE * 0.85 / .pt, family = FONT_FAMILY,
               color = "grey25") +
      inset_element(pie, left = 0.02, bottom = 0.06,
                    right = 0.66, top = 1.00, align_to = "panel") +
      inset_element(wrap_elements(full = pie_leg),
                    left = 0.62, bottom = 0.22,
                    right = 1.00, top = 0.86, align_to = "panel"))
}))
