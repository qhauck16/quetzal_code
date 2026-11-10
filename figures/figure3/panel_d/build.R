# Figure 3, panel d. Sourced by ../compose.R. Must leave `panel_d`.
#
# Rebuild just this panel:   Rscript compose.R d
#
# d. GSP #13 phenotype captures consistent displacement of canonical
#    3'ss downstream. Histogram of cryptic-3'ss position relative to
#    the annotated 3'ss, with an event-type pie inset.
#
# Input: df_final.tsv from the leafcutter_analyses/diff_splicing/
# unified_results pipeline, nkap 20v60 dpsi0.01_fc2.5 instance (2806
# rows); the histogram and pie follow its scripts/plot_3ss_barplots.R
# (zoom branch) and scripts/plot_proportions.R.

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
  # cryptic-relative-to-annotated, and the axis runs the other way
  # ("position of cryptic 3'ss relative to annotated").
  d_hist <- df %>%
    filter(type == "alt 3'ss", !is.na(distance_to_cryptic)) %>%
    mutate(pos = -distance_to_cryptic) %>%
    filter(pos >= D_ZOOM[1], pos <= D_ZOOM[2])

  # ---- pie ---------------------------------------------------------
  # Over ALL rows, not just alt 3'ss.
  d_pie <- df %>%
    count(type) %>%
    mutate(prop  = n / sum(n),
           label = sprintf("%.2f%%", prop * 100))

  d_lvl <- sort(unique(d_pie$type))
  d_pal <- setNames(RColorBrewer::brewer.pal(5, "Set2"), d_lvl)

  # Shared helper: cartesian wedges, every percentage ON its own wedge;
  # small wedges get boxed chips placed in 2-D so they fan out without
  # covering other wedges. See ../scripts/pie_chart.R.
  #
  # D_PIE_R_IN is the drawn pie radius (in). It cannot be computed from
  # the inset's width: coord_equal makes the pie HEIGHT-limited and the
  # drawn height is only settled by patchwork at layout time, so it is a
  # measured value. Re-measure if D_PIE_BOX, ang_tol_deg or ROW_H changes.
  D_PIE_BOX  <- c(left = 0.04, bottom = 0.16, right = 0.46, top = 0.99)   # clear of the y axis
  D_PIE_R_IN <- 0.6450
  # Angular span of every wedge that gets a boxed chip (all but alt 3'ss).
  D_SMALL_DEG <- 360 * sum(d_pie$prop[d_pie$type != "alt 3'ss"])
  pie <- pie_labels_on_wedges(d_pie, "type", d_pal,
                              levels_order = d_lvl,
                              pie_r_in = D_PIE_R_IN,
                              label_pt = SIZE_BASE * 0.75,
                              # Chips sit ON the pie, filled with their
                              # wedge's colour; wedges start at 12 o'clock.
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
    # Both titles are wrapped to fit the cell width (the axis title is
    # centred on the panel, not the plot).
    labs(x = "Position of cryptic 3'ss\nrelative to annotated 3'ss",
         y = "Count",
         title = paste0("GSP #13 phenotype captures\n",
                        "consistent displacement of\n",
                        "canonical 3'ss downstream")) +
    theme_publication() +
    theme(plot.title  = element_text(hjust = 0.5),
          axis.title  = element_text(size = SIZE_BASE),
          # Zero top margin keeps the title level with panels e and f.
          plot.margin = margin(0, 6, 4, 4, "pt")) +
    angled_x_labels(size = SIZE_BASE * 0.8)

  # Pie and its legend are inset into the histogram's empty left half
  # (all the mass is at pos >= 3). The legend is placed separately so
  # it centres at its natural size instead of stealing pie width; it
  # sits in the headroom above the bars.
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
