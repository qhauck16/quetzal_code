# ------------------------------------------------------------------
# Figure 4, panel g -- "YLPM1 LoF confers widespread gene expression
# downregulation through mis-splicing".
#
# Empirical CDF of mean scaled RPM across the top-20 YLPM1 (GSP #30)
# samples for the 604 genes carrying a Leaf2-unproductive junction,
# against a 100-iteration stratified-random-sample bootstrap of the
# SAME gene set (mean + 95% band).
#
# All of the arithmetic lives in ../scripts/build_cache_panel_g.R; this
# file only reads cache/panel_g_*.{tsv,rds} and lays the panel out.
# See panel_g/build.md for the analysis and the layout reasoning.
#
# Sourced by ../compose.R, which supplies theme.R's tokens, the
# tidyverse, ggrepel, CACHE and panel_cached().
# ------------------------------------------------------------------

panel_g <- panel_cached("g", local({

  # ---------------- constants -------------------------------------
  G_XLIM <- c(-1.30, 1.30)      # view only -- NEVER a scale limit
  # 1.10 -> 1.16: the three-line legend in the top-left came down to the
  # 1.00 plateau, "with GSP #30 < 1e-4" nearly touching the ribbon (QA).
  G_YLIM <- c(-0.04, 1.16)
  G_LW   <- 0.70
  G_SEED <- 20260922L

  # Panel geometry, MEASURED off the composed render at 400 dpi (row 4
  # right cell): 1.174 in per x-unit, 1.723 in per y-unit. Everything
  # positional below is derived from these, so re-measure if compose.R's
  # ROW_H / COL_W change.
  G_IN_PER_X <- 1.125
  G_IN_PER_Y <- 1.610
  G_EXPAND   <- 0.05            # ggplot's default continuous mult
  g_pan_x <- G_XLIM + c(-1, 1) * G_EXPAND * diff(G_XLIM)
  g_pan_y <- G_YLIM + c(-1, 1) * G_EXPAND * diff(G_YLIM)

  # Series names double as the legend text, so the wrapping lives here
  # and nothing has to re-map labels later. The bootstrap entry is
  # wrapped to two lines: unwrapped it is 2.4 in wide and there is no
  # 2.4 in x 0.3 in empty pocket anywhere in this panel.
  G_TOP  <- "Top 20 Samples on GSP #30"
  # Was "...from samples with GSP #30 < 1e-4": the bootstrap never applies that filter (it draws cohort-matched
  # samples from all of TCGA), so the label just names the curve (Quinn, 2026-11-09).
  G_BOOT <- "Mean bootstrapped CDF"
  G_PAL  <- c("#B22222", "grey30");  names(G_PAL) <- c(G_TOP, G_BOOT)
  G_LTY  <- c("solid", "dashed");    names(G_LTY) <- c(G_TOP, G_BOOT)

  # Legend block size, measured with ggplotGrob at SIZE_LEGEND (in).
  G_LEG_W <- 1.400 / G_IN_PER_X
  G_LEG_H <- 0.375 / G_IN_PER_Y
  G_LEG_NPC_X <- 0.005          # must match legend.position.inside
  g_leg <- list(xmin = g_pan_x[1] + G_LEG_NPC_X * diff(g_pan_x),
                ymax = g_pan_y[2])
  g_leg$xmax <- g_leg$xmin + G_LEG_W
  g_leg$ymin <- g_leg$ymax - G_LEG_H

  # ---------------- data ------------------------------------------
  red    <- data.table::fread(file.path(CACHE, "panel_g_red.tsv"))
  ribbon <- data.table::fread(file.path(CACHE, "panel_g_ribbon.tsv"))
  hl     <- data.table::fread(file.path(CACHE, "panel_g_highlights.tsv"))
  meta   <- readRDS(file.path(CACHE, "panel_g_meta.rds"))

  # The series name has to travel as a DATA COLUMN, not as a scalar
  # captured from this local() -- panel_cached() round-trips the plot
  # through saveRDS and tidy-eval would then fail to find it, rendering
  # the whole figure blank. Same reason `p_lab` becomes a column below.
  red$grp    <- G_TOP
  ribbon$grp <- G_BOOT

  cat("[panel g] n target =", meta$n_target, " ", meta$p_label,
      "| labelled", meta$n_highlight, "of 10\n")

  # ---------------- repel obstacles -------------------------------
  # ggrepel only avoids its OWN layer's data, so neither the step curve,
  # the dashed mean nor the grey band exists as far as it is concerned.
  # Feed them all in as zero-text rows so a gene name cannot be parked
  # on top of them.
  red_curve <- red |>
    dplyr::filter(is.finite(x)) |>
    dplyr::arrange(x) |>
    dplyr::mutate(ecdf_y = dplyr::row_number() / dplyr::n()) |>
    dplyr::filter(x > G_XLIM[1], x < G_XLIM[2])
  obstacles <- dplyr::bind_rows(
    red_curve |>
      dplyr::slice(unique(round(seq(1, dplyr::n(), length.out = 150)))) |>
      dplyr::transmute(gene_name = "", x, ecdf_y),
    ribbon |>
      dplyr::slice(unique(round(seq(1, dplyr::n(), length.out = 90)))) |>
      dplyr::select(x, lower, mean_cdf, upper) |>
      tidyr::pivot_longer(-x, values_to = "ecdf_y") |>
      dplyr::transmute(gene_name = "", x, ecdf_y))

  # The legend is invisible to ggrepel too; tile its rectangle.
  obstacles <- dplyr::bind_rows(obstacles, expand.grid(
    x      = seq(g_leg$xmin, g_leg$xmax, length.out = 11),
    ecdf_y = seq(g_leg$ymin, g_leg$ymax, length.out = 4)) |>
    dplyr::transmute(gene_name = "", x, ecdf_y))

  # ---------------- label bands -----------------------------------
  # One column of names hugging the LEFT flank of the red curve, as in
  # the reference. The curve's lower tail runs out to x = -1.49, so the
  # free corridor is widest in the middle of the panel and closes up
  # near the bottom -- PBRM1 (y = 0.02) therefore gets its own band down
  # in the corner rather than competing for the corridor.
  #
  # TWO repel layers with DISJOINT y ranges so they cannot collide with
  # each other; the corridor layer is capped below the legend block.
  G_CEIL   <- g_leg$ymin - 0.04
  hl_low   <- hl |> dplyr::filter(ecdf_y <  0.10)
  hl_mid   <- hl |> dplyr::filter(ecdf_y >= 0.10)

  ann <- data.frame(x = G_XLIM[2], y = 0.78, label = meta$p_label)

  # ---------------- plot ------------------------------------------
  p <- ggplot() +
    geom_ribbon(data = ribbon, aes(x = x, ymin = lower, ymax = upper),
                fill = "grey65", alpha = 0.45) +
    geom_line(data = ribbon,
              aes(x = x, y = mean_cdf, color = grp, linetype = grp),
              linewidth = G_LW) +
    stat_ecdf(data = red, aes(x = x, color = grp, linetype = grp),
              linewidth = G_LW, na.rm = TRUE) +
    geom_point(data = hl, aes(x = x, y = ecdf_y),
               size = 0.45, color = "grey15", inherit.aes = FALSE) +
    ggrepel::geom_text_repel(
      data = dplyr::bind_rows(hl_mid, obstacles),
      aes(x = x, y = ecdf_y, label = gene_name),
      inherit.aes = FALSE, color = "grey15", fontface = "italic",
      family = FONT_FAMILY, size = SIZE_BASE / .pt,
      nudge_x = -0.30, direction = "y", hjust = 1,
      xlim = c(-1.27, -0.18), ylim = c(0.13, G_CEIL),
      box.padding = 0.10, point.padding = 0.02,
      min.segment.length = 0, segment.size = 0.25,
      segment.color = "grey55", segment.alpha = 0.9,
      max.overlaps = Inf, max.time = 1e6, max.iter = 30000,
      seed = G_SEED, show.legend = FALSE) +
    ggrepel::geom_text_repel(
      data = dplyr::bind_rows(hl_low, obstacles),
      aes(x = x, y = ecdf_y, label = gene_name),
      inherit.aes = FALSE, color = "grey15", fontface = "italic",
      family = FONT_FAMILY, size = SIZE_BASE / .pt,
      nudge_x = -0.28, nudge_y = 0.05, hjust = 1,
      xlim = c(-1.27, -0.82), ylim = c(0.02, 0.11),
      box.padding = 0.10, point.padding = 0.02,
      min.segment.length = 0, segment.size = 0.25,
      segment.color = "grey55", segment.alpha = 0.9,
      max.overlaps = Inf, max.time = 1e6, max.iter = 30000,
      seed = G_SEED, show.legend = FALSE) +
    geom_text(data = ann, aes(x = x, y = y, label = label),
              inherit.aes = FALSE, hjust = 1, vjust = 0.5,
              family = FONT_FAMILY, size = SIZE_BASE / .pt,
              color = "grey15") +
    scale_color_manual(values = G_PAL, breaks = names(G_PAL)) +
    scale_linetype_manual(values = G_LTY, breaks = names(G_LTY)) +
    scale_x_continuous(breaks = seq(-1, 1, 0.5)) +
    scale_y_continuous(breaks = seq(0, 1, 0.25),
                       labels = scales::number_format(accuracy = 0.01)) +
    coord_cartesian(xlim = G_XLIM, ylim = G_YLIM, clip = "on") +
    # The x title is wrapped: on one line it is 4.0 in wide against a
    # 4.07 in cell, so it ran off the page edge and into panel f.
    labs(x = paste0("Mean gene expression over genes\n",
                    "with unproductive junctions"),
         y = "Empirical CDF", color = NULL, linetype = NULL,
         title = gene_md(paste0("YLPM1 LoF confers widespread gene expression\n",
                        "downregulation through mis-splicing"))) +
    theme_publication() +
    theme(
      plot.title             = md_title(),
      legend.position        = "inside",
      legend.position.inside = c(G_LEG_NPC_X, 1.00),
      legend.justification   = c(0, 1),
      legend.background      = element_blank(),
      legend.key             = element_blank(),
      legend.text            = element_text(size = SIZE_LEGEND,
                                            family = FONT_FAMILY,
                                            lineheight = 0.85),
      legend.key.height      = unit(2.8, "mm"),
      legend.key.width       = unit(4.6, "mm"),
      legend.key.spacing.y   = unit(0.4, "mm"),
      legend.margin          = margin(0, 0, 0, 0)
    )

  wrap_elements(full = p)
}))
