# Figure 2c (built here as panel_b). Sourced by ../compose.R, which supplies the
# shared setup (theme tokens, FLASHIER/CACHE paths, ROW_H, panel_cached,
# missing_panel, sf3b1_info). Must leave a `panel_b` object defined.
#
# Rebuild just this panel:   Rscript compose.R b

# ---- GSP #1 neuronal-junction ECDF + WDR11 sashimi inset ---------
if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
if (!have_cache("panel_b_running_psi.tsv")) {
  panel_b <- missing_panel("b  GSP #1 neuronal ECDF",
                           "build_cache_panel_b.R")
} else if (!identical(REBUILD, "ALL") && !("b" %in% REBUILD) &&
           file.exists(file.path(PANEL_DIR, "b.rds"))) {
  cat("[panel] b (cached)\n"); panel_b <- readRDS(file.path(PANEL_DIR, "b.rds"))
} else {
  cat("[panel] b building ...\n")
  b_psi <- fread(file.path(CACHE, "panel_b_running_psi.tsv"),
                  data.table = FALSE)

  # Column names in the cache are TCGA_brain_psi etc.; map to the
  # display names used in the legend and PAL$cohort.
  B_MAP <- c(GTEx_brain_psi = "GTEx-brain", TCGA_brain_psi = "TCGA-brain",
             GTEx_other_psi = "GTEx-other", TCGA_other_psi = "TCGA-other")
  b_long <- b_psi %>%
    pivot_longer(cols = any_of(names(B_MAP)), names_to = "category",
                 values_to = "psi") %>%
    filter(!is.na(psi)) %>%
    mutate(category = factor(unname(B_MAP[category]),
                             levels = names(PAL$cohort)))

  # Sashimi inset, drawn natively from the pyGenomeTracks working
  # directory's bigwigs + .links (same helpers as figure 1) rather
  # than insetting its rendered PDF -- keeps the text vector and
  # consistent with the rest of the figure.
  WD <- file.path(QROOT, "coverage_tracks/pygenome/neuronal_gsp/WDR11")
  # Framed on the cryptic-junction cluster. The .links also contain two
  # flanking ~99.8% junctions (120885938-120886688, 120889184-120889894)
  # that would dominate the inset. The left edge includes the whole
  # upstream exon (coverage 120886687-120886843) while staying right of
  # 120885938, so that flanking link is dropped by sashimi_track_gg's
  # `s1 >= start` test.
  b_region <- c(120886650, 120889260)
  B_K23_JUNC <- c(120886836, 120887209)       # WDR11-k23 unique junction (links coordinates)
  B_K23_COL  <- "#D62728"                     # red: orange (TCGA brain) is already the coverage colour
  # Inset geometry (panel npc). The inset sits lower-right, where every
  # ECDF curve is already near 1, so it does not cover the WDR11 point
  # on the TCGA-brain curve. The wedge, bounding rect and connectors all
  # derive from B_INS. Drawn sizes below feed the PSI-label clearance in
  # sashimi_track_gg.
  B_INS <- c(left = 0.50, bottom = 0.12, right = 0.965, top = 0.66)
  PANEL_B_W_IN <- (COL_W$bc[["ecdf"]] / sum(COL_W$bc)) * PAGE_W * 0.82
  PANEL_B_H_IN <- ROW_H[2] * 0.80 * PANEL_H_SCALE[["b"]]   # drawn height (compose.R v_shrink)
  B_INS_W_IN <- (B_INS["right"] - B_INS["left"])  * PANEL_B_W_IN
  B_INS_H_IN <- (B_INS["top"]   - B_INS["bottom"]) * PANEL_B_H_IN
  # Row heights: tcga, matched-control, gene body. The control track is
  # shorter because it carries one arc where the TCGA track carries four;
  # y_max scales with drawn height (below), so this removes headroom
  # without shrinking the coverage peaks.
  B_ROWS     <- c(1, 0.58, 0.38)
  # Gutter below each coverage track, derived from the label point size.
  B_LAB_PT  <- SIZE_BASE * 0.8                         # track/PSI label size
  # The cohort caption hangs under each coverage track, so the gutter
  # holds a whole line of it plus air.
  B_GAP_PT  <- B_LAB_PT * 1.55
  # The matched controls are non-brain TCGA samples, so they carry the
  # TCGA-other colour from the ECDF legend.
  B_TOP_COL <- unname(PAL$cohort["TCGA-brain"])
  B_BOT_COL <- unname(PAL$cohort["TCGA-other"])
  # Drawn height of each coverage track = its share of the inset minus
  # its gutter. Both feed sashimi_track_gg's point-size -> data-unit
  # conversion, so the label clearances follow the real drawn size.
  B_TRACK_H <- unname(B_INS_H_IN * B_ROWS[1:2] / sum(B_ROWS)) - B_GAP_PT / 72

  # Coverage scale shared between the tracks (their maxima are within 2%),
  # with y_max scaled by drawn height so one read is the same number of
  # points in both. 3.2x leaves headroom for four stacked arc tiers.
  b_bws     <- c("High_k2-3-..bw", "Matched_Controls-3-..bw")
  b_cov_max <- max(vapply(b_bws, function(f)
    max(bigwig_binned(file.path(WD, f), "chr10", b_region[1], b_region[2]
                      )$coverage, na.rm = TRUE), numeric(1)))
  b_y_max <- b_cov_max * 3.2 * B_TRACK_H / B_TRACK_H[1]

  b_track_theme <- theme(plot.margin = margin(0, 4, B_GAP_PT, 0, "pt"))

  # ---- PSI-label placement -------------------------------------
  # Explicit PSI-label placement: pins a label a given number of label
  # heights above/below a named arc, using the arc geometry inside the
  # helper's returned object. Defined but not applied below -- the
  # labels use the helper's inline style.
  b_place_psi <- function(p, spec, y_max, track_h_in, psi_pt) {
    arcs <- Filter(function(l) inherits(l$geom, "GeomPath"), p$layers)
    ti <- which(vapply(p$layers, function(l)
      inherits(l$geom, "GeomText") && is.data.frame(l$data) &&
        !is.null(l$data$label) && all(grepl("%$", l$data$label)),
      logical(1)))[1]
    d <- p$layers[[ti]]$data
    lab_h <- (psi_pt * 0.55 + 1.5) * y_max / (track_h_in * 72)
    for (i in seq_len(nrow(d))) {
      s <- spec[[d$label[i]]]
      if (is.null(s)) next
      own <- arcs[[i]]$data
      # `at`: "apex" or a fraction along this arc's own x-span.
      xi <- if (identical(s$at, "apex")) own$x[which.max(own$y)]
            else min(own$x) + s$at * diff(range(own$x))
      # `rel`: which arc to clear (defaults to the label's own arc).
      ref <- if (is.null(s$rel)) own else arcs[[match(s$rel, d$label)]]$data
      d$x[i] <- xi
      d$y[i] <- approx(ref$x, ref$y, xout = xi, rule = 2)$y + s$dy * lab_h
    }
    p$layers[[ti]]$data <- d
    p
  }
  # Top track, four arcs: 7.5% cryptic-exon inclusion; 1.3% / 4.9% donor
  # twins 3 bp apart; 71% the thick annotated junction.
  B_PSI_POS <- list(
    "7.5%" = list(at = "apex", dy =  1.15),
    "1.3%" = list(at = "apex", dy =  1.15),
    "71%"  = list(at = 0.32,   dy =  1.30),
    "4.9%" = list(at = "apex", dy = -1.30, rel = "71%")
  )
  b_sash_top <- sashimi_track_gg(
    bw_path = file.path(WD, "High_k2-3-..bw"),
    chrom = "chr10", start = b_region[1], end = b_region[2],
    track_color = B_TOP_COL,
    track_label = "Top 100 TCGA on GSP #1",
    links_path = file.path(WD, "MeanPSI.High_k2-3-..links"),
    # WDR11-k23's unique (driving) junction, Snaptron chr10:120886837-120887209 = links 120886836-120887209
    # (cache/panel_b_top_junc.tsv), drawn in red.
    highlight = list(B_K23_JUNC), highlight_color = B_K23_COL,
    min_psi = 1, y_max = b_y_max[1], label_size = B_LAB_PT / .pt,
    # Caption hangs under this track's coverage, left-aligned, in the
    # track's own colour.
    label_hjust = 0, label_color = B_TOP_COL, label_below = TRUE,
    # Arc endpoints sit at ~0.30 of y_max; keep the rise above them modest.
    # Leave min_arc_rise at 0 (see the control track below).
    apex_range = c(0.78, 0.36),
    # Per-arc peaks: the 1.3%/4.9% twins share their acceptor with 71%,
    # so no span-ranked tier order leaves room for three inline numbers.
    # Twins high, 71% low, 7.5% (short cryptic exon) moderate so it does
    # not read as a spike.
    apex_override = c("1.3%" = 0.84, "4.9%" = 0.66, "71%" = 0.42,
                      "7.5%" = 0.60),
    psi_size = B_LAB_PT / .pt, label_family = FONT_FAMILY,
    track_height_in = B_TRACK_H[1], track_width_in = unname(B_INS_W_IN),
    lw_scale = 0.45)
  # PSI numbers use the shared helper's inline style (each on its own arc,
  # in a gap cut out of it), as in every sashimi plot in the paper.
  b_sash_top <- b_sash_top + b_track_theme
  b_sash_bot <- sashimi_track_gg(
    bw_path = file.path(WD, "Matched_Controls-3-..bw"),
    chrom = "chr10", start = b_region[1], end = b_region[2],
    track_color = B_BOT_COL,
    track_label = "Matched controls (400)",
    links_path = file.path(WD, "MeanPSI.Matched_Controls-3-..links"),
    highlight = list(B_K23_JUNC), highlight_color = B_K23_COL,
    min_psi = 1, y_max = b_y_max[2], label_size = B_LAB_PT / .pt,
    label_hjust = 0, label_color = B_BOT_COL, label_below = TRUE,
    # Single arc -> apex at mean(apex_range). Its endpoints sit at ~0.54
    # of y_max (the coverage peaks), so the apex must sit well above that
    # or the arc draws as a straight chord.
    # min_arc_rise stays at its 0 default: as a floor,
    # `max(apex, max(y1,y2) + rise)` clamps a small value up to the higher
    # endpoint and draws every arc as a straight chord.
    apex_range = c(0.78, 0.74),
    psi_size = B_LAB_PT / .pt, label_family = FONT_FAMILY,
    track_height_in = B_TRACK_H[2], track_width_in = unname(B_INS_W_IN),
    lw_scale = 0.45) + b_track_theme
  # WDR11 MANE gene body under the tracks (exons flanking the cryptic
  # exon; the cryptic exon itself is unannotated, which is the point).
  b_gene <- gene_body_gg(readRDS(file.path(CACHE, "wdr11_mane_exons.rds")),
                         "chr10", b_region[1], b_region[2],
                         gene_label = NULL, exon_h = 0.45, n_arrows = 6,
                         label_family = FONT_FAMILY) +
    # The gene is named under its gene body; the box's bottom edge sits
    # below this label (B_BOX_BOT_PAD).
    annotate("text", x = b_region[1] + diff(b_region) * 0.5, y = 1,
             label = "WDR11", fontface = "italic", vjust = 1.9,
             size = B_LAB_PT / .pt, family = FONT_FAMILY,
             color = "grey15") +
    theme(plot.margin = margin(0, 4, 0, 0, "pt"))
  # Legend row for the highlighted arc: an arc glyph + label, under the gene body.
  b_leg <- ggplot() +
    # Sits in the LOWER part of its row: the italic WDR11 gene label above hangs (vjust 1.9) into the top.
    geom_path(data = data.frame(x = 0.0 + seq(0, 1, length.out = 50) * 0.06,
                                y = 0.08 + sin(seq(0, pi, length.out = 50)) * 0.30),
              aes(x, y), colour = B_K23_COL, linewidth = 0.8) +
    annotate("text", x = 0.085, y = 0.16, label = "WDR11-k23 unique junction", hjust = 0, vjust = 0.5,
             fontface = "plain", size = B_LAB_PT / .pt, family = FONT_FAMILY, colour = "grey15") +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    coord_cartesian(clip = "off") + theme_void()
  b_sashimi <- (b_sash_top / b_sash_bot / b_gene / b_leg) +
    plot_layout(heights = c(B_ROWS, 0.30 * sum(B_ROWS)))

  # Connector: WDR11's mean PSI in the top-100 GSP #1 samples, on the
  # TCGA-brain ECDF, to the inset's top-left corner (in data units).
  b_wdr11_x <- b_psi$TCGA_brain_psi[b_psi$gene_name == "WDR11"]
  b_wdr11_y <- ecdf(b_long$psi[b_long$category == "TCGA-brain"])(b_wdr11_x)

  # inset_element()'s left/bottom/right/top are panel npc, but the
  # connector is drawn in data units, so npc is mapped through the
  # scale's expansion. Must stay in step with the scale_*_continuous()
  # calls below.
  B_XLIM <- c(0, 0.11); B_YLIM <- c(0, 1); B_EXP <- 0.01
  b_npc_to_data <- function(f, lim) {
    pad <- diff(lim) * B_EXP
    (lim[1] - pad) + f * (diff(lim) + 2 * pad)
  }
  # ---- callout wedge -----------------------------------------------
  # A light bounding box around the inset, with two leaders from the
  # WDR11 point to the box's two top corners (zoom-callout idiom).
  # sashimi_track_gg draws its y-axis numbers outside its own panel
  # (clip off), so the box is widened on the left to enclose them.
  B_AXIS_PAD <- 0.115
  b_ins_x0 <- b_npc_to_data(unname(B_INS["left"]) - B_AXIS_PAD, B_XLIM)
  b_ins_x1 <- b_npc_to_data(unname(B_INS["right"]),  B_XLIM)
  # Bottom edge sits below the inset's lower edge by one WDR11 label
  # height, so the gene name under the gene body is inside the box.
  B_BOX_BOT_PAD <- 0.045
  b_ins_y0 <- b_npc_to_data(unname(B_INS["bottom"]) - B_BOX_BOT_PAD, B_YLIM)
  # Box top sits a little above the inset's top edge, clearing the top
  # track's y-axis label centred on that edge.
  B_BOX_TOP_PAD <- 0.035
  b_ins_y1 <- b_npc_to_data(unname(B_INS["top"]) + B_BOX_TOP_PAD, B_YLIM)

  # Both leaders as one data frame so they share a single layer.
  b_wedge <- data.frame(
    x    = b_wdr11_x,        y    = b_wdr11_y,
    xend = c(b_ins_x0, b_ins_x1), yend = b_ins_y1)

  panel_b <- wrap_elements(full =
    (ggplot(b_long, aes(x = psi, color = category)) +
       stat_ecdf(linewidth = 0.6) +
       # Light bounding box around the inset; the leaders end on its two
       # top corners.
       annotate("rect", xmin = b_ins_x0, xmax = b_ins_x1,
                ymin = b_ins_y0, ymax = b_ins_y1,
                fill = NA, color = "grey75", linewidth = 0.3) +
       geom_segment(data = b_wedge, inherit.aes = FALSE,
                    aes(x = x, y = y, xend = xend, yend = yend),
                    color = "grey60", linewidth = 0.3,
                    linetype = "dashed") +
       annotate("point", x = b_wdr11_x, y = b_wdr11_y, size = 1.6,
                color = "black", fill = unname(PAL$cohort["TCGA-brain"]),
                shape = 21, stroke = 0.5) +
       scale_color_manual(values = PAL$cohort, name = NULL) +
       # The view is set by coord, not scale limits: a scale limit would
       # drop junctions with PSI > 0.11 before stat_ecdf runs, and
       # b_wdr11_y comes from ecdf() on the full vector.
       scale_x_continuous(expand = c(0.01, 0)) +
       scale_y_continuous(expand = c(0.01, 0)) +
       # clip must stay "on": the ECDF extends past x = 0.11 and would
       # otherwise be drawn across the rest of the figure row.
       coord_cartesian(xlim = B_XLIM, ylim = B_YLIM, clip = "on") +
       labs(x = "Mean junction PSI", y = "Cumulative distribution",
            title = "GSP #1 represents neuronal junctions") +
       theme_publication() +
       theme(plot.title      = element_text(hjust = 0.5),
             plot.margin     = margin(4, 6, 4, 14, "pt"),
             # Legend lower-left: the inset owns the lower-right.
             legend.position.inside = c(0.37, 0.02),
             legend.position = "inside",
             legend.justification = c(1, 0),
             legend.key.size = unit(3, "mm"),
             legend.background = element_blank())) +
      inset_element(b_sashimi,
                    left = B_INS["left"], bottom = B_INS["bottom"],
                    right = B_INS["right"], top = B_INS["top"],
                    align_to = "panel"))
  saveRDS(panel_b, file.path(PANEL_DIR, "b.rds"))
}

