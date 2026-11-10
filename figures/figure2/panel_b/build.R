# Figure 2, panel b. Sourced by ../compose.R, which supplies the
# shared setup (theme tokens, FLASHIER/CACHE paths, ROW_H, panel_cached,
# missing_panel, sf3b1_info). Must leave a `panel_b` object defined.
#
# Rebuild just this panel:   Rscript compose.R b

# ---- b. GSP #1 neuronal-junction ECDF + WDR11 sashimi inset ------
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
  # Framed on the cryptic-junction cluster the panel is about. The
  # .links also contain two flanking ~99.8% junctions
  # (120885938-120886688 and 120889184-120889894); including them
  # would dominate the inset and is not what the reference shows.
  # Left edge is pulled back to 120886650 so the WHOLE upstream exon
  # (coverage runs 120886687-120886843) is inside the frame: cutting it
  # at 120886700 sliced 13 bp off the rising flank and the peak read as
  # clipped. Still right of 120885938, so the flanking 99.8% link stays
  # filtered out by sashimi_track_gg's `s1 >= start` test.
  b_region <- c(120886650, 120889260)
  B_K23_JUNC <- c(120886836, 120887209)       # WDR11-k23 unique junction (links coordinates)
  B_K23_COL  <- "#D62728"                     # red: panel b already uses orange (TCGA brain) for coverage
  # Inset geometry. The inset sits lower-right so it does NOT cover the
  # WDR11 point on the TCGA-brain curve; a connector runs from that point
  # to the inset's top-left corner. Drawn sizes below feed the PSI-label
  # clearance in sashimi_track_gg (see [[derive-layout-offsets]]).
  # Widened/deepened vs. the first pass: every ECDF curve is above
  # y = 0.84 for x > 0.04, so the whole lower-right quadrant is free
  # and the arcs need the horizontal room to separate their labels.
  # Raised: every ECDF curve is above y = 0.86 for x > 0.044, so the
  # inset can sit much closer to the WDR11 point than the old top =
  # 0.74 did -- which left a long connector dangling through empty
  # space. Top = 0.84 puts the inset's top-left corner just under the
  # point, so the connector is a short drop the eye follows.
  # Shrunk (was left 0.40 / top 0.82, i.e. 0.60 x 0.74 npc): at that
  # size the inset owned most of panel b and read as the subject of
  # the panel rather than a callout on it. 0.50 x 0.60 npc is ~33%
  # less area. The wedge, the bounding rect and the connector all
  # derive from B_INS, so they follow automatically.
  # right 1.00 -> 0.965: with no box drawn, a leader ending ON the panel's
  # right edge read as clipped (QA review). Ending it just inside fixes that.
  B_INS <- c(left = 0.50, bottom = 0.12, right = 0.965, top = 0.66)
  PANEL_B_W_IN <- (COL_W$bc[["ecdf"]] / sum(COL_W$bc)) * PAGE_W * 0.82
  PANEL_B_H_IN <- ROW_H[2] * 0.80 * PANEL_H_SCALE[["b"]]   # drawn height (compose.R v_shrink)
  B_INS_W_IN <- (B_INS["right"] - B_INS["left"])  * PANEL_B_W_IN
  B_INS_H_IN <- (B_INS["top"]   - B_INS["bottom"]) * PANEL_B_H_IN
  # tcga, matched-control, gene body. The control track is deliberately
  # SHORTER: it carries one arc where the TCGA track carries four, so at
  # equal heights everything above its single arc was dead air and its
  # own "Matched controls (400)" label floated in the middle of it.
  # The control track is SHORTER still: it carries one arc where the
  # TCGA track carries four, and y_max scales with drawn height (below)
  # so shrinking it removes headroom without shrinking the coverage
  # peaks. That, plus the tighter gutter, closes the gap between the
  # two tracks without letting them touch.
  B_ROWS     <- c(1, 0.58, 0.38)
  # Gutter below each coverage track. Without it the lower track's
  # top-right label ("Matched controls (400)") sits directly on the
  # upper track's coverage baseline -- a rule touching text. One label
  # height of air is the smallest gap that reads as a separation, so
  # derive it from the label's point size rather than eyeballing pt.
  B_LAB_PT  <- SIZE_BASE * 0.8                         # track/PSI label size
  # The cohort caption now hangs UNDER each coverage track, so the
  # gutter has to hold a whole line of it plus air, not just separate
  # the tracks.
  B_GAP_PT  <- B_LAB_PT * 1.55
  # The matched controls are the non-brain TCGA samples, so they carry
  # the TCGA-other colour from the ECDF legend -- not GTEx-brain blue,
  # which named a cohort that is not in this track.
  B_TOP_COL <- unname(PAL$cohort["TCGA-brain"])
  B_BOT_COL <- unname(PAL$cohort["TCGA-other"])
  # Drawn height of each coverage track = its share of the inset minus
  # its gutter. Both feed sashimi_track_gg's point-size -> data-unit
  # conversion, so the label clearances follow the real drawn size.
  B_TRACK_H <- unname(B_INS_H_IN * B_ROWS[1:2] / sum(B_ROWS)) - B_GAP_PT / 72

  # Coverage scale shared between the tracks. Their maxima are within 2%
  # of each other, so autoscaling each one separately would imply a
  # depth difference that isn't there; the claim of this inset lives in
  # the ARCS, not the peak heights. Sharing the scale across tracks of
  # UNEQUAL height means scaling y_max with the drawn height -- then one
  # read equals the same number of points in both. 3.2x leaves headroom
  # for four stacked arc tiers plus a label above the highest and puts
  # the exon peaks at ~31% of the TCGA track (at the helper's 3.8x
  # default they read as stubs).
  b_bws     <- c("High_k2-3-..bw", "Matched_Controls-3-..bw")
  b_cov_max <- max(vapply(b_bws, function(f)
    max(bigwig_binned(file.path(WD, f), "chr10", b_region[1], b_region[2]
                      )$coverage, na.rm = TRUE), numeric(1)))
  b_y_max <- b_cov_max * 3.2 * B_TRACK_H / B_TRACK_H[1]

  b_track_theme <- theme(plot.margin = margin(0, 4, B_GAP_PT, 0, "pt"))

  # ---- PSI-label placement -------------------------------------
  # sashimi_track_gg's solver anchors every PSI label at its arc's
  # MIDPOINT and then pushes it up/down until it clears the arcs. That
  # is fine when arcs are spread out, but three of the four junctions
  # here share the 120889077 acceptor and so share a midpoint to within
  # a few hundred bp: the solver can only stack them, and at this track
  # height two of the labels end up printing on top of each other.
  # Worse, the 4.89% and 1.27% junctions differ by 3 bp at the donor,
  # so their arcs are ~parallel with only ~7 pt between apexes -- no
  # label can ever be threaded between them.
  #
  # So place them explicitly. The arcs' true geometry is already inside
  # the object the helper returns (one GeomPath layer per arc, in the
  # same row order as the PSI-label layer), so a label can be pinned a
  # given number of LABEL HEIGHTS above/below any named arc rather than
  # by a hand-tuned data-unit constant. Label height uses the same
  # point-size -> data-unit conversion the helper does, so it tracks
  # the drawn track size automatically.
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
  # Top track, four arcs (tiers run shortest-highest -> longest-lowest):
  #   7.5%  cryptic-exon inclusion, the tall narrow arc on the left
  #   1.3% / 4.9%  the 3-bp-apart donor twins, nested
  #   71%   the thick annotated junction, lowest tier
  # 7.5% and 1.3% each sit above their own apex, where nothing else
  # reaches. 71% is pinned a third of the way along its own arc, left of
  # where the twins climb across it -- there it is the topmost line, so
  # the label reads onto the thick arc directly beneath it. 4.9% cannot
  # go between the twins, so it drops into the empty intronic band,
  # cleared off the *71%* arc (the lowest thing above it) and aligned on
  # its own apex x so the eye runs straight up the stack to it.
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
    # (cache/panel_b_top_junc.tsv), drawn in red (Quinn, 2026-11-09; orange would clash with the TCGA-brain coverage).
    highlight = list(B_K23_JUNC), highlight_color = B_K23_COL,
    min_psi = 1, y_max = b_y_max[1], label_size = B_LAB_PT / .pt,
    # Caption hangs UNDER this track's coverage, left-aligned, in the
    # track's own colour.
    label_hjust = 0, label_color = B_TOP_COL, label_below = TRUE,
    # FLATTER (was c(0.84, 0.40)). The 7.5% junction has much the
    # shortest span, so it takes the highest tier -- and with its right
    # endpoint down on the intronic baseline it drew as a near-vertical
    # spike that dominated the inset. What matters visually is each
    # arc's RISE ABOVE ITS OWN ENDPOINTS, not the raw apex fraction, so
    # these are back-solved for ~30% less rise: the endpoints here sit
    # at ~0.30 of y_max, so 0.84 -> 0.30 + 0.70*(0.84-0.30) = 0.68 and
    # 0.40 -> 0.30 + 0.70*(0.40-0.30) = 0.37.
    # Do NOT reach for min_arc_rise to compensate -- see the comment on
    # the control track below and in ../scripts/coverage_plot.R.
    # 0.68 -> 0.70 top / 0.37 bottom. (0.60/0.42 packed the tiers too close for
    # INLINE numbers: 71%, 4.9% and 1.3% jammed together.) On the top tier (Quinn: "the 7.5% PSI curve is too
    # high"). 7.5% has much the shortest span, so it takes the top
    # tier, and a short arc with a big rise reads as a SPIKE. (0.52 was
    # tried: 71%, 1.3% and 4.9% then merged into one line at the right
    # exon.) 0.60 gives
    # it ~0.22 of y_max above its left endpoint instead of ~0.38.
    apex_range = c(0.78, 0.36),
    # Per-arc peaks: the 1.3%/4.9% donor twins share their acceptor with
    # 71% and ran together with it along the right half, so no span-
    # ranked tier order left room for three inline numbers. Twins high,
    # 71% low, 7.5% (short, cryptic exon) moderate so it does not spike.
    apex_override = c("1.3%" = 0.84, "4.9%" = 0.66, "71%" = 0.42,
                      "7.5%" = 0.60),
    psi_size = B_LAB_PT / .pt, label_family = FONT_FAMILY,
    track_height_in = B_TRACK_H[1], track_width_in = unname(B_INS_W_IN),
    lw_scale = 0.45)
  # Numbers are now placed by the shared helper's INLINE style (each on
  # its own arc, in a gap cut out of it) -- one label style for every
  # sashimi plot in the paper. b_place_psi()/B_PSI_POS are kept above
  # for reference but no longer applied: moving a label after the arc
  # has been cut would leave the gap in the wrong place.
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
    # Single arc -> apex at mean(apex_range). This track is short, so
    # y_max is only ~1.7x its own coverage peak; at the default 0.545
    # the apex landed AT the peak height and the 98% arc drew as a
    # straight line. 0.85 put it ~40% above the peaks so it read as an
    # arc again.
    # FLATTENED with the top track, but by RISE, not by apex fraction:
    # this arc's endpoints sit at ~0.54 of y_max (the peaks), so a flat
    # 30% cut of 0.85 would land at 0.60 -- essentially ON the
    # endpoints, i.e. a straight chord, the same failure the default
    # 0.545 produced. 0.54 + 0.70*(0.85-0.54) = 0.76 takes 30% off the
    # visible rise and keeps the curvature.
    # min_arc_rise stays at its 0 default: it is a FLOOR, and
    # `max(apex, max(y1,y2) + rise)` clamps the apex up to the higher
    # endpoint the moment it is switched on with a small value, which
    # draws every arc as a straight chord.
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
    # The gene is named HERE, under its gene body, not at the ECDF point
    # (Quinn, 2026-10-05). The box's bottom edge is dropped below this
    # label (B_BOX_BOT_PAD) -- it used to cut through it.
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

  # inset_element()'s left/bottom/right/top are PANEL NPC, but the
  # connector is drawn in DATA units -- so the npc corner has to be
  # mapped through the scale's expansion or the line stops short of the
  # inset. These must stay in step with the two scale_*_continuous()
  # calls below.
  B_XLIM <- c(0, 0.11); B_YLIM <- c(0, 1); B_EXP <- 0.01
  b_npc_to_data <- function(f, lim) {
    pad <- diff(lim) * B_EXP
    (lim[1] - pad) + f * (diff(lim) + 2 * pad)
  }
  # ---- callout wedge -----------------------------------------------
  # A light bounding box around the inset, with TWO leaders running
  # from the WDR11 point to the box's two TOP corners. This is the
  # standard zoom-callout idiom and it removes the problem a single
  # leader had here: with the inset raised there is only ~0.6 in
  # between the point and the arcs, so one line either stopped in
  # mid-air or ran nearly horizontally along the GTEx-brain curve. The
  # wedge instead reads as "this point expands into this box", and the
  # box edge is what the lines terminate on.
  # `sashimi_track_gg` draws its y-axis (tick + "38.7"/"0") at
  # `start - n_pt` in DATA units with clip off, i.e. OUTSIDE its own
  # panel and therefore outside the npc region inset_element gives it.
  # The callout box has to be widened by that much on the left or it
  # cuts the axis numbers off outside the frame. Measured at 400 dpi:
  # the numbers overhang the inset's left edge by 65 px against the
  # inset's 1400 px per npc unit, so 0.055 npc covers it with air.
  B_AXIS_PAD <- 0.115   # 0.055 -> 0.085: the box edge still ran through "38.7" / "19.4"
  b_ins_x0 <- b_npc_to_data(unname(B_INS["left"]) - B_AXIS_PAD, B_XLIM)
  b_ins_x1 <- b_npc_to_data(unname(B_INS["right"]),  B_XLIM)
  # Bottom edge sits below the inset's lower edge by one WDR11 label
  # height, so the gene name under the gene body is inside the box.
  B_BOX_BOT_PAD <- 0.045
  b_ins_y0 <- b_npc_to_data(unname(B_INS["bottom"]) - B_BOX_BOT_PAD, B_YLIM)
  # Box top sits a little ABOVE the inset's top edge: the top track's
  # "38.7" is centred on its top tick, i.e. on that edge, so a box line
  # (and the left leader ending at the box corner) ran through it.
  B_BOX_TOP_PAD <- 0.035
  b_ins_y1 <- b_npc_to_data(unname(B_INS["top"]) + B_BOX_TOP_PAD, B_YLIM)

  # Both leaders as one data frame so they share a single layer.
  b_wedge <- data.frame(
    x    = b_wdr11_x,        y    = b_wdr11_y,
    xend = c(b_ins_x0, b_ins_x1), yend = b_ins_y1)

  panel_b <- wrap_elements(full =
    (ggplot(b_long, aes(x = psi, color = category)) +
       stat_ecdf(linewidth = 0.6) +
       # Light bounding box around the inset (Quinn, 2026-09-23: "Figure
       # 2B should still have a bounding box around the WDR11 inset" --
       # supersedes the earlier "no box needed"). The leaders end on its
       # two top corners.
       annotate("rect", xmin = b_ins_x0, xmax = b_ins_x1,
                ymin = b_ins_y0, ymax = b_ins_y1,
                fill = NA, color = "grey75", linewidth = 0.3) +
       geom_segment(data = b_wedge, inherit.aes = FALSE,
                    aes(x = x, y = y, xend = xend, yend = yend),
                    color = "grey60", linewidth = 0.3,
                    linetype = "dashed") +   # Quinn, 2026-10-05
       annotate("point", x = b_wdr11_x, y = b_wdr11_y, size = 1.6,
                color = "black", fill = unname(PAL$cohort["TCGA-brain"]),
                shape = 21, stroke = 0.5) +
       scale_color_manual(values = PAL$cohort, name = NULL) +
       # The view is set by coord, NOT by scale limits. A scale limit
       # DROPS the 141 junctions with PSI > 0.11 before stat_ecdf ever
       # runs, so the drawn curve was the ECDF of a truncated sample --
       # it reached 1.0 at x = 0.11 and sat ABOVE the true ECDF
       # everywhere. b_wdr11_y comes from ecdf() on the full vector, so
       # the WDR11 point was placed off its own curve. Same trap as the
       # sashimi tracks (see ../scripts/coverage_plot.R).
       scale_x_continuous(expand = c(0.01, 0)) +
       scale_y_continuous(expand = c(0.01, 0)) +
       # clip MUST stay "on": the ECDF now extends past x = 0.11 (that
       # is the whole point of using coord rather than scale limits),
       # and with clip off those tails were drawn across the rest of
       # the figure row.
       coord_cartesian(xlim = B_XLIM, ylim = B_YLIM, clip = "on") +
       labs(x = "Mean junction PSI", y = "Cumulative distribution",
            title = "GSP #1 represents neuronal junctions") +
       theme_publication() +
       theme(plot.title      = element_text(hjust = 0.5),
             plot.margin     = margin(4, 6, 4, 14, "pt"),
             # Legend lower-left: the inset now owns the lower-right.
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

