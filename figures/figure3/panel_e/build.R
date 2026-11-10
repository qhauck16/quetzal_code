# ------------------------------------------------------------------
# Figure 3, panel e -- "Mis-splicing leads to widespread effects on
# gene expression".
#
# Four empirical CDFs of mean scaled RPM across the top-20 NKAP (GSP
# #13) samples:
#   frameshifting cryptic 3'ss genes, in-frame cryptic 3'ss genes, and
#   an RPM-matched control gene set for each.
#
# Transcribed from diff_splicing/combined_expression_analysis.Rmd
# (build_set_pair / match_by_rpm / plot_single_panel_cdf). See
# panel_e/build.md for the analysis rules and the layout reasoning.
#
# Sourced by ../compose.R, which supplies theme.R's tokens, the
# tidyverse, ggrepel, CACHE and panel_cached().
# ------------------------------------------------------------------

if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
panel_e <- panel_cached("e", local({

  # ---------------- constants -------------------------------------
  E_XLIM   <- c(-1.30, 0.80)      # view only -- never a scale limit
  E_YLIM   <- c(-0.04, 1.08)      # headroom for labels above the CDF
  E_LW     <- 0.70                # ECDF stroke
  E_SEED   <- 20260922L
  E_LAB_PT <- SIZE_BASE           # gene callouts
  E_LEG_PT <- SIZE_LEGEND         # legend text
                                  # Both stay at the shared tokens.
                                  # They were 7 pt while row 3 was
                                  # 2.45 in; row 3 has since been grown
                                  # and then NARROWED (cell 3.40 ->
                                  # 3.23 in, row 2.95 -> 3.05 in) when
                                  # figure 3 went from four rows to
                                  # three. The narrowing was paid for
                                  # by re-wrapping the legend and
                                  # re-cutting the tiers, NOT by
                                  # shrinking the type -- see build.md.
                                  # E_LAB_PT is the LAST constant to
                                  # touch, never the first.
  E_BOXPAD <- 0.20                # ggrepel box.padding (DDB1 only)

  # Panel geometry, MEASURED off the composed render at 400 dpi (row 3
  # middle cell, ROW_H[3] = 3.05 in, COL_W$def[2] = 0.38 of 8.5 in):
  # x ticks at -1.0/-0.5/0.0/0.5 land 219.3 px apart and the y ticks
  # 184.0 px per 0.25, i.e. 438.7 px (1.0967 in) per x-unit and
  # 736.0 px (1.8400 in) per y-unit; the drawn panel is 1014 x 907 px
  # = 2.535 x 2.267 in.
  #
  # RE-MEASURE THESE FIRST whenever compose.R's ROW_H or COL_W change,
  # then re-derive E_PITCH, the bands and the tier columns, and only
  # then consider the point sizes. The history of this pair:
  #   ROW_H[3] 2.45, cell 3.40 -> 1.4718 / 1.3937
  #   ROW_H[3] 2.95, cell 3.40 -> 1.4718 / 1.7600
  #   ROW_H[3] 3.05, cell 3.23 -> 1.0967 / 1.8400
  #   ROW_H[3] 2.75, cell 4.50 -> 1.6688 / 1.6025   (current, 2026-10-05:
  #     d|e two columns, f full-width row 4; x from axis length 3.855 in)
  # The x figure fell by 25%: everything measured in x-units -- the
  # legend block, the label widths, the tier columns -- got 34% WIDER
  # in data units at a stroke, which is what broke the old layout.
  E_IN_PER_X <- 1.6688
  E_IN_PER_Y <- 1.6025
  E_EXPAND   <- 0.05              # ggplot's default continuous mult
  e_pan_x <- E_XLIM + c(-1, 1) * E_EXPAND * diff(E_XLIM)
  e_pan_y <- E_YLIM + c(-1, 1) * E_EXPAND * diff(E_YLIM)

  # Drawn width of an italic gene name at E_LAB_PT, in x-units.
  # 0.630 em/char is MEASURED (MAD1L1 renders 168 px wide at 400 dpi =
  # 0.420 in for 6 characters at 8 pt), not the 0.56 used for upright
  # text elsewhere in this figure.
  E_CHAR_EM  <- 0.630
  lab_w <- function(g) nchar(g) * E_CHAR_EM * E_LAB_PT / 72 / E_IN_PER_X

  # Legend block, measured as an ink bounding box off a 400 dpi render
  # of this legend at this cell size (not with ggplotGrob -- the grob
  # misses ~0.15 in of the width). E_LEG_TOP is the drop from the
  # guide box's top (= the panel top) to the first line's ink.
  E_LEG_W   <- 1.032 / E_IN_PER_X
  E_LEG_TOP <- 0.068 / E_IN_PER_Y
  E_LEG_H   <- 0.983 / E_IN_PER_Y
  E_LEG_NPC_X <- 0.005            # must match legend.position.inside
  e_leg <- list(
    xmin = e_pan_x[1] + E_LEG_NPC_X * diff(e_pan_x),
    ymax = e_pan_y[2] - E_LEG_TOP)
  e_leg$xmax <- e_leg$xmin + E_LEG_W
  e_leg$ymin <- e_leg$ymax - E_LEG_H

  # Series colours. Names double as the legend text, so the wrapping
  # lives here and nothing has to re-map labels later.
  #
  # THE WORDING IS FIXED; only the wrapping is a layout choice, and in
  # a 2.535 in panel the wrapping is load-bearing. A legend is as wide
  # as its widest LINE across ALL entries, so the previous wrap was
  # 1.515 in wide because of "Frameshifted (n = 971)" -- re-wrapping
  # only the two long entries does nothing at all. Measured on the
  # render at this cell size:
  #   2-line RPM entries, 1-line n entries   1.515 x 0.642 in
  #   all four entries wrapped (this one)    1.032 x 0.983 in
  #   n entries wrapped, RPM 2-line          1.407 x 0.812 in
  #   n entries 1-line, RPM 3-line           1.515 x 0.812 in
  # Only the narrow one works. 1.515 in is 1.381 x-units here, so that
  # block reaches x = -0.013 and there is no room at all between it
  # and the frameshift curve for the top tier of names -- which is
  # exactly the overlap this rebuild had to fix. 1.032 in is 0.941
  # x-units, so the block ends at x -0.452 and the top tier clears it.
  # The extra height (0.983 vs 0.642 in) is paid for out of the tier
  # split: see E_SPLIT.

  HIGHLIGHT <- c("SMAD4", "ARID1A", "PALB2", "FANCD2", "DEPDC5", "DDB1",
                 "BRD2", "BRCA1", "RBBP8", "XRCC5", "NHEJ1", "CHFR",
                 "MAD1L1", "MAD2L2", "TOPBP1", "CTC1", "EZH2", "BRD9")

  # ---------------- data ------------------------------------------
  # 20-column slices of the two 13 GB expression matrices, built once
  # into cache/ (see build.md). header=TRUE is required: fread's
  # header sniffing gets these wrong and returns V1..V21.
  read_mat <- function(f) {
    data.table::fread(file.path(CACHE, f), header = TRUE) |>
      as.data.frame() |>
      tibble::column_to_rownames("gene_name") |>
      as.matrix()
  }
  expr_mat     <- read_mat("panel_e_scaled.tsv")   # scaled RPM  (curves)
  expr_mat_raw <- read_mat("panel_e_raw.tsv")      # raw RPM     (matching)
  top_ids      <- readLines(file.path(CACHE, "panel_e_top_ids.txt"))

  df_final <- data.table::fread(file.path(QROOT, "leafcutter_analyses/diff_splicing/unified_results/nkap/20v60/dpsi0.01_fc2.5/df_final.tsv"))

  # -- target / complement -----------------------------------------
  # alt-3'ss hits whose cryptic site sits strictly inside 20 nt
  # UPSTREAM of the annotated 3'ss; the shift is frame-preserving iff
  # it is a multiple of 3. A gene with any frameshifting junction is
  # removed from the in-frame set (target takes precedence).
  hits <- df_final |>
    dplyr::filter(deltapsi > 0.01, test / control > 2.5, type == "alt 3'ss")
  base_hits <- hits |>
    dplyr::filter(!is.na(distance_to_cryptic),
                  distance_to_cryptic > -20, distance_to_cryptic < 0)

  target <- base_hits |>
    dplyr::filter(distance_to_cryptic %% 3 != 0) |>
    dplyr::pull(gene_name) |> unique() |> intersect(rownames(expr_mat))
  complem <- base_hits |>
    dplyr::filter(distance_to_cryptic %% 3 == 0) |>
    dplyr::pull(gene_name) |> unique() |> intersect(rownames(expr_mat)) |>
    setdiff(target)

  if (length(target) != 971L || length(complem) != 866L)
    message("[panel e] group sizes differ from the published 971 / 866: ",
            length(target), " / ", length(complem))
  # Legend keys carry the n, so they are built from the data.
  G_FS  <- sprintf("Frameshifted\n(n = %d)", length(target))
  G_IF  <- sprintf("In-frame\n(n = %d)", length(complem))
  G_CFS <- "RPM-matched\ncontrol genes\nfor frameshifts"
  G_CIF <- "RPM-matched\ncontrol genes\nfor in-frame"
  E_PAL <- c("#B22222", "#4682B4", "#EE5C42", "#9AC0CD")
  names(E_PAL) <- c(G_FS, G_IF, G_CFS, G_CIF)

  # -- RPM-matched controls ----------------------------------------
  # 1:1 nearest-neighbour on log10(median RAW RPM in the top-20
  # samples). matchit's tie-breaking is stochastic -> seed it.
  match_by_rpm <- function(treated, top_samples) {
    universe <- intersect(rownames(expr_mat), rownames(expr_mat_raw))
    med_rpm  <- matrixStats::rowMedians(
      expr_mat_raw[universe, top_samples, drop = FALSE], na.rm = TRUE)
    names(med_rpm) <- universe
    med_rpm <- med_rpm[is.finite(med_rpm) & med_rpm > 0]
    treated <- intersect(treated, names(med_rpm))
    pool    <- setdiff(names(med_rpm), treated)
    d <- data.frame(
      gene_name = c(treated, pool),
      treat     = c(rep(1L, length(treated)), rep(0L, length(pool))),
      log_rpm   = log10(c(med_rpm[treated], med_rpm[pool])))
    m  <- MatchIt::matchit(treat ~ log_rpm, data = d, method = "nearest",
                           distance = "glm", replace = FALSE)
    md <- MatchIt::match.data(m)
    md$gene_name[md$treat == 0]
  }
  set.seed(E_SEED)
  ctrl_fs <- match_by_rpm(target,  top_ids)
  ctrl_if <- match_by_rpm(complem, top_ids)

  # -- curve values: per-gene mean over the 20 top samples ----------
  curve_of <- function(genes, label) {
    m <- apply(expr_mat[genes, top_ids, drop = FALSE], 1, mean, na.rm = TRUE)
    tibble::tibble(gene_name = names(m), x = as.numeric(m), group = label)
  }
  curves <- dplyr::bind_rows(
    curve_of(target,  G_FS),
    curve_of(complem, G_IF),
    curve_of(ctrl_fs, G_CFS),
    curve_of(ctrl_if, G_CIF)) |>
    dplyr::mutate(group = factor(group, levels = names(E_PAL)))

  # -- one-sided Wilcoxon, frameshift < in-frame -------------------
  tx <- curves$x[curves$group == G_FS];  tx <- tx[is.finite(tx)]
  cx <- curves$x[curves$group == G_IF];  cx <- cx[is.finite(cx)]
  pv <- wilcox.test(tx, cx, alternative = "less")$p.value
  p_lab <- sprintf("p = %.2g", pv)
  cat("[panel e] n target =", length(target), " n complement =",
      length(complem), " ", p_lab, "\n")

  # -- gene labels sit at their percentile ON THE FRAMESHIFT CURVE --
  # SMAD4 and DEPDC5 fall in the complement set, so they are (and must
  # be) dropped here -- 16 of the 18 requested names survive.
  fs_curve <- curves |>
    dplyr::filter(group == G_FS) |>
    dplyr::mutate(ecdf_y = rank(x, ties.method = "first") / dplyr::n())
  hl <- fs_curve |>
    dplyr::filter(gene_name %in% HIGHLIGHT, is.finite(x)) |>
    dplyr::select(gene_name, x, ecdf_y)
  cat("[panel e] labelled:", nrow(hl), "of", length(HIGHLIGHT),
      "| suppressed:", paste(setdiff(HIGHLIGHT, hl$gene_name),
                             collapse = ", "), "\n")
  # The percentiles the tier splits are cut at. Printed because
  # E_SPLIT has to be re-cut whenever the tier CAPACITIES change, and
  # the natural gaps in this list are where to cut.
  cat("[panel e] label percentiles:\n",
      paste0(sprintf("%s %.3f", hl$gene_name[order(hl$ecdf_y)],
                     sort(hl$ecdf_y)), collapse = " | "), "\n")

  # ---------------- repel obstacles -------------------------------
  # Only DDB1 is still placed by ggrepel (see the ladder block below),
  # but ggrepel only avoids its OWN layer's data, so the four step
  # curves and the legend are invisible to it. Feed them in as
  # zero-text rows; DDB1 then cannot be pushed onto a curve.
  obstacles <- curves |>
    dplyr::filter(is.finite(x)) |>
    dplyr::group_by(group) |>
    dplyr::arrange(x, .by_group = TRUE) |>
    dplyr::mutate(ecdf_y = dplyr::row_number() / dplyr::n()) |>
    dplyr::ungroup() |>
    dplyr::filter(x > -1.30, x < 0.85) |>
    dplyr::group_by(group) |>
    dplyr::slice(unique(round(seq(1, dplyr::n(), length.out = 120)))) |>
    dplyr::ungroup() |>
    dplyr::transmute(gene_name = "", x, ecdf_y)

  # The legend is invisible to ggrepel too. Tile its rectangle with
  # zero-text rows so a label cannot be parked underneath it even if a
  # layer's ylim would formally allow it.
  leg_block <- expand.grid(
    x      = seq(e_leg$xmin, e_leg$xmax, length.out = 9),
    ecdf_y = seq(e_leg$ymin, e_leg$ymax, length.out = 6)) |>
    dplyr::transmute(gene_name = "", x, ecdf_y)
  obstacles <- dplyr::bind_rows(obstacles, leg_block)

  # ---------------- label ladders ---------------------------------
  # Same idiom as figure 4 panel g: names in a tidy right-aligned
  # column on the free flank of the curve, short leaders running
  # RIGHTWARDS to the point. Because every label in a column sits at
  # the SAME x and every leader leaves that column heading right, no
  # leader can ever cross a label in its own column.
  #
  # Panel g gets all ten of its names into one column; this panel has
  # sixteen in 2.267 in of drawn height, with a legend in the same
  # column. At SIZE_BASE a name plus a readable gutter is 0.110 in
  # tall, so one column tops out at ~20 names with no legend -- but
  # the legend eats the top 0.983 in, so the region left of it holds
  # only ~8.
  #
  # The three tiers step RIGHT as they climb, following the curve's
  # left flank, which is what keeps the leaders short; they do NOT buy
  # vertical room, because they are y-disjoint. What the narrow cell
  # changed is the SPLIT. The legend now reaches down to y 0.565, and
  # tiers A and B sit in its x-shadow (their left edges are at -0.79
  # and -0.45, the legend ends at -0.452), so A and B have to fit
  # BELOW it: y 0.070..0.548, which is 0.478 y-units and holds 7
  # names, not 9. Tier C's column is right of the legend, so it takes
  # the other 8.
  #
  # The splits therefore move from (0.35, 0.70) to (0.35, 0.55), which
  # are still natural gaps in the label percentiles -- the build
  # prints them; BRD9 | XRCC5 and EZH2 | MAD1L1 are the two widest.
  E_SPLIT  <- c(0.35, 0.55)       # 5 / 2 / 8
  E_COL_A  <- -0.47               # right edge, bottom tier
  E_COL_B  <- -0.13               # right edge, middle tier
  E_COL_C  <-  0.00               # right edge, top tier
  E_COL_DD <- -0.95               # right edge of DDB1's corner label

  # Each band is sized by its LABEL COUNT (n x E_PITCH), not by where
  # its points sit, and the slack is spent as the two inter-tier
  # gutters. A band at exactly capacity pins its stack to the band
  # instead of to its points, so keep ~0.01 y-units of slack.
  #   A 5 names -> 0.298   B 2 names -> 0.119   C 8 names -> 0.477
  E_BAND_A <- c(0.070, 0.385)     # 5 names, slack 0.017
  E_BAND_B <- c(0.415, 0.548)     # 2 names, slack 0.014
  E_BAND_C <- c(0.600, 1.110)     # 8 names, slack 0.033
  # Tiers A and B must clear the legend block's bottom edge; tier C is
  # to the RIGHT of it and so is unconstrained here.
  e_a_top  <- min(E_BAND_A[2], e_leg$ymin - 0.030)
  hl_corner <- hl |> dplyr::filter(ecdf_y <  0.05)
  hl_a <- hl |> dplyr::filter(ecdf_y >= 0.05, ecdf_y < E_SPLIT[1])
  hl_b <- hl |> dplyr::filter(ecdf_y >= E_SPLIT[1], ecdf_y < E_SPLIT[2])
  hl_c <- hl |> dplyr::filter(ecdf_y >= E_SPLIT[2])

  # ggrepel pins x when direction = "y", so a per-row nudge is what
  # turns the ladder into a true column. The nudge has to be a plain
  # numeric vector computed HERE (it is stored as a layer parameter and
  # survives saveRDS); the obstacle rows get a zero nudge. Only DDB1
  # still uses this -- see below.
  col_nudge <- function(d, n_obs, x_col) c(x_col - d$x, rep(0, n_obs))

  # The three ladders are placed BY HAND rather than by ggrepel.
  # ggrepel's equilibrium for a column that is packed against its band
  # is "boxes just touching", and box.padding is the only lever on the
  # gutter -- but pushing padding far enough to open a readable gutter
  # makes the stack taller than the band, at which point repel silently
  # lets the boxes overlap instead (measured: padding 0.06/0.11/0.20 ->
  # gutters of 2/4/2-with-merges px at 400 dpi). Even pitch, order
  # preserved, stack centred on the tier's own points and clamped into
  # the band gives an exact, reproducible ladder instead.
  # DERIVED, not measured: one line of cap height plus the gutter the
  # approved build used (0.030 in = 12 px at 400 dpi), converted
  # through E_IN_PER_Y. Stating it as a y-unit constant is what went
  # stale when the row changed shape -- 0.0625 y-units was 0.110 in at
  # E_IN_PER_Y 1.76 and would be 0.115 in now.
  E_GUTTER_IN <- 0.030
  E_PITCH   <- (E_LAB_PT * 0.717 / 72 + E_GUTTER_IN) / E_IN_PER_Y
  E_SEG_GAP <- 0.022 * (1.4718 / E_IN_PER_X)   # ~0.022 in of air
  cat(sprintf("[panel e] pitch %.4f y-units = %.4f in (cap %.4f + gutter %.4f)\n",
              E_PITCH, E_PITCH * E_IN_PER_Y, E_LAB_PT * 0.717 / 72, E_GUTTER_IN))
  cat(sprintf("[panel e] legend block x %.3f..%.3f  y %.3f..%.3f\n",
              e_leg$xmin, e_leg$xmax, e_leg$ymin, e_leg$ymax))
  ladder <- function(d, x_col, band) {
    d <- d[order(d$ecdf_y), , drop = FALSE]
    n     <- nrow(d)
    pitch <- min(E_PITCH, diff(band) / n)
    span  <- pitch * (n - 1)
    ctr   <- min(max(mean(d$ecdf_y), band[1] + pitch / 2 + span / 2),
                 band[2] - pitch / 2 - span / 2)
    d$lab_x <- x_col
    d$lab_y <- ctr - span / 2 + (seq_len(n) - 1) * pitch
    d$seg_x <- x_col + E_SEG_GAP
    d
  }
  lad <- dplyr::bind_rows(
    ladder(hl_a, E_COL_A, c(E_BAND_A[1], e_a_top)),
    ladder(hl_b, E_COL_B, E_BAND_B),
    ladder(hl_c, E_COL_C, E_BAND_C))

  # Clearance audit: the frameshift CDF is the leftmost curve, so a
  # label is clear of ALL four iff the FS curve is right of the column
  # at the label's row. The curve climbs left-to-right, so the tightest
  # row is the box's BOTTOM edge, not its centre -- checking the centre
  # is what let MAD2L2's last glyph sit 0.007 in off the red curve.
  # Printed every build; under ~0.05 units (0.07 in) a tier has to move.
  E_LAB_HALF_H <- (E_LAB_PT * 0.717 / 72) / 2 / E_IN_PER_Y
  fs_x_at <- approxfun(fs_curve$ecdf_y, fs_curve$x, rule = 2)
  lad$clear <- fs_x_at(lad$lab_y - E_LAB_HALF_H) - lad$lab_x
  # ... and clear of the legend block on the LEFT. A label whose box
  # overlaps the legend's rows must start right of its right edge.
  # This is the check the narrow cell broke: six names printed inside
  # the legend text before the rebuild.
  lad$lab_left <- lad$lab_x - lab_w(lad$gene_name)
  over_leg <- lad$lab_y + E_LAB_HALF_H > e_leg$ymin &
              lad$lab_y - E_LAB_HALF_H < e_leg$ymax
  lad$leg_gap <- ifelse(over_leg, lad$lab_left - e_leg$xmax, NA_real_)
  cat("[panel e] ladder clearance (x-units; FS curve at the box bottom,",
      "legend at the box left):\n")
  cat(paste0("   ", sprintf("%-7s y %.3f  curve %.3f  legend %s",
                            lad$gene_name, lad$lab_y, lad$clear,
                            ifelse(is.na(lad$leg_gap), "  -   ",
                                   sprintf("%.3f", lad$leg_gap))),
             collapse = "\n"), "\n")
  stopifnot(all(lad$clear > 0.04))
  stopifnot(all(is.na(lad$leg_gap) | lad$leg_gap > 0.03))

  # p sits in the empty wedge right of the curves, right-aligned at the
  # x limit exactly as in figure 4 panel g. Its ROW had to come down
  # from 0.78 to 0.40 when the cell narrowed, and that is measured, not
  # taste: the string is 0.685 in = 0.625 x-units wide here, so parked
  # at x = 0.80 it reaches back to x = 0.175 -- and at y = 0.78 the
  # rightmost of the four curves is already at x = 0.302 (measured off
  # the 400 dpi render), with RBBP8's leader running to its point at
  # x = 0.20. There is no y near 0.78 where the string fits: clearing
  # the curves would need a right edge at x = 0.955, past the panel.
  # At y = 0.40 the curves reach only x = 0.047, so the string clears
  # them by 0.128 x-units (0.140 in) and clears every leader (tier B's
  # longest ends at x = -0.02).
  ann <- data.frame(x = 0.80, y = 0.40, label = p_lab)

  # ---------------- plot ------------------------------------------
  p <- ggplot() +
    stat_ecdf(data = curves, aes(x = x, color = group),
              linewidth = E_LW, na.rm = TRUE) +
    geom_segment(data = lad,
                 aes(x = seg_x, y = lab_y, xend = x, yend = ecdf_y),
                 inherit.aes = FALSE, linewidth = 0.25,
                 color = "grey55", alpha = 0.9) +
    geom_point(data = hl, aes(x = x, y = ecdf_y),
               size = 0.45, color = "grey15", inherit.aes = FALSE) +
    geom_text(data = lad,
              aes(x = lab_x, y = lab_y, label = gene_name),
              inherit.aes = FALSE, hjust = 1, vjust = 0.5,
              color = "grey15", fontface = "italic",
              family = FONT_FAMILY, size = E_LAB_PT / .pt) +
    ggrepel::geom_text_repel(
      data = dplyr::bind_rows(hl_corner, obstacles),
      aes(x = x, y = ecdf_y, label = gene_name),
      inherit.aes = FALSE, color = "grey15", fontface = "italic",
      family = FONT_FAMILY, size = E_LAB_PT / .pt,
      nudge_x = col_nudge(hl_corner, nrow(obstacles), E_COL_DD),
      nudge_y = 0.055,
      direction = "y", hjust = 1,
      xlim = c(-1.29, E_COL_DD), ylim = c(0.04, 0.09),
      box.padding = E_BOXPAD, min.segment.length = 0, segment.size = 0.25,
      segment.color = "grey55", segment.alpha = 0.9,
      max.overlaps = Inf, max.time = 1e6, seed = E_SEED,
      show.legend = FALSE) +
    geom_text(data = ann, aes(x = x, y = y, label = label),
              inherit.aes = FALSE, hjust = 1, vjust = 0.5,
              family = FONT_FAMILY, size = SIZE_BASE / .pt,
              color = "grey15") +
    scale_color_manual(values = E_PAL, breaks = names(E_PAL),
                       drop = FALSE) +
    scale_x_continuous(breaks = seq(-1, 0.5, 0.5)) +
    scale_y_continuous(breaks = seq(0, 1, 0.25),
                       labels = scales::number_format(accuracy = 0.01)) +
    coord_cartesian(xlim = E_XLIM, ylim = E_YLIM, clip = "on") +
    labs(x = "Mean scaled RPM in top-20 samples",
         y = "Empirical CDF",
         color = NULL,
         title = paste0("Mis-splicing leads to widespread effects\n",
                        "on gene expression")) +
    theme_publication() +
    theme(
      plot.title              = element_text(hjust = 0.5),
      legend.position         = "inside",
      legend.position.inside  = c(E_LEG_NPC_X, 1.00),
      legend.justification    = c(0, 1),
      legend.background       = element_blank(),
      legend.key              = element_blank(),
      legend.text             = element_text(size = E_LEG_PT,
                                             family = FONT_FAMILY,
                                             lineheight = 0.85),
      legend.key.height       = unit(2.8, "mm"),
      legend.key.width        = unit(4.6, "mm"),
      legend.key.spacing.y    = unit(0.4, "mm"),
      legend.margin           = margin(0, 0, 0, 0)
    )

  wrap_elements(full = p)
}))
