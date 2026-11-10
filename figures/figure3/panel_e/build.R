# ------------------------------------------------------------------
# Figure 3, panel e -- "Mis-splicing leads to widespread effects on
# gene expression".
#
# Four empirical CDFs of mean scaled RPM across the top-20 NKAP (GSP
# #13) samples:
#   frameshifting cryptic 3'ss genes, in-frame cryptic 3'ss genes, and
#   an RPM-matched control gene set for each.
#
# Gene sets follow diff_splicing/unified_results df_final.tsv; controls
# are 1:1 RPM-matched with MatchIt.
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
                                  # Shared tokens; grow the cell
                                  # before shrinking E_LAB_PT.
  E_BOXPAD <- 0.20                # ggrepel box.padding (DDB1 only)

  # Panel geometry: inches per data unit in the drawn panel at the
  # current ROW_H / COL_W. Re-measure whenever those change, then
  # re-derive E_PITCH, the bands and the tier columns.
  E_IN_PER_X <- 1.6688
  E_IN_PER_Y <- 1.6025
  E_EXPAND   <- 0.05              # ggplot's default continuous mult
  e_pan_x <- E_XLIM + c(-1, 1) * E_EXPAND * diff(E_XLIM)
  e_pan_y <- E_YLIM + c(-1, 1) * E_EXPAND * diff(E_YLIM)

  # Drawn width of an italic gene name at E_LAB_PT, in x-units.
  # 0.630 em/char is measured for italic text (wider than upright).
  E_CHAR_EM  <- 0.630
  lab_w <- function(g) nchar(g) * E_CHAR_EM * E_LAB_PT / 72 / E_IN_PER_X

  # Legend block size in inches (ink bounding box; ggplotGrob
  # underestimates the width). E_LEG_TOP is the drop from the guide
  # box's top (= the panel top) to the first line's ink.
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
  # lives here and nothing has to re-map labels later. All four entries
  # are wrapped: a legend is as wide as its widest line, and the narrow
  # block leaves room for the top tier of names beside the frameshift
  # curve.

  HIGHLIGHT <- c("SMAD4", "ARID1A", "PALB2", "FANCD2", "DEPDC5", "DDB1",
                 "BRD2", "BRCA1", "RBBP8", "XRCC5", "NHEJ1", "CHFR",
                 "MAD1L1", "MAD2L2", "TOPBP1", "CTC1", "EZH2", "BRD9")

  # ---------------- data ------------------------------------------
  # 20-column slices of the two expression matrices, built once into
  # cache/ by scripts/build_cache_panel_e.sh. header=TRUE is required:
  # fread's header sniffing gets these wrong and returns V1..V21.
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
  # SMAD4 and DEPDC5 fall in the complement set, so they are dropped
  # here -- 16 of the 18 requested names survive.
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
  # Only DDB1 is placed by ggrepel (see the ladder block below), but
  # ggrepel only avoids its OWN layer's data, so the four step curves
  # and the legend are invisible to it. Feed them in as zero-text rows;
  # DDB1 then cannot be pushed onto a curve.
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
  # Names in tidy right-aligned columns on the free flank of the curve,
  # short leaders running RIGHTWARDS to the point, so no leader crosses
  # a label in its own column.
  #
  # Sixteen names do not fit one column beside the legend, so they are
  # split into three tiers that step RIGHT as they climb, following the
  # curve's left flank to keep leaders short. Tiers A and B sit in the
  # legend's x-shadow and must fit below it; tier C is right of it.
  # E_SPLIT is cut at natural gaps in the label percentiles (printed
  # above).
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
  # uses this -- see below.
  col_nudge <- function(d, n_obs, x_col) c(x_col - d$x, rep(0, n_obs))

  # The three ladders are placed BY HAND rather than by ggrepel: in a
  # packed band ggrepel either lets boxes touch or, with enough padding
  # for a gutter, silently overlaps them. Even pitch, order preserved,
  # stack centred on the tier's own points and clamped into the band
  # gives an exact, reproducible ladder.
  # Pitch = one line of cap height plus a 0.030 in gutter, converted
  # through E_IN_PER_Y so it tracks the panel's drawn height.
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
  # row is the box's BOTTOM edge, not its centre.
  # Printed every build; under ~0.05 units (0.07 in) a tier has to move.
  E_LAB_HALF_H <- (E_LAB_PT * 0.717 / 72) / 2 / E_IN_PER_Y
  fs_x_at <- approxfun(fs_curve$ecdf_y, fs_curve$x, rule = 2)
  lad$clear <- fs_x_at(lad$lab_y - E_LAB_HALF_H) - lad$lab_x
  # ... and clear of the legend block on the LEFT. A label whose box
  # overlaps the legend's rows must start right of its right edge.
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
  # x limit. At y = 0.40 the string clears all four curves and every
  # leader; higher up the curves reach too far right for it to fit.
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
