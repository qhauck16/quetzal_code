# Coverage + sashimi plot helpers built directly from raw bigwigs (via
# rtracklayer) and MeanPSI .links files, so every track is a native ggplot
# rather than a rasterised pyGenomeTracks crop.

suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(ggplot2)
  library(patchwork)
  library(data.table)
})

# ------------------------------------------------------------------
# Read a bigwig over a region, average into `n_bins` equal-width bins,
# and return a data.frame with columns bin_start, bin_end, coverage.
# ------------------------------------------------------------------
bigwig_binned <- function(bw_path, chrom, start, end, n_bins = 2000) {
  stopifnot(file.exists(bw_path))
  gr <- GRanges(seqnames = chrom, ranges = IRanges(start, end))
  bw <- rtracklayer::import(bw_path, which = gr, as = "NumericList")[[1]]
  # bw is a numeric vector of length (end - start + 1) with per-bp coverage
  if (length(bw) == 0) return(data.frame(bin_center = numeric(0),
                                          coverage   = numeric(0)))
  # Each bin is INTEGER bp wide. Use floor so the bins exactly tile
  # the region (drop the last up-to-(bin_size-1) bp), then place each
  # bin center at its true genomic position, so coverage stays aligned
  # with the junction arcs.
  n_bp <- length(bw)
  bin_size <- max(1L, as.integer(floor(n_bp / n_bins)))
  n_bins_actual <- as.integer(floor(n_bp / bin_size))
  bw <- bw[seq_len(n_bins_actual * bin_size)]
  bin_mat <- matrix(bw, nrow = bin_size)
  bin_mean <- colMeans(bin_mat)
  bin_center <- start + (seq_len(n_bins_actual) - 1L) * bin_size +
    (bin_size - 1) / 2
  data.frame(bin_center = bin_center,
             coverage   = bin_mean)
}

# ------------------------------------------------------------------
# Build one coverage-track ggplot (single cohort). Returns a ggplot
# that is intentionally minimal: colored area + right-aligned label,
# no axes, no gridlines, no background.
# ------------------------------------------------------------------
coverage_track_gg <- function(bw_path, chrom, start, end,
                              track_color, track_label,
                              y_max = NULL,
                              show_x_axis = FALSE,
                              label_size = 3.2,
                              label_family = "Nimbus Sans",
                              # A real y axis -- rule, a tick at 0 and at
                              # y_max, the numbers at their ticks -- in
                              # the same style as sashimi_track_gg's, so
                              # every coverage track in the paper reads
                              # alike. Off by default.
                              y_axis = FALSE,
                              y_axis_pt = 6) {
  df <- bigwig_binned(bw_path, chrom, start, end)
  if (is.null(y_max)) y_max <- max(df$coverage, na.rm = TRUE) * 1.05
  p <- ggplot(df, aes(x = bin_center, y = coverage)) +
    geom_area(fill = track_color, color = track_color, linewidth = 0.15) +
    scale_x_continuous(limits = c(start, end), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, y_max), expand = c(0, 0)) +
    # Plain face (not bold) and sized to match the structure plot's
    # axis labels, so cohort names read consistently across panel c.
    annotate("text", x = end, y = y_max * 0.65,
              label = track_label, hjust = -0.05, vjust = 0.5,
              size = label_size, color = "grey15",
              family = label_family, fontface = "plain") +
    coord_cartesian(clip = "off") +
    theme_void() +
    # Right gutter holds the overflowing cohort label (hjust = -0.05
    # places it outside the panel). Sized for the largest label at
    # the default label_size.
    theme(plot.margin      = margin(0, 44, 0, 0, "pt"),
          plot.background  = element_rect(fill = NA, color = NA),
          panel.background = element_rect(fill = NA, color = NA))
  if (show_x_axis) {
    p <- p + theme(axis.text.x = element_text(size = 7))
  }
  if (y_axis) {
    p <- p +
      scale_y_continuous(limits = c(0, y_max), expand = c(0, 0),
                         breaks = c(0, y_max),
                         labels = c("0", sprintf("%.1f", y_max))) +
      theme(axis.line.y  = element_line(colour = "grey15", linewidth = 0.3),
            axis.ticks.y = element_line(colour = "grey15", linewidth = 0.3),
            axis.ticks.length.y = grid::unit(2.5, "pt"),
            axis.text.y  = element_text(size = y_axis_pt, colour = "grey15",
                                        family = label_family,
                                        margin = margin(r = 1.5, unit = "pt")),
            # Stacked tracks butt flush; this bottom gutter keeps one
            # track's "0" off the next track's top number.
            plot.margin  = margin(0, 44, 1.3 * y_axis_pt, 0, "pt"))
  }
  p
}

# ------------------------------------------------------------------
# Stack N coverage tracks vertically with matched aspect + optional
# junction arcs drawn above the top track. Returns a ggplot built with
# patchwork; caller can save it via ggsave.
# ------------------------------------------------------------------
coverage_stack <- function(bw_paths, chrom, start, end,
                            track_colors, track_labels,
                            share_ylim = TRUE,
                            arcs_df = NULL,           # data.frame(x_start, x_end, color, width)
                            arc_height_scale = 0.5,
                            label_size = 3.2,
                            label_family = "Nimbus Sans") {
  stopifnot(length(bw_paths) == length(track_colors))
  stopifnot(length(bw_paths) == length(track_labels))
  dfs <- lapply(bw_paths, bigwig_binned,
                 chrom = chrom, start = start, end = end)
  y_max <- if (share_ylim) max(sapply(dfs, function(d) max(d$coverage, na.rm = TRUE))) * 1.05
           else               NA
  tracks <- vector("list", length(bw_paths))
  for (i in seq_along(bw_paths)) {
    tracks[[i]] <- coverage_track_gg(
      bw_paths[i], chrom, start, end,
      track_color = track_colors[i],
      track_label = track_labels[i],
      y_max       = if (share_ylim) y_max else NULL,
      show_x_axis = FALSE,
      label_size  = label_size,
      label_family = label_family
    )
  }
  # optional arc row above the top track
  if (!is.null(arcs_df) && nrow(arcs_df) > 0) {
    arc_plot <- ggplot() +
      scale_x_continuous(limits = c(start, end), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
      coord_cartesian(clip = "off") +
      theme_void() +
      theme(plot.margin = margin(0, 30, 0, 0, "pt"))
    max_w <- max(arcs_df$width, na.rm = TRUE)
    for (j in seq_len(nrow(arcs_df))) {
      x1 <- arcs_df$x_start[j]; x2 <- arcs_df$x_end[j]
      xs <- seq(x1, x2, length.out = 120)
      arc <- sin((xs - x1) / (x2 - x1) * pi)
      thickness <- 0.05 * (arcs_df$width[j] / max_w) * arc_height_scale
      y_top    <- arc * arc_height_scale        + thickness / 2
      y_bottom <- arc * arc_height_scale        - thickness / 2
      poly_df <- data.frame(x = c(xs, rev(xs)),
                             y = c(y_top, rev(y_bottom)))
      arc_plot <- arc_plot + geom_polygon(
        data = poly_df, aes(x = x, y = y),
        fill = arcs_df$color[j], alpha = 1, color = NA)
    }
    tracks <- c(list(arc_plot), tracks)
    heights <- c(0.35, rep(1, length(bw_paths)))
  } else {
    heights <- rep(1, length(bw_paths))
  }
  Reduce(`/`, tracks) + plot_layout(heights = heights)
}

# ------------------------------------------------------------------
# Hand-rolled gene body plot (no ggbio, no ggplot viewport surprises).
# Takes a list of transcripts, each a GRanges of exons, and draws:
#   - one horizontal intron line per transcript at its row y
#   - one filled rectangle per exon
#   - directional arrowheads along the intron line (strand-aware)
#   - italic gene label at the right edge if supplied
# Returns a ggplot with theme_void + transparent backgrounds so it
# stacks cleanly with the coverage tracks above.
# ------------------------------------------------------------------
gene_body_gg <- function(t_models, chrom, start, end,
                          gene_label = NULL,
                          exon_h = 0.35,
                          n_arrows = 15,
                          arrow_size = 0.08,     # arrow length in npc-ish units
                          label_size = 3.4,
                          label_family = "Nimbus Sans") {
  # t_models is a GRangesList (one entry per transcript). Stack them
  # vertically -- one row per transcript.
  n_tx <- length(t_models)
  exon_dfs <- vector("list", n_tx)
  for (i in seq_len(n_tx)) {
    ex <- t_models[[i]]
    ex_start <- pmax(start(ex), start)
    ex_end   <- pmin(end(ex),   end)
    keep <- ex_start < ex_end
    if (!any(keep)) next
    exon_dfs[[i]] <- data.frame(
      x_start = ex_start[keep],
      x_end   = ex_end[keep],
      tx_y    = n_tx - i + 1,    # row from top down
      strand  = as.character(strand(ex))[1]
    )
  }
  exon_df <- do.call(rbind, exon_dfs)
  if (is.null(exon_df) || nrow(exon_df) == 0) {
    return(ggplot() + theme_void())
  }
  # intron backbone lines (one per transcript row)
  tx_lines <- do.call(rbind, lapply(seq_len(n_tx), function(i) {
    data.frame(x_start = start, x_end = end,
                tx_y = n_tx - i + 1,
                strand = as.character(strand(t_models[[i]]))[1])
  }))
  # arrow positions along each backbone
  arrow_df <- do.call(rbind, lapply(seq_len(nrow(tx_lines)), function(i) {
    row <- tx_lines[i, ]
    xs <- seq(row$x_start, row$x_end,
               length.out = n_arrows + 2)[-c(1, n_arrows + 2)]
    data.frame(x = xs, tx_y = row$tx_y,
                strand = row$strand)
  }))

  p <- ggplot() +
    # backbone
    geom_segment(data = tx_lines,
                 aes(x = x_start, xend = x_end,
                      y = tx_y, yend = tx_y),
                 color = "black", linewidth = 0.35) +
    # exon boxes
    geom_rect(data = exon_df,
              aes(xmin = x_start, xmax = x_end,
                   ymin = tx_y - exon_h / 2,
                   ymax = tx_y + exon_h / 2),
              fill = "grey20", color = "grey20") +
    # directional arrows along the intron backbone
    geom_segment(
      data = arrow_df,
      aes(x    = x - ifelse(strand == "-", -1, 1) * (end - start) * arrow_size / 200,
           xend = x + ifelse(strand == "-", -1, 1) * (end - start) * arrow_size / 200,
           y    = tx_y, yend = tx_y),
      arrow = arrow(length = unit(1.6, "mm"), type = "open", ends = "last"),
      color = "grey20", linewidth = 0.35) +
    scale_x_continuous(limits = c(start, end), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.3, n_tx + 0.7), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(plot.margin      = margin(2, 44, 2, 0, "pt"),
          plot.background  = element_rect(fill = NA, color = NA),
          panel.background = element_rect(fill = NA, color = NA))
  if (!is.null(gene_label)) {
    p <- p + annotate("text", x = end, y = (n_tx + 1) / 2,
                       label = gene_label, hjust = -0.05, vjust = 0.5,
                       size = label_size, fontface = "italic",
                       family = label_family, color = "grey15")
  }
  p
}

# ------------------------------------------------------------------
# Custom sashimi: coverage + arcs with PSI-value labels above them.
# ------------------------------------------------------------------
sashimi_track_gg <- function(bw_path, chrom, start, end,
                              track_color, track_label,
                              links_path,
                              min_psi = 1.0,
                              y_max = NULL,
                              label_size = 3.2,
                              psi_size = 2.6,
                              label_family = "Nimbus Sans",
                              # Cohort label placement. Default is the
                              # top-RIGHT in grey, INSIDE the track.
                              # label_below = TRUE instead hangs it under
                              # the coverage baseline, which is the only
                              # placement that is unambiguous when tracks
                              # are stacked: a label inside the track
                              # floats in the arc headroom, nearer to the
                              # track above it than to its own coverage.
                              # The caller must leave room for it in
                              # plot.margin's bottom.
                              label_hjust = 1,
                              label_color = "grey15",
                              label_below = FALSE,
                              # Apex heights of the arc tiers, as
                              # fractions of y_max (shortest arc highest).
                              # Raise the first element to use up dead
                              # air above the arcs.
                              apex_range = c(0.72, 0.37),
                              # Minimum rise of an arc's apex above its
                              # higher endpoint, as a fraction of y_max.
                              # A FLOOR on apex_range that rescues arcs
                              # which would otherwise sag below their own
                              # endpoints (an exon peak taller than the
                              # arc's tier makes `apex - max(y1,y2)`
                              # negative and the curve bows DOWNWARD).
                              #
                              # DEFAULT 0 = off, deliberately. A caller
                              # that back-solves apex_range from chosen
                              # midpoint heights -- so it can place PSI
                              # labels at known y's, as figure 3 panel c
                              # does -- must not have those apexes
                              # silently raised, or its labels end up
                              # struck through by the arcs they belong
                              # to. Only set this when apex_range is a
                              # plain tier spec.
                              min_arc_rise = 0,
                              # Drawn size of THIS track on the page.
                              # Needed to convert the PSI label's point
                              # size into data units for the collision
                              # clearance -- hardcoded data-unit offsets
                              # silently break whenever the track's
                              # drawn size changes.
                              track_height_in = 0.85,
                              track_width_in = 6.8,
                              # Arc stroke = 0.35 + lw_scale * sqrt(PSI/100).
                              # Small insets want ~0.5; full-width tracks 1.0.
                              lw_scale = 1.0,
                              # How far back from a splice site to look
                              # for the "exon top" an arc end is pinned
                              # to (see peak_y below). 200 bp suits
                              # figures 1 and 2; a small value (~15)
                              # instead pins the arc to the coverage AT
                              # the junction, which is what you want when
                              # the exon next door is much taller than
                              # the splice site itself and would
                              # otherwise flatten the arc.
                              # X is never affected -- only the height.
                              anchor_bp = 200,
                              # How the PSI number is rendered. NULL gives
                              # the "96%" / "1.4%" form; pass a
                              # function(psi) -> character for anything
                              # else (figure 3 panel c prints bare
                              # 3-significant-digit values, no % sign).
                              psi_format = NULL,
                              # FALSE draws the arcs but no PSI numbers,
                              # for callers that place the numbers
                              # themselves. The built-in solver only ever
                              # moves a label in Y, at its own arc's
                              # midpoint; two junctions that share an
                              # endpoint put both labels at the same x,
                              # and in a track only a few text-heights
                              # tall there is then no vertical slot that
                              # clears both arcs AND the other label.
                              show_psi = TRUE,
                              # "inline" (default): each PSI number sits ON
                              # its own arc, near the arc's middle, in a gap
                              # cut out of the arc. The number slides along
                              # its arc only if another label or another
                              # arc would touch it.
                              # "offset": number above or below the apex.
                              psi_style = "inline",
                              # Extra boxes an inline label must not touch,
                              # as a list of c(x0, x1, y0, y1) in data units
                              # (e.g. a zoom source box, a "-11 nt" note, a
                              # dashed zoom leader's corridor).
                              psi_obstacles = list(),
                              # Per-arc peak heights, overriding the span-
                              # ranked tiers: a named numeric vector of
                              # fractions of y_max, named by the PSI label
                              # AS DISPLAYED (e.g. c("71%" = 0.40)). For a
                              # track where two arcs share an end and run
                              # together, so no tier order separates them
                              # (figure 2 panel b).
                              apex_override = NULL,
                              # ---- y axis (0 .. y_max) -------------------
                              # A real axis: a vertical rule at the panel's
                              # LEFT edge spanning 0 .. y_max, a tick at each
                              # end, and the two numbers sitting AT their
                              # ticks (vjust 0.5), so the "0" is level with
                              # the coverage baseline instead of floating
                              # above it. Everything is placed in the left
                              # gutter via clip = "off", so the caller must
                              # leave room in plot.margin's LEFT -- and,
                              # where tracks are stacked flush, a few pt in
                              # the BOTTOM too, or this track's "0" and the
                              # next track's y_max number land on the same
                              # page row.
                              y_axis = TRUE,
                              # NULL -> psi_size (the smallest type already
                              # in the track). Pass explicitly to match a
                              # figure's own axis-text token.
                              y_axis_size = NULL,
                              y_axis_fmt = function(v) sprintf("%.1f", v),
                              # Junctions to draw in highlight_color instead
                              # of black: a list of c(start, end) pairs in the
                              # links' coordinates (matched within 2 bp).
                              highlight = NULL,
                              highlight_color = "#D62728",
                              y_axis_color = "grey15",
                              y_axis_linewidth = 0.3,
                              y_axis_tick_pt = 2.5,   # tick length, pt
                              y_axis_gap_pt  = 1.5    # tick -> number, pt
                              ) {
  df <- bigwig_binned(bw_path, chrom, start, end)
  cov_max <- max(df$coverage, na.rm = TRUE)
  # Room above coverage for arcs: coverage in [0, cov_max], arcs above.
  if (is.null(y_max)) y_max <- cov_max * 3.8
  links <- read.table(links_path, header = FALSE, sep = "\t",
                       col.names = c("c1","s1","e1","c2","s2","e2","psi"))
  # keep only links fully inside the shown region and above min_psi
  links <- links[links$s1 >= start & links$e2 <= end &
                  links$psi >= min_psi, , drop = FALSE]
  p <- ggplot(df, aes(x = bin_center, y = coverage)) +
    geom_area(fill = track_color, color = track_color, linewidth = 0.15) +
    # The view is set by coord in BOTH directions, NOT by scale limits:
    # a scale limit DROPS out-of-range values, which would silently
    # delete a label placed below the baseline by label_below -- or the
    # y-axis ticks and numbers, which live LEFT of `start`, in the
    # gutter.
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(start, end), ylim = c(0, y_max),
                    clip = "off") +
    theme_void() +
    theme(plot.margin      = margin(0, 30, 0, 0, "pt"),
          plot.background  = element_rect(fill = NA, color = NA),
          panel.background = element_rect(fill = NA, color = NA))
  if (nrow(links) > 0) {
    # Sort SHORT -> LONG and assign apex heights MONOTONICALLY:
    # shortest arc highest, longest arc lowest, so a long arc never
    # sweeps across shorter arcs and their labels.
    ord <- order(links$e2 - links$s1, decreasing = FALSE)
    links <- links[ord, , drop = FALSE]
    n_arc <- nrow(links)
    # Apex band stops short of y_max so a label placed ABOVE the
    # highest arc still has headroom and cannot reach the track
    # stacked above.
    # Fractions of y_max, leaving headroom ABOVE the top apex so a
    # label can sit above it (narrow arcs have no room inside).
    apex_rows <- if (n_arc == 1) mean(apex_range) * y_max
                 else seq(apex_range[1], apex_range[2],
                          length.out = n_arc) * y_max
    if (!is.null(apex_override)) {
      .lab0 <- if (is.null(psi_format))
                 ifelse(links$psi < 10, sprintf("%.1f%%", links$psi),
                        sprintf("%.0f%%", links$psi))
               else psi_format(links$psi)
      .hit <- .lab0 %in% names(apex_override)
      apex_rows[.hit] <- unname(apex_override[.lab0[.hit]]) * y_max
    }
    # Y-only lookup: search the exon-facing side of the splice site
    # for the peak coverage height, so the arc's HEIGHT matches the
    # local exon top. X is always the raw junction coordinate --
    # never shifted.
    peak_y <- function(x, side) {
      lo <- if (side == "left")  x - anchor_bp else x
      hi <- if (side == "left")  x             else x + anchor_bp
      w <- df$bin_center >= lo & df$bin_center <= hi
      if (!any(w)) return(df$coverage[which.min(abs(df$bin_center - x))])
      max(df$coverage[w])
    }
    # Thickness ~ sqrt(PSI).
    # Build every arc curve FIRST, then place labels -- label
    # positions need to know about all arcs to avoid landing on one.
    curves <- vector("list", n_arc)
    for (i in seq_len(n_arc)) {
      x1 <- links$s1[i]
      x2 <- links$e2[i]
      y1 <- peak_y(x1, "left"); y2 <- peak_y(x2, "right")
      # The apex is a fraction of y_max, but the ENDPOINTS sit on the
      # local coverage tops -- so whenever an exon peak is taller than
      # its arc's tier, `apex_y - max(y1, y2)` goes negative and the
      # "arc" is drawn sagging BELOW its own endpoints. Guarantee a
      # visible rise above the higher endpoint instead. Purely a floor:
      # arcs whose tier already clears their endpoints are untouched.
      # NB the floor must be skipped entirely when off: `max(apex,
      # max(y1,y2) + 0)` still clamps the apex up to the higher
      # endpoint, which draws every arc as a STRAIGHT CHORD.
      apex_y <- if (min_arc_rise > 0)
                  max(apex_rows[i], max(y1, y2) + min_arc_rise * y_max)
                else apex_rows[i]
      # Parametric sine interpolating (x1,y1) -> apex -> (x2,y2) so
      # the curve ends exactly on the coverage tops.
      xs <- seq(x1, x2, length.out = 200)
      t  <- (xs - x1) / (x2 - x1)
      arc_y <- (1 - t) * y1 + t * y2 + sin(t * pi) * (apex_y - max(y1, y2))
      hl <- any(vapply(highlight, function(h)
        abs(links$s1[i] - h[1]) <= 2 && abs(links$e2[i] - h[2]) <= 2, logical(1)))
      curves[[i]] <- list(xs = xs, ys = arc_y, apex = apex_y,
                           mid = (x1 + x2) / 2,
                           lw = 0.35 + lw_scale * sqrt(links$psi[i] / 100),
                           col = if (hl) highlight_color else "black")
    }
    # Arcs are DRAWN after the labels are placed (below), so an inline
    # label can have a gap cut out of its own arc.
    draw_arc <- function(p, cv, cut = NULL) {
      keep <- rep(TRUE, length(cv$xs))
      for (ct in cut)
        keep <- keep & !(abs(cv$xs - ct$x) <= ct$hw & abs(cv$ys - ct$y) <= ct$hh)
      grp <- cumsum(c(TRUE, diff(keep) != 0))
      d   <- data.frame(x = cv$xs, y = cv$ys, g = grp)[keep, , drop = FALSE]
      if (!nrow(d)) return(p)
      p + geom_path(data = d, aes(x = x, y = y, group = g),
                    inherit.aes = FALSE, color = cv$col,
                    linewidth = cv$lw, lineend = "round")
    }

    # Resolve each PSI label to a y that clears EVERY arc. Try
    # candidate offsets (below the apex first, then above, stepping
    # outward) and take the first with no arc within `clearance`.
    # Clearance DERIVED from the label's real page size rather than
    # guessed in data units: convert half the text's point height
    # (plus the arc's stroke) into data units via the track's drawn
    # height. psi_size is in mm, so *.pt gives points.
    psi_pt    <- psi_size * .pt
    track_pt  <- track_height_in * 72
    clearance <- (psi_pt * 0.55 + 1.5) * y_max / track_pt
    # Half-width of the label in DATA-x units. A label has width, so a
    # slot inside a narrow arc's bowl is unusable even when the
    # vertical gap at the label's centre looks fine -- the limbs close
    # in sideways. Checking across the label's x-span catches that.
    x_per_pt  <- (end - start) / (track_width_in * 72)
    half_w_x  <- function(lab) nchar(lab) * 0.52 * psi_pt / 2 * x_per_pt
    # Try just above the apex first, then just below, then further out.
    steps     <- clearance * c(1.15, -1.15, 1.9, -1.9, 2.6, -2.6)
    # Hard ceiling so no label can ride up into the adjacent track.
    y_ceiling <- y_max - clearance
    label_ys  <- numeric(n_arc)
    # One decimal below 10%: rounding a low-PSI junction (e.g. the
    # cryptic-exon inclusion at 3.9%) to an integer throws away the
    # precision that matters most at that end of the range.
    labels    <- if (is.null(psi_format))
                   ifelse(links$psi < 10,
                           sprintf("%.1f%%", links$psi),
                           sprintf("%.0f%%", links$psi))
                 else psi_format(links$psi)
    # Labels already placed, as (x, half-width, y). Near-duplicate
    # junctions (same acceptor, donors a few bp apart) put two labels at
    # almost the same x, so a label must clear PREVIOUS LABELS as well
    # as the arcs -- checking arcs alone let those two overprint.
    placed <- list()
    if (show_psi && identical(psi_style, "inline")) {
      y_per_pt <- y_max / track_pt
      # REAL string widths, measured on cairo (the device the figures are
      # written with); an nchar()-based estimate undershoots this font, so
      # the gap cut into each arc would be too narrow.
      lab_w_pt <- local({
        old <- grDevices::dev.cur(); tf <- tempfile(fileext = ".pdf")
        grDevices::cairo_pdf(tf, width = 4, height = 4)
        on.exit({ grDevices::dev.off(); unlink(tf)
                  if (old > 1) try(grDevices::dev.set(old), silent = TRUE) })
        vapply(labels, function(l) grid::convertWidth(grid::grobWidth(
          grid::textGrob(l, gp = grid::gpar(fontsize = psi_pt,
                                            fontfamily = label_family))),
          "pt", valueOnly = TRUE), numeric(1))
      })
      # Half the INK height of a digits-only label (cap height is ~0.72
      # em, so half is ~0.36 em) plus 0.4 pt of air. Sized to the ink,
      # not the em box: nested arcs in these tracks run only ~5 pt apart,
      # and an em-box clearance made every slot between them look taken.
      hh_lab   <- (psi_pt * 0.36 + 1.2) * y_per_pt
      cuts     <- replicate(n_arc, list(), simplify = FALSE)
      # Longest arcs first: they have the most room to slide.
      ord_i <- order(-(links$e2 - links$s1))
      # Slide up to 0.37 of the span either way: in figure 3 panel c the
      # dashed zoom leaders cross the middle of the High track, and the
      # only clear slots for 43.6 / 70.8 lie further out.
      ts <- 0.5 + c(0, rbind(seq(0.03, 0.37, 0.03), -seq(0.03, 0.37, 0.03)))
      for (i in ord_i) {
        cv <- curves[[i]]
        # half the measured width + 2.5 pt air + the stroke's round cap
        hw <- (lab_w_pt[i] / 2 + 2.5 + curves[[i]]$lw * .pt / 2) * x_per_pt
        best <- NULL
        for (tt in ts) {
          lx <- cv$xs[1] + tt * (cv$xs[length(cv$xs)] - cv$xs[1])
          ly <- approx(cv$xs, cv$ys, xout = lx)$y
          if (ly + hh_lab > y_max * 0.995) next
          hit_lab <- any(vapply(placed, function(q)
            abs(q$x - lx) < q$hw + hw && abs(q$y - ly) < 2 * hh_lab, logical(1)))
          hit_arc <- any(vapply(seq_len(n_arc)[-i], function(j) {
            cj <- curves[[j]]
            any(abs(cj$xs - lx) <= hw & abs(cj$ys - ly) <= hh_lab)
          }, logical(1)))
          hit_obs <- any(vapply(psi_obstacles, function(o)
            lx + hw > o[1] && lx - hw < o[2] &&
            ly + hh_lab > o[3] && ly - hh_lab < o[4], logical(1)))
          if (!hit_lab && !hit_arc && !hit_obs) { best <- c(lx, ly); break }
        }
        notch <- FALSE
        if (is.null(best)) {
          # Second tier (figure 4's rule): clear of every LABEL and
          # obstacle, but allowed to sit where a neighbouring arc passes;
          # that arc is then NOTCHED behind the number, so no line ever
          # runs through the digits.
          for (tt in ts) {
            lx <- cv$xs[1] + tt * (cv$xs[length(cv$xs)] - cv$xs[1])
            ly <- approx(cv$xs, cv$ys, xout = lx)$y
            if (ly + hh_lab > y_max * 0.995) next
            hit_lab <- any(vapply(placed, function(q)
              abs(q$x - lx) < q$hw + hw && abs(q$y - ly) < 2 * hh_lab, logical(1)))
            hit_obs <- any(vapply(psi_obstacles, function(o)
              lx + hw > o[1] && lx - hw < o[2] &&
              ly + hh_lab > o[3] && ly - hh_lab < o[4], logical(1)))
            if (!hit_lab && !hit_obs) { best <- c(lx, ly); notch <- TRUE; break }
          }
        }
        if (is.null(best)) {
          message(sprintf("[sashimi] inline label '%s': no clear slot on its arc, fell back to the midpoint", labels[i]))
          lx <- cv$mid; best <- c(lx, approx(cv$xs, cv$ys, xout = lx)$y)
          notch <- TRUE
        }
        box <- list(x = best[1], y = best[2], hw = hw, hh = hh_lab)
        cuts[[i]] <- c(cuts[[i]], list(box))
        if (notch) for (j in seq_len(n_arc)[-i])
          cuts[[j]] <- c(cuts[[j]], list(box))
        placed[[length(placed) + 1L]] <- list(x = best[1], hw = hw, y = best[2])
        label_ys[i] <- best[2]
        curves[[i]]$mid <- best[1]
      }
      for (i in seq_len(n_arc)) p <- draw_arc(p, curves[[i]], cuts[[i]])
    } else {
      for (i in seq_len(n_arc)) p <- draw_arc(p, curves[[i]], list())
    }
    if (!(show_psi && identical(psi_style, "inline")))
    for (i in seq_len(n_arc)) {
      lx <- curves[[i]]$mid
      hw <- half_w_x(labels[i])
      # Sample EVERY arc (including this one) across the label's whole
      # x-footprint, not just its centre.
      probe <- seq(lx - hw, lx + hw, length.out = 9)
      blockers <- unlist(lapply(seq_len(n_arc), function(j) {
        cj <- curves[[j]]
        px <- probe[probe >= min(cj$xs) & probe <= max(cj$xs)]
        if (length(px) == 0) return(numeric(0))
        approx(cj$xs, cj$ys, xout = px)$y
      }))
      # y of any already-placed label whose x-footprint overlaps ours.
      lab_blockers <- unlist(lapply(placed, function(q) {
        if (abs(q$x - lx) > (q$hw + hw)) NULL else q$y
      }))
      chosen <- NA_real_
      for (s in steps) {
        cand <- curves[[i]]$apex + s
        if (cand < min(df$coverage) || cand > y_ceiling) next
        if (all(abs(blockers - cand) >= clearance) &&
            all(abs(lab_blockers - cand) >= clearance)) {
          chosen <- cand; break
        }
      }
      # Fall back to the point furthest from any arc or placed label.
      if (is.na(chosen)) {
        cands <- curves[[i]]$apex + steps
        cands <- cands[cands > 0 & cands <= y_ceiling]
        allb  <- c(blockers, lab_blockers)
        chosen <- cands[which.max(vapply(cands,
                    function(cd) min(abs(allb - cd)), numeric(1)))]
      }
      label_ys[i] <- chosen
      placed[[length(placed) + 1L]] <- list(x = lx, hw = hw, y = chosen)
    }
    label_df <- data.frame(x = vapply(curves, `[[`, numeric(1), "mid"),
                            y = label_ys,
                            label = labels)
    # PSI labels: match figure font family, larger size, plain-face.
    if (show_psi)
      p <- p + geom_text(data = label_df, aes(x = x, y = y, label = label),
                          inherit.aes = FALSE, size = psi_size,
                          family = label_family, fontface = "plain",
                          color = "grey20")
  }
  # ---- y axis: rule + two ticks + the two numbers -----------------
  # Every value travels as a DATA COLUMN (or a plain layer parameter):
  # an aes() closing over a local scalar renders blank after a cached
  # panel's saveRDS/readRDS round trip.
  ax_size <- if (is.null(y_axis_size)) psi_size else y_axis_size
  if (y_axis) {
    x_per_pt <- (end - start) / (track_width_in * 72)
    tick_x   <- start - y_axis_tick_pt * x_per_pt
    lab_x    <- tick_x - y_axis_gap_pt * x_per_pt
    p <- p +
      # the rule, at the panel's own left edge
      geom_segment(data = data.frame(x = start, xend = start,
                                      y = 0, yend = y_max),
                   aes(x = x, xend = xend, y = y, yend = yend),
                   inherit.aes = FALSE, color = y_axis_color,
                   linewidth = y_axis_linewidth, lineend = "square") +
      # one tick at 0, one at y_max
      geom_segment(data = data.frame(x = c(tick_x, tick_x),
                                      xend = c(start, start),
                                      y = c(0, y_max), yend = c(0, y_max)),
                   aes(x = x, xend = xend, y = y, yend = yend),
                   inherit.aes = FALSE, color = y_axis_color,
                   linewidth = y_axis_linewidth, lineend = "butt") +
      # the numbers, CENTRED ON their ticks (vjust 0.5): the "0" is
      # then level with the coverage baseline, half of it below.
      geom_text(data = data.frame(x = c(lab_x, lab_x),
                                   y = c(y_max, 0),
                                   lab = c(y_axis_fmt(y_max), "0")),
                aes(x = x, y = y, label = lab),
                inherit.aes = FALSE, hjust = 1, vjust = 0.5,
                size = ax_size, family = label_family,
                color = y_axis_color)
  }
  # Track label sits INSIDE the panel at the top-right, so it doesn't
  # rely on plot.margin room and never gets clipped.
  # Sits inside its OWN track, a touch below the very top, so it
  # doesn't hug the baseline of the track stacked above it and read
  # as belonging to that track.
  # Default bottom margin reserves ONE HALF-LINE of axis text under the
  # baseline. Tracks are normally stacked flush, and the "0" is
  # centred ON the baseline while the next track's y_max number is
  # centred on its top tick -- without this gap those two land on the
  # same page row and overprint. A caller that sets plot.margin itself
  # (figure 3 panel c, figure 2 panel b) must leave the room instead.
  p +
    theme(plot.margin = margin(0, 4,
                               if (y_axis) 0.8 * ax_size * .pt else 0,
                               0, "pt")) +
    annotate("text", x = if (label_hjust >= 0.5) end else start,
              # Below the baseline: hang the text's TOP a quarter of a
              # line under y = 0 (vjust = 1 anchors the top, so the gap
              # does not depend on the string). Inside: a touch below
              # the very top of this track.
              y = if (label_below)
                    -0.25 * (label_size * .pt) / 72 * y_max / track_height_in
                  else y_max * 0.99,
              label = track_label, hjust = label_hjust, vjust = 1,
              size = label_size, color = label_color,
              family = label_family, fontface = "plain")
}

# ------------------------------------------------------------------
# Single-transcript "zoom-out" gene body: draws intron backbone, exons,
# strand-aware directional arrows over the FULL transcript range, plus
# an optional translucent highlight rect over the region shown in the
# coverage panel above. Returns a ggplot with theme_void.
# ------------------------------------------------------------------
transcript_zoom_gg <- function(exons, x_lo, x_hi,
                                highlight = NULL,    # c(lo, hi) in genomic bp
                                highlight_fill = "#B0C4FF",
                                highlight_alpha = 0.45,
                                tx_label = NULL,
                                exon_h = 0.35,
                                n_arrows = 20,
                                label_size = 3.2,
                                label_family = "Nimbus Sans",
                                # --- optional genomic-axis ticks -------------
                                # Accepted but IGNORED: the gene model
                                # draws no coordinate ticks (they read as
                                # thin exons). Kept for caller compatibility.
                                x_ticks = NULL,
                                tick_h = 0.30,
                                tick_color = "grey25",
                                tick_linewidth = 0.3,
                                # --- which end the gene label hangs off ------
                                # "right" (default) puts it over the
                                # transcript's 3'-most edge.
                                # "left" is for genes whose terminal exon at
                                # that end is a fat UTR block the label would
                                # otherwise sit on (e.g. SENP2).
                                label_side = c("right", "left")) {
  label_side <- match.arg(label_side)
  strand_s <- as.character(strand(exons))[1]
  ex_df <- data.frame(x_start = start(exons), x_end = end(exons))
  # Backbone and directional arrows ONLY span the transcript's own
  # extent, so a short isoform (e.g. NF1-216) does not paint arrows or
  # a backbone into the region past its edge exons.
  tx_lo <- min(start(exons))
  tx_hi <- max(end(exons))
  intron_df <- data.frame(x_start = tx_lo, x_end = tx_hi)
  n_arrows_tx <- max(2, round(n_arrows * (tx_hi - tx_lo) / (x_hi - x_lo)))
  arrow_xs <- seq(tx_lo, tx_hi,
                   length.out = n_arrows_tx + 2)[-c(1, n_arrows_tx + 2)]
  arrow_df <- data.frame(x = arrow_xs, strand = strand_s)
  p <- ggplot() +
    (if (!is.null(highlight))
       annotate("rect", xmin = highlight[1], xmax = highlight[2],
                 ymin = 0, ymax = 2, fill = highlight_fill,
                 alpha = highlight_alpha, color = NA)
     else NULL) +
    geom_segment(data = intron_df,
                 aes(x = x_start, xend = x_end, y = 1, yend = 1),
                 color = "black", linewidth = 0.35) +
    geom_rect(data = ex_df,
              aes(xmin = x_start, xmax = x_end,
                   ymin = 1 - exon_h / 2, ymax = 1 + exon_h / 2),
              fill = "grey20", color = "grey20") +
    geom_segment(
      data = arrow_df,
      aes(x    = x - ifelse(strand == "-", -1, 1) * (x_hi - x_lo) / 200,
           xend = x + ifelse(strand == "-", -1, 1) * (x_hi - x_lo) / 200,
           y    = 1, yend = 1),
      arrow = arrow(length = unit(1.4, "mm"), type = "open", ends = "last"),
      color = "grey20", linewidth = 0.35) +
    # Tick y-extents are carried as DATA COLUMNS, never as scalars closed
    # over by aes(): a cached panel is restored with readRDS into a
    # different environment, where tidy evaluation can no longer find
    # them, and the whole panel renders blank.
    # No coordinate ticks on the gene model: any vertical rule near the
    # backbone reads as a thin exon. x_ticks is accepted and ignored; the
    # coordinate numbers are drawn by the callers' axis rows.
    scale_x_continuous(limits = c(x_lo, x_hi), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 2), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_void() +
    # Zero top/bottom margins so the highlight rects flow continuously
    # between the trapezoid above and consecutive transcript rows.
    theme(plot.margin      = margin(0, 4, 0, 0, "pt"),
          plot.background  = element_rect(fill = NA, color = NA),
          panel.background = element_rect(fill = NA, color = NA))
  if (!is.null(tx_label)) {
    # Label sits ON TOP of this transcript, right-aligned on the
    # transcript's OWN right end (tx_hi, not the panel edge) -- so a
    # short isoform's name stays over the isoform it names.
    #
    # vjust = 0 anchors the text's BOTTOM edge at y, which is what
    # guarantees clearance above the exon boxes (top edge at
    # 1 + exon_h/2). Centring with vjust = 0.5 would make the gap depend
    # on the rendered text height, which is unknown here.
    p <- p + annotate("text",
                       x     = if (label_side == "right") tx_hi else tx_lo,
                       y     = 1 + exon_h / 2 + 0.16,
                       hjust = if (label_side == "right") 1 else 0,
                       vjust = 0,
                       label = tx_label, size = label_size,
                       family = label_family, fontface = "italic",
                       color = "grey15")
  }
  p
}

# ------------------------------------------------------------------
# Trapezoid "zoom connector": wide at top (full transcript range),
# narrows to the highlighted window at bottom. Draw ABOVE the
# transcript row(s); its x-scale is the WIDE transcript range so
# patchwork aligns its bottom corners with the transcript highlight.
# ------------------------------------------------------------------
trapezoid_gg <- function(x_lo, x_hi, win_lo, win_hi,
                          fill = "#B0C4FF", alpha = 0.45) {
  ggplot() +
    annotate("polygon",
             x = c(x_lo, x_hi, win_hi, win_lo),
             y = c(1, 1, 0, 0),
             fill = fill, alpha = alpha, color = NA) +
    scale_x_continuous(limits = c(x_lo, x_hi), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_void() +
    # Zero top/bottom margin so the trapezoid butts flush against the
    # transcript row's highlight rect immediately below -- no visible
    # gap between the tapering shape and the highlight it points to.
    theme(plot.margin      = margin(0, 4, 0, 0, "pt"),
          plot.background  = element_rect(fill = NA, color = NA),
          panel.background = element_rect(fill = NA, color = NA))
}

# ------------------------------------------------------------------
# ZOOM CALLOUT ROW  (figure 3 panel c)
# ------------------------------------------------------------------
# One ggplot that sits directly ABOVE a stack of coverage/sashimi
# tracks and holds N boxed "magnifier" insets. Each inset re-reads its
# cohort's bigwig over a NARROW sub-region and draws that coverage at
# nucleotide scale inside a bordered box; two straight leader lines run
# from the box's bottom corners down to the top corners of a small
# source box drawn ON the corresponding wide track.
#
# Why it is one plot and not several
# ----------------------------------
# The leaders have to cross from the callout row into (and past) the
# track rows. patchwork gives each row its own panel, so the only way
# to draw a line that spans rows is to draw it in ONE plot with
# `clip = "off"` and let it spill. That fixes the coordinate system:
#
#   x  = the WIDE region's genomic coordinates (identical to the
#        tracks below, so patchwork aligns them and the source box's
#        x-extent is DERIVED from the zoom coordinates -- never
#        hand-placed).
#   y  = 0 at the callout row's panel BOTTOM, 1 at its TOP. A callout's
#        `src_drop` is therefore measured in units of the callout row's
#        own height: 0 is the top of the first track below, and
#        h_track/h_callout is the top of the second, etc. Expressing
#        the drop as a RATIO of patchwork heights (rather than in
#        inches) is what makes it survive the figure being rescaled --
#        patchwork's `heights` are relative, so only ratios are stable.
#
# Arguments
# ---------
# callouts   list of lists, one per inset, each with:
#              bw        path to the bigwig to render inside the box
#              color     fill for that coverage
#              src_drop  y (in callout-row heights, POSITIVE = down)
#                        of the top edge of this inset's source box
#              box       c(lo, hi) as FRACTIONS of the wide region,
#                        giving where the box sits horizontally
# zoom_start/zoom_end   the magnified sub-region (bp)
# nt_label   e.g. "+9 nt"; drawn inside every box with a leader curve
# nt_site    bp the leader curve points at (use the midpoint of the
#            canonical and cryptic splice sites)
# cov_frac   fraction of the box interior the coverage may use; the
#            rest is reserved headroom for nt_label, so the label can
#            never land on the coverage.
#
# Returns a ggplot; give it zero top/bottom plot.margin and the same
# left/right margins as the tracks so the x scales line up.
# ------------------------------------------------------------------
zoom_callout_gg <- function(callouts,
                             chrom, wide_start, wide_end,
                             zoom_start, zoom_end,
                             nt_label = NULL,
                             nt_site  = NULL,
                             box_y    = c(0.28, 0.99),
                             cov_frac = 0.66,
                             box_color = "black",
                             box_fill  = "white",
                             box_linewidth = 0.45,
                             leader_color = "grey35",
                             leader_linetype = "22",
                             leader_linewidth = 0.32,
                             nt_curve_curvature = 0.30,
                             label_size = 2.4,
                             label_family = "Nimbus Sans",
                             label_color = "grey10",
                             n_bins = 600) {
  wide_w <- wide_end - wide_start
  # Source box: the zoom window mapped into the wide window, CLAMPED
  # (PRKDC's zoom starts 7 bp before its wide region does).
  src_lo <- max(zoom_start, wide_start)
  src_hi <- min(zoom_end,   wide_end)

  rects <- list(); covs <- list(); leads <- list()
  labs_df <- list(); curves <- list()
  for (i in seq_along(callouts)) {
    co    <- callouts[[i]]
    bx_lo <- wide_start + co$box[1] * wide_w
    bx_hi <- wide_start + co$box[2] * wide_w
    by_lo <- box_y[1]; by_hi <- box_y[2]
    rects[[i]] <- data.frame(xmin = bx_lo, xmax = bx_hi,
                              ymin = by_lo, ymax = by_hi)

    d <- bigwig_binned(co$bw, chrom, zoom_start, zoom_end, n_bins = n_bins)
    # Map bp -> box x, coverage -> box y. The coverage is scaled to its
    # OWN maximum inside the box: an inset exists to show shape at
    # nucleotide resolution, and sharing the wide track's y_max would
    # flatten it.
    cmax  <- max(d$coverage, na.rm = TRUE)
    if (!is.finite(cmax) || cmax <= 0) cmax <- 1
    inner_lo <- by_lo + 0.035 * (by_hi - by_lo)
    inner_hi <- by_lo + cov_frac * (by_hi - by_lo)
    px <- bx_lo + (d$bin_center - zoom_start) / (zoom_end - zoom_start) *
          (bx_hi - bx_lo)
    py <- inner_lo + d$coverage / cmax * (inner_hi - inner_lo)
    covs[[i]] <- data.frame(
      x    = c(px[1], px, px[length(px)]),
      y    = c(inner_lo, py, inner_lo),
      fill = co$color, grp = i)

    leads[[i]] <- data.frame(
      x    = c(bx_lo, bx_hi),
      y    = c(by_lo, by_lo),
      xend = c(src_lo, src_hi),
      yend = c(-co$src_drop, -co$src_drop))

    if (!is.null(nt_label) && !is.null(nt_site)) {
      sx <- bx_lo + (nt_site - zoom_start) / (zoom_end - zoom_start) *
            (bx_hi - bx_lo)
      ly <- by_lo + 0.855 * (by_hi - by_lo)
      # Label goes on the roomier side of the splice site and the
      # curve reaches back to it, so the text never sits on the notch
      # it is describing.
      on_left <- (sx - bx_lo) > (bx_hi - sx)
      lx <- if (on_left) sx - 0.10 * (bx_hi - bx_lo)
            else         sx + 0.10 * (bx_hi - bx_lo)
      labs_df[[i]] <- data.frame(x = lx, y = ly, label = nt_label,
                                  hjust = if (on_left) 1 else 0)
      curves[[i]] <- data.frame(
        x = lx + (if (on_left) -0.012 else 0.012) * (bx_hi - bx_lo),
        y = ly,
        xend = sx,
        yend = inner_hi * 0.5 + inner_lo * 0.5)
    }
  }
  rect_df <- do.call(rbind, rects)
  cov_df  <- do.call(rbind, covs)
  lead_df <- do.call(rbind, leads)

  p <- ggplot() +
    geom_segment(data = lead_df,
                 aes(x = x, y = y, xend = xend, yend = yend),
                 color = leader_color, linewidth = leader_linewidth,
                 linetype = leader_linetype) +
    geom_rect(data = rect_df,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = box_fill, color = NA) +
    geom_polygon(data = cov_df,
                 aes(x = x, y = y, group = grp, fill = fill),
                 color = NA) +
    scale_fill_identity() +
    geom_rect(data = rect_df,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = NA, color = box_color, linewidth = box_linewidth)
  if (length(labs_df) > 0) {
    lab_df <- do.call(rbind, labs_df); crv_df <- do.call(rbind, curves)
    p <- p +
      geom_curve(data = crv_df,
                 aes(x = x, y = y, xend = xend, yend = yend),
                 curvature = nt_curve_curvature, color = label_color,
                 linewidth = 0.3, inherit.aes = FALSE) +
      geom_text(data = lab_df,
                aes(x = x, y = y, label = label, hjust = hjust),
                vjust = 0.5, size = label_size, family = label_family,
                color = label_color, inherit.aes = FALSE)
  }
  p +
    scale_x_continuous(limits = c(wide_start, wide_end), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, 1), clip = "off") +
    theme_void() +
    theme(plot.margin      = margin(0, 4, 0, 0, "pt"),
          plot.background  = element_rect(fill = NA, color = NA),
          panel.background = element_rect(fill = NA, color = NA))
}

# ------------------------------------------------------------------
# The other half of a zoom callout: the small box drawn ON a wide
# track marking the sub-region the callout magnifies. Add to a
# `sashimi_track_gg()` with `+`. `y_max` must be the track's own y_max
# so the box spans the full track height and its TOP corners coincide
# with the panel's top edge -- which is where `zoom_callout_gg()`'s
# leaders are aimed.
# ------------------------------------------------------------------
zoom_source_box <- function(zoom_start, zoom_end, wide_start, wide_end,
                             y_max, color = "black", linewidth = 0.45) {
  geom_rect(data = data.frame(xmin = max(zoom_start, wide_start),
                               xmax = min(zoom_end,   wide_end),
                               ymin = 0, ymax = y_max),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = NA, color = color, linewidth = linewidth,
            inherit.aes = FALSE)
}
