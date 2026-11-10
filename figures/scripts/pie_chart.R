# Pie chart with EVERY percentage drawn at its own wedge. Shared by
# figure 3 panel d and figure 4 panel d (figure 2 panel e has its own
# tuned copy).
#
# Percentages are never in the legend and never on leader lines. Every
# label sits on or beside its own wedge; a label too big for its wedge
# gets a filled rectangle behind it so it stays legible.
#
# WHY CARTESIAN, NOT coord_polar: polar coords give no control over
# where a label lands. Drawing the wedges as polygons under coord_equal()
# makes every label position an ordinary x/y whose bounding box can be
# compared against its neighbours'.
#
# Wedge order is clockwise from 12 o'clock in REVERSE level order,
# which reproduces coord_polar(theta = "y") + position_stack (first
# level ends up last clockwise).

# Clockwise angle from 12 o'clock -> cartesian. The math angle is
# a = pi/2 - theta.
cw2xy <- function(theta, r) cbind(x = r * cos(pi / 2 - theta),
                                  y = r * sin(pi / 2 - theta))

# Drawn width and cap height of each string, in inches, at the size it
# will actually be set. Guessing from nchar() is not good enough: the
# labels are digits plus "." plus "%", whose advance widths differ by
# 3x in Helvetica, and the fit test below turns a few percent of width
# error into a wrongly boxed (or unboxed) label.
#
# MEASURE ON cairo_pdf, which is the device the figures are written
# with. The base pdf() device uses the built-in Helvetica AFM metrics,
# but cairo resolves "Helvetica" through fontconfig to a substitute that
# is 12-18% wider.
.pie_str_size_in <- function(labels, pt, family) {
  old <- grDevices::dev.cur()
  tf  <- tempfile(fileext = ".pdf")
  grDevices::cairo_pdf(tf, width = 7, height = 7)
  on.exit({
    grDevices::dev.off()
    if (old > 1) try(grDevices::dev.set(old), silent = TRUE)
    unlink(tf)
  }, add = TRUE)
  g <- lapply(labels, function(s)
    grid::textGrob(s, gp = grid::gpar(fontsize = pt, fontfamily = family)))
  list(w = vapply(g, function(x)
         grid::convertWidth(grid::grobWidth(x), "in", valueOnly = TRUE), 0),
       h = vapply(g, function(x)
         grid::convertHeight(grid::grobHeight(x), "in", valueOnly = TRUE), 0))
}

# df must have: a category column, a `prop` column summing to 1, and a
# `label` column (the text to draw). `levels_order` fixes the fill
# order; wedges are drawn in reverse of it, clockwise from 12 o'clock.
#
#   pie_r_in       the pie's DRAWN radius in inches. MEASURE it off the
#                  render and pass it -- it is what converts the label's
#                  point size into pie radii, and every fit/collision
#                  decision below depends on that conversion. It cannot
#                  be computed here: coord_equal makes the pie
#                  height-limited and the drawn height is only settled
#                  by patchwork at layout time.
#   r_inside       preferred label radius when the label fits inside
#   r_box          radius a boxed label's search STARTS at
#   r_box_max      how far out the collision search may push a box
#   ang_tol_deg    how far off its own mid-ray a boxed label may sit.
#                  Widened per wedge to that wedge's half-angle. This
#                  is what buys the horizontal space; with 0 the
#                  solver only slides along the radius.
#   x_max          optional hard right-hand bound for a boxed label's
#                  right EDGE, in pie radii from the centre. For a pie
#                  inset beside something tall (fig 3 panel d's 650-count
#                  bar): chips that cannot go right go UP instead.
#   start_deg      rotates the whole pie clockwise by this many degrees
#                  (negative = anticlockwise). Wedge ORDER is unchanged.
#                  Use it to centre a cluster of small wedges on 12
#                  o'clock so their chips fan symmetrically upward.
#   disc_pad       clear space between a boxed label and the pie's rim.
#                  A boxed label is never allowed to touch the disc, so
#                  it can never cover a wedge it does not name.
#   box_pad_pt     padding between the text's ink and the box edge
#   gap_pt         clear space required between two label boxes
#   box_colour     border colour for boxed labels. NULL (the default)
#                  draws each box in ITS OWN WEDGE'S colour -- see the
#                  note at the geom_rect below for why that matters.
pie_labels_on_wedges <- function(df, cat_col, palette,
                                 pie_r_in     = 0.70,
                                 label_pt     = 6.8,
                                 family       = "Nimbus Sans",
                                 levels_order = NULL,
                                 r_inside     = 0.58,
                                 r_box        = 1.00,
                                 r_box_max    = 2.60,
                                 r_step       = 0.01,
                                 ang_tol_deg  = 90,
                                 disc_pad     = 0.015,
                                 on_disc      = TRUE,
                                 r_pref       = 0.80,
                                 r_min_box    = 0.30,
                                 chip_r_pt    = 1.6,
                                 rim_cost     = 0,
                                 r_on_max     = 0.92,
                                 x_max        = Inf,
                                 start_deg    = 0,
                                 box_pad_pt   = 1.3,
                                 gap_pt       = 2.4,
                                 box_fill     = "white",
                                 box_colour   = NULL,
                                 box_lw       = 0.6) {
  stopifnot(all(c("prop", "label") %in% names(df)))
  df$.cat <- as.character(df[[cat_col]])
  if (is.null(levels_order)) levels_order <- sort(unique(df$.cat))
  df  <- df[match(levels_order, df$.cat), ]
  df  <- df[!is.na(df$.cat), ]
  ord <- df[rev(seq_len(nrow(df))), ]          # clockwise draw order
  n   <- nrow(ord)

  ord$th0 <- c(0, head(cumsum(ord$prop), -1)) * 2 * pi + start_deg * pi / 180
  ord$th1 <- cumsum(ord$prop) * 2 * pi + start_deg * pi / 180
  ord$thm <- (ord$th0 + ord$th1) / 2           # the wedge's mid-ray
  ord$ha  <- (ord$th1 - ord$th0) / 2           # half-angle

  wedge <- do.call(rbind, lapply(seq_len(n), function(i) {
    th <- seq(ord$th0[i], ord$th1[i], length.out = 120)
    xy <- cw2xy(th, 1)
    data.frame(.cat = ord$.cat[i], x = c(0, xy[, "x"]), y = c(0, xy[, "y"]))
  }))

  # ---- label half-extents, in pie radii ---------------------------
  pad   <- (box_pad_pt / 72) / pie_r_in
  gap   <- (gap_pt     / 72) / pie_r_in
  sz    <- .pie_str_size_in(ord$label, label_pt, family)
  # grobHeight is the string's INK height, not the em box: the labels
  # are digits/./% with no descender, so this is the cap height. Using
  # the em box would inflate every box by ~40% and start a chain
  # reaction in the collision search.
  ord$hw <- (sz$w / pie_r_in) / 2 + pad       # half-width  incl. padding
  ord$hh <- (sz$h / pie_r_in) / 2 + pad       # half-height incl. padding

  ux <- sin(ord$thm); uy <- cos(ord$thm)      # unit mid-ray (cw from 12)

  # ---- does the label fit INSIDE its wedge? -----------------------
  # Two constraints on the label-box centre, at radius r on the mid-ray:
  #   edges: the perpendicular distance to each straight edge is
  #          r*sin(halfangle); the box's support in that direction is
  #          hw*|cos th| + hh*|sin th| for an edge at clockwise angle
  #          th. -> r >= support / sin(halfangle).
  #   rim:   the outermost corner must stay inside radius 1, i.e.
  #          (|ux|r + hw)^2 + (|uy|r + hh)^2 <= 1, which is a quadratic
  #          in r with |u| = 1.
  # A wedge wider than a half circle is not convex at the apex, so the
  # edge constraint does not apply -- only the rim binds there.
  supp <- function(th, hw, hh) hw * abs(cos(th)) + hh * abs(sin(th))
  r_min <- pmax(supp(ord$th0, ord$hw, ord$hh),
                supp(ord$th1, ord$hw, ord$hh)) / pmax(sin(ord$ha), 1e-9)
  r_min[ord$ha >= pi / 2] <- 0
  b     <- abs(ux) * ord$hw + abs(uy) * ord$hh
  disc  <- b^2 - (ord$hw^2 + ord$hh^2 - 1)
  r_max <- ifelse(disc > 0, -b + sqrt(pmax(disc, 0)), -1)
  ord$fits <- r_max > 0 & r_min <= r_max

  ord$r <- NA_real_
  ord$r[ord$fits] <- pmin(pmax(r_inside, r_min[ord$fits]), r_max[ord$fits])

  # ---- boxed labels: a 2-D search, not a slide along the ray ------
  # Pinning every box to its own mid-ray leaves the radius as the only
  # free parameter; with several adjacent narrow wedges that stacks the
  # boxes into a tight fan and parks inner ones on neighbouring wedges.
  #
  # So each boxed label searches a polar NEIGHBOURHOOD of its wedge
  # -- an angular offset as well as a radius -- under two hard
  # constraints and one soft one:
  #   HARD  the box may not intersect the pie disc at all, so a chip
  #         can never cover a wedge that is not its own;
  #   HARD  no box may touch another box (or an inside label);
  #   SOFT  the centre stays within `ang_tol` of its own mid-ray, so
  #         the chip still reads as belonging to that wedge -- widened
  #         to the wedge's own half-angle when that is larger.
  # Cost prefers small radius then small angular deviation, and labels
  # are placed in ANGULAR ORDER, which makes the result a clean fan
  # that tracks the wedge order instead of an arbitrary stack.
  boxed <- which(!ord$fits)
  if (length(boxed)) {
    # Distance from the origin to the nearest point of an axis-aligned
    # box. Zero when the origin is inside it.
    box_r0 <- function(cx, cy, hw, hh)
      sqrt(max(0, abs(cx) - hw)^2 + max(0, abs(cy) - hh)^2)

    obs <- lapply(which(ord$fits), function(i)
      c(ord$r[i] * ux[i], ord$r[i] * uy[i], ord$hw[i], ord$hh[i]))

    ang_tol <- ang_tol_deg * pi / 180
    das <- c(0, outer(seq_len(ceiling(ang_tol / (2 * pi / 180))),
                      c(-1, 1))) * (2 * pi / 180)

    # Points along the OUTER part of every boxed (small) wedge's mid-ray.
    # Each sliver should stay visible at the rim with its chip beside it
    # or further in -- so a candidate chip pays for every one of these
    # points it would cover, its own wedge's included.
    rim <- do.call(rbind, lapply(boxed, function(k)
      cbind(x = seq(0.70, 1.0, by = 0.05) * ux[k],
            y = seq(0.70, 1.0, by = 0.05) * uy[k])))
    px <- py <- rep(NA_real_, nrow(ord))
    for (i in boxed[order(ord$thm[boxed])]) {
      tol  <- max(ang_tol, ord$ha[i])
      best <- NULL; best_cost <- Inf
      for (da in das[abs(das) <= tol + 1e-9]) {
        th <- ord$thm[i] + da
        vx <- sin(th); vy <- cos(th)
        rr <- if (on_disc) r_min_box else r_box
        while (rr <= r_box_max) {
          cx <- rr * vx; cy <- rr * vy
          if ((if (on_disc) rr <= r_on_max
               else box_r0(cx, cy, ord$hw[i], ord$hh[i]) >= 1 + disc_pad) &&
              cx + ord$hw[i] <= x_max) {
            ok <- TRUE
            for (o in obs)
              if (abs(cx - o[1]) < ord$hw[i] + o[3] + gap &&
                  abs(cy - o[2]) < ord$hh[i] + o[4] + gap) {
                ok <- FALSE; break
              }
            if (ok) {
              # The chip sits ON the pie, as close as possible to
              # its own wedge without overlapping another label. The
              # cost is the chip centre's distance from its OWN wedge
              # (zero anywhere inside the wedge), plus a light pull
              # towards r_pref so a chip inside its wedge sits nearer the
              # rim than the apex.
              dang <- max(0, abs(da) - ord$ha[i])
              cost <- (if (on_disc) rr * sin(min(dang, pi / 2)) +
                         0.15 * abs(rr - r_pref)
                       else abs(rr - r_box) + 1.4 * abs(da))
              if (on_disc)
                cost <- cost + rim_cost *
                  sum(abs(rim[, 1] - cx) < ord$hw[i] &
                      abs(rim[, 2] - cy) < ord$hh[i])
              if (cost < best_cost) { best_cost <- cost; best <- c(cx, cy) }
              if (!on_disc) break   # off-disc: radius only grows
            }
          }
          rr <- rr + r_step
        }
      }
      if (is.null(best)) {
        # Nothing legal inside the tolerance: fall back to
        # straight out along the ray rather than dropping
        # the number, and say so.
        best <- c(r_box_max * ux[i], r_box_max * uy[i])
        warning(sprintf("pie_labels_on_wedges: no clear slot for '%s'",
                        ord$label[i]), call. = FALSE)
      }
      px[i] <- best[1]; py[i] <- best[2]
      obs[[length(obs) + 1]] <- c(best[1], best[2], ord$hw[i], ord$hh[i])
    }
    ord$px <- px; ord$py <- py
  } else {
    ord$px <- ord$py <- NA_real_
  }
  # Inside labels keep their polar placement; boxed ones use the solved
  # 2-D point.
  ord$lx <- ifelse(is.na(ord$px), ord$r * ux, ord$px)
  ord$ly <- ifelse(is.na(ord$py), ord$r * uy, ord$py)

  # Everything is a COLUMN of the layer data -- an aes() closing over a
  # local scalar renders blank after the panel cache's readRDS.
  # A box pushed out along its mid-ray no longer touches its wedge, and
  # adjacent slivers' boxes become ambiguous. Colouring each box's
  # BORDER with its own wedge's fill identifies it without a leader
  # line. White fill stays, so the number is still black-on-white.
  lab <- data.frame(
    label = ord$label,
    x     = ord$lx,
    y     = ord$ly,
    hw    = ord$hw,
    hh    = ord$hh,
    boxed = !ord$fits,
    bcol  = if (is.null(box_colour)) unname(palette[ord$.cat])
            else rep(box_colour, nrow(ord)),
    bfill = unname(palette[ord$.cat]),
    stringsAsFactors = FALSE)
  lab$xmin <- lab$x - lab$hw; lab$xmax <- lab$x + lab$hw
  lab$ymin <- lab$y - lab$hh; lab$ymax <- lab$y + lab$hh
  box <- lab[lab$boxed, , drop = FALSE]

  # Limits = the pie plus whatever the labels actually reach, computed
  # rather than padded by a guess, so the pie is drawn as large as the
  # inset box allows. clip = "off" stays on as a belt-and-braces.
  e   <- 0.02
  xl  <- c(min(-1, lab$xmin) - e, max(1, lab$xmax) + e)
  yl  <- c(min(-1, lab$ymin) - e, max(1, lab$ymax) + e)

  p <- ggplot() +
    geom_polygon(data = wedge, aes(x = x, y = y, fill = .cat, group = .cat),
                 colour = "white", linewidth = 0.35)
  # Labels that fit inside their wedge: plain text. Labels that do not:
  # a rounded CHIP FILLED WITH ITS OWN WEDGE'S COLOUR, no border, black
  # text. The fill is what ties a chip to its wedge, so a chip may sit
  # across neighbouring wedges.
  inl <- lab[!lab$boxed, , drop = FALSE]
  if (nrow(inl))
    p <- p + geom_text(data = inl, aes(x = x, y = y, label = label),
                       size = label_pt / .pt, family = family,
                       hjust = 0.5, vjust = 0.5)
  if (nrow(box))
    p <- p + geom_label(data = box, aes(x = x, y = y, label = label,
                                        fill = I(bfill)),
                        size = label_pt / .pt, family = family,
                        hjust = 0.5, vjust = 0.5, colour = "black",
                        label.size = 0,
                        label.padding = grid::unit(box_pad_pt, "pt"),
                        label.r = grid::unit(chip_r_pt, "pt"))

  p +
    scale_fill_manual(values = palette, breaks = levels_order,
                      limits = levels_order, guide = "none") +
    coord_equal(xlim = xl, ylim = yl, expand = FALSE, clip = "off") +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
}

# The pie's legend, built by ggplot from a throwaway bar so that key
# metrics and spacing stay ggplot's, then extracted as its own grob to
# be placed independently. CATEGORY NAMES ONLY -- percentages live on
# the wedges, never here.
pie_legend_grob <- function(df, cat_col, palette, levels_order = NULL,
                            label_pt = 6.8, family = "Nimbus Sans",
                            key_mm = 3.0) {
  if (is.null(levels_order)) levels_order <- sort(unique(df[[cat_col]]))
  src <- ggplot(df, aes(x = 1, y = prop, fill = .data[[cat_col]])) +
    geom_col() +
    scale_fill_manual(values = palette, breaks = levels_order,
                      limits = levels_order, name = NULL) +
    theme_void() +
    theme(legend.text          = element_text(size = label_pt, family = family),
          legend.key.size      = unit(key_mm, "mm"),
          legend.key.spacing.y = unit(0.6, "mm"),
          legend.margin        = margin(0, 0, 0, 0),
          legend.box.spacing   = unit(0, "mm")) +
    guides(fill = guide_legend(ncol = 1))
  cowplot::get_plot_component(src, "guide-box-right", return_all = TRUE)
}
