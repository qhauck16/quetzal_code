# figure3/panel_c/build.R -- "Mutations are associated with recurrent
# mis-splicing of cancer-relevant genes".
#
# A 2x2 grid of loci (SENP2 / PRKDC on top, EZH2 / BRD9 below). Every
# quadrant is a vertical patchwork stack drawn ENTIRELY from raw
# bigwigs + MeanPSI .links with ../scripts/coverage_plot.R -- no
# pygenometracks, no ggbio, no raster crops:
#
#   [top row only] zoom callout row   zoom_callout_gg()
#   High GSP #13 coverage + arcs      sashimi_track_gg()
#   Matched Controls coverage + arcs  sashimi_track_gg()
#   lavender taper                    trapezoid_gg()
#   full MANE gene model              transcript_zoom_gg()
#   genomic coordinate labels         bespoke row (same x scale)
#

if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
panel_c <- panel_cached("c", local({

PY   <- file.path(QROOT, "coverage_tracks/pygenome/NKAP")
BLUE <- "#3498db"
GREY <- "#969696"
LAV  <- "#ccccff"

HIGH_LAB <- "High GSP #13"
CTRL_LAB <- "Matched Controls"
HIGH_N   <- "(15)"
CTRL_N   <- "(150)"

BW_HIGH <- "High_k13-3-..bw"            # k13 high-loading cohort
BW_CTRL <- "Matched_Controls-3-..bw"
LK_HIGH <- "MeanPSI.High_k13-3-..links"
LK_CTRL <- "MeanPSI.Matched_Controls-3-..links"

# ------------------------------------------------------------------
# Loci. y_max is the pipeline's own ini max_value for that region
# (shared by both tracks so the two cohorts are comparable). PRKDC's
# published "5.4" is a hand-retype error; 9.433 is the ini value.
# ------------------------------------------------------------------
LOCI <- list(
  SENP2 = list(gene = "SENP2", chrom = "chr3",
               wide = c(185609022, 185611878),
               zoom = c(185611603, 185611781),
               y_max = 11.6660145759583, label_side = "left"),
  PRKDC = list(gene = "PRKDC", chrom = "chr8",
               wide = c(47858866, 47859738),
               zoom = c(47858859, 47858999),
               y_max = 9.43334007263184, label_side = "right"),
  EZH2  = list(gene = "EZH2",  chrom = "chr7",
               wide = c(148846440, 148847314),
               zoom = NULL,
               y_max = 5.70414457321167, label_side = "right"),
  BRD9  = list(gene = "BRD9",  chrom = "chr5",
               wide = c(889033, 889686),
               zoom = NULL,
               y_max = 10.5822075843811, label_side = "right")
)
region_dir <- function(L, which = "wide") {
  r <- if (which == "wide") L$wide else L$zoom
  file.path(PY, L$gene, sprintf("%s_%d-%d", L$chrom, r[1], r[2]))
}

# ------------------------------------------------------------------
# Row heights. patchwork heights are RELATIVE, so only the ratios
# here are stable. Everything that needs
# a DRAWN size (the PSI-label clearance solver, the callout's leader
# drop) is derived from these, never restated.
# ------------------------------------------------------------------
# Stretch the plot rather than shrink the type. The
# two PSI numbers per track are what sets the floor on H_TRK -- they
# are anchored APEX_GAP apart, so the track has to be tall enough for
# APEX_GAP * H_TRK to exceed one line of SZ_PSI text. Everything else
# in the stack was squeezed to pay for it.
H_CALL <- 0.38     # zoom-callout row (top quadrants only)
H_TRK  <- 0.36     # one coverage/sashimi track
H_TRAP <- 0.09     # lavender taper
H_GENE <- 0.15     # MANE gene model
H_AXIS <- 0.12     # genomic coordinate labels
TOP_H  <- H_CALL + 2 * H_TRK + H_TRAP + H_GENE + H_AXIS
BOT_H  <-          2 * H_TRK + H_TRAP + H_GENE + H_AXIS

# ROW_H[["c"]] is panel c's whole row; the title band eats ~0.34 in of it.
# Scale the nominal stack above onto what is actually left so
# sashimi_track_gg's point-size -> data-unit conversion is right.
C_TITLE_H  <- 0.34
# Gap under each coverage track, in pt. The tracks would otherwise be
# stacked flush, and this track's y-axis "0" (now centred ON the
# baseline, see sashimi_track_gg's y_axis) and the next track's y_max
# number would then be drawn on the same page row and overprint. One
# line of SZ_AXIS text is ~4.5 pt tall, so 5 pt buys both their half
# line.
# 9 pt between rows; the tracks flatten to pay for it because the gaps come off the budget first (below).
TRK_GAP_PT <- 9
# Whitespace between the top (SENP2 / PRKDC) and bottom (EZH2 / BRD9)
# locus pairs, so the upper axis labels don't sit on the next track.
QUAD_GAP_PT <- 8
# The four gaps (two per stacked quadrant) are plot MARGINS, and
# patchwork's relative heights govern the PANELS -- so the gaps come
# off the budget BEFORE the rows are scaled. Get this wrong and every
# point-size -> data-unit conversion below is off by the same factor.
C_SCALE    <- (ROW_H[["c"]] - C_TITLE_H - (4 * TRK_GAP_PT + QUAD_GAP_PT) / 72) /
              (TOP_H + BOT_H)
# DRAWN (panel) height of one track.
TRACK_H_IN <- H_TRK * C_SCALE

# Gutters, in pt. Only the RIGHT column carries the cohort captions
# so the two columns need different margins to
# end up with comparable track widths.
MARG <- list(left  = c(30, 40),      # c(left, right) in pt
             right = c(22, 62))
gutter <- function(side) margin(0, MARG[[side]][2], 0, MARG[[side]][1], "pt")
# Same gutter plus the inter-track gap; only the coverage tracks use it.
trk_gutter <- function(side)
  margin(0, MARG[[side]][2], TRK_GAP_PT, MARG[[side]][1], "pt")
TRACK_W_IN <- list(left  = (PAGE_W / 2 - sum(MARG$left)  / 72),
                   right = (PAGE_W / 2 - sum(MARG$right) / 72))

SZ_AXIS   <- SIZE_BASE * 0.80 / .pt   # "11.7" / "0" and coordinate labels
SZ_COHORT <- SIZE_BASE * 0.90 / .pt
SZ_PSI    <- SIZE_BASE * 0.78 / .pt
SZ_NT     <- SIZE_BASE * 0.85 / .pt
SZ_GENE   <- SIZE_BASE * 0.95 / .pt

# ------------------------------------------------------------------
# .links helpers. sashimi_track_gg filters s1 >= start & e2 <= end;
# reproduce that here so the arcs we reason about are exactly the
# arcs that get drawn.
# ------------------------------------------------------------------
read_links <- function(path, start, end, min_psi = 1) {
  lk <- read.table(path, header = FALSE, sep = "\t",
                   col.names = c("c1","s1","e1","c2","s2","e2","psi"))
  lk[lk$s1 >= start & lk$e2 <= end & lk$psi >= min_psi, , drop = FALSE]
}

# Derive the canonical/cryptic splice-site pair and the nt offset
# from the coordinates themselves -- never hardcoded. Among the drawn
# HIGH junctions, find the pair that shares one endpoint; the shared
# end is the anchor, the other end is the splice site that moves. The
# CANONICAL site is whichever of the two the matched controls use
# (highest control PSI); offset = cryptic - canonical.
derive_offset <- function(hi, ct) {
  cand <- NULL
  for (i in seq_len(nrow(hi))) for (j in seq_len(nrow(hi))) {
    if (i >= j) next
    if (hi$s1[i] == hi$s1[j]) cand <- c(hi$e2[i], hi$e2[j])
    if (hi$e2[i] == hi$e2[j]) cand <- c(hi$s1[i], hi$s1[j])
  }
  stopifnot(!is.null(cand))
  psi_of <- function(x) {
    m <- ct$psi[ct$s1 == x | ct$e2 == x]
    if (length(m) == 0) 0 else max(m)
  }
  p <- vapply(cand, psi_of, numeric(1))
  canonical <- cand[which.max(p)]
  cryptic   <- cand[which.min(p)]
  list(canonical = canonical, cryptic = cryptic,
       offset = cryptic - canonical,
       # every cryptic site here is downstream of the canonical 3'ss in the
       # transcript, so label the size of the displacement, not its genomic sign
       label  = sprintf("+%d nt", abs(cryptic - canonical)))
}

# ------------------------------------------------------------------
# One coverage/sashimi track.
# ------------------------------------------------------------------
# ---- Arc geometry, computed HERE rather than left to the helper ----
#
# Two things forced this. (1) Arc ends: with sashimi_track_gg's 200 bp
# default lookback, SENP2's arc end was pinned to the top of the exon
# next door (0.70 of y_max) instead of to the coverage at the splice
# site itself, which left no room for the arc. ANCHOR_BP = 0 pins each
# end to the coverage AT the junction coordinate -- the bins are 1 bp
# wide in every one of these windows, so that is the exact value.
# (2) PSI labels: the helper's solver only moves a label vertically at
# its own arc's midpoint, and these junctions share an endpoint, so
# both labels want the same x. In a track five text-heights tall there
# is no vertical slot that clears both arcs and the other label, and
# the two numbers overprinted. `show_psi = FALSE`; place_psi() below
# uses the horizontal axis, which is the one with room to spare.
#
# The lever is that sashimi_track_gg's curve is
#     arc_y(t) = (1-t) y1 + t y2 + sin(t pi) (apex - max(y1, y2))
# so `apex - max(y1, y2)` is the sine bump's amplitude. Back-solving it
# lets each arc's PEAK be chosen, which is what guarantees the arcs are
# visibly separated and never sag.  <-- keep in step with the helper.
ANCHOR_BP <- 0

# ---- Arc HEIGHTS -------------------------------------------------
# Every number below is a fraction of y_max.
#
# The arcs are drawn from the coverage top AT each splice site, and in
# these windows those tops are high: the cryptic site sits inside the
# exon plateau, so its end of the arc starts as high as 0.81 of y_max
# (PRKDC controls). An arc height expressed as an ABSOLUTE fraction of
# y_max would therefore be below its own
# endpoints, and sin(t pi) * (apex - max(y1,y2)) goes NEGATIVE: the
# "arc" would sag down to the baseline and back up. Heights are
# relative to each arc's own higher endpoint.
#
# ARC_RISE is the bulge of the LOWEST arc above its own higher
# endpoint; each arc above it is lifted by at least ARC_SEP more so the
# two stay apart. ARC_CEIL keeps the drawn peak (which is NOT at the
# midpoint when the chord is steep -- it is solved for numerically)
# inside the track. When the ceiling bites, the arcs BELOW are pulled
# down with it so the separation survives, but never below
# ARC_RISE_MIN, which is the floor on "still reads as an arc".
ARC_RISE     <- 0.20
ARC_SEP      <- 0.16
ARC_CEIL     <- 0.93
ARC_RISE_MIN <- 0.09

# Peak of (1-t) y1 + t y2 + sin(t pi) * a, over t in [0, 1].
arc_peak <- function(y1, y2, a, tt = seq(0, 1, length.out = 401))
  max((1 - tt) * y1 + tt * y2 + sin(tt * pi) * a)

arc_geom <- function(bw, lk_path, L) {
  d  <- bigwig_binned(bw, L$chrom, L$wide[1], L$wide[2])
  lk <- read_links(lk_path, L$wide[1], L$wide[2])
  # sashimi_track_gg's own ordering: SHORT -> LONG. It INTERPOLATES
  # apex_range linearly over this order, so handing back a 2-element
  # apex_range only reproduces the apexes we solved for while n <= 2.
  lk <- lk[order(lk$e2 - lk$s1, decreasing = FALSE), , drop = FALSE]
  n  <- nrow(lk)
  if (n > 2)
    warning("panel c: ", L$gene, " has ", n, " arcs; sashimi_track_gg ",
            "interpolates apex_range linearly, so the solved peaks ",
            "will not be honoured. Pass per-arc apexes instead.")
  covat <- function(x) {
    w <- d$bin_center >= x - ANCHOR_BP & d$bin_center <= x + ANCHOR_BP
    if (!any(w)) d$coverage[which.min(abs(d$bin_center - x))] else max(d$coverage[w])
  }
  y1 <- vapply(lk$s1, covat, numeric(1)) / L$y_max
  y2 <- vapply(lk$e2, covat, numeric(1)) / L$y_max
  E  <- pmax(y1, y2)

  # Stack the arcs by their own endpoint height, NOT by junction
  # length. "Shortest junction highest" is the helper's convention and
  # is right when arcs spring from the baseline, but here the LONGER
  # junction of each pair is the one anchored deeper inside the exon,
  # so it starts higher at one end and its chord dominates the other's
  # everywhere. Forcing the short arc on top would mean lifting it
  # clear over the long arc's endpoint -- out through the ceiling -- or
  # letting the two cross, which is what the old constants did.
  ord  <- order(E, (y1 + y2) / 2)
  peak <- numeric(n)
  for (k in seq_len(n)) {
    i <- ord[k]
    peak[i] <- E[i] + ARC_RISE
    if (k > 1) peak[i] <- max(peak[i], peak[ord[k - 1]] + ARC_SEP)
  }
  for (k in rev(seq_len(n))) {                     # apply the ceiling downwards
    i   <- ord[k]
    cap <- if (k == n) ARC_CEIL else peak[ord[k + 1]] - ARC_SEP
    peak[i] <- min(peak[i], max(cap, E[i] + ARC_RISE_MIN))
  }

  # Back-solve the sine amplitude that puts the DRAWN peak there.
  amp <- vapply(seq_len(n), function(i) {
    f <- function(a) arc_peak(y1[i], y2[i], a) - peak[i]
    if (f(0) >= 0) 0 else uniroot(f, c(0, 4), tol = 1e-7)$root
  }, numeric(1))
  apex <- amp + E                # helper: bump = apex - max(y1, y2)

  # The exact curves the helper will draw, in DATA units, so the PSI
  # labels can be solved against them rather than against an estimate.
  curves <- lapply(seq_len(n), function(i) {
    xs <- seq(lk$s1[i], lk$e2[i], length.out = 200)
    tt <- (xs - lk$s1[i]) / (lk$e2[i] - lk$s1[i])
    list(x = xs,
         y = ((1 - tt) * y1[i] + tt * y2[i] + sin(tt * pi) * amp[i]) * L$y_max)
  })
  list(lk = lk, n = n, d = d, curves = curves,
       apex_range = c(apex[1], apex[n]), peak = peak, E = E)
}

# SIGNED vertical gap between the arcs, taken in peak order: negative
# anywhere means the higher-peaked arc dips under the lower one, i.e.
# a crossing, which is what the old absolute tiers produced. Reported
# by make_track so a regression shows up in the build log.
arc_min_gap <- function(curves, peak) {
  n <- length(curves)
  if (n < 2) return(NA_real_)
  o <- order(peak)
  g <- Inf
  for (k in 2:n) {
    up <- curves[[o[k]]]; dn <- curves[[o[k - 1]]]
    lo <- max(min(up$x), min(dn$x)); hi <- min(max(up$x), max(dn$x))
    if (hi <= lo) next
    # Ignore the very ends: these junction pairs SHARE an endpoint, so
    # the gap there is 0 by construction and would mask a crossing.
    m  <- 0.08 * (hi - lo)
    xs <- seq(lo + m, hi - m, length.out = 300)
    g  <- min(g, min(approx(up$x, up$y, xs)$y - approx(dn$x, dn$y, xs)$y))
  }
  g
}

# Bare 3-significant-digit values.
# format() on the whole VECTOR pads to a common width ("92.40"), so
# each value is formatted on its own.
psi_fmt <- function(p)
  vapply(p, function(x) format(signif(x, 3), trim = TRUE), character(1))

# ------------------------------------------------------------------
# PSI labels: each one rides ON ITS OWN ARC.
# ------------------------------------------------------------------
# The old scheme stacked both numbers in a column above every arc. That
# only works if there is a band of clear air above the highest arc, and
# there is not: one line of SZ_PSI text is ~0.20 of a track's height
# here, the arcs already start as high as 0.81, and the two junctions
# of a pair share an endpoint so both numbers want the same x.
#
# Instead each label is offset one clearance from its own curve --
# above it where that fits, otherwise below it into the empty bowl --
# and slid ALONG the arc until it clears everything else. Horizontal
# room is the cheap axis here: an arc spans hundreds of bp and a label
# is ~50 of them, so two labels at different t never interact at all,
# which is what makes the crowded tracks solvable. Labels are placed
# highest-arc-first, because that is the one with the least room.
#
# Preference order: centre of the arc first, then outward in t,
# trying above the curve before below at each stop.
#
# NB the candidates are ONE flat list, deliberately. As a nested
# `for (t) for (side)` an R `break` leaves only the inner loop, so the
# search silently returned the LAST acceptable slot instead of the
# first and every number drifted out to t = 0.78.
PSI_CAND <- local({
  ts <- c(0.50, 0.44, 0.56, 0.38, 0.62, 0.32, 0.68, 0.27, 0.73, 0.22, 0.78)
  data.frame(t = rep(ts, each = 2), sgn = rep(c(1, -1), times = length(ts)))
})

place_psi <- function(g, L, side, obstacles = list()) {
  y_max <- L$y_max
  s <- L$wide[1]; e <- L$wide[2]
  # text metrics -> data units, via the track's DRAWN size
  hh   <- 0.5 * (SZ_PSI * .pt) / 72 * y_max / TRACK_H_IN   # half line height
  pad  <- 1.8 / 72 * y_max / TRACK_H_IN                    # breathing room
  off  <- hh + pad                                         # curve -> label centre
  bh   <- hh + 0.35 * pad                                  # collision half-height
  xpp  <- (e - s) / (TRACK_W_IN[[side]] * 72)              # data-x per pt
  labs <- psi_fmt(g$lk$psi)
  bw_  <- (nchar(labs) * 0.56 * (SZ_PSI * .pt) / 2 + 1.5) * xpp

  # Highest label CENTRE that keeps the glyphs inside the track. Using
  # the collision half-height here instead would double-count the pad
  # and cost ~0.06 of y_max -- enough that the top label in four of the
  # eight tracks found no slot at all and fell through to the fallback.
  ceil_y <- y_max - (hh + 0.5 * pad)
  cov_at <- function(x0, x1) {
    w <- g$d$bin_center >= x0 & g$d$bin_center <= x1
    if (any(w)) max(g$d$coverage[w]) else 0
  }
  curve_at <- function(cv, xs) {
    keep <- xs >= min(cv$x) & xs <= max(cv$x)
    if (!any(keep)) return(numeric(0))
    approx(cv$x, cv$y, xs[keep])$y
  }
  hits <- function(box) {                    # box = c(x0, x1, y0, y1)
    xs <- seq(box[1], box[2], length.out = 25)
    for (cv in g$curves) {
      yy <- curve_at(cv, xs)
      if (length(yy) && any(yy > box[3] & yy < box[4])) return(TRUE)
    }
    if (cov_at(box[1], box[2]) > box[3]) return(TRUE)
    for (o in obstacles)
      if (box[1] < o[2] && box[2] > o[1] && box[3] < o[4] && box[4] > o[3])
        return(TRUE)
    FALSE
  }

  placed <- list()
  out    <- data.frame(x = numeric(0), y = numeric(0), lab = character(0))
  for (i in order(g$peak, decreasing = TRUE)) {
    cv   <- g$curves[[i]]
    span <- range(cv$x)
    got  <- NULL
    # Pass 1 also demands the label be clear of every other label in X.
    # Two numbers at the same x, one over the top arc and one under the
    # bottom one, are legible but they re-read as the old stacked
    # column; separating them horizontally is what makes each one
    # obviously the property of the arc it touches. Pass 2 drops the
    # requirement so a crowded track still gets a placement.
    for (strict in c(TRUE, FALSE)) {
      for (ci in seq_len(nrow(PSI_CAND))) {
        t <- PSI_CAND$t[ci]; sgn <- PSI_CAND$sgn[ci]
        xc <- span[1] + t * diff(span)
        # Clear the arc over the label's WHOLE width, not just under
        # its centre: on a steep stretch the limb closes in sideways
        # and a centre-only offset lets the curve clip the text.
        yb <- curve_at(cv, seq(xc - bw_[i], xc + bw_[i], length.out = 25))
        yl <- if (sgn > 0) max(yb) + off else min(yb) - off
        box <- c(xc - bw_[i], xc + bw_[i], yl - bh, yl + bh)
        if (yl > ceil_y || yl - hh < 0) next
        if (box[1] < s || box[2] > e) next
        if (hits(box)) next
        if (length(placed) &&
            any(vapply(placed, function(p)
                  box[1] < p[2] + strict * bw_[i] &&
                  box[2] > p[1] - strict * bw_[i] &&
                  (strict || (box[3] < p[4] && box[4] > p[3])),
                  logical(1)))) next
        got <- list(x = xc, y = yl, box = box)
        break
      }
      if (!is.null(got)) break
    }
    if (is.null(got)) {                      # never leave a number unplaced
      xc <- span[1] + 0.5 * diff(span)
      yc <- approx(cv$x, cv$y, xc)$y
      yl <- min(yc + off, ceil_y)
      got <- list(x = xc, y = yl,
                  box = c(xc - bw_[i], xc + bw_[i], yl - bh, yl + bh))
      cat(sprintf("[panel c]   !! no clear slot for PSI %s (%s)\n",
                  labs[i], L$gene))
    }
    placed[[length(placed) + 1]] <- got$box
    out <- rbind(out, data.frame(x = got$x, y = got$y, lab = labs[i]))
  }
  out
}

# ------------------------------------------------------------------
# The inline "-11 nt" / "-5 nt" label, for the loci with no callout.
# ------------------------------------------------------------------
# It hangs just LEFT of the leftmost splice site
# -- over the exon's descending shoulder -- at a height derived from
# the coverage under its own width, so it clears the exon below it.
# The one place in the track guaranteed to be arc-free: every arc here
# starts at a splice site and runs rightwards.
#
# It is drawn on the HIGH track only, once per locus: there is no room
# above the grey exon on the control tracks, because y_max is the pipeline's ini value and is driven by the CONTROLS, so
# a control exon plateau fills 0.91 of y_max while half a line of
# SZ_NT is 0.11 of it. Placed there the text printed ON the coverage.
# One label per locus, on the cohort the cryptic junction belongs to,
# is what both fits and reads.
nt_place <- function(g, L, side, nt) {
  y_max <- L$y_max; s <- L$wide[1]; e <- L$wide[2]
  xpp  <- (e - s) / (TRACK_W_IN[[side]] * 72)
  hh   <- 0.5 * (SZ_NT * .pt) / 72 * y_max / TRACK_H_IN
  # 1 pt of air, not more: EZH2's High exon tops out at 0.72 of y_max
  # and the label itself is 0.22 of it, so the whole budget above the
  # coverage is 0.28 and a 1.8 pt gap overran the track by a hair.
  clr  <- hh + 1.0 / 72 * y_max / TRACK_H_IN
  wid  <- nchar(nt$label) * 0.56 * (SZ_NT * .pt) * xpp
  site <- min(nt$canonical, nt$cryptic)
  x1   <- site - 0.006 * (e - s)
  x0   <- x1 - wid
  k    <- g$d$bin_center >= x0 - 2 & g$d$bin_center <= x1 + 2
  cv   <- if (any(k)) max(g$d$coverage[k]) else 0
  ok   <- cv + clr + hh <= y_max
  list(x0 = x0, x1 = x1, xpp = xpp, clr = clr, hh = hh, cov = cv,
       y = if (ok) cv + clr else y_max - hh, ok = ok)
}

# ------------------------------------------------------------------
# Zoom callouts: the EDGE OF THE JUNCTION ARC, at nucleotide scale.
# ------------------------------------------------------------------
# `zoom_callout_gg`'s built-in
# annotation is a decorative geom_curve leader from the text down to
# the coverage; it is switched off here (`nt_label = NULL`) and
# replaced by the real thing:
#
#   * the FOOT of each junction's arc, anchored at that junction's true
#     splice-site coordinate on the zoomed coverage and rising away in
#     the direction the junction actually runs, drawn in the same black
#     stroke as the main track's arcs and with the same PSI-derived
#     linewidth (0.35 + LW_SCALE * sqrt(psi/100)). Both junctions of
#     the pair get a foot, so the box IS the main track's arc pair
#     magnified -- and the 9 nt (SENP2) / 6 nt (PRKDC) offset between
#     the two feet is the thing the label names.
#   * a PLAIN, straight, horizontal "+9 nt" / "-6 nt" set in the
#     headroom immediately beside the feet. No leader is drawn: the
#     text abuts the feet, which are themselves anchored on the sites.
#
# Shape. sashimi_track_gg's arc is (1-t) y1 + t y2 + sin(t pi) a, so
# the quarter running from a foot up to the apex is exactly
# sin(s pi/2) in the foot's own coordinates. The visible stretch stops
# at s = ZOOM_ARC_U (< 1) so the stroke is still CLIMBING where it
# meets the frame -- it reads as an arc that continues outside the box,
# which is true: the other end of each junction is hundreds of bp away.
#
# Both feet end on the same box edge, and their END heights are fanned
# out by ZOOM_ARC_SEP in anchor order. Without the fan they arrive
# separated only by the coverage step they sprang from, and in the
# SENP2 controls that step is 0.08 of the box -- the two strokes read
# as one doubled line. Ordering the ends the same way as the anchors
# also keeps the pair from crossing; `zoom_feet` measures the signed
# minimum gap between the drawn curves every build and stops if it
# ever goes negative, so this is checked rather than assumed.
#
# COUPLING: the box interior geometry is computed inside
# zoom_callout_gg and is only partly parameterised. ZOOM_BOX_Y,
# ZOOM_COVF and ZOOM_NBINS are passed through so they cannot drift;
# ZOOM_INNER0 mirrors the helper's hardcoded 0.035 baseline inset and
# MUST be kept in step with it by hand.
ZOOM_BOX_Y  <- c(0.28, 0.99)   # box top/bottom, callout-row units
ZOOM_COVF   <- 0.66            # fraction of the box the coverage may use
ZOOM_NBINS  <- 600
ZOOM_INNER0 <- 0.035           # <- coverage_plot.R zoom_callout_gg
ZOOM_ARC_U  <- 0.80            # stop short of the apex (see above)
ZOOM_ARC_TOP<- 0.80            # box height the highest foot climbs to
ZOOM_ARC_SEP<- 0.13            # box heights the feet's ENDS are fanned by
ZOOM_PROBE  <- 3               # bp into the exon the anchor height is read
ZOOM_LAB_GAP<- 0.060           # box widths of air between feet and text
LW_SCALE    <- 0.45            # identical to the main tracks' lw_scale

# Drawn height of the callout row, in inches -- the point-size -> y-unit
# conversion for the nt label. H_CALL is in stack units and the callout
# row's y axis runs 0..1 over its own height.
CALL_H_IN <- H_CALL * C_SCALE

zoom_feet <- function(L, nt, spec, side) {
  zs <- L$zoom[1]; ze <- L$zoom[2]
  ww    <- L$wide[2] - L$wide[1]
  bx_lo <- L$wide[1] + spec$box[1] * ww
  bx_hi <- L$wide[1] + spec$box[2] * ww
  bxw   <- bx_hi - bx_lo
  by_lo <- ZOOM_BOX_Y[1]; by_hi <- ZOOM_BOX_Y[2]; byh <- by_hi - by_lo
  inner_lo <- by_lo + ZOOM_INNER0 * byh
  inner_hi <- by_lo + ZOOM_COVF   * byh

  d    <- bigwig_binned(spec$bw, L$chrom, zs, ze, n_bins = ZOOM_NBINS)
  cmax <- max(d$coverage, na.rm = TRUE)
  if (!is.finite(cmax) || cmax <= 0) cmax <- 1
  px  <- function(bp) bx_lo + (bp - zs) / (ze - zs) * bxw
  bpx <- function(x)  zs + (x - bx_lo) / bxw * (ze - zs)   # inverse
  py  <- function(cv) inner_lo + cv / cmax * (inner_hi - inner_lo)
  cov_bp <- function(lo, hi) {
    w <- d$bin_center >= lo & d$bin_center <= hi
    if (any(w)) max(d$coverage[w]) else 0
  }

  lk <- spec$links
  # Which end of the junction moves? derive_offset() already found the
  # pair; whichever column carries both candidate coordinates is the
  # splice site being magnified. An e2 site has its exon to the RIGHT,
  # so its arc runs LEFT, and vice versa for an s1 site.
  is_e2    <- all(c(nt$canonical, nt$cryptic) %in% lk$e2)
  exon_dir <- if (is_e2) 1 else -1
  arc_dir  <- -exon_dir
  site     <- if (is_e2) lk$e2 else lk$s1
  keep     <- site >= zs & site <= ze
  lk <- lk[keep, , drop = FALSE]; site <- site[keep]
  if (!nrow(lk)) return(NULL)

  anchor <- vapply(site, function(s0)
    py(if (exon_dir > 0) cov_bp(s0, s0 + ZOOM_PROBE)
       else              cov_bp(s0 - ZOOM_PROBE, s0)), numeric(1))
  x_site <- px(site)
  x_end  <- if (arc_dir < 0) bx_lo else bx_hi
  kk     <- sin(ZOOM_ARC_U * pi / 2)
  ord    <- order(anchor, decreasing = TRUE)
  end_y  <- numeric(length(anchor))
  end_y[ord] <- (by_lo + ZOOM_ARC_TOP * byh) -
                (seq_along(ord) - 1) * ZOOM_ARC_SEP * byh
  A      <- end_y - anchor

  feet <- lapply(seq_along(site), function(i) {
    W  <- abs(x_end - x_site[i])
    ss <- seq(0, 1, length.out = 120)
    data.frame(x  = x_site[i] + arc_dir * ss * W,
               y  = anchor[i] + A[i] * sin(ss * ZOOM_ARC_U * pi / 2) / kk,
               lw = 0.35 + LW_SCALE * sqrt(lk$psi[i] / 100),
               grp = i)
  })

  # Signed minimum gap between the pair over their common x range,
  # measured on the drawn curves. Negative = the feet cross.
  fgap <- NA_real_
  if (length(feet) == 2) {
    lo <- max(min(feet[[1]]$x), min(feet[[2]]$x))
    hi <- min(max(feet[[1]]$x), max(feet[[2]]$x))
    if (hi > lo) {
      xs <- seq(lo, hi, length.out = 200)
      y1 <- approx(feet[[1]]$x, feet[[1]]$y, xs)$y
      y2 <- approx(feet[[2]]$x, feet[[2]]$y, xs)$y
      fgap <- min(abs(y1 - y2)) / byh
      if (min((y1 - y2) * sign(anchor[1] - anchor[2])) < 0)
        stop("panel c: zoom arc feet cross in ", L$gene, " ", spec$cohort)
    }
  }

  # ---- the nt label: plain, horizontal, beside the feet ----------
  # It goes on the side of the sites the arcs do NOT sweep across, one
  # gap clear of the outermost foot, at a height derived from the
  # coverage under its own width. Widths/heights come from the point
  # size through the DRAWN size of this row, never hardcoded.
  nt_pt <- SZ_NT * .pt
  xpp   <- ww / (TRACK_W_IN[[side]] * 72)              # wide-x per pt
  wid   <- nchar(nt$label) * 0.56 * nt_pt * xpp
  hh    <- 0.36 * nt_pt / 72 / CALL_H_IN               # cap half-height
  pad   <- 1.0 / 72 / CALL_H_IN
  gap   <- ZOOM_LAB_GAP * bxw
  if (arc_dir < 0) { lx <- max(x_site) + gap; hj <- 0; sp <- c(lx, lx + wid) }
  else             { lx <- min(x_site) - gap; hj <- 1; sp <- c(lx - wid, lx) }
  sp   <- pmin(pmax(sp, bx_lo), bx_hi)
  cvy  <- py(cov_bp(bpx(sp[1]), bpx(sp[2])))
  ceil <- by_hi - hh - pad
  want <- cvy + hh + pad
  # `want` is this box's own floor for the label. The caller raises
  # every box of the locus to the tallest of them so the two "+9 nt"
  # do not sit at visibly different heights side by side.
  list(feet = feet, lab_x = lx, lab_hj = hj, lab_want = want,
       lab_ceil = ceil, lab_half = hh, label = nt$label,
       bx = c(bx_lo, bx_hi), by = c(by_lo, by_hi),
       diag = sprintf("%s feet x=%s anchor=%s end=%s gap=%.3f | '%s' x=%.2f",
                      spec$cohort,
                      paste(sprintf("%.3f", (x_site - bx_lo) / bxw), collapse = "/"),
                      paste(sprintf("%.3f", (anchor - by_lo) / byh), collapse = "/"),
                      paste(sprintf("%.3f", (end_y  - by_lo) / byh), collapse = "/"),
                      fgap, nt$label, (lx - bx_lo) / bxw))
}

# Turn the solved feet + labels into layers on the callout plot. Every
# value the layers need travels as a data-frame COLUMN or as a plain
# numeric layer parameter -- an aes() closing over a local scalar
# renders blank after panel_cached()'s saveRDS/readRDS round trip.
zoom_arc_layers <- function(L, nt, specs, side) {
  zz <- Filter(Negate(is.null), lapply(specs, function(s)
    zoom_feet(L, nt, s, side)))
  if (!length(zz)) return(list())
  ly <- min(max(vapply(zz, function(z) z$lab_want, numeric(1))),
            min(vapply(zz, function(z) z$lab_ceil, numeric(1))))
  out <- list()
  for (z in zz) {
    cat(sprintf("[panel c]     zoom %s y=%.3f%s\n", z$diag,
                (ly - z$by[1]) / diff(z$by),
                if (ly + 1e-9 < z$lab_want) "  !! label capped" else ""))
    for (f in z$feet)
      out[[length(out) + 1]] <- geom_path(
        data = f, aes(x = x, y = y, group = grp), inherit.aes = FALSE,
        color = "black", linewidth = f$lw[1], lineend = "butt")
    out[[length(out) + 1]] <- geom_text(
      data = data.frame(x = z$lab_x, y = ly, label = z$label, hj = z$lab_hj),
      aes(x = x, y = y, label = label, hjust = hj),
      vjust = 0.5, size = SZ_NT, family = FONT_FAMILY,
      color = "grey10", inherit.aes = FALSE)
  }
  out
}

# Controls-callout leaders, as obstacles for the HIGH track's PSI
# numbers. Those two dashed lines run from the grey inset's bottom
# corners down to the grey source box, so they cross the whole High
# track diagonally -- and a PSI number solved without knowing about
# them gets struck through (PRKDC's "22.1" was). The obstacle test is
# rectangle-vs-rectangle, so each leader is approximated as a stack of
# thin boxes, one per y band, holding the leader's x-extent in that
# band plus a 2 pt margin.
#
# Geometry mirrors zoom_callout_gg's `leads`: corner (box edge, by_lo)
# -> (source box edge, -src_drop), in callout-PANEL units, with the
# High track's panel spanning callout y 0 (its top, = y_max) down to
# -H_TRK / H_CALL (its baseline, = 0).
leader_obst <- function(L, box, side, src_drop, n_band = 10) {
  if (is.null(box) || is.null(L$zoom)) return(list())
  s <- L$wide[1]; e <- L$wide[2]; ww <- e - s
  bx  <- s + box * ww
  src <- c(max(L$zoom[1], s), min(L$zoom[2], e))
  by_lo <- ZOOM_BOX_Y[1]
  pad <- 7 * ww / (TRACK_W_IN[[side]] * 72)       # 7 pt, in data-x: at 2 pt an inline "43.6" still sat on the dashed line
  ys  <- seq(0, L$y_max, length.out = n_band + 1)
  cy  <- function(y) -(1 - y / L$y_max) * (H_TRK / H_CALL)
  out <- list()
  for (k in seq_len(n_band)) {
    f <- (by_lo - cy(ys[k:(k + 1)])) / (by_lo + src_drop)
    for (j in 1:2) {
      xs <- bx[j] + (src[j] - bx[j]) * f
      out[[length(out) + 1]] <- c(min(xs) - pad, max(xs) + pad,
                                  ys[k], ys[k + 1])
    }
  }
  out
}

make_track <- function(L, cohort, col, side,
                       caption = NULL, caption_n = NULL,
                       src_box = FALSE, nt = NULL,
                       extra_obst = list()) {
  wd    <- region_dir(L, "wide")
  bw    <- file.path(wd, if (cohort == "high") BW_HIGH else BW_CTRL)
  lk    <- file.path(wd, if (cohort == "high") LK_HIGH else LK_CTRL)
  y_max <- L$y_max
  s <- L$wide[1]; e <- L$wide[2]

  g <- arc_geom(bw, lk, L)
  cat(sprintf("[panel c]   %-5s %-4s psi=%-14s end=%s peak=%s gap=%.2f\n",
              L$gene, cohort, paste(g$lk$psi, collapse = ","),
              paste(sprintf("%.2f", g$E),    collapse = "/"),
              paste(sprintf("%.2f", g$peak), collapse = "/"),
              arc_min_gap(g$curves, g$peak) / y_max))

  p <- sashimi_track_gg(
        bw_path = bw, chrom = L$chrom, start = s, end = e,
        track_color = col, track_label = "",     # caption drawn below
        links_path = lk, min_psi = 1, y_max = y_max,
        psi_size = SZ_PSI, label_family = FONT_FAMILY,
        apex_range = g$apex_range,
        track_height_in = TRACK_H_IN,
        track_width_in  = TRACK_W_IN[[side]],
        lw_scale = 0.45, anchor_bp = ANCHOR_BP,
        show_psi = FALSE, y_axis_size = SZ_AXIS)
  # NB the arcs + numbers are re-drawn below by a SECOND call once the
  # obstacles (nt note, zoom source box, leader corridors) are known;
  # this first call only fixes geometry the obstacles are derived from.

  # ---- inline "-11 nt" / "-5 nt", for the loci with no callout ----
  # Solved BEFORE the PSI numbers so it can be handed to them as an
  # obstacle. Placed to the LEFT of the leftmost splice site, which is
  # the one place in the track guaranteed to be arc-free (every arc
  # here starts at a splice site and runs rightwards), at a height
  # derived from the coverage immediately left of it so it clears the
  # exon below.
  nt_lab <- NULL; obst <- list()
  if (!is.null(nt)) {
    f <- nt_place(g, L, side, nt)
    cat(sprintf("[panel c]     nt '%s' cov=%.2f y=%.2f top=%.2f%s\n",
                nt$label, f$cov / y_max, f$y / y_max,
                (f$y + f$hh) / y_max,
                if (f$ok) "" else "  !! does not fit in the track"))
    nt_lab <- data.frame(x = f$x1, y = f$y, lab = nt$label)
    obst[[length(obst) + 1]] <- c(f$x0 - 2 * f$xpp, f$x1 + 2 * f$xpp,
                                  f$y - f$clr, f$y + f$clr)
  }
  # The zoom source box is a full-height rectangle; keep the numbers
  # out of its whole x span rather than just off its outline.
  if (src_box)
    obst[[length(obst) + 1]] <- c(max(L$zoom[1], s), min(L$zoom[2], e),
                                  -Inf, Inf)
  obst <- c(obst, extra_obst)

  # ---- PSI numbers: INLINE, on their own arcs ----------------------
  # One sashimi label style for the whole paper: the number
  # sits in the middle of its own arc, in a gap cut out of the arc (the
  # figure 4 look). The shared helper does exactly that, so the track is
  # rebuilt with show_psi = TRUE and every obstacle this panel knows
  # about. place_psi() (above-/below-the-arc numbers) is no longer used.
  p <- sashimi_track_gg(
        bw_path = bw, chrom = L$chrom, start = s, end = e,
        track_color = col, track_label = "",
        links_path = lk, min_psi = 1, y_max = y_max,
        psi_size = SZ_PSI, label_family = FONT_FAMILY,
        apex_range = g$apex_range,
        track_height_in = TRACK_H_IN,
        track_width_in  = TRACK_W_IN[[side]],
        lw_scale = 0.45, anchor_bp = ANCHOR_BP,
        show_psi = TRUE, psi_style = "inline", psi_format = psi_fmt,
        # SENP2 High: its two arcs share both feet and ran ~4 pt apart
        # end to end, so 43.6 had no slot clear of 47.3's arc and the
        # dashed zoom leaders. Its peak is dropped to open a lane.
        apex_override = if (L$gene == "SENP2" && cohort == "high")
                          c("43.6" = 0.50) else NULL,
        psi_obstacles = obst, y_axis_size = SZ_AXIS)

  if (!is.null(nt_lab))
    p <- p + geom_text(data = nt_lab, aes(x = x, y = y, label = lab),
                       hjust = 1, vjust = 0.5, size = SZ_NT,
                       family = FONT_FAMILY, color = "grey10",
                       inherit.aes = FALSE)

  # The y axis (rule + ticks at 0 / y_max + the two numbers) is drawn
  # by sashimi_track_gg itself, in the left gutter -- see y_axis there.
  if (!is.null(caption)) {
    # hjust = 0 plus a FIXED 3 pt gap, not hjust = -0.06: grid applies
    # hjust per LINE as a fraction of that line's own width, so the
    # shorter second line ("(15)") sat left of the first (QA review).
    cap <- data.frame(x = e + 3 * (e - s) / (TRACK_W_IN[[side]] * 72), y = y_max * 0.5,
                      lab = paste0(caption, "\n", caption_n))
    p <- p + geom_text(data = cap, aes(x = x, y = y, label = lab),
                       hjust = 0, vjust = 0.5, size = SZ_COHORT,
                       lineheight = 1.25,
                       family = FONT_FAMILY, color = "grey15",
                       inherit.aes = FALSE)
  }

  if (src_box)
    p <- p + zoom_source_box(L$zoom[1], L$zoom[2], s, e, y_max)

  p + theme(plot.margin = trk_gutter(side))
}

# ------------------------------------------------------------------
# Gene model + coordinate row.
# ------------------------------------------------------------------
EXON_RDS <- file.path(CACHE, "panel_c_exons.rds")
if (!file.exists(EXON_RDS)) {
  cat("[panel c] building", EXON_RDS, "(hg38_granges.RDS is 26 MB) ...\n")
  gencode <- readRDS(file.path(QROOT, "inputs/gencode_v46_granges.RDS"))
  mane <- data.table::fread(
    file.path(QROOT, "inputs/MANE.GRCh38.v1.4.summary.txt"),
    data.table = FALSE)
  ex_list <- lapply(names(LOCI), function(g) {
    ex <- gencode[gencode$gene_name == g & gencode$type == "exon"]
    tx <- mane$Ensembl_nuc[mane$symbol == g]
    e2 <- ex[ex$transcript_id %in% tx]
    if (length(e2) == 0) e2 else e2
  })
  names(ex_list) <- names(LOCI)
  dir.create(CACHE, showWarnings = FALSE)
  saveRDS(ex_list, EXON_RDS)
}
EXONS <- readRDS(EXON_RDS)

make_gene_rows <- function(L, side) {
  ex    <- EXONS[[L$gene]]
  tx_lo <- min(start(ex)); tx_hi <- max(end(ex))
  span  <- tx_hi - tx_lo
  brks  <- scales::breaks_pretty(5)(c(tx_lo, tx_hi))
  # Drop breaks that would collide with the chromosome label parked at
  # the left edge, or run off the right-hand end of the gene model.
  # Both cutoffs are DERIVED from the labels' real drawn width at this
  # point size rather than being eyeballed offsets.
  lab_w <- 9   * SZ_AXIS * .pt * 0.56 / 72        # inches, "185590000"
  chr_w <- 4.5 * SZ_AXIS * .pt * 0.56 / 72        # inches, "chr3"
  tw    <- TRACK_W_IN[[side]]
  at    <- (brks - tx_lo) / span * tw
  # The axis row's right gutter is empty (the cohort captions live on
  # the track rows), so the last label may overhang a little.
  keep  <- at - lab_w / 2 > chr_w + 0.05 & at + lab_w / 2 < tw + 0.22 &
           brks >= tx_lo & brks <= tx_hi
  brks  <- brks[keep]

  trap <- trapezoid_gg(tx_lo, tx_hi, L$wide[1], L$wide[2],
                       fill = LAV, alpha = 1) +
    theme(plot.margin = gutter(side))
  gene <- transcript_zoom_gg(
            ex, tx_lo, tx_hi,
            highlight = L$wide, highlight_fill = LAV, highlight_alpha = 1,
            tx_label = L$gene, exon_h = 0.85, n_arrows = 26,
            label_size = SZ_GENE, label_family = FONT_FAMILY,
            x_ticks = brks, tick_h = 0.55, label_side = L$label_side) +
    theme(plot.margin = gutter(side))
  axis <- ggplot() +
    geom_text(data = data.frame(x = brks, y = 1,
                                lab = format(brks, scientific = FALSE,
                                             trim = TRUE)),
              aes(x = x, y = y, label = lab), hjust = 0.5, vjust = 1,
              size = SZ_AXIS, family = FONT_FAMILY, color = "grey15") +
    geom_text(data = data.frame(x = tx_lo, y = 1, lab = L$chrom),
              aes(x = x, y = y, label = lab), hjust = 0, vjust = 1,
              size = SZ_AXIS, family = FONT_FAMILY, color = "grey15") +
    scale_x_continuous(limits = c(tx_lo, tx_hi), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, 1), clip = "off") +
    theme_void() +
    theme(plot.margin      = gutter(side),
          plot.background  = element_rect(fill = NA, color = NA),
          panel.background = element_rect(fill = NA, color = NA))
  list(trap = trap, gene = gene, axis = axis)
}

# ------------------------------------------------------------------
# Assemble one quadrant.
# ------------------------------------------------------------------
build_quadrant <- function(key, side, box_hi = NULL, box_ct = NULL) {
  L  <- LOCI[[key]]
  wd <- region_dir(L, "wide")
  hi <- read_links(file.path(wd, LK_HIGH), L$wide[1], L$wide[2])
  ct <- read_links(file.path(wd, LK_CTRL), L$wide[1], L$wide[2])
  nt <- derive_offset(hi, ct)
  cat(sprintf("[panel c] %-5s  high PSI %s | ctrl PSI %s | %s\n", key,
              paste(sort(hi$psi, decreasing = TRUE), collapse = ", "),
              paste(sort(ct$psi, decreasing = TRUE), collapse = ", "),
              nt$label))

  has_zoom <- !is.null(L$zoom)
  cap      <- if (side == "right") HIGH_LAB else NULL

  # The Controls inset's leaders drop past the High track (see the
  # zoom_callout_gg call below); keep the High track's PSI numbers out
  # of their way.
  ctrl_drop <- H_TRK / H_CALL + (TRK_GAP_PT / 72) / CALL_H_IN
  trk_hi <- make_track(L, "high", BLUE, side,
                       caption = cap, caption_n = HIGH_N,
                       src_box = has_zoom,
                       nt = if (has_zoom) NULL else nt,
                       extra_obst = if (has_zoom)
                         leader_obst(L, box_ct, side, ctrl_drop)
                       else list())
  # nt on the High track only -- see nt_place.
  trk_ct <- make_track(L, "ctrl", GREY, side,
                       caption = if (side == "right") CTRL_LAB else NULL,
                       caption_n = CTRL_N, src_box = has_zoom)
  g <- make_gene_rows(L, side)

  if (has_zoom) {
    # src_drop is in units of the CALLOUT row's own height, so it
    # survives the figure being rescaled (patchwork heights are
    # relative). EACH callout leads from ITS OWN track's source box:
    # the High (blue) inset drops 0 -- the top of the High track, i.e.
    # the callout panel's own bottom edge -- and the Controls (grey)
    # inset drops one whole track row, H_TRK / H_CALL, which is the top
    # edge of the grey source box. The insets are ordered along x so the two pairs of
    # leaders cannot cross: the inset whose leaders stop HIGHER (High)
    # is the one nearer the source box -- see box_hi / box_ct below.
    specs <- list(
      list(bw = file.path(region_dir(L, "zoom"), BW_HIGH), links = hi,
           color = BLUE, cohort = "high", box = box_hi),
      list(bw = file.path(region_dir(L, "zoom"), BW_CTRL), links = ct,
           color = GREY, cohort = "ctrl", box = box_ct))
    # nt_label = NULL switches OFF the helper's decorative curved
    # leader; zoom_arc_layers() draws the real arc edge plus a plain
    # straight label instead (see the block above make_track).
    call <- zoom_callout_gg(
      callouts = lapply(specs, function(s)
        list(bw = s$bw, color = s$color,
             # One whole track PANEL plus the gap under it -- the drop
             # is measured in callout-PANEL heights, and the gap is a
             # margin, so it is not part of H_TRK / H_CALL.
             src_drop = if (s$cohort == "high") 0 else ctrl_drop,
             box = s$box)),
      chrom = L$chrom, wide_start = L$wide[1], wide_end = L$wide[2],
      zoom_start = L$zoom[1], zoom_end = L$zoom[2],
      box_y = ZOOM_BOX_Y, cov_frac = ZOOM_COVF, n_bins = ZOOM_NBINS,
      nt_label = NULL, nt_site = NULL,
      label_size = SZ_NT, label_family = FONT_FAMILY) +
      zoom_arc_layers(L, nt, specs, side) +
      theme(plot.margin = gutter(side))
    (call / trk_hi / trk_ct / g$trap / g$gene / g$axis) +
      plot_layout(heights = c(H_CALL, H_TRK, H_TRK, H_TRAP, H_GENE, H_AXIS))
  } else {
    (trk_hi / trk_ct / g$trap / g$gene / g$axis) +
      plot_layout(heights = c(H_TRK, H_TRK, H_TRAP, H_GENE, H_AXIS))
  }
}

# Callout box x-positions, as fractions of the wide region. The source
# box sits at the FAR RIGHT for SENP2 (90-97%) and the FAR LEFT for
# PRKDC (0-15%), so the pair of boxes is anchored to that same side to
# keep the leaders short and un-crossed.
#
# Ordering rule now that each inset leads to its OWN source box: the
# HIGH inset's leaders stop one track higher than the Controls', so
# the High inset must be the one NEARER the source box in x. Otherwise
# its leaders sweep across the Controls' leaders on their way to the
# source. That is why SENP2 (source at the right) puts High on the
# right, and PRKDC (source at the left) keeps High on the left.
q_senp2 <- build_quadrant("SENP2", "left",
                          box_hi = c(0.58, 0.98), box_ct = c(0.14, 0.54))
q_prkdc <- build_quadrant("PRKDC", "right",
                          box_hi = c(0.02, 0.42), box_ct = c(0.46, 0.86))
q_ezh2  <- build_quadrant("EZH2",  "left")
q_brd9  <- build_quadrant("BRD9",  "right")

top <- (q_senp2 | q_prkdc)
bot <- (q_ezh2  | q_brd9)

wrap_elements(full =
  (top / plot_spacer() / bot) +
    plot_layout(heights = c(TOP_H, QUAD_GAP_PT / 72 / C_SCALE, BOT_H)) +
  plot_annotation(
    title = paste0("Mutations are associated with recurrent mis-splicing\n",
                   "of cancer-relevant genes"),
    theme = theme(plot.title = element_text(
                    size = SIZE_TITLE, face = "bold", hjust = 0.5,
                    family = FONT_FAMILY, margin = margin(0, 0, 3, 0)),
                  plot.margin = margin(0, 0, 0, 0))))
}))
