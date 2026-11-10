# figure4/panel_e/build.R -- "Tumor suppressor genes and frequently
# mutated oncogenes are mis-spliced in samples with high loadings upon
# GSP #30".
#
# A 2x2 grid of loci (TP53BP1 / PBRM1 on top, MSH6 / ASXL1 below).
# Every quadrant is a vertical patchwork stack drawn ENTIRELY from raw
# bigwigs + MeanPSI .links with ../scripts/coverage_plot.R:
#
#   High GSP #30 coverage + arcs      sashimi_track_gg() + own arc layer
#   Matched Controls coverage + arcs  sashimi_track_gg() + own arc layer
#   lavender taper                    trapezoid_gg()
#   full MANE gene model              transcript_zoom_gg()
#   genomic coordinate labels         bespoke row (same x scale)
#
# The arc/label machinery is a COPY of panel_b/build.R's rather than a
# shared helper, so the two panels stay independently editable. Keep
# them in step by hand.

if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
panel_e <- panel_cached("e", local({

PY <- file.path(QROOT, "coverage_tracks/pygenome/YLPM1")

RED  <- "#C0392B"
GREY <- "#969696"
LAV  <- "#ccccff"        # same taper colour as figure 3 panel c

# Track labels: 15 top GSP #30 samples vs 150 matched controls.
HIGH_LAB <- "High GSP"
CTRL_LAB <- "Matched"
HIGH_N   <- "#30 (15)"
CTRL_N   <- "Controls\n(150)"

BW_HIGH <- "High_k30-3-..bw"
BW_CTRL <- "Matched_Controls-3-..bw"
LK_HIGH <- "MeanPSI.High_k30-3-..links"
LK_CTRL <- "MeanPSI.Matched_Controls-3-..links"

# ------------------------------------------------------------------
# Loci. `dir` IS the region (pygenome names the directory after it);
# y_max is that region's own ini max_value, shared by both tracks so
# the two cohorts are comparable.
# ------------------------------------------------------------------
LOCI <- list(
  TP53BP1 = list(gene = "TP53BP1", chrom = "chr15",
                 wide = c(43469677, 43474976),
                 y_max = 6.93896007537842),
  PBRM1   = list(gene = "PBRM1",   chrom = "chr3",
                 wide = c(52633846, 52642490),
                 y_max = 4.9629430770874),
  MSH6    = list(gene = "MSH6",    chrom = "chr2",
                 wide = c(47800733, 47803847),
                 y_max = 27.2977384567261),
  ASXL1   = list(gene = "ASXL1",   chrom = "chr20",
                 wide = c(32358130, 32366548),
                 y_max = 10.4825664520264)
)
# NB: the italic gene name's side is NOT a property of the locus and is
# not listed here -- it is solved in make_gene_rows() from where the
# zoom window falls on the transcript, so it never lands on the taper.
region_dir <- function(L)
  file.path(PY, L$gene, sprintf("%s_%d-%d", L$chrom, L$wide[1], L$wide[2]))

# ------------------------------------------------------------------
# Row heights. patchwork heights are RELATIVE, so only the ratios here
# are stable.
# ------------------------------------------------------------------
H_TRK  <- 0.36     # one coverage/sashimi track
H_TRAP <- 0.09     # lavender taper
H_GENE <- 0.15     # MANE gene model
H_AXIS <- 0.12     # genomic coordinate labels
STACK_H <- 2 * H_TRK + H_TRAP + H_GENE + H_AXIS

# Title band height, in inches. It does not change the layout (heights
# are relative); it makes TRACK_H_IN the track's true drawn height,
# which is what the PSI solver converts point sizes through.
E_TITLE_H  <- 0.29
# The high track's "0" is centred on its baseline, so it hangs half a
# text height into the control track below, where that track's y_max
# number sits. AX_BOT_PT gives each high/control boundary a gutter,
# taken off the content first and handed back as extra row height, so
# every track keeps the same drawn coverage height (the two cohorts of
# a locus share a y_max). Same as panel b.
AX_BOT_PT  <- 4.5
E_SCALE    <- (ROW_H[["e"]] - E_TITLE_H - 2 * AX_BOT_PT / 72) / (2 * STACK_H)
H_AXGAP    <- (AX_BOT_PT / 72) / E_SCALE        # the gutter, in row units
TRACK_H_IN <- H_TRK * E_SCALE

# Gutters, in pt. Only the RIGHT column carries the cohort captions
#, so the two columns need different margins to
# end up with comparable track widths.
MARG <- list(left  = c(30, 40),      # c(left, right) in pt
             right = c(22, 62))
gutter <- function(side) margin(0, MARG[[side]][2], 0, MARG[[side]][1], "pt")
TRACK_W_IN <- list(left  = PAGE_W / 2 - sum(MARG$left)  / 72,
                   right = PAGE_W / 2 - sum(MARG$right) / 72)

SZ_AXIS   <- SIZE_BASE * 0.80 / .pt   # "27.3" / "0" and coordinate labels
SZ_COHORT <- SIZE_BASE * 0.90 / .pt
SZ_PSI    <- SIZE_BASE * 0.80 / .pt
SZ_GENE   <- SIZE_BASE * 0.95 / .pt
AX_IN     <- 0.40                     # y-extreme labels, lines inward

# Cap height + descender of a string, in inches, at a ggplot `size`.
txt_h_in <- function(sz) sz * .pt * 0.72 / 72

# ------------------------------------------------------------------
# .links helpers. sashimi_track_gg filters s1 >= start & e2 <= end;
# reproduce that here so the arcs we reason about are exactly the arcs
# that get drawn.
# ------------------------------------------------------------------
read_links <- function(path, start, end, min_psi = 1) {
  lk <- read.table(path, header = FALSE, sep = "\t",
                   col.names = c("c1","s1","e1","c2","s2","e2","psi"))
  lk[lk$s1 >= start & lk$e2 <= end & lk$psi >= min_psi, , drop = FALSE]
}

# PSI as printed: the .links value at its own
# precision, but never bare-integer ("74" -> "74.0", "99" -> "99.0").
psi_fmt <- function(p) {
  s <- vapply(p, function(x) sprintf("%g", x), character(1))
  ifelse(grepl(".", s, fixed = TRUE), s, paste0(s, ".0"))
}

# ------------------------------------------------------------------
# MANE exons, cached (the GENCODE GRanges is 26 MB).
# ------------------------------------------------------------------
EXON_RDS <- file.path(CACHE, "panel_e_exons.rds")
if (!file.exists(EXON_RDS)) {
  cat("[panel e] building", EXON_RDS, "(hg38_granges.RDS is 26 MB) ...\n")
  gencode <- readRDS(file.path(QROOT, "inputs/gencode_v46_granges.RDS"))
  mane <- data.table::fread(
    file.path(QROOT, "inputs/MANE.GRCh38.v1.4.summary.txt"),
    data.table = FALSE)
  ex_list <- lapply(names(LOCI), function(g) {
    ex <- gencode[gencode$gene_name == g & gencode$type == "exon"]
    tx <- mane$Ensembl_nuc[mane$symbol == g]
    me <- ex[ex$transcript_id %in% tx]
    stopifnot(length(me) > 0)
    me
  })
  names(ex_list) <- names(LOCI)
  dir.create(CACHE, showWarnings = FALSE)
  saveRDS(ex_list, EXON_RDS)
}
EXONS <- readRDS(EXON_RDS)

# ------------------------------------------------------------------
# Arc geometry (same as panel_b/build.R).
#
# EACH FOOT SITS ON THE COVERAGE AT ITS OWN JUNCTION COORDINATE. The
# curve is the same parametric form sashimi_track_gg() uses:
#
#     y(t) = (1 - t) y1 + t y2 + sin(t pi)^ARC_P * amp
#
# with y1 / y2 the coverage AT s1 and AT e2. Every foot is drawn at the
# raw .links coordinate.
#
# ANCHOR_BP = 0: the foot takes the coverage of the BIN the junction
# coordinate falls in, which is a bin that IS drawn, so the foot sits
# ON the coverage silhouette. A wider window would walk up into the
# exon plateau and the foot would hang over the fill.
#
# Heights are solved in PEAK space (the drawn maximum of the curve),
# never as an absolute fraction of y_max: with the feet at coverage
# height an absolute apex can be BELOW its own feet and the arc would
# droop instead of bulge.
#
# The heights are SOLVED, arc by arc, bottom-up in the CONTAINMENT
# order (Kahn's algorithm, tie-broken on endpoint height then span).
# Arc i's amplitude is the smallest that satisfies all of:
#   * peak_i >= max(y1, y2) + ARC_RISE            (it reads as an arc)
#   * peak_i >= peak_j + ARC_SEP                  for every already
#     placed j that CONTAINS i or is contained by it
#   * min over their shared x of (y_i - y_j) >= ARC_GAP, same j
#     (dominance: i never passes under j)
# Each is monotone in the amplitude, so the binding one is found by
# bisection, and the amplitude for a chosen drawn peak is back-solved
# with uniroot (with a steep chord the maximum sits a long way from
# t = 0.5 -- TP53BP1's 12.9 runs 0.55 -> 0.01 and peaks at t = 0.24).
# DISJOINT pairs get no constraints at all, which is what keeps the
# budget affordable: TP53BP1's 12.9 costs 7.81 and 1.54 nothing.
# An infeasible dominance constraint is SKIPPED, not driven to the
# ceiling -- parking an arc on the ceiling makes everything above it
# cross.
#
# ---- Parameter values ------------------------------------------------
# ARC_RISE / ARC_SEP / ARC_GAP are VISUAL quantities (bulge, tier
# separation, clearance), so they are fixed in DRAWN INCHES -- the same
# inches as panel b -- and converted through this panel's track height,
# so both sashimi panels of the figure match.
#
# ARC_TRIM and ARC_P are dimensionless. Most arcs here share an endpoint
# with another (ASXL1's four all land on one acceptor), and at a shared
# foot the gap is 0 by construction, so the trim keeps the dominance
# tests satisfiable; 0.10 leaves the smallest residual dip.
# ------------------------------------------------------------------
ANCHOR_BP <- 0
ARC_P     <- 0.50     # bump shape, sin(t pi)^ARC_P (~ a half-ellipse)
ARC_TRIM  <- 0.10     # of the shared span ignored at each end
ARC_RISE  <- 0.0265 / TRACK_H_IN   # an arc's bulge above its higher foot
ARC_SEP   <- 0.0281 / TRACK_H_IN   # peak-to-peak separation, nested pairs
ARC_GAP   <- 0.0099 / TRACK_H_IN   # drawn clearance, nested pairs

# ARC_CEIL is DERIVED, not chosen: it is exactly place_psi()'s own
# ceiling for a label CENTRE, so the top arc's apex is always a legal
# place for that arc's number. Both come from this one expression, so a
# change to SZ_PSI or to the row heights moves them together.
lab_ceil_frac <- function(trk_h_in)
  1 - (txt_h_in(SZ_PSI) / 2 + 0.45 / 72) / trk_h_in
ARC_CEIL  <- lab_ceil_frac(TRACK_H_IN)

arc_bump <- function(t) sin(t * pi)^ARC_P

# Drawn maximum of (1-t) y1 + t y2 + bump(t) * amp over t in [0, 1].
arc_peak_val <- function(y1, y2, amp, tt = seq(0, 1, length.out = 1201))
  max((1 - tt) * y1 + tt * y2 + arc_bump(tt) * amp)

# Amplitude whose drawn maximum is exactly `p`.
amp_for_peak <- function(y1, y2, p) {
  f <- function(a) arc_peak_val(y1, y2, a) - p
  if (f(0) >= 0) 0 else uniroot(f, c(0, 6), tol = 1e-9)$root
}

# 600 points, not 200: with ARC_P = 0.5 the limbs turn hard within a
# few percent of t and a coarse grid visibly chamfers the feet.
arc_curve <- function(a, n = 600) {
  t <- seq(0, 1, length.out = n)
  data.frame(x = a$x1 + t * (a$x2 - a$x1),
             y = ((1 - t) * a$y1 + t * a$y2 + arc_bump(t) * a$amp) * a$y_max)
}
arc_y_at <- function(a, x) {
  t  <- (x - a$x1) / (a$x2 - a$x1)
  tc <- pmin(pmax(t, 0), 1)
  y  <- ((1 - tc) * a$y1 + tc * a$y2 + arc_bump(tc) * a$amp) * a$y_max
  ifelse(t < 0 | t > 1, NA_real_, y)
}
# y in FRACTION-of-y_max units, for the solver (which works there).
arc_frac_at <- function(y1, y2, amp, x1, x2, x) {
  t <- (x - x1) / (x2 - x1)
  (1 - t) * y1 + t * y2 + arc_bump(t) * amp
}

# The x at which two arcs are compared; the outer ARC_TRIM of the span
# they share is dropped (see above).
shared_xs <- function(a1, a2, b1, b2, n = 240) {
  lo <- max(a1, b1); hi <- min(a2, b2)
  if (hi <= lo) return(numeric(0))
  m <- ARC_TRIM * (hi - lo)
  seq(lo + m, hi - m, length.out = n)
}

nests <- function(a1, a2, b1, b2)
  (a1 >= b1 && a2 <= b2) || (b1 >= a1 && b2 <= a2)

# Signed vertical gap, in fractions of y_max, between every pair of
# arcs that is NOT a forced crossing: (upper) - (lower) over the x they
# share. Negative anywhere means an arc passes UNDER one it is tiered
# above. Printed on every build.
arc_min_gap <- function(arcs) {
  n <- nrow(arcs)
  out <- list(gap = Inf, worst = "", forced = character(0))
  if (n < 2) { out$gap <- NA_real_; return(out) }
  for (i in seq_len(n - 1)) for (j in (i + 1):n) {
    xs <- shared_xs(arcs$x1[i], arcs$x2[i], arcs$x1[j], arcs$x2[j])
    if (!length(xs)) next                                  # disjoint
    if (!nests(arcs$x1[i], arcs$x2[i], arcs$x1[j], arcs$x2[j])) {
      out$forced <- c(out$forced,
                      sprintf("%s/%s [partial overlap]", arcs$lab[i],
                              arcs$lab[j]))
      next
    }
    u <- if (arcs$peak[i] > arcs$peak[j]) i else j
    l <- if (arcs$peak[i] > arcs$peak[j]) j else i
    g <- min(arc_frac_at(arcs$y1[u], arcs$y2[u], arcs$amp[u],
                         arcs$x1[u], arcs$x2[u], xs) -
             arc_frac_at(arcs$y1[l], arcs$y2[l], arcs$amp[l],
                         arcs$x1[l], arcs$x2[l], xs))
    if (g <= 0)
      out$forced <- c(out$forced, sprintf("%s under %s [%+.3f, ceiling]",
                                          arcs$lab[u], arcs$lab[l], g))
    else if (g < out$gap) {
      out$gap <- g; out$worst <- sprintf("%s over %s", arcs$lab[u], arcs$lab[l])
    }
  }
  if (is.infinite(out$gap)) out$gap <- NA_real_
  out
}

arc_frame <- function(lk, d, y_max) {
  n <- nrow(lk)
  if (n == 0) return(NULL)
  lk <- lk[order(lk$e2 - lk$s1), , drop = FALSE]

  covat <- function(x) {
    w <- d$bin_center >= x - ANCHOR_BP & d$bin_center <= x + ANCHOR_BP
    if (!any(w)) d$coverage[which.min(abs(d$bin_center - x))]
    else         max(d$coverage[w])
  }
  y1   <- vapply(lk$s1, covat, numeric(1)) / y_max
  y2   <- vapply(lk$e2, covat, numeric(1)) / y_max
  E    <- pmax(y1, y2)
  span <- lk$e2 - lk$s1

  # nested_in[a, b]: arc a lies inside arc b, so a is solved (and sits)
  # BELOW it. Kahn's algorithm over that partial order, taking the
  # lowest (E, span) among the arcs whose nested children are all done.
  nested_in <- matrix(FALSE, n, n)
  for (a in seq_len(n)) for (b in seq_len(n))
    nested_in[a, b] <- a != b && lk$s1[a] >= lk$s1[b] &&
                       lk$e2[a] <= lk$e2[b] && span[a] < span[b]
  ord <- integer(0); rem <- rep(TRUE, n)
  for (k in seq_len(n)) {
    ready <- which(rem & vapply(seq_len(n),
                     function(b) !any(rem & nested_in[, b]), logical(1)))
    if (!length(ready)) ready <- which(rem)   # unreachable: strict order
    pick <- ready[order(E[ready], span[ready])][1]
    ord  <- c(ord, pick); rem[pick] <- FALSE
  }

  amp <- rep(NA_real_, n)
  for (k in seq_len(n)) {
    i    <- ord[k]
    cap  <- amp_for_peak(y1[i], y2[i], ARC_CEIL)
    need <- amp_for_peak(y1[i], y2[i], E[i] + ARC_RISE)
    for (j in ord[seq_len(k - 1)]) {
      if (!nests(lk$s1[i], lk$e2[i], lk$s1[j], lk$e2[j])) next
      xs <- shared_xs(lk$s1[i], lk$e2[i], lk$s1[j], lk$e2[j])
      if (!length(xs)) next
      # (a) stay ARC_SEP above j's peak
      need <- max(need, amp_for_peak(y1[i], y2[i],
                          arc_peak_val(y1[j], y2[j], amp[j]) + ARC_SEP))
      # (b) dominate j by ARC_GAP everywhere they share x. Monotone in
      #     the amplitude (sin(t pi) > 0 in the interior), so bisect.
      yj <- arc_frac_at(y1[j], y2[j], amp[j], lk$s1[j], lk$e2[j], xs)
      gp <- function(a) min(arc_frac_at(y1[i], y2[i], a,
                                        lk$s1[i], lk$e2[i], xs) - yj)
      if (gp(cap) >= ARC_GAP) {
        lo <- need; hi <- cap
        if (gp(lo) < ARC_GAP) {
          for (it in seq_len(60)) {
            mid <- (lo + hi) / 2
            if (gp(mid) >= ARC_GAP) hi <- mid else lo <- mid
          }
          need <- max(need, hi)
        }
      }
      # else: not reachable inside the ceiling. Leave `need` alone --
      # spending the budget on it would not clear the pair and would
      # push every arc above this one through the top of the track.
      # arc_min_gap() reports the pair as a crossing.
    }
    amp[i] <- min(need, cap)
  }
  peak <- vapply(seq_len(n),
                 function(i) arc_peak_val(y1[i], y2[i], amp[i]), numeric(1))

  data.frame(x1 = lk$s1, x2 = lk$e2, psi = lk$psi,
             y1 = y1, y2 = y2, amp = amp, peak = peak, y_max = y_max,
             lab = psi_fmt(lk$psi), stringsAsFactors = FALSE)
}

# ------------------------------------------------------------------
# PSI labels: CONTOUR-STYLE, i.e. the number sits in a break in the arc
# it names, exactly the way a contour map labels an isoline. Same as
# panel_b/build.R. Stacking numbers above the arcs cannot work here:
# tiers are much closer together than a digit is tall, so a stacked
# number would not say which arc it names.
#
# So the arcs are CUT instead. Every label gets a box; each arc's curve
# is split wherever it enters one, and only the surviving runs are
# drawn. Nothing white is painted, the PDF stays vector-clean, and the
# invariant holds exactly: no line ever touches a glyph. The number
# lands ON the arc it belongs to, so attribution needs no leader.
#
# The placement search is scored, not first-fit. Hard constraints:
# inside the panel, under the ceiling, clear of the coverage area,
# clear of every protected landing, and no overlap with a label already
# placed. Among what is left it minimises
#   CUT_W * (number of OTHER arcs the box would cut)
#     + 3.0 * (label pushed off its own arc)
#     + 3.0 * min(on-page slope of its own arc under the box, 2)
#     + 0.9 * |t - 0.5|
# The SLOPE term is what keeps the number readable as part of its arc:
# a break in a near-horizontal line reads as a labelled line, a break
# in a steep one reads as a line that stopped.
# ------------------------------------------------------------------
LAB_TS    <- seq(0.10, 0.90, by = 0.02)
CUT_W     <- 1.5
FOOT_KEEP <- 5.0   # pt of arc that must survive leading into every foot

# Helvetica advance widths: digits are 0.556 em, the period 0.278 em.
# The width sizes the cut in the arc, so it is computed exactly.
str_w_em <- function(s)
  vapply(strsplit(s, "", fixed = TRUE),
         function(cs) sum(ifelse(cs == ".", 0.278, 0.556)), numeric(1))

place_psi <- function(arcs, d, y_max, trk_h_in, trk_w_in, s, e) {
  psi_pt  <- SZ_PSI * .pt
  txt_h   <- txt_h_in(SZ_PSI)
  x_pt    <- (e - s)  / (trk_w_in * 72)     # data-x per pt
  y_pt    <- y_max    / (trk_h_in * 72)     # data-y per pt
  n       <- nrow(arcs)
  hw      <- str_w_em(arcs$lab) / 2 * psi_pt * x_pt
  hh      <- txt_h / 2 * y_max / trk_h_in
  # TWO break sizes, which is the whole trick for reading a crowded
  # track. The arc the number NAMES is cut generously (bw/bh), so the
  # number visibly sits in a gap in that line. Every OTHER arc the box
  # touches is cut as tightly as "no line may touch a glyph" allows
  # (bw2/bh2) -- a hairline notch that the eye reads straight through.
  # That is the contour-map convention: only the labelled isoline
  # breaks, so in a bundle the number still says which arc it names.
  bw      <- hw + 1.7 * x_pt
  bh      <- hh + 0.9 * y_pt
  bw2     <- hw + 0.5 * x_pt
  bh2     <- hh + 0.5 * y_pt
  # Extra air demanded BETWEEN two numbers, on top of their boxes:
  # boxes that merely fail to overlap still print as "4.8:1.28".
  air_x   <- 6.0 * x_pt
  air_y   <- 1.0 * y_pt
  # How far a label may be pushed off its own arc when it cannot sit on
  # it (a low arc whose apex is inside the coverage).
  off     <- hh + 1.6 * y_pt
  # Highest label CENTRE whose GLYPHS stay in the track. bh carries the
  # box's air as well as the glyphs, so testing `ly + bh > y_max - bh`
  # would charge the pad twice.
  ceil_y  <- lab_ceil_frac(trk_h_in) * y_max     # == ARC_CEIL * y_max

  # The protected run leading into each foot, as drawn points.
  prot <- lapply(seq_len(n), function(j) {
    cv  <- arc_curve(arcs[j, ], n = 200)
    dl  <- c(0, cumsum(sqrt((diff(cv$x) / x_pt)^2 + (diff(cv$y) / y_pt)^2)))
    tot <- dl[length(dl)]
    cv[dl < FOOT_KEEP | (tot - dl) < FOOT_KEEP, c("x", "y"), drop = FALSE]
  })

  px <- py <- rep(NA_real_, n)
  # Leader endpoints for labels that had to sit beside their arc; NA
  # means the number is ON its arc and needs no leader.
  lead_x <- lead_y <- rep(NA_real_, n)
  # Most-constrained first: increasing (arc span / label width). The
  # number that barely fits on its own arc chooses before the ones with
  # a whole track to slide along.
  for (i in order((arcs$x2 - arcs$x1) / (2 * bw))) {
    a     <- arcs[i, ]
    done  <- which(!is.na(py))
    best  <- NULL; best_s <- Inf
    for (sgn in c(0, 1, -1)) {
      for (tt in LAB_TS) {
        lx <- a$x1 + tt * (a$x2 - a$x1)
        if (lx - bw[i] < s || lx + bw[i] > e) next
        ly <- arc_y_at(a, lx) + sgn * off
        if (ly > ceil_y || ly - bh < 0) next
        inw  <- d$bin_center >= lx - bw[i] & d$bin_center <= lx + bw[i]
        covb <- if (any(inw)) max(d$coverage[inw]) else 0
        if (ly - bh < covb + 0.6 * y_pt) next
        if (length(done) &&
            any(abs(px[done] - lx) < (bw[done] + bw[i] + air_x) &
                abs(py[done] - ly) < (bh + bh + air_y))) next
        # The break must never eat an arc's LANDING: a number can be
        # wider than a short arc and would otherwise cut a neighbour
        # that shares its splice site short of the exon. So the last
        # FOOT_KEEP points of every arc, at BOTH ends, are off limits
        # to every box.
        if (any(vapply(seq_len(n), function(j) {
              hwj <- if (j == i) bw[i] else bw2[i]
              hhj <- if (j == i) bh   else bh2
              any(abs(prot[[j]]$x - lx) < hwj & abs(prot[[j]]$y - ly) < hhj)
            }, logical(1)))) next
        # How many OTHER arcs would this box cut?
        probe <- seq(lx - bw[i], lx + bw[i], length.out = 11)
        cuts  <- sum(vapply(setdiff(seq_len(n), i), function(j) {
                   yj <- arc_y_at(arcs[j, ], probe)
                   any(!is.na(yj) & abs(yj - ly) < bh)
                 }, logical(1)))
        # On-page slope of the OWN arc across the label's width.
        rise  <- abs(diff(range(arc_y_at(a, probe), na.rm = TRUE))) /
                 y_max * trk_h_in
        run   <- 2 * bw[i] / (e - s) * trk_w_in
        slope <- if (run > 0) rise / run else 0
        sc <- CUT_W * cuts + 3.0 * abs(sgn) +
              3.0 * min(slope, 2) + 0.9 * abs(tt - 0.5)
        if (sc < best_s) { best_s <- sc; best <- c(lx, ly) }
      }
    }
    if (is.null(best)) {
      # No slot ON the arc. This happens when an arc is both SHORT and
      # LOW and shares a foot with a tall neighbour (MSH6's 11.5), so
      # there is nowhere on or just above it that is not inside the
      # neighbour.
      #
      # Put it BESIDE the arc instead, level with its apex, on whichever
      # flank has room, and run a short leader back to the apex. Scan
      # outward so the label sits as close to its arc as it can.
      apex_x <- (a$x1 + a$x2) / 2
      apex_y <- arc_y_at(a, apex_x)
      step   <- 0.35 * bw[i]
      cand   <- NULL
      # Also try heights ABOVE the apex, so a label blocked by intron
      # coverage at apex height can still sit right next to its arc.
      for (k in 1:40) for (dir in c(-1, 1)) for (dy in c(0, 1, 2, 3, 4)) {
        # R's `break` only leaves the INNERMOST loop; this stops the
        # scan at the nearest clear slot.
        if (!is.null(cand)) break
        lx <- if (dir < 0) a$x1 - k * step - bw[i] else a$x2 + k * step + bw[i]
        if (lx - bw[i] < s || lx + bw[i] > e) next
        ly <- apex_y + dy * 0.9 * bh
        if (ly > ceil_y || ly - bh < 0) next
        inw  <- d$bin_center >= lx - bw[i] & d$bin_center <= lx + bw[i]
        covb <- if (any(inw)) max(d$coverage[inw]) else 0
        if (ly - bh < covb + 0.6 * y_pt) next
        # clear of EVERY arc, including its own
        probe <- seq(lx - bw[i], lx + bw[i], length.out = 11)
        if (any(vapply(seq_len(n), function(j) {
              yj <- arc_y_at(arcs[j, ], probe)
              any(!is.na(yj) & abs(yj - ly) < bh + 0.4 * y_pt)
            }, logical(1)))) next
        # clear of every protected landing and of labels already placed
        if (any(vapply(seq_len(n), function(j)
              any(abs(prot[[j]]$x - lx) < bw2[i] &
                  abs(prot[[j]]$y - ly) < bh2), logical(1)))) next
        if (length(done) &&
            any(abs(px[done] - lx) < (bw[done] + bw[i] + air_x) &
                abs(py[done] - ly) < (bh + bh + air_y))) next
        cand <- c(lx, ly); break
      }
      if (!is.null(cand)) {
        best <- cand
        lead_x[i] <- apex_x; lead_y[i] <- apex_y
      } else {
        best <- c(apex_x, min(ceil_y, apex_y + off))
        warning(sprintf("place_psi: no clear slot for PSI label %s", a$lab),
                call. = FALSE)
      }
    }
    px[i] <- best[1]; py[i] <- best[2]
  }
  # Every column travels as DATA, never as an aes() reference to a
  # scalar in this environment -- a panel restored from
  # panel_cached()'s saveRDS renders blank otherwise.
  data.frame(x = px, y = py, lab = arcs$lab,
             arc = seq_len(n), bw = bw, bh = bh, bw2 = bw2, bh2 = bh2,
             lead_x = lead_x, lead_y = lead_y,
             stringsAsFactors = FALSE)
}

# ------------------------------------------------------------------
# Draw one arc, split wherever it runs into a PSI label's box. NA-ing
# the points would make geom_path drop rows with a warning and join
# across the hole, so the surviving points are grouped into contiguous
# runs and the group travels as a data column.
# ------------------------------------------------------------------
arc_layer <- function(a, boxes, color, linewidth, x_pt, y_pt, min_pt = 3,
                      own = NA_integer_) {
  cv  <- arc_curve(a)
  hit <- rep(FALSE, nrow(cv))
  if (!is.null(boxes) && nrow(boxes) > 0)
    for (j in seq_len(nrow(boxes))) {
      # Wide break for this arc's own number, hairline notch for the rest.
      mine <- !is.na(own) && identical(boxes$arc[j], own)
      hw_j <- if (mine) boxes$bw[j] else boxes$bw2[j]
      hh_j <- if (mine) boxes$bh[j] else boxes$bh2[j]
      hit  <- hit | (abs(cv$x - boxes$x[j]) <= hw_j &
                     abs(cv$y - boxes$y[j]) <= hh_j)
    }
  cv$grp <- cumsum(c(1L, abs(diff(hit))))
  cv     <- cv[!hit, , drop = FALSE]
  if (nrow(cv) == 0) return(NULL)
  # A run only a few points long is the sliver of arc left between two
  # neighbouring numbers. Drawn, it reads as a dash -- "74.0 - 4.8" --
  # not as a line, so drop it and let the two labels share one break.
  len_pt <- vapply(split(cv, cv$grp), function(g)
              sum(sqrt((diff(g$x) / x_pt)^2 + (diff(g$y) / y_pt)^2)),
              numeric(1))
  cv <- cv[as.character(cv$grp) %in% names(len_pt)[len_pt >= min_pt], ,
           drop = FALSE]
  if (nrow(cv) == 0) return(NULL)
  geom_path(data = cv, aes(x = x, y = y, group = grp), inherit.aes = FALSE,
            color = color, linewidth = linewidth, lineend = "round")
}

# ------------------------------------------------------------------
# One coverage / sashimi track.
#
# sashimi_track_gg draws the coverage (and owns the theming); its own
# arc engine is switched OFF with min_psi = Inf because it cannot
# express per-arc solved peaks. The FORM of the curve is identical to
# the helper's -- see arc_frame.
#
# The arcs are built HERE, not by the caller, because their feet are
# sampled from this track's own coverage.
# ------------------------------------------------------------------
make_track <- function(L, cohort, col, side, lk_df,
                       caption = NULL, caption_n = NULL) {
  wd    <- region_dir(L)
  bw    <- file.path(wd, if (cohort == "high") BW_HIGH else BW_CTRL)
  lk    <- file.path(wd, if (cohort == "high") LK_HIGH else LK_CTRL)
  y_max <- L$y_max
  s <- L$wide[1]; e <- L$wide[2]
  trk_w <- TRACK_W_IN[[side]]

  p <- sashimi_track_gg(
        bw_path = bw, chrom = L$chrom, start = s, end = e,
        track_color = col, track_label = "",
        links_path = lk, min_psi = Inf, y_max = y_max,
        label_family = FONT_FAMILY,
        # The shared axis is OFF: this panel draws its own below, in the
        # same chrome as panel b, because its numbers must be pushed
        # inward (AX_IN) where tracks butt flush.
        y_axis = FALSE,
        track_height_in = TRACK_H_IN, track_width_in = trk_w)

  d    <- bigwig_binned(bw, L$chrom, s, e)
  arcs <- arc_frame(lk_df, d, y_max)

  if (!is.null(arcs) && nrow(arcs) > 0) {
    # Regression check against "junctions go under other junctions".
    ag <- arc_min_gap(arcs)
    cat(sprintf("[panel e]   %-7s %-4s feet/peaks (of y_max): %s\n",
                L$gene, cohort,
                paste(sprintf("%s %.3f->%.3f ^%.3f", arcs$lab, arcs$y1,
                              arcs$y2, arcs$peak), collapse = " | ")))
    cat(sprintf("[panel e]   %-7s %-4s min signed gap %s   %s\n",
                L$gene, cohort,
                if (is.na(ag$gap)) "n/a"
                else sprintf("%+.4f of y_max (%s)", ag$gap, ag$worst),
                if (length(ag$forced))
                  paste0("| crossings: ", paste(ag$forced, collapse = "; "))
                else "| no crossings"))
    if (!is.na(ag$gap) && ag$gap <= 0)
      stop(sprintf(paste("[panel e] %s %s: a SOLVABLE pair still crosses",
                         "(min signed gap %+.4f)."), L$gene, cohort, ag$gap))

    # Labels are solved FIRST: the arcs are then drawn with a gap at
    # every label box, so nothing overprints anything.
    lab_df <- place_psi(arcs, d, y_max, TRACK_H_IN, trk_w, s, e)
    # Every number must be sitting IN its own arc. "off" is how far its
    # centre had to leave the curve, in fractions of y_max; anything
    # non-zero means the contour idiom broke down for that label.
    cat(sprintf("[panel e]   %-7s %-4s labels: %s\n", L$gene, cohort,
                paste(sprintf("%s t=%.2f off=%+.3f", lab_df$lab,
                        (lab_df$x - arcs$x1) / (arcs$x2 - arcs$x1),
                        (lab_df$y - vapply(seq_len(nrow(arcs)), function(i)
                           arc_y_at(arcs[i, ], lab_df$x[i]), numeric(1))) /
                          y_max),
                      collapse = " | ")))
    for (i in seq_len(nrow(arcs)))
      p <- p + arc_layer(arcs[i, ], lab_df, "black",
                         0.25 + 0.60 * sqrt(arcs$psi[i] / 100),
                         x_pt = (e - s) / (trk_w * 72),
                         y_pt = y_max  / (TRACK_H_IN * 72), own = i)
    # A number that had to sit BESIDE its arc gets a short leader back
    # to that arc's apex, so it still reads as belonging to it. Labels
    # sitting ON their arc have lead_x = NA and draw nothing.
    ld <- lab_df[!is.na(lab_df$lead_x), , drop = FALSE]
    if (nrow(ld) > 0) {
      # Stop the leader at the label's box edge, not its centre, or it
      # runs under the digits.
      ld$xe <- ld$x + sign(ld$lead_x - ld$x) * ld$bw
      p <- p + geom_segment(
        data = ld, aes(x = xe, y = y, xend = lead_x, yend = lead_y),
        colour = "grey45", linewidth = 0.22, inherit.aes = FALSE)
    }
    p <- p + geom_text(
      data = lab_df,
      aes(x = x, y = y, label = lab), hjust = 0.5, vjust = 0.5,
      size = SZ_PSI, family = FONT_FAMILY, color = "grey15",
      inherit.aes = FALSE)
  }

  # y-axis extremes, hard against the left edge (clip = "off" puts them
  # in the left gutter). One decimal.
  #
  # Both labels are pushed INWARD by AX_IN lines, so the gap between
  # this track's "0" and the next track's "27.3" is 2 * AX_IN lines
  # (the rows butt flush), enough that they do not read as one pair.
  #
  # Axis chrome: a rule over 0..y_max, a tick at each end, and each
  # number hung off its tick -- the same as panel b and as sashimi_track_gg's own
  # axis in figures 1-3 (grey15, 0.3 pt rule, 2.5 pt ticks, 1.5 pt air).
  # The 0 is CENTRED on the baseline; the top number hangs just inside
  # its tick. Every value is a DATA COLUMN (panel_cached's saveRDS).
  x_pt_ax <- (e - s) / (trk_w * 72)
  p <- p +
    geom_segment(data = data.frame(x = s, xend = s, y = 0, yend = y_max),
                 aes(x = x, xend = xend, y = y, yend = yend),
                 inherit.aes = FALSE, color = "grey15", linewidth = 0.3,
                 lineend = "square") +
    geom_segment(data = data.frame(x = s, xend = s - 2.5 * x_pt_ax,
                                   y = c(0, y_max), yend = c(0, y_max)),
                 aes(x = x, xend = xend, y = y, yend = yend),
                 inherit.aes = FALSE, color = "grey15", linewidth = 0.3,
                 lineend = "square") +
    geom_text(data = data.frame(x = s - 4.0 * x_pt_ax, y = c(y_max, 0),
                                lab = c(sprintf("%.1f", y_max), "0"),
                                vj = c(1, 0.5)),
              aes(x = x, y = y, label = lab, vjust = vj),
              hjust = 1, size = SZ_AXIS, family = FONT_FAMILY,
              color = "grey15", inherit.aes = FALSE)

  if (!is.null(caption))
    p <- p + geom_text(
      data = data.frame(x = e + 3 * (e - s) / (trk_w * 72), y = y_max * 0.5,
                        lab = paste0(caption, "\n", caption_n)),
      # hjust 0 + fixed 3 pt gap (hjust -0.06 shifts each line by its own width)
      aes(x = x, y = y, label = lab), hjust = 0, vjust = 0.5,
      size = SZ_COHORT, lineheight = 1.15, family = FONT_FAMILY,
      color = "grey15", inherit.aes = FALSE)

  mg <- gutter(side)
  if (cohort == "high") mg[3] <- grid::unit(AX_BOT_PT, "pt")
  p + theme(plot.margin = mg)
}

# ------------------------------------------------------------------
# Gene model + coordinate rows.
# ------------------------------------------------------------------
make_gene_rows <- function(L, side) {
  ex    <- EXONS[[L$gene]]
  tx_lo <- min(GenomicRanges::start(ex)); tx_hi <- max(GenomicRanges::end(ex))
  span  <- tx_hi - tx_lo
  # ASXL1's pygenome region starts 201 bp BEFORE its MANE transcript
  # does; clamp so the taper's feet and the highlight rect stay on the
  # gene model instead of spilling into the left gutter.
  win   <- c(max(L$wide[1], tx_lo), min(L$wide[2], tx_hi))
  brks  <- scales::breaks_pretty(5)(c(tx_lo, tx_hi))
  # Drop breaks that would collide with the chromosome label parked at
  # the left edge, or run off the right-hand end of the gene model.
  # Both cutoffs are DERIVED from the labels' real drawn width at this
  # point size rather than being eyeballed offsets.
  lab_w <- 8   * SZ_AXIS * .pt * 0.56 / 72
  chr_w <- 5.5 * SZ_AXIS * .pt * 0.56 / 72
  tw    <- TRACK_W_IN[[side]]
  at    <- (brks - tx_lo) / span * tw
  keep  <- at - lab_w / 2 > chr_w + 0.05 & at + lab_w / 2 < tw + 0.22 &
           brks >= tx_lo & brks <= tx_hi
  brks  <- brks[keep]

  # ---- which end the italic gene name hangs off, DERIVED ----------
  # transcript_zoom_gg anchors the name on the transcript's own tx_lo /
  # tx_hi edge and draws it INSIDE the gene row -- the same row the
  # lavender highlight rect spans from y = 0 to y = 2. So the name
  # lands on the taper's foot whenever the zoom window reaches that end
  # of the gene, which is exactly what happened to ASXL1 (its window
  # starts AT tx_lo).
  #
  # The rule, applied identically to all four loci: put the name on an
  # end the window does not reach, testing with the label's REAL drawn
  # width converted to bp; prefer LEFT when both ends are free, so the
  # column reads consistently and the only exception is forced by the
  # data rather than chosen. Helvetica caps/digits run ~0.62 em -- the
  # same estimator this file already uses for `lab_w` above.
  lab_in   <- nchar(L$gene) * 0.62 * SZ_GENE * .pt / 72
  lab_bp   <- lab_in / tw * span
  gap_bp   <- 0.07  / tw * span          # visible clearance, inches
  free_l   <- tx_lo + lab_bp + gap_bp < win[1]
  free_r   <- tx_hi - lab_bp - gap_bp > win[2]
  lab_side <- if (free_l) "left" else if (free_r) "right" else
              if ((win[1] - tx_lo) >= (tx_hi - win[2])) "left" else "right"
  cat(sprintf("[panel e]   %-7s gene name -> %-5s (free L/R: %s/%s)\n",
              L$gene, lab_side, free_l, free_r))

  trap <- trapezoid_gg(tx_lo, tx_hi, win[1], win[2], fill = LAV, alpha = 1) +
    theme(plot.margin = gutter(side))
  gene <- transcript_zoom_gg(
            ex, tx_lo, tx_hi,
            highlight = win, highlight_fill = LAV, highlight_alpha = 1,
            tx_label = L$gene, exon_h = 0.85, n_arrows = 26,
            label_size = SZ_GENE, label_family = FONT_FAMILY,
            x_ticks = brks, tick_h = 0.55, label_side = lab_side) +
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
build_quadrant <- function(key, side) {
  L  <- LOCI[[key]]
  wd <- region_dir(L)
  stopifnot(dir.exists(wd))
  hi <- read_links(file.path(wd, LK_HIGH), L$wide[1], L$wide[2])
  ct <- read_links(file.path(wd, LK_CTRL), L$wide[1], L$wide[2])
  cat(sprintf("[panel e] %-7s %s:%d-%d  y_max %.4f\n",
              key, L$chrom, L$wide[1], L$wide[2], L$y_max))
  cat(sprintf("[panel e]   high PSI %s | ctrl PSI %s\n",
              paste(sort(hi$psi, decreasing = TRUE), collapse = ", "),
              paste(sort(ct$psi, decreasing = TRUE), collapse = ", ")))

  trk_hi <- make_track(L, "high", RED,  side, hi,
                       caption   = if (side == "right") HIGH_LAB else NULL,
                       caption_n = HIGH_N)
  trk_ct <- make_track(L, "ctrl", GREY, side, ct,
                       caption   = if (side == "right") CTRL_LAB else NULL,
                       caption_n = CTRL_N)
  g <- make_gene_rows(L, side)

  (trk_hi / trk_ct / g$trap / g$gene / g$axis) +
    plot_layout(heights = c(H_TRK + H_AXGAP, H_TRK, H_TRAP, H_GENE, H_AXIS))
}

q_tp53bp1 <- build_quadrant("TP53BP1", "left")
q_pbrm1   <- build_quadrant("PBRM1",   "right")
q_msh6    <- build_quadrant("MSH6",    "left")
q_asxl1   <- build_quadrant("ASXL1",   "right")

top <- (q_tp53bp1 | q_pbrm1)
bot <- (q_msh6    | q_asxl1)

wrap_elements(full =
  (top / bot) + plot_layout(heights = c(STACK_H + H_AXGAP, STACK_H + H_AXGAP)) +
  plot_annotation(
    title = paste0("Tumor suppressor genes and frequently mutated ",
                   "oncogenes are mis-spliced in\n",
                   "samples with high loadings upon GSP #30"),
    theme = theme(plot.title = element_text(
                    size = SIZE_TITLE, face = "bold", hjust = 0.5,
                    family = FONT_FAMILY, margin = margin(0, 0, 3, 0)),
                  plot.margin = margin(0, 0, 0, 0))))
}))
