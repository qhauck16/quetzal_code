# Figure 3, panel a. Sourced by ../compose.R, which supplies the shared
# setup (theme tokens, ROW_H, COL_W, panel_cached).
# Must leave `panel_a` defined.
#
# Rebuild just this panel:   Rscript compose.R a
#
# a. Somatic mutations in NKAP, PRPF8, FAM50A and SLU7 are enriched
#    upon GSP #13. Left: full CDF with the zoom region boxed.
#    Right: the zoom, with each sample's protein change labelled.
#
# Labels are placed by `cdf_place_labels()` below rather than ggrepel:
# 47 labels crowd the 85-100% band, and ggrepel only avoids its OWN
# layer's data, so it cannot see the markers or the step lines. The
# solver knows every obstacle and searches a cost-ordered candidate grid
# (slide out in x AND step in y), spilling overflow into the >100%
# headroom and the sparser 68-85% band. The same function places
# figure 4 panel a.

source("scripts/nkap_variants.R")

A_INFO <- vt_path("NKAP_full_plus_variants")

# ==================================================================
#  Shared "CDF + labelled zoom" label machinery.
#  KEEP IN SYNC with ../../figure4/panel_a/build.R -- the two copies
#  are deliberately identical.
# ==================================================================

# Exact Helvetica metrics. Widths must come from a *cairo* device:
# the figure is written with cairo_pdf, and cairo resolves "Helvetica"
# through fontconfig, so grDevices::pdf()'s built-in AFM metrics are
# ~6% narrow -- enough to let two labels print flush.
cdf_text_metrics <- function(txt, lab_pt, family) {
  dev0 <- grDevices::dev.cur()
  tf   <- tempfile(fileext = ".pdf")
  grDevices::cairo_pdf(tf, width = 4, height = 4)
  on.exit({
    grDevices::dev.off()
    unlink(tf)
    if (dev0 > 1) grDevices::dev.set(dev0)
  }, add = TRUE)
  g <- grid::gpar(fontfamily = family, fontsize = lab_pt)
  list(
    w = unname(vapply(txt, function(s) grid::convertWidth(
          grid::grobWidth(grid::textGrob(s, gp = g)), "pt", TRUE),
          numeric(1))),
    # "Rgp" spans cap height to descender: the full ink box a
    # neighbouring label has to clear.
    h = grid::convertHeight(
          grid::grobHeight(grid::textGrob("Rgp", gp = g)), "pt", TRUE))
}

# geom_step(direction = "hv") as explicit axis-aligned segments, one
# staircase per group. These are obstacles the labels must clear;
# ggrepel is structurally unable to see them.
cdf_step_segments <- function(df, xv, yv, gv) {
  out <- lapply(split(df, df[[gv]]), function(z) {
    z <- z[order(z[[xv]], z[[yv]]), ]
    k <- nrow(z)
    if (k < 2) return(NULL)
    x <- z[[xv]]; y <- z[[yv]]
    rbind(data.frame(x0 = x[-k], y0 = y[-k], x1 = x[-1], y1 = y[-k]),
          data.frame(x0 = x[-1], y0 = y[-k], x1 = x[-1], y1 = y[-1]))
  })
  do.call(rbind, out[!vapply(out, is.null, logical(1))])
}

# ---- geometry primitives, vectorised over the obstacle list -------
# Axis-aligned overlap of one box against many boxes.
.cdf_hit <- function(b, R) {
  if (!NROW(R)) return(logical(0))
  (b[1] < R[, 2]) & (b[2] > R[, 1]) & (b[3] < R[, 4]) & (b[4] > R[, 3])
}
# Liang-Barsky: does one segment cut any of many boxes?
.cdf_seg_boxes <- function(sg, B) {
  if (!NROW(B)) return(logical(0))
  m  <- nrow(B)
  dx <- sg[3] - sg[1]; dy <- sg[4] - sg[2]
  t0 <- rep(0, m); t1 <- rep(1, m); ok <- rep(TRUE, m)
  f <- function(p, q) {
    if (p == 0) { ok <<- ok & (q >= 0); return(invisible(NULL)) }
    r <- q / p
    if (p < 0) t0 <<- pmax(t0, r) else t1 <<- pmin(t1, r)
  }
  f(-dx, sg[1] - B[, 1]); f(dx, B[, 2] - sg[1])
  f(-dy, sg[2] - B[, 3]); f(dy, B[, 4] - sg[2])
  ok & (t0 <= t1)
}
# ... and the transpose: do many segments cut one box?
.cdf_segs_box <- function(S, b) {
  if (!NROW(S)) return(logical(0))
  m  <- nrow(S)
  dx <- S[, 3] - S[, 1]; dy <- S[, 4] - S[, 2]
  t0 <- rep(0, m); t1 <- rep(1, m); ok <- rep(TRUE, m)
  f <- function(p, q) {
    z <- p == 0
    ok[z] <<- ok[z] & (q[z] >= 0)
    r <- q / ifelse(z, 1, p)
    neg <- !z & p < 0; pos <- !z & p > 0
    t0[neg] <<- pmax(t0[neg], r[neg]); t1[pos] <<- pmin(t1[pos], r[pos])
  }
  f(-dx, S[, 1] - b[1]); f(dx, b[2] - S[, 1])
  f(-dy, S[, 2] - b[3]); f(dy, b[4] - S[, 2])
  ok & (t0 <= t1)
}
# Proper crossing of one segment against many.
.cdf_cross <- function(sg, S) {
  if (!NROW(S)) return(logical(0))
  cw <- function(ax, ay, bx, by, cx, cy)
    (cy - ay) * (bx - ax) - (by - ay) * (cx - ax)
  d1 <- cw(S[, 1], S[, 2], S[, 3], S[, 4], sg[1], sg[2])
  d2 <- cw(S[, 1], S[, 2], S[, 3], S[, 4], sg[3], sg[4])
  d3 <- cw(sg[1], sg[2], sg[3], sg[4], S[, 1], S[, 2])
  d4 <- cw(sg[1], sg[2], sg[3], sg[4], S[, 3], S[, 4])
  (d1 * d2 < 0) & (d3 * d4 < 0)
}

# ------------------------------------------------------------------
# The placer.
#
#   labs  x, y, text, pref (-1 left / +1 right, a SOFT preference)
#   pts   every marker drawn in the panel
#   segs  x0,y0,x1,y1 from cdf_step_segments()
#
# Returns labs with lx/ly (text anchor = the box's inner edge), hj,
# and sx,sy,ex,ey for the leader (NA when the label is adjacent).
#
# Candidates are (side, x-offset, y-offset) triples scored by
#   w_off * extra x-offset + w_dy * |y-offset| + w_flip * wrong side
# and tried in increasing cost order, so a label only leaves its
# point's side when every near slot is genuinely taken. Hardest
# labels (most neighbours in y, then longest string) go first.
# ------------------------------------------------------------------
cdf_place_labels <- function(labs, pts, segs, xlim, ylim, w_in, h_in,
                             lab_pt, marker_pt, family,
                             pad_x = 0.95, pad_y = 0.40, lw_pt = 0.7,
                             mr_extra = 0.85, lead_pad = 1.3,
                             off_pt = c(2.2, 3.6, 5.2, 7, 9.2, 12, 15.5,
                                        19.5, 24, 29, 35, 42, 50, 59, 69),
                             voff_pt = c(1.4, 2.6, 4.2, 6.2, 8.6, 11.5,
                                         15, 19, 24, 30, 37),
                             dy_pt = 2.6, dy_n = 18,
                             w_off = 1.0, w_dy = 1.0, w_flip = 8,
                             w_cen = 5, w_cross = 8, w_graze = 30,
                             lead_max_pt = 40,
                             lead_min_pt = 5.5, rounds = 26,
                             w_far = 20, polish = 80, shuffle = 8,
                             tag = "", verbose = TRUE) {
  ux <- diff(xlim) / (w_in * 72)     # data units per pt, x
  uy <- diff(ylim) / (h_in * 72)     # data units per pt, y

  m  <- cdf_text_metrics(labs$text, lab_pt, family)
  n  <- nrow(labs)
  # The box has to be as tall as the FONT, not as the measured ink.
  # geom_text(vjust = 0.5) centres each string on its OWN ink box, so
  # "R366Q" (caps only) and "N1881_V1883dup" (underscore below the
  # baseline) sit differently inside the same nominal row; ink-sized
  # boxes let one string's underscore print through its neighbour.
  # 0.94 em spans cap height + descender with a little slack.
  th <- max(m$h, lab_pt * 0.92)
  hw <- (m$w / 2 + pad_x) * ux
  hh <- (th / 2 + pad_y) * uy
  mr <- marker_pt * .pt / 2 + mr_extra  # marker radius, pt

  PR <- cbind(pts$x - mr * ux, pts$x + mr * ux,
              pts$y - mr * uy, pts$y + mr * uy)
  SR <- cbind(pmin(segs$x0, segs$x1) - lw_pt * ux,
              pmax(segs$x0, segs$x1) + lw_pt * ux,
              pmin(segs$y0, segs$y1) - lw_pt * uy,
              pmax(segs$y0, segs$y1) + lw_pt * uy)

  # Candidate grid. Three anchors, because two are not enough: a flank
  # label needs a free slot the width of the whole string, whereas a
  # CENTRED one (mode 0, sitting squarely above or below its own point)
  # only needs the height -- and along the NKAP diagonal and in the
  # >100% headroom that is often the only slot left. `o` is the
  # clearance from the point to the box's NEAR edge; `d` is the
  # y-offset for the flanks and the up/down sign for the centred ones.
  dys <- c(0, as.vector(rbind(seq_len(dy_n), -seq_len(dy_n)))) * dy_pt
  G   <- rbind(
    cbind(expand.grid(o = off_pt,  d = dys),        mode = -1),
    cbind(expand.grid(o = off_pt,  d = dys),        mode =  1),
    cbind(expand.grid(o = voff_pt, d = c(1, -1)),   mode =  0))
  # Displacement of the label from its point, in points, plus a
  # QUADRATIC far term: without it the solver happily parks one label
  # 50 pt away to save two labels 25 pt each -- the same linear cost,
  # but one leader flying across the panel instead of two short ones.
  # `gvis` is the part that is VISIBLE (leader length) and is what the
  # polish loop minimises; the side/centred preferences only steer the
  # per-label search and must not enter the objective, or the solver
  # trades real distance for a preference nobody can see.
  gdis <- ifelse(G$mode == 0,
                 w_dy * (G$o + th / 2),
                 w_off * (G$o - off_pt[1]) + w_dy * abs(G$d))
  gvis  <- gdis + w_far * (gdis / 25)^2
  gcost <- gvis + w_cen * (G$mode == 0)

  # Difficulty: linear text demand within a few text-rows in y,
  # weighted by how little of it the panel can hold. Long strings in a
  # busy band must choose before short ones in quiet ones, or they
  # find nothing left. Failures are promoted to the front and the
  # whole thing is re-run (see the round loop below) -- a static
  # ordering is never right for every label at once.
  dem <- vapply(seq_len(n), function(i)
           sum(m$w * pmax(0, 1 - abs(labs$y - labs$y[i]) / (5 * hh))),
           numeric(1))

  # Candidate k for label i, in data units. `ax/ay` is where the
  # leader lands: the box's inner edge for a flank label, the middle
  # of its near edge for a centred one.
  geom_k <- function(i, k) {
    md <- G$mode[k]
    if (md != 0) {
      ax <- labs$x[i] + md * G$o[k] * ux
      ay <- labs$y[i] + G$d[k] * uy
      cx <- ax + md * hw[i]; cy <- ay
      hv <- if (md < 0) 1 else 0
      # The text is anchored pad_x INSIDE the box edge, so the gap
      # sits between the string and its own marker.
      lx <- ax + md * pad_x * ux
    } else {
      cx <- labs$x[i]
      cy <- labs$y[i] + G$d[k] * (G$o[k] * uy + hh)
      ax <- cx; ay <- cy - G$d[k] * hh
      hv <- 0.5; lx <- cx
    }
    dxp <- (ax - labs$x[i]) / ux
    dyp <- (ay - labs$y[i]) / uy
    list(b = c(cx - hw[i], cx + hw[i], cy - hh, cy + hh),
         ax = ax, ay = ay, cy = cy, hv = hv, lx = lx,
         dxp = dxp, dyp = dyp, len = sqrt(dxp^2 + dyp^2))
  }

  # Hard constraints vs soft penalties.
  #
  # HARD at level 0: the box must be inside the panel, off every
  # marker, off every step line, off every other box, and no leader
  # may run clean THROUGH another string.
  #
  # SOFT: a leader crossing another leader (w_cross) and a leader
  # passing within lead_pad of a string without entering it
  # (w_graze -- a hairline 1 pt past the final "I" of "R290I" prints
  # as "R290Il"). As hard constraints these take slots away from the
  # greedy too early and force later labels onto the relaxation ladder.
  #
  # The relaxation ladder is the last resort, one step at a time:
  # level 1 lets a leader cross a string, level 2 lets a step line
  # cross it. Level 2 exists only so a label is never silently dropped.
  place_once <- function(ord, relax) {
    lx <- ly <- sx <- sy <- ex <- ey <- rep(NA_real_, n)
    hj   <- rep(0, n)
    BOX  <- matrix(NA_real_, 0, 4)
    BOXL <- matrix(NA_real_, 0, 4)   # ... inflated by lead_pad in x
    LEAD <- matrix(NA_real_, 0, 4)
    bad  <- integer(0)
    tot  <- 0
    for (i in ord) {
      rl  <- relax[i]
      cst <- gcost + w_flip * (labs$pref[i] != 0 & G$mode != 0 &
                                G$mode != labs$pref[i])
      kk  <- order(cst)
      pick <- NA_integer_; pick_pen <- 0; best_c <- Inf
      for (k in kk) {
        if (cst[k] >= best_c) break     # cst is sorted: nothing better left
        g <- geom_k(i, k)
        b <- g$b
        if (b[1] < xlim[1] || b[2] > xlim[2] ||
            b[3] < ylim[1] || b[4] > ylim[2]) next
        if (g$len > lead_max_pt)   next
        if (any(.cdf_hit(b, PR)))  next
        if (rl < 2 && any(.cdf_hit(b, SR)))  next
        if (any(.cdf_hit(b, BOX))) next
        sg <- c(labs$x[i], labs$y[i], g$ax, g$ay)
        if (rl < 1 && (any(.cdf_seg_boxes(sg, BOX)) ||
                       any(.cdf_segs_box(LEAD, b)))) next
        bl  <- b + c(-lead_pad * ux, lead_pad * ux,
                     -lead_pad * uy, lead_pad * uy)
        pen <- w_graze * (any(.cdf_seg_boxes(sg, BOXL)) ||
                          any(.cdf_segs_box(LEAD, bl))) +
               w_cross * any(.cdf_cross(sg, LEAD))
        if (cst[k] + pen < best_c) {
          best_c <- cst[k] + pen; pick <- k; pick_pen <- pen
        }
        if (pen == 0) break
      }
      if (is.na(pick)) { bad <- c(bad, i); next }
      tot <- tot + gvis[pick] + pick_pen
      g <- geom_k(i, pick)
      lx[i] <- g$lx; ly[i] <- g$cy; hj[i] <- g$hv
      BOX   <- rbind(BOX, g$b)
      BOXL  <- rbind(BOXL, g$b + c(-lead_pad * ux, lead_pad * ux,
                                   -lead_pad * uy, lead_pad * uy))
      if (g$len > mr + lead_min_pt) {
        sx[i] <- labs$x[i] + g$dxp / g$len * mr * ux
        sy[i] <- labs$y[i] + g$dyp / g$len * mr * uy
        ex[i] <- g$ax - g$dxp / g$len * 0.5 * ux
        ey[i] <- g$ay - g$dyp / g$len * 0.5 * uy
        LEAD  <- rbind(LEAD, c(sx[i], sy[i], ex[i], ey[i]))
      }
    }
    list(lx = lx, ly = ly, hj = hj, sx = sx, sy = sy, ex = ex, ey = ey,
         bad = bad, tot = tot + 1e6 * length(bad))
  }

  # Greedy placement is order-dependent, and "hardest first" is not
  # always the best order -- a label that squeezes into a near slot
  # can push a neighbour thirty points away, which reads as a long
  # leader flying across the panel. So: promote failures to the front
  # until nothing fails (the relax ladder backs that up), then keep
  # shuffling the order and KEEP THE CHEAPEST layout seen, cost being
  # the summed displacement of every label. That is what turns "all
  # 47 placed" into "all 47 placed next to their own points".
  ord   <- order(-dem, -m$w)
  relax <- integer(n)
  best  <- NULL
  keep  <- function(res) if (is.null(best) || res$tot < best$tot) best <<- res
  for (lvl in 0:4) {
    for (r in seq_len(rounds)) {
      res <- place_once(ord, relax)
      keep(res)
      if (!length(res$bad)) break
      ord <- c(res$bad, setdiff(ord, res$bad))
    }
    if (!length(best$bad)) break
    # Escalate the offenders by ONE step each, not to the loop's
    # level: a label that only started failing because its neighbours
    # were relaxed must still try level 1 before it is allowed to
    # print a leader through someone's string.
    relax[best$bad] <- relax[best$bad] + 1
  }
  if (!length(best$bad) && polish > 0) {
    # Deterministic: own RNG stream, restored on the way out, so the
    # panel is byte-reproducible and nothing else in compose.R moves.
    old <- if (exists(".Random.seed", .GlobalEnv))
             get(".Random.seed", .GlobalEnv) else NULL
    set.seed(20260923)
    base_ord <- order(-dem, -m$w)
    for (r in seq_len(polish)) {
      o   <- order(order(base_ord) + runif(n, 0, shuffle))
      res <- place_once(o, relax)
      if (!length(res$bad)) keep(res)
    }
    if (is.null(old)) rm(".Random.seed", envir = .GlobalEnv) else
      assign(".Random.seed", old, envir = .GlobalEnv)
  }
  if (length(best$bad))
    warning(sprintf("cdf_place_labels%s: no clear slot for %s", tag,
                    paste(labs$text[best$bad], collapse = ", ")),
            call. = FALSE)

  # Audit the layout that is actually going to be drawn, rather than
  # trusting the search. Anything printed here is a visible defect.
  ok <- which(!seq_len(n) %in% best$bad)
  # lx is the TEXT anchor, pad_x inside the box edge (geom_k).
  bax <- best$lx + ifelse(best$hj == 1, pad_x * ux,
                   ifelse(best$hj == 0, -pad_x * ux, 0))
  BB <- cbind(ifelse(best$hj == 1, bax - 2 * hw,
               ifelse(best$hj == 0, bax, bax - hw)),
              ifelse(best$hj == 1, bax,
               ifelse(best$hj == 0, bax + 2 * hw, bax + hw)),
              best$ly - hh, best$ly + hh)
  LL   <- cbind(best$sx, best$sy, best$ex, best$ey)
  lown <- which(!is.na(best$sx))          # which label each leader is for
  LL   <- LL[lown, , drop = FALSE]
  viol <- c(text = 0, marker = 0, step = 0, leader = 0, graze = 0)
  for (i in ok) {
    b  <- BB[i, ]
    # A label's OWN leader stops 0.5 pt short of its box and would
    # otherwise register against the padded box every single time.
    others <- which(lown != i)
    bl <- b + c(-lead_pad * ux, lead_pad * ux,
                -lead_pad * uy, lead_pad * uy)
    viol["text"]   <- viol["text"] +
      sum(.cdf_hit(b, BB[setdiff(ok, i), , drop = FALSE]))
    viol["marker"] <- viol["marker"] + any(.cdf_hit(b, PR))
    viol["step"]   <- viol["step"]   + any(.cdf_hit(b, SR))
    viol["leader"] <- viol["leader"] +
      any(.cdf_segs_box(LL[others, , drop = FALSE], b))
    viol["graze"]  <- viol["graze"] +
      any(.cdf_segs_box(LL[others, , drop = FALSE], bl))
  }
  viol["text"]  <- viol["text"] / 2
  viol["graze"] <- viol["graze"] - viol["leader"]
  if (any(viol[c("text", "marker", "step", "leader")] > 0))
    warning(sprintf(paste("cdf_place_labels%s: %g label-label,",
                          "%g label-marker, %g label-stepline,",
                          "%g label-leader collisions remain"), tag,
                    viol["text"], viol["marker"], viol["step"],
                    viol["leader"]), call. = FALSE)
  if (verbose)
    cat(sprintf(paste("[label solver%s] %d labels, %d leaders,",
                      "%d grazed, %d unplaced%s\n"),
                tag, n, sum(!is.na(best$sx)),
                unname(viol["graze"]), length(best$bad),
                if (length(best$bad))
                  paste0(": ", paste(labs$text[best$bad], collapse = ", "))
                else ""))

  labs$lx <- best$lx; labs$ly <- best$ly; labs$hj <- best$hj
  labs$sx <- best$sx; labs$sy <- best$sy
  labs$ex <- best$ex; labs$ey <- best$ey
  labs$placed <- !seq_len(n) %in% best$bad
  labs
}

# Symmetric gutters, in POINTS, around the point cloud. Solving needs
# xlim and xlim needs to know how far labels reach, so the circle is
# cut by fixing the gutters up front and letting long labels overhang
# the cloud (at a y the solver has proved is clear) instead.
cdf_xlim <- function(x, w_in, gut_l_pt, gut_r_pt) {
  wpt  <- w_in * 72
  span <- diff(range(x)) / (1 - (gut_l_pt + gut_r_pt) / wpt)
  c(min(x) - span * gut_l_pt / wpt, max(x) + span * gut_r_pt / wpt)
}

# ==================================================================

panel_a <- panel_cached("a", local({
  d <- fread(A_INFO, data.table = FALSE)
  stopifnot("k13" %in% names(d), "logk13" %in% names(d))
  d <- classify_variants(d)

  a_cum <- d %>%
    filter(!is.na(logk13)) %>%
    mutate(variant_class = factor(variant_class, levels = VARIANT_LEVELS)) %>%
    arrange(logk13) %>%
    group_by(variant_class) %>%
    mutate(cumulative = row_number() / n() * 100) %>%
    ungroup()

  A_ZOOM_FROM <- -4        # the zoom panel's left edge, boxed on the CDF
  A_YLIM      <- c(0, 108) # the CDF's own range

  # Where the two zoom-callout lines END, in the CDF's y units. The
  # zoom panel is taller than the CDF panel (the CDF gives a third of
  # its own cell to the legend below it), so both ends sit OUTSIDE
  # A_YLIM -- the lines splay to reach the zoom's corners. Measured with
  # ../../scripts/measure_zoom_connectors.py (box corners -> zoom panel's
  # left axis line); re-run it after any layout change to this row.
  A_CALL_DX   <- 1.5268
  A_CALL_TOP  <- 107.450
  A_CALL_BOT  <- -94.369

  # ---- left: full CDF --------------------------------------------
  panel_a_cdf <- ggplot(a_cum, aes(x = logk13, y = cumulative,
                                    color = variant_class)) +
    geom_step(linewidth = 0.55) +
    annotate("rect", xmin = A_ZOOM_FROM, xmax = max(a_cum$logk13) + 0.25,
             ymin = 0, ymax = 102, fill = NA, color = "black",
             linewidth = 0.4) +
    scale_color_manual(values = VARIANT_COL, drop = FALSE, name = NULL,
                       labels = function(x) gene_pm(x)) +  # italic genes (plotmath; see theme.R)
    scale_y_continuous(labels = function(x) paste0(x, "%"),
                        breaks = c(0, 25, 50, 75, 100)) +
    # xlim is pinned as well as ylim because the callout segments added
    # below reach BEYOND the data on both axes, and a layer's data
    # drives the panel range -- without this the panel stretches to fit
    # them and drags the axis labels with it.
    coord_cartesian(xlim = range(a_cum$logk13), ylim = A_YLIM,
                     expand = FALSE, clip = "off") +
    labs(x = "log(GSP #13)", y = "Cumulative percentage") +
    theme_publication() +
    theme(legend.position = "none",
          axis.title  = element_text(size = SIZE_BASE),
          plot.margin = margin(4, 6, 4, 4, "pt"))

  # ---- right: the zoom, labelled ---------------------------------
  a_zoom <- filter(a_cum, logk13 >= A_ZOOM_FROM)
  a_lab  <- filter(a_zoom, nzchar(proteinChange))

  # The zoom gets its OWN y range, not the CDF's. Every class is
  # already past 60% by log = -4, so sharing c(0, 108) would squeeze all
  # 47 labelled points into the top third. The headroom above 100% is
  # where the solver parks the overflow from the crowded 85-100% band;
  # 114 is the smallest top that holds four text rows there.
  A_ZLIM <- c(floor(min(a_zoom$cumulative)) - 2, 110)

  # Drawn size of the zoom's PANEL (not its cell), in inches. Every
  # point->data conversion in the solver hangs off these two numbers,
  # so they are DERIVED from the layout rather than typed in, and the
  # two constants below are the only measured quantities:
  #
  #   A_FURN_W  the cell width that is NOT panel -- the CDF's y-axis
  #             title + tick labels + the three plot margins.
  #   A_FURN_H  the row height that is NOT panel -- the plot_annotation
  #             title band, the zoom's x tick labels and axis title.
  #
  # Both are measured on the render (zoom x-axis line and left border
  # against the cell rectangle). Re-measure after any change to
  # ROW_H[["ab"]], COL_W$ab, A_SPLIT or the axis furniture. The left
  # cell's CDF width is pinned by the legend's longest entry, so
  # widening COL_W$ab[1] goes almost entirely to the zoom.
  A_FURN_W  <- 0.9475
  A_FURN_H  <- 0.8075
  A_SPLIT   <- c(1, 1.62)   # CDF panel : zoom panel (patchwork widths)
  A_CELL_W  <- PAGE_W * (COL_W$ab[1] / sum(COL_W$ab))
  A_ZOOM_W_IN <- (A_CELL_W - A_FURN_W) * A_SPLIT[2] / sum(A_SPLIT)
  A_ZOOM_H_IN <- ROW_H[["ab"]] - A_FURN_H

  A_LAB_PT  <- SIZE_BASE * 0.62
  A_PT_SIZE <- 0.9
  # Legend text size, as a fraction of SIZE_BASE. The long legend entry
  # has to fit the left column; A_SPLIT buys the room rather than
  # shrinking the type.
  A_LEG_F <- 1

  # Gutters in points. The long strings overhang the point cloud
  # rather than pushing the axis out.
  # 48/20 minimises leader count and length in a sweep over the solver.
  # Raising `polish` above ~80 rounds does not help: the solver minimises
  # total displacement, not leader count.
  a_lim <- cdf_xlim(a_zoom$logk13, A_ZOOM_W_IN, 48, 20)

  # Side preference, SOFT (w_flip) and per class. Only NKAP asks for
  # one: its 20 points are the long lower diagonal with room on both
  # sides, and alternating is what keeps it reading as a DIAGONAL
  # rather than a column. The other three choose freely: PRPF8's stack
  # sits ~10 pt from the panel's left edge, so pinning it to a flank
  # would cost long vertical travel and a leader across the panel.
  A_SIDE <- c("NKAP variant"   = 9,    # 9 = alternate
              "PRPF8 variant"  = 0,    # 0 = let the solver choose
              "FAM50A variant" = 0,
              "SLU7 variant"   = 0)
  a_lab <- a_lab %>%
    arrange(variant_class, cumulative, logk13) %>%
    group_by(variant_class) %>%
    mutate(.pin = A_SIDE[as.character(variant_class)],
           pref = ifelse(.pin != 9, .pin,
                         ifelse(row_number() %% 2 == 1, -1, 1))) %>%
    ungroup() %>%
    select(-.pin) %>%
    mutate(x = logk13, y = cumulative, text = proteinChange)

  a_segs <- cdf_step_segments(as.data.frame(a_zoom),
                              "logk13", "cumulative", "variant_class")

  a_lab <- cdf_place_labels(
    as.data.frame(a_lab),
    pts    = data.frame(x = a_zoom$logk13, y = a_zoom$cumulative),
    segs   = a_segs,
    xlim   = a_lim, ylim = A_ZLIM,
    w_in   = A_ZOOM_W_IN, h_in = A_ZOOM_H_IN,
    lab_pt = A_LAB_PT, marker_pt = A_PT_SIZE, family = FONT_FAMILY,
    lead_max_pt = 40,
    tag = " fig3a")

  a_lead <- filter(a_lab, !is.na(sx))

  panel_a_zoom <- ggplot(a_zoom, aes(x = logk13, y = cumulative,
                                      color = variant_class)) +
    geom_step(linewidth = 0.35, alpha = 0.5) +
    geom_point(size = A_PT_SIZE) +
    geom_segment(data = a_lead,
                 aes(x = sx, y = sy, xend = ex, yend = ey),
                 inherit.aes = FALSE, linewidth = 0.12,
                 color = "grey55", alpha = 0.9) +
    geom_text(data = a_lab,
              aes(x = lx, y = ly, label = text, hjust = hj,
                  color = variant_class),
              inherit.aes = FALSE, vjust = 0.5,
              size = A_LAB_PT / .pt, family = FONT_FAMILY,
              show.legend = FALSE) +
    scale_color_manual(values = VARIANT_COL, drop = FALSE, name = NULL,
                       labels = function(x) gene_pm(x)) +  # italic genes (plotmath; see theme.R)
    # Ticks only over the DATA. The x limits carry a label gutter on
    # each side (cdf_xlim), and a tick at -6 in a region the zoom does
    # not contain reads as "the zoom starts at -6", contradicting the
    # box drawn on the CDF at -4.
    scale_x_continuous(breaks = seq(ceiling(A_ZOOM_FROM), 0, 1)) +
    scale_y_continuous(labels = function(x) paste0(x, "%"),
                        breaks = c(0, 25, 50, 75, 100)) +
    coord_cartesian(xlim = a_lim, ylim = A_ZLIM, expand = FALSE,
                     clip = "off") +
    labs(x = "log(GSP #13)", y = NULL) +
    theme_publication() +
    theme(legend.position = "none",
          axis.title  = element_text(size = SIZE_BASE),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin = margin(6, 6, 4, 4, "pt"))

  # Shared legend, drawn once under the CDF.
  a_legend <- with_cairo_metrics(
  cowplot::get_plot_component(
    panel_a_cdf +
      theme(legend.position   = "bottom",
            legend.key.size   = unit(2.3, "mm"),
            legend.text       = element_text(size = SIZE_BASE * A_LEG_F, family = FONT_FAMILY),
            legend.box.spacing = unit(0, "mm")) +
      guides(color = guide_legend(ncol = 1,
                                   override.aes = list(linewidth = 1.2))),
    "guide-box-bottom", return_all = TRUE))

  # ---- zoom callout ----------------------------------------------
  # Lines from the CDF's zoom box to the zoom panel mark it as a
  # zoom-in. They leave the box's two right corners and fan out towards
  # the zoom panel's two left corners, which sit outside the CDF's own
  # y range because the zoom is full-cell height while the CDF gives a
  # third of its cell to the legend.
  # The x reach is the physical gap pushed through the CDF panel's
  # drawn width, so it tracks COL_W$ab rather than being typed in.
  A_CDF_W_IN <- (A_CELL_W - A_FURN_W) * A_SPLIT[1] / sum(A_SPLIT)
  A_GAP_IN   <- 0.139             # CDF's 6 pt right + zoom's 4 pt left
  a_xr       <- range(a_cum$logk13)
  a_call <- data.frame(
    x    = max(a_cum$logk13) + 0.25,
    y    = c(102, 0),
    xend = max(a_cum$logk13) + 0.25 + A_CALL_DX,   # measured, see A_CALL_*
    yend = c(A_CALL_TOP, A_CALL_BOT))
  panel_a_cdf <- panel_a_cdf +
    geom_segment(data = a_call, aes(x = x, y = y, xend = xend, yend = yend),
                 inherit.aes = FALSE, linewidth = 0.35, colour = "black")

  a_left <- (panel_a_cdf / wrap_elements(full = a_legend)) +
    plot_layout(heights = c(1, 0.5))

  wrap_elements(full =
    ((a_left | panel_a_zoom) + plot_layout(widths = c(1, 1.62))) +
      plot_annotation(
        title = gene_md(paste0("Somatic mutations in NKAP, PRPF8, FAM50A\n",
                       "and SLU7 are enriched upon GSP #13")),
        # Zero top margin on both: patchwork stacks a plot_annotation
        # title on its own plot margin, which the labs()-titled panels
        # do not have.
        theme = theme(plot.title = ggtext::element_markdown(halign = 0.5, lineheight = 0.9,
          size = SIZE_TITLE, face = "bold", hjust = 0.5,
          # 2 pt top: panel a is titled via plot_annotation() and panel
          # b via labs(), which patchwork spaces differently; this levels
          # the two titles.
          family = FONT_FAMILY, margin = margin(2, 0, 2, 0)),
          plot.margin = margin(0, 0, 0, 0))))
}))
