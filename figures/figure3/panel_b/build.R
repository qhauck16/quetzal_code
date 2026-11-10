# Figure 3, panel b. Sourced by ../compose.R, which supplies the shared
# setup (theme tokens, FLASHIER, PAGE_W, ROW_H, COL_W, panel_cached).
# Must leave `panel_b` defined.
#
# Rebuild just this panel:   Rscript compose.R b
#
# b. "Mutations cluster in specific regions of NKAP and PRPF8".
#    Two stacked protein cartoons (NKAP over PRPF8), x = amino-acid
#    position on a BROKEN axis, lollipops (no stems) stacked on the
#    body's top edge, one circle per carrier sample, coloured by TCGA
#    project.
#
# The reference was drawn by hand in Illustrator. Everything here is
# derived: stack membership from the TSV, label text by prefix-collapse,
# x positions by a piecewise aa->inch map plus a collision sweep, and
# every font/circle size from one auto-fit scale factor (see B_FIT).
#
# See panel_b/build.md for the top-60 rule and the broken-axis rationale.

source("scripts/nkap_variants.R")

B_INFO <- vt_path("NKAP_full_plus_variants")

panel_b <- panel_cached("b", local({

  # ================================================================
  # 0. Typography / geometry tokens.
  # ================================================================
  # All of these are BASE values; the auto-fit in section 5 scales them
  # by one factor so the two cartoons exactly fill the cell's height.
  # 62 -> 68. Swept at the CURRENT cell width (a|b = 0.58/0.42) with
  # PRPF8's new segment fractions: 62 -> S 1.002, 68 -> 1.074,
  # 72 -> 1.066, 76 -> 1.060. 62 was the optimum at the OLD, wider cell
  # and went stale with it. The mechanism: a shallower angle spends the
  # label's length on WIDTH, and width is what binds here (the auto-fit
  # reports xover = 0 at the binding S, i.e. it stops because any more
  # type would overrun, not because it ran out of height). Steepening
  # trades that width for height, and there WAS 0.22 in of unused
  # height to trade into -- at 68 h_req 2.563 against h_avail 2.576, so
  # the panel now fills its cell and the dead band is gone.
  # RE-SWEEP after any change to COL_W$ab or the segment fractions.
  B_ANGLE     <- 68                    # rotated stack labels, degrees
  B_LAB_PT0   <- SIZE_BASE * 0.58      # mutation labels
  B_TICK_PT0  <- SIZE_BASE * 0.55      # aa tick numbers
  B_TICK_LEN  <- 0.20                  # tick mark length below the body, in circle diameters
  B_TICK_GAP  <- 0.08                  # tick end -> top of its number
  B_GENE_PT0  <- SIZE_BASE * 0.85      # italic gene name
  B_DOM_PT0   <- SIZE_BASE * 0.62      # domain names inside the boxes
  B_NOTE_PT0  <- SIZE_BASE * 0.60      # the two prose annotations
  B_CIRC_D0   <- 0.076                 # lollipop diameter, inches
  B_LAB_MAXW0 <- 0.62                  # max label LINE length, inches
  B_LAB_LH    <- 0.88                  # drawn lineheight inside a label
  # Legend type is NOT auto-fit (it does not compete for the cartoons'
  # height -- it lives in dead space ABOVE NKAP's label staircase), so
  # its size is solved in section 7 against the width that staircase
  # actually leaves free. These are only the ceiling, the swatch size in
  # lines of that solved type, and the block's shape.
  B_LEG_PT_MAX <- SIZE_LEGEND * 0.80
  B_LEG_D_LINE <- 1.25                 # swatch diameter, in legend lines
  B_LEG_NCOL   <- 4                    # 16 projects as 4 x 4, not 1 x 16
  B_LEG_NROW   <- 4
  B_LEG_ROWF   <- 1.34                 # row pitch, in legend lines
  B_LEG_GAPC   <- 0.35                 # column gutter, in swatch diameters
  B_LEG_TOP    <- 0.02                 # inset from the drawing's top edge
  B_LEG_CLR    <- 0.10                 # gutter to the nearest label ink

  # The two leaders from the below-body circles down to their caption.
  # Both are the SAME length and the SAME slope magnitude, mirrored
  # about the caption's centre -- neither is parallel to a page edge and
  # the pair reads as one treatment. Length is in circle diameters, so
  # it scales with the auto-fit like everything else.
  B_STUB_ANG <- 30                     # degrees off VERTICAL
  B_STUB_LEN <- 1.15                   # leader length, in circle diameters
  B_STUB_PAD <- 0.30                   # leader-to-caption gap, same units

  # The two prose annotations. Held here because their DRAWN width is
  # part of the layout solve (section 5), not just of the drawing.
  B_CAP1_TXT <- "Mutations in the top 60\nsamples on GSP #13"
  B_CAP2_TXT <- "TCGA samples carrying a plotted\nmutation with GSP #13 < 1e-4"
  # ONE line now: the footnote used to be broken over two because it sat
  # in a 0.62 in legend COLUMN. Under the 4 x 4 block it has ~2.9 in.
  B_FOOT_TXT <- "*Called via WXS or RNA-seq"

  B_SIN <- sin(B_ANGLE * pi / 180)
  B_COS <- cos(B_ANGLE * pi / 180)

  # TCGA project colours, as published in the reference's legend.
  B_PROJ_COL <- c(
    "TCGA-UCEC" = "#9fd4fa", "TCGA-BLCA" = "#bce4a1", "TCGA-BRCA" = "#ffc98e",
    "TCGA-LIHC" = "#ff8e92", "TCGA-ACC"  = "#ffbdd5", "TCGA-READ" = "#b8c7d0",
    "TCGA-HNSC" = "#b296c8", "TCGA-LUSC" = "#7fb39b", "TCGA-KIRC" = "#8c999c",
    "TCGA-SKCM" = "#ffff7f", "TCGA-LUAD" = "#af9b89", "TCGA-CESC" = "#f8f8f8",
    "TCGA-UCS"  = "#e2cdb6", "TCGA-COAD" = "#ff7fff", "TCGA-STAD" = "#7fffff",
    "TCGA-PRAD" = "#b47f83")
  B_STROKE <- "#7f7f7f"

  # Cell geometry. The panel is drawn in INCHES: x and y data units are
  # both inches, and coord_cartesian pins the view to exactly the drawn
  # size, so a circle of radius r data-units is r inches on the page.
  #
  # The cell-minus-margins arithmetic below is only an ESTIMATE of the
  # drawn panel -- patchwork adds its own margins around a
  # wrap_elements() and the title band is not exactly 2 * lineheight.
  # B_CAL_W / B_CAL_H correct it BY MEASUREMENT off the 400 dpi render.
  # Re-measure with the [panel b] cal line printed below whenever the
  # row height, the column split or the title changes; they ARE only
  # ratios, but patchwork's own margins are absolute, so a ratio fitted
  # at one cell width is wrong at another. (Both went stale when the
  # a|b split moved 0.52/0.48: width was out by +0.65%, height by
  # -1.4%.) Method: crop the drawing area, find the rows carrying a
  # >1000 px horizontal rule -- those are the two protein bodies' top
  # and bottom edges -- and check that
  #   the body's drawn width       == B_BODY_W
  #   NKAP top - PRPF8 top (drawn) == yN_body_t - yP_body_t
  # both printed by the cal line. Scale the factor by the ratio and
  # repeat. If they drift, lollipops stop touching and every derived
  # text offset goes stale. Same practice as figure 2's title
  # alignment -- see ../figure2/compose.md.
  B_CELL_W <- PAGE_W * COL_W$ab[2] / sum(COL_W$ab)
  B_MAR    <- c(t = 0, r = 2, b = 2, l = 2)          # plot.margin, pt
  B_TITLE_H <- (2 * SIZE_TITLE * 1.15 + 2) / 72       # 2 title lines + gap
  # Re-measured 2026-09-23 at the a|b split 0.58/0.42 (cell 3.57 in).
  B_CAL_W  <- 0.9852
  B_CAL_H  <- 0.9882
  B_W_TOT  <- (B_CELL_W - (B_MAR[["r"]] + B_MAR[["l"]]) / 72) * B_CAL_W
  B_H_TOT  <- (ROW_H[["ab"]] - B_TITLE_H -
               (B_MAR[["t"]] + B_MAR[["b"]]) / 72) * B_CAL_H

  # The legend no longer takes a column off the right-hand edge: a
  # single 1 x 16 column starved both cartoons of 0.67 in (17% of the
  # cell) and forced the broken axis to spend 30% of NKAP's width on a
  # region with no mutations. It is now a 4 x 4 block in the dead space
  # above NKAP's label staircase, so BOTH proteins get the full width.
  B_BODY_W  <- B_W_TOT                    # protein cartoon width
  # Nothing -- circle, tick number or rotated label -- may reach past
  # this; it is the drawing area's own right edge less a hairline.
  B_X_RIGHT <- B_W_TOT - 0.01

  # ================================================================
  # 1. Text metrics. Every offset below comes from a real measurement
  #    of the drawn string, never from a hardcoded data-unit guess.
  # ================================================================
  .dev <- tempfile(fileext = ".pdf")
  grDevices::cairo_pdf(.dev, width = 8, height = 8)
  on.exit({ grDevices::dev.off(); unlink(.dev) }, add = TRUE)

  # cairo_pdf reports string widths rounded to WHOLE POINTS, which at
  # 4 pt type quantises every measurement to ~15% -- enough to make the
  # auto-fit below land on a wrap discontinuity. Measure once at 100 pt
  # (where the rounding is 1 part in 1000) and scale: font metrics are
  # linear in size. Cached, because the bisection re-measures ~25x.
  .wc <- new.env(parent = emptyenv())
  str_w <- function(s, pt, face = "plain") {
    vapply(seq_along(s), function(i) {
      key <- paste0(face, "\r", s[i])
      v <- .wc[[key]]
      if (is.null(v)) {
        v <- grid::convertWidth(grid::grobWidth(grid::textGrob(
          s[i], gp = grid::gpar(fontsize = 100, fontfamily = FONT_FAMILY,
                                fontface = face))), "in", valueOnly = TRUE)
        assign(key, v, envir = .wc)
      }
      v * pt / 100
    }, numeric(1))
  }

  # ================================================================
  # 2. Data -> stacks.
  # ================================================================
  d <- fread(B_INFO, data.table = FALSE)
  stopifnot(all(c("rail_id", "k13", "project_id", "gene_mutation") %in% names(d)))
  d <- d[order(-d$k13), ]
  d$rank <- seq_len(nrow(d))

  B_WIN <- list(NKAP = c(300, 413), PRPF8 = c(1570, 1900))
  B_TOPN <- 60

  # Long form: one row per (sample, gene, protein change) inside the
  # gene's functional window, across all three mutation sources.
  long <- do.call(rbind, lapply(names(B_WIN), function(g) {
    lo <- B_WIN[[g]][1]; hi <- B_WIN[[g]][2]
    do.call(rbind, lapply(seq_len(nrow(d)), function(i) {
      row <- d[i, ]
      ch  <- gene_changes(row, g)
      if (!length(ch)) return(NULL)
      pos <- suppressWarnings(as.numeric(regmatches(ch, regexpr("\\d+", ch))))
      k   <- which(!is.na(pos) & pos >= lo & pos <= hi)
      if (!length(k)) return(NULL)
      # A circle is starred when the same residue was ALSO called from
      # the WXS or RNA-seq pipeline, not only from the MAF/AlphaMissense
      # column -- that is what the reference's footnote means.
      alt <- c(extract_all_positions(row[["wxs_mutation"]],    g),
               extract_all_positions(row[["rnaseq_mutation"]], g))
      data.frame(rail_id = row$rail_id, rank = row$rank, k13 = row$k13,
                 project = row$project_id, gene = g,
                 change = ch[k], pos = pos[k], star = pos[k] %in% alt,
                 stringsAsFactors = FALSE)
    }))
  }))

  # The panel's circles ABOVE the body are the carriers among the top
  # `B_TOPN` samples by GSP #13. A sample carrying two changes in the
  # same window (866663 = R402Q + L313M, 111261 = L1872P + S1887L)
  # contributes ONE circle: the one at the more recurrent residue. That
  # is what removes L313M / S1887L, which the reference does not draw.
  top <- long[long$rank <= B_TOPN, ]
  rec <- as.data.frame(table(paste(top$gene, top$pos)), stringsAsFactors = FALSE)
  top$nrec <- rec$Freq[match(paste(top$gene, top$pos), rec$Var1)]
  top <- top[order(top$rail_id, top$gene, -top$nrec, top$pos), ]
  top <- top[!duplicated(paste(top$rail_id, top$gene)), ]

  # Stacks: one column per set of changes whose RESIDUE SPANS overlap.
  # "P1880_N1881del", "N1881_V1883dup" and "N1881insKN" all touch
  # residue 1881, so they share a column (as the reference draws them);
  # N371K and Q372H are 1 aa apart but disjoint, so they stay separate.
  # Position numbers are taken only from `<AA-letter><digits>` so the
  # "13" in "R402Sfs*13" is not mistaken for a residue.
  chg_span <- function(s) {
    v <- as.numeric(sub("^[A-Z]", "",
           regmatches(s, gregexpr("[A-Z][0-9]+", s))[[1]]))
    if (!length(v)) c(NA_real_, NA_real_) else range(v)
  }
  spans <- t(vapply(top$change, chg_span, numeric(2), USE.NAMES = FALSE))
  top$lo <- ifelse(is.na(spans[, 1]), top$pos, spans[, 1])
  top$hi <- ifelse(is.na(spans[, 2]), top$pos, spans[, 2])

  top <- top[order(top$gene, top$lo, top$hi, top$rank), ]
  top$stack <- NA_character_
  for (g in names(B_WIN)) {
    ix <- which(top$gene == g)
    lo <- top$lo[ix]; hi <- top$hi[ix]
    grp <- integer(length(ix)); reach <- -Inf; k <- 0L
    for (j in seq_along(ix)) {
      if (lo[j] > reach) { k <- k + 1L; reach <- hi[j] } else reach <- max(reach, hi[j])
      grp[j] <- k
    }
    top$stack[ix] <- paste0(g, "_", ave(lo, grp, FUN = min))
  }

  # Circles BELOW the body: samples outside the top 60 that carry one of
  # the plotted changes and have GSP #13 < 1e-4 (three as of 2026-10-08).
  plotted <- unique(top$change)
  # Only carriers with NO detectable phenotype (GSP #13 < 1e-4, the control
  # threshold used throughout), so the contrast with the top 60 is clean.
  # (Quinn, 2026-10-08: drops PRPF8 R1832H at rank 130, GSP #13 = 3.1e-3.)
  B_BELOW_MAX <- 1e-4
  below <- long[long$rank > B_TOPN & long$k13 < B_BELOW_MAX & long$change %in% plotted, ]
  below <- below[order(below$gene, below$pos, below$rank), ]
  below$stack <- top$stack[match(below$change, top$change)]
  stopifnot(nrow(below) >= 1L, all(below$k13 < B_BELOW_MAX), !anyNA(below$stack))
  cat(sprintf("[panel b] %d circles below the bodies (GSP #13 < %g): %s\n", nrow(below),
              B_BELOW_MAX, paste(below$gene, below$change, collapse = ", ")))

  # ----------------------------------------------------------------
  # Stack label: collapse the shared prefix progressively, exactly as
  # the reference does. R402* + R402Q + R402Sfs*13 -> "R402*/Q/Sfs*13".
  # ----------------------------------------------------------------
  collapse_changes <- function(ch, pos) {
    u <- unique(ch[order(pos, ch)])
    if (length(u) == 1) return(u)
    out <- u[1]
    for (i in 2:length(u)) {
      a <- u[i - 1]; b <- u[i]
      n <- 0L
      m <- min(nchar(a), nchar(b))
      while (n < m && substr(a, n + 1, n + 1) == substr(b, n + 1, n + 1)) n <- n + 1L
      # Only collapse a prefix long enough to be recognised as such, and
      # never to nothing.
      if (n >= 3L && n < nchar(b)) b <- substring(b, n + 1)
      out <- c(out, b)
    }
    paste(out, collapse = "/")
  }

  stacks <- do.call(rbind, lapply(split(top, top$stack), function(s) {
    data.frame(stack = s$stack[1], gene = s$gene[1], pos = min(s$lo),
               n = nrow(s), label = collapse_changes(s$change, s$lo),
               stringsAsFactors = FALSE)
  }))
  stacks <- stacks[order(stacks$gene, stacks$pos), ]
  rownames(stacks) <- NULL

  # Stack order, bottom -> top: descending GSP #13, so the strongest
  # carrier sits on the body. (The reference's own order is not
  # reproducible from the current loadings -- see build.md.)
  top <- top[order(top$stack, top$rank), ]
  top$slot <- ave(top$rank, top$stack, FUN = seq_along) - 1L
  below$slot <- ave(below$rank, below$stack, FUN = seq_along) - 1L

  # ================================================================
  # 3. Domain cartoons + the broken aa -> inch map.
  # ================================================================
  # Each protein is drawn as GROUPS of abutting boxes joined by short
  # connectors. Ticks fall on group edges, so they can never disagree
  # with the boxes.
  B_PROTEIN <- list(
    NKAP = list(
      len = 415, ticks = c(0, 272, 415),
      brk = c(0, 272, 415), frac = c(0.24, 0.76),
      groups = list(
        list(c(0, 140, 272), c("RS Domain", "Basic Domain")),
        list(c(285, 415),    c("Domain of Unknown Function")))),
    PRPF8 = list(
      len = 2305, ticks = c(0, 1579, 1754, 2000, 2305),
      # frac[1] 0.15 -> 0.22, paid for by frac[3] 0.47 -> 0.40.
      # Quinn: "3b the mutations are falling off the edge of the RNase
      # H-domain, which, you can make smaller so as to allow the text
      # for PRPF8's earlier domains to actually be readable".
      # The four N-terminal boxes (NTD/HB/RT/Thumb/X) all live in
      # segment 1, and the auto-fit sizes a domain name from its BOX
      # width, so at 0.15 both NTD and Thumb/X were pinned to the
      # 3.2 pt floor and "Thumb/X" ran its X into the box border.
      # Segment 3 (aa 1754-2000) holds the RNase H-Domain and the
      # mutation-dense run; it can afford the loss because the width
      # constraint has slack here (the auto-fit reports xover = 0, i.e.
      # it is HEIGHT-bound, not width-bound). Re-check that xover stays
      # 0 after any further change.
      brk = c(0, 1579, 1754, 2000, 2305), frac = c(0.22, 0.24, 0.40, 0.14),
      groups = list(
        list(c(0, 289, 558, 905, 1406), c("NTD", "HB", "RT", "Thumb/X")),
        list(c(1579, 1754, 2000), c("Endonuclease-like Domain", "RNase H-Domain")),
        list(c(2090, 2305), c("Jab1/MPN")))))

  aa2in <- function(g) {
    p <- B_PROTEIN[[g]]
    cum <- c(0, cumsum(p$frac)) * B_BODY_W
    function(aa) {
      aa <- pmin(pmax(aa, p$brk[1]), p$brk[length(p$brk)])
      i <- pmax(1, pmin(length(p$frac),
                        findInterval(aa, p$brk, rightmost.closed = TRUE)))
      cum[i] + (aa - p$brk[i]) / (p$brk[i + 1] - p$brk[i]) * (cum[i + 1] - cum[i])
    }
  }
  B_MAP <- lapply(names(B_PROTEIN), aa2in); names(B_MAP) <- names(B_PROTEIN)

  # ================================================================
  # 4. Label wrapping. Breaks only at "/" and "_", which are the real
  #    token boundaries in an HGVS-ish protein change, and packs
  #    greedily to a maximum drawn LINE width.
  # ================================================================
  wrap_change <- function(s, max_w, pt) {
    split_keep <- function(x, ch) {
      parts <- strsplit(x, ch, fixed = TRUE)[[1]]
      if (length(parts) <= 1) return(x)
      out <- paste0(parts, ch)
      if (!endsWith(x, ch)) out[length(out)] <- parts[length(parts)]
      out[nzchar(out)]
    }
    # "/" separates the merged changes and is the only break a reader
    # can re-join unambiguously, so it is always tried first; "_" (a
    # residue RANGE inside one change) is only used when a single
    # "/"-chunk is still too wide on its own.
    toks <- split_keep(s, "/")
    toks <- unlist(lapply(toks, function(t)
      if (str_w(t, pt) > max_w) split_keep(t, "_") else t))
    lines <- character(0); cur <- ""
    for (t in toks) {
      cand <- paste0(cur, t)
      if (nzchar(cur) && str_w(cand, pt) > max_w) { lines <- c(lines, cur); cur <- t }
      else cur <- cand
    }
    paste(c(lines, cur), collapse = "\n")
  }

  # ================================================================
  # 5. Layout. One function of the fit scale S; bisected in section 6
  #    so the content fills the cell height exactly.
  # ================================================================
  geom_for <- function(S) {
    lab_pt  <- B_LAB_PT0  * S; tick_pt <- B_TICK_PT0 * S
    gene_pt <- B_GENE_PT0 * S; note_pt <- B_NOTE_PT0 * S
    dom_pt  <- B_DOM_PT0  * S
    cd      <- B_CIRC_D0  * S
    # lh is the line pitch ASSUMED by the collision model. The drawn
    # lineheight (B_LAB_LH, below 1) is deliberately tighter, so the
    # lines of one label sit closer together than any two labels do --
    # that is what makes a 3-line label read as one object.
    lh      <- 1.00 * lab_pt / 72
    maxw    <- B_LAB_MAXW0 * S

    st <- stacks
    st$wrapped <- vapply(st$label, wrap_change, character(1),
                         max_w = maxw, pt = lab_pt, USE.NAMES = FALSE)
    st$nlines  <- lengths(strsplit(st$wrapped, "\n", fixed = TRUE))
    st$linew   <- str_w(st$wrapped, lab_pt)
    st$below   <- vapply(st$stack, function(k) sum(below$stack == k), integer(1))

    # --- x: nominal position, then two collision sweeps -------------
    st$x <- NA_real_
    xover <- 0                      # ink left OUTSIDE the cell after (c)
    for (g in names(B_MAP)) {
      ix <- which(st$gene == g)
      x  <- B_MAP[[g]](st$pos[ix])
      yA <- st$n[ix] * cd            # anchor height above the body top
      # (a) circles must not overlap, and (b) two labels are parallel
      #     lines that collide unless their PERPENDICULAR separation
      #     clears half a text block each (perp = dx*sin - dy*cos).
      #
      # ONE pass, taking the max of both. They used to be two passes,
      # and that silently broke (a): pass (b) moves x[i-1] right AFTER
      # pass (a) has spaced x[i] off its old position, so wherever the
      # perpendicular term came out below `cd` -- which happens whenever
      # the stack gets SHORTER (dy < 0, e.g. R368Q's 3 circles then
      # N371K's 1) -- the two stacks' circles ended up overlapping.
      gapmin <- rep(0, length(x))       # required x[i] - x[i-1]
      for (i in seq_along(x)[-1]) {
        # (n/2)*lh is a block's perpendicular ink half-extent; the extra
        # 0.85*lh is white space, and it is LARGER than the drawn line
        # pitch inside a label -- without that ordering two adjacent
        # labels read as one string ("insQ" + "P1880_..." run together).
        need <- (st$nlines[ix][i - 1] + st$nlines[ix][i]) / 2 * lh + 0.85 * lh
        dy   <- yA[i] - yA[i - 1]
        gapmin[i] <- max(cd, (need + dy * B_COS) / B_SIN)
        x[i] <- max(x[i], x[i - 1] + gapmin[i])
      }
      # (d) A stack may not leave the drawn stretch its residue lives in.
      #     The forward sweep only pushes RIGHT, and it pushed PRPF8's
      #     P1880/N1881 column (aa 1880-1883, inside the RNase H domain,
      #     which ends at 2000) out over the LINKER to Jab1/MPN -- a
      #     circle floating over a region its residue is not in. So cap
      #     each stack at the right edge of its own segment (a domain
      #     group, or the unannotated gap between two groups -- e.g.
      #     PRPF8 aa 1406-1579, where Q1575H genuinely sits) and sweep
      #     BACKWARD, pulling earlier stacks left to keep their spacing.
      seg_edges <- sort(unique(unlist(lapply(B_PROTEIN[[g]]$groups, function(gr)
                          range(gr[[1]])))))
      seg_edges <- sort(unique(c(0, seg_edges, B_PROTEIN[[g]]$len)))
      k  <- findInterval(st$pos[ix], seg_edges, rightmost.closed = TRUE)
      # Keep a stack ~2 circle diameters inside its segment's right edge,
      # not half of one: at cd/2 the aa-1880 column sat right on the
      # "2000" tick and read as if it were at residue 2000 (Quinn).
      ub <- B_MAP[[g]](seg_edges[pmin(k + 1, length(seg_edges))]) - 2.2 * cd
      lb <- B_MAP[[g]](seg_edges[k]) + cd / 2
      n  <- length(x)
      x[n] <- min(x[n], ub[n])
      if (n > 1) for (i in (n - 1):1)
        x[i] <- min(x[i], ub[i], x[i + 1] - gapmin[i + 1])
      # ...and forward again against the LEFT edges, in case the pull-back
      # dragged a stack out of its own segment on that side. If a segment
      # is genuinely too short for its stacks this pass wins (overlap is
      # preferable to a circle over the wrong region), and is reported.
      for (i in seq_len(n)) {
        lo <- if (i > 1) max(lb[i], x[i - 1] + gapmin[i]) else lb[i]
        if (x[i] < lo) x[i] <- lo
      }
      bad <- which(x > ub + 1e-9)
      if (length(bad) && !isTRUE(getOption(".pb_warned_seg")))
        options(.pb_warned_seg = TRUE)
      if (length(bad) && identical(getOption(".pb_warned_seg_n"), NULL) &&
          { options(.pb_warned_seg_n = 1); TRUE })
        cat(sprintf("[panel b] %s: %d stack(s) could not fit inside their segment: %s\n",
                    g, length(bad), paste(st$label[ix][bad], collapse = ", ")))
      # (c) both sweeps only ever push RIGHT, so the last label can run
      #     into the legend column. Slide the whole gene's set back by
      #     the overshoot, but never left of the functional window's
      #     own start.
      # The slide is bounded by how far the first stack may travel left,
      # so it can FAIL to clear the edge. What it cannot absorb is
      # returned as `xover` and the auto-fit treats it exactly like a
      # height overflow -- otherwise the tallest label silently printed
      # outside the cell (it used to land on the legend).
      reach <- x + st$nlines[ix] / 2 * lh * B_SIN + st$linew[ix] * B_COS
      over  <- max(reach) - B_X_RIGHT
      if (over > 0) {
        shift <- min(over, max(0, x[1] - B_MAP[[g]](B_WIN[[g]][1])))
        x     <- x - shift
        xover <- max(xover, over - shift)
      }
      st$x[ix] <- x
    }

    # --- y: build each protein block from the bottom up -------------
    body_h  <- max(0.070 * S, 1.35 * dom_pt / 72)
    # tick numbers and the italic gene name share one row, both hung
    # from the row's top edge; the row must hold a full line of the
    # LARGER of the two. Neither digits nor NKAP/PRPF8 have a descender,
    # so 1.30 lineheights was 0.03 in of nothing, twice over.
    # B_TICK_LEN + B_TICK_GAP: the numbers start below the tick marks (they overlapped, Quinn 2026-11-09).
    tick_h  <- 1.18 * max(tick_pt, gene_pt) / 72 + (B_TICK_LEN + B_TICK_GAP) * cd
    lab_pad <- 0.28 * cd

    # Per-stack label INK extent. The label is anchored at
    # BODY_T + n*cd + pad with hjust = 0, vjust = 0.5 and rotated by
    # B_ANGLE, so its ink runs from that anchor up-and-right and spreads
    # (nlines*lh/2) perpendicular to the text direction -- which is
    # +-(nlines*lh/2)*sin in x and +-(nlines*lh/2)*cos in y.
    #
    # Reserving `max(n)*cd + max(ink)` (what this used to do) over-
    # reserves whenever the tallest STACK and the longest LABEL are
    # different stacks, and it gives no per-stack extents at all -- and
    # those extents are what let the "outside the top 60" caption slide
    # under PRPF8's label band below instead of costing a full row.
    # The perpendicular half-extent is nlines*lh/2 -- the same half the
    # x extents below and the collision sweep above both use. `ink` used
    # the FULL nlines*lh, which reserved half a label block of empty air
    # over the top of each staircase, twice, for 0.10 in of the cell.
    # (nlines*lh/2 still carries ~20% margin over the real ink: the
    # block's true half-extent is ((nlines-1)*B_LAB_LH + 0.72)/2 lines.)
    st$ink    <- st$linew * B_SIN + st$nlines * lh / 2 * B_COS
    st$ink_t  <- st$n * cd + lab_pad + st$ink
    st$ink_x1 <- st$x - st$nlines * lh / 2 * B_SIN
    st$ink_x2 <- st$x + st$linew * B_COS + st$nlines * lh / 2 * B_SIN

    # The band UNDER a body carries two things: the circles for carriers
    # outside the top 60 (two or three x positions) and the aa-axis row
    # (tick numbers on the GROUP EDGES, italic gene name hard left).
    # They only have to be stacked where their ink actually shares an x
    # span -- which, because the ticks sit on group edges, it never
    # does here. Reserving the SUM unconditionally cost 0.21 in of a
    # 2.58 in cell for a clash that does not happen; reserve the MAX and
    # start the axis row at the body's own bottom edge. `off` is how far
    # below the body that row begins, so the test drives the drawing too.
    foot_of <- function(g) {
      ix <- which(st$gene == g)
      nb <- st$below[ix]
      bh <- max(nb) * cd
      if (bh <= 0) return(c(foot = tick_h, off = 0))
      cx  <- st$x[ix][nb > 0]
      p   <- B_PROTEIN[[g]]
      tx  <- B_MAP[[g]](p$ticks)
      tw  <- str_w(as.character(p$ticks), tick_pt)
      hj  <- c(0, rep(0.5, length(p$ticks) - 2), 1)
      row <- cbind(tx - hj * tw, tx + (1 - hj) * tw)
      gx  <- tw[1] + 0.55 * cd
      row <- rbind(row, c(gx, gx + str_w(g, gene_pt, "italic")))
      pad <- 0.25 * cd
      hit <- any(outer(cx - cd / 2, row[, 2] + pad, "<") &
                 outer(cx + cd / 2, row[, 1] - pad, ">"))
      if (hit) c(foot = bh + tick_h, off = bh)
      else     c(foot = max(bh, tick_h), off = 0)
    }
    blk <- function(g) {
      ix <- st$gene == g
      f  <- foot_of(g)
      c(above = max(st$ink_t[ix]), body = body_h,
        belowb = max(st$below[ix]) * cd,
        foot = f[["foot"]], tickoff = f[["off"]])
    }
    bN <- blk("NKAP"); bP <- blk("PRPF8")

    # --- where the "outside the top 60" caption goes ----------------
    # It hangs from the two NKAP circles it annotates on two SHORT
    # leaders of equal length and equal slope, mirrored about its own
    # centre (section 7). That geometry only closes if the caption is
    # CENTRED on those two circles, so its x is no longer slid sideways
    # to dodge PRPF8's tall labels -- it is pinned, and what the
    # clearance costs is paid vertically, in `drop`.
    #
    # The earlier solve slid it left until it was free, which bought the
    # height back but left one leader vertical and the other a 0.6 in
    # near-horizontal diagonal reaching for a circle 0.5 in outside the
    # caption. That pair is what Quinn read as inconsistent. The height
    # this costs is roughly what `foot_of` above just gave back.
    note_ln <- strsplit(B_CAP2_TXT, "\n", fixed = TRUE)[[1]]
    note_w <- max(str_w(note_ln, note_pt))
    note_h  <- length(note_ln) * 1.15 * note_pt / 72
    foot_w <- str_w(B_FOOT_TXT, note_pt)
    ixP    <- which(st$gene == "PRPF8")
    nb_x   <- sort(st$x[st$gene == "NKAP" &
                        st$stack %in% below$stack[below$gene == "NKAP"]])
    # Clamped only so the caption cannot run off the drawing area or
    # into the footnote that shares the band, hard left.
    capx1  <- min(max(mean(range(nb_x)) - note_w / 2, foot_w + 0.14),
                  B_X_RIGHT - note_w)
    capx2  <- capx1 + note_w
    # How far the caption's ink hangs BELOW the band's own top edge: the
    # circles' own depth, then the leader, then the two text lines.
    stub_dy <- B_STUB_LEN * cd * cos(B_STUB_ANG * pi / 180)
    # Clear air between where a leader STOPS and where the caption's
    # first line starts. Without it the two coincide, and the left-hand
    # leader drove its stroke straight through the "A" of "TCGA"
    # (visible at 1200 dpi). In circle diameters, so it scales with the
    # auto-fit like every other offset here.
    stub_pad <- B_STUB_PAD * cd
    cap_dep <- max(0, bN[["belowb"]] + stub_dy + stub_pad + note_h -
                      bN[["foot"]])
    ov      <- ixP[st$ink_x2[ixP] > capx1 - lab_pad &
                   st$ink_x1[ixP] < capx2 + lab_pad]
    drop    <- max(bP[["above"]] + lab_pad,
                   cap_dep + lab_pad +
                     if (length(ov)) max(st$ink_t[ov]) else 0)

    list(S = S, lab_pt = lab_pt, tick_pt = tick_pt, gene_pt = gene_pt,
         note_pt = note_pt, dom_pt = dom_pt, cd = cd, lh = lh,
         st = st, bN = bN, bP = bP, body_h = body_h, note_h = note_h,
         lab_pad = lab_pad, capx1 = capx1, capx2 = capx2, drop = drop,
         note_w = note_w, foot_w = foot_w, xover = xover, maxw = maxw,
         stub_dy = stub_dy, stub_pad = stub_pad,
         stub_dx = B_STUB_LEN * cd * sin(B_STUB_ANG * pi / 180),
         h_req = bN[["above"]] + bN[["body"]] + bN[["foot"]] + drop +
                 bP[["body"]] + bP[["foot"]])
  }

  # ================================================================
  # 6. Auto-fit: one scalar drives every size, bisected so the drawn
  #    height matches the cell. Nothing below is hand-tuned.
  # ================================================================
  lo <- 0.35; hi <- 1.60
  fits <- function(S) { g <- geom_for(S); g$h_req <= B_H_TOT && g$xover <= 1e-4 }
  for (.i in 1:26) {
    mid <- (lo + hi) / 2
    if (fits(mid)) lo <- mid else hi <- mid
  }
  G  <- geom_for(lo * 0.995)
  st <- G$st; cd <- G$cd
  cat(sprintf("[panel b] fit S=%.3f  h_req=%.3f  h_avail=%.3f in  xover=%.4f\n",
              G$S, G$h_req, B_H_TOT, G$xover))

  slack <- B_H_TOT - G$h_req

  # y anchors, top of the panel downwards. Every band is exactly the ink
  # it holds, so the two cartoons plus the caption fill the cell with
  # nothing left over -- `slack` is what the bisection could not spend
  # and is ~1e-3 in.
  yN_lab_top <- B_H_TOT
  yN_body_t  <- yN_lab_top - G$bN[["above"]]
  yN_body_b  <- yN_body_t - G$body_h
  # The aa-axis row starts `tickoff` below the body -- 0 wherever the
  # below-body circles and the axis row do not share an x span, which is
  # both proteins here (see foot_of). `foot` is the whole band.
  yN_tick_t  <- yN_body_b - G$bN[["tickoff"]]
  yN_note_t  <- yN_body_b - G$bN[["foot"]]     # top edge of the middle band
  yP_body_t  <- yN_note_t - G$drop - slack
  yP_body_b  <- yP_body_t - G$body_h
  yP_tick_t  <- yP_body_b - G$bP[["tickoff"]]

  # Calibration handle: at 400 dpi the NKAP and PRPF8 body TOP edges
  # must land (crop_top + (B_H_TOT - y) * 400) px apart and the NKAP
  # body must span x = 0 .. B_BODY_W. Measure those two edges in the
  # render and reset B_CAL_W / B_CAL_H if they disagree.
  cat(sprintf("[panel b] cal: body tops y = %.4f / %.4f in, body_w = %.4f in\n",
              yN_body_t, yP_body_t, B_BODY_W))

  BODY_T <- c(NKAP = yN_body_t, PRPF8 = yP_body_t)
  BODY_B <- c(NKAP = yN_body_b, PRPF8 = yP_body_b)
  TICK_T <- c(NKAP = yN_tick_t, PRPF8 = yP_tick_t)

  # ================================================================
  # 7. Drawing data. Everything an aes() touches is a COLUMN -- a
  #    cached panel is restored by readRDS into an empty environment.
  # ================================================================
  # ---- protein bodies --------------------------------------------
  boxes <- do.call(rbind, lapply(names(B_PROTEIN), function(g) {
    p <- B_PROTEIN[[g]]; m <- B_MAP[[g]]
    do.call(rbind, lapply(p$groups, function(gr) {
      e <- gr[[1]]; nm <- gr[[2]]
      data.frame(gene = g, x1 = m(e[-length(e)]), x2 = m(e[-1]),
                 lab = nm, stringsAsFactors = FALSE)
    }))
  }))
  boxes$ymin <- BODY_B[boxes$gene]; boxes$ymax <- BODY_T[boxes$gene]
  boxes$xm   <- (boxes$x1 + boxes$x2) / 2
  boxes$ym   <- (boxes$ymin + boxes$ymax) / 2
  # Domain names are shrunk per box only where the box is too narrow
  # ("Endonuclease-like Domain"); everything else keeps one size.
  boxes$pt <- pmax(G$dom_pt * pmin(1, (boxes$x2 - boxes$x1) * 0.92 /
                                     str_w(boxes$lab, G$dom_pt)), 3.2)
  boxes$size <- boxes$pt / .pt

  conn <- do.call(rbind, lapply(names(B_PROTEIN), function(g) {
    p <- B_PROTEIN[[g]]; m <- B_MAP[[g]]
    if (length(p$groups) < 2) return(NULL)
    ends <- vapply(p$groups, function(gr) max(gr[[1]]), numeric(1))
    begs <- vapply(p$groups, function(gr) min(gr[[1]]), numeric(1))
    data.frame(gene = g, x = m(ends[-length(ends)]), xend = m(begs[-1]),
               y = unname((BODY_T[g] + BODY_B[g]) / 2),
               row.names = NULL, stringsAsFactors = FALSE)
  }))
  conn$yend <- conn$y

  ticks <- do.call(rbind, lapply(names(B_PROTEIN), function(g) {
    p <- B_PROTEIN[[g]]
    n <- length(p$ticks)
    # The END tick numbers are aligned INWARD, not centred. A centred
    # "0" hangs half its width outside the cell (the body starts at
    # x = 0) and a centred "2305"/"415" runs into the legend column,
    # where it printed on top of the "*Called via WXS" footnote.
    data.frame(gene = g, x = B_MAP[[g]](p$ticks), lab = as.character(p$ticks),
               hj = c(0, rep(0.5, n - 2), 1),
               stringsAsFactors = FALSE)
  }))
  ticks$y    <- BODY_B[ticks$gene]
  ticks$yend <- BODY_B[ticks$gene] - B_TICK_LEN * G$cd
  ticks$ytx  <- TICK_T[ticks$gene] - (B_TICK_LEN + B_TICK_GAP) * G$cd

  # Gene name shares the tick row; it starts clear of the "0" tick
  # label, whose half-width is measured, not guessed.
  # Dropped 0.30 cd lower than the tick numbers (Quinn, 2026-11-09): hung at the tick row's top edge the
  # italic names touched the bottom of the domain maps.
  gname <- data.frame(
    gene = names(B_PROTEIN),
    y    = TICK_T[names(B_PROTEIN)] - (0.34 + B_TICK_LEN + B_TICK_GAP - 0.04) * G$cd,
    lab  = names(B_PROTEIN), stringsAsFactors = FALSE)
  gname$x <- str_w(ticks$lab[match(gname$gene, ticks$gene)], G$tick_pt) +
             0.55 * G$cd
  gname$size <- G$gene_pt / .pt

  # ---- lollipops --------------------------------------------------
  circ <- data.frame(
    x0 = st$x[match(top$stack, st$stack)],
    y0 = BODY_T[top$gene] + (top$slot + 0.5) * cd,
    r  = cd / 2 * 0.94,
    fill = unname(B_PROJ_COL[top$project]),
    star = top$star, stringsAsFactors = FALSE)
  circ_b <- data.frame(
    x0 = st$x[match(below$stack, st$stack)],
    y0 = BODY_B[below$gene] - (below$slot + 0.5) * cd,
    r  = cd / 2 * 0.94,
    fill = unname(B_PROJ_COL[below$project]),
    star = below$star, stringsAsFactors = FALSE)
  circ_all <- rbind(circ, circ_b)
  stars <- circ_all[circ_all$star, c("x0", "y0")]
  names(stars) <- c("x", "y")

  # ---- rotated stack labels ---------------------------------------
  labs_df <- data.frame(
    x = st$x, y = BODY_T[st$gene] + st$n * cd + G$lab_pad,
    lab = st$wrapped, angle = B_ANGLE, size = G$lab_pt / .pt,
    stringsAsFactors = FALSE)

  # ---- legend: a 4 x 4 block in the dead upper-left ----------------
  cnt   <- sort(table(c(top$project, below$project)), decreasing = TRUE)
  projs <- names(cnt)[order(-as.integer(cnt), names(cnt))]

  # How far right the block may run is not a guess. Each rotated label
  # is a strip anchored at (x, y0) running up-right at B_ANGLE, so
  # inside a y band its leftmost ink is
  #   x + max(0, y - y0)/tan(angle) - (nlines*lh/2)*sin(angle).
  # The minimum of that over NKAP's labels IS the free width.
  lab_y0 <- BODY_T[st$gene] + st$n * cd + G$lab_pad
  free_x <- function(yb, yt) {
    ix   <- which(st$gene == "NKAP")
    half <- st$nlines[ix] * G$lh / 2
    ylo  <- lab_y0[ix] - half * B_COS
    yhi  <- lab_y0[ix] + st$linew[ix] * B_SIN + half * B_COS
    ok   <- yhi > yb & ylo < yt
    if (!any(ok)) return(B_W_TOT)
    yy <- pmax(yb, ylo[ok])
    min(st$x[ix][ok] + pmax(0, yy - lab_y0[ix][ok]) / B_SIN * B_COS -
        half[ok] * B_SIN)
  }

  # Solve the legend's type size against that free width, exactly as the
  # old right-hand column solved it against its 0.62 in. It has to be a
  # BISECTION, not a one-way shrink: the block's width grows with the
  # type while the free width SHRINKS with it (a taller block reaches
  # down into the staircase), so a proportional step overshoots and a
  # shrink-only loop then exits at the first size that happens to fit --
  # it settled on 3.9 pt where 5.2 pt fits.
  leg_need <- function(pt) {
    dsw <- B_LEG_D_LINE * pt / 72
    ent <- dsw * (1.32 + B_LEG_GAPC) + max(str_w(projs, pt))
    c(ent = ent, dsw = dsw, need = B_LEG_NCOL * ent - B_LEG_GAPC * dsw)
  }
  leg_fits <- function(pt) {
    blkh <- B_LEG_NROW * B_LEG_ROWF * pt / 72
    leg_need(pt)[["need"]] <=
      free_x(B_H_TOT - B_LEG_TOP - blkh, B_H_TOT - B_LEG_TOP) - B_LEG_CLR
  }
  # The ceiling is the PROSE size, not SIZE_LEGEND: a key set larger
  # than the mutation labels it explains inverts the hierarchy, and in a
  # 2.5 in block the free width alone would have allowed 5.6 pt.
  B_LEG_CEIL <- min(B_LEG_PT_MAX, G$note_pt)
  if (leg_fits(B_LEG_CEIL)) {
    B_LEG_PT <- B_LEG_CEIL
  } else {
    plo <- 2.5; phi <- B_LEG_CEIL
    for (.i in 1:30) {
      pm <- (plo + phi) / 2
      if (leg_fits(pm)) plo <- pm else phi <- pm
    }
    B_LEG_PT <- plo
  }
  rowh   <- B_LEG_ROWF * B_LEG_PT / 72
  .ld    <- leg_need(B_LEG_PT)
  entw   <- .ld[["ent"]]; B_LEG_D <- .ld[["dsw"]]
  # Columns are JUSTIFIED across the free width rather than packed hard
  # left: when the size is capped by the prose rather than by the space,
  # packing leaves a wedge of white between the block and the tall
  # labels. Capped at 1.5x the packed pitch so it never reads as gappy.
  leg_fw  <- free_x(B_H_TOT - B_LEG_TOP - B_LEG_NROW * rowh,
                    B_H_TOT - B_LEG_TOP) - B_LEG_CLR
  ent_ink <- B_LEG_D * 1.32 + max(str_w(projs, B_LEG_PT))
  pitch_c <- min(1.5 * entw,
                 max(entw, (leg_fw - ent_ink) / (B_LEG_NCOL - 1)))
  ii     <- seq_along(projs) - 1L
  leg <- data.frame(
    x0 = 0.01 + (ii %% B_LEG_NCOL) * pitch_c + B_LEG_D / 2,
    y0 = B_H_TOT - B_LEG_TOP - rowh / 2 - (ii %/% B_LEG_NCOL) * rowh,
    r  = B_LEG_D / 2 * 0.94,
    fill = unname(B_PROJ_COL[projs]), stringsAsFactors = FALSE)
  leg_tx <- data.frame(x = leg$x0 + B_LEG_D * 0.82, y = leg$y0,
                       lab = projs, size = B_LEG_PT / .pt,
                       stringsAsFactors = FALSE)
  leg_bot <- B_H_TOT - B_LEG_TOP - B_LEG_NROW * rowh
  cat(sprintf("[panel b] legend: %.2f pt, block %.2f x %.2f in (free %.2f), bottom y = %.3f\n",
              B_LEG_PT, (B_LEG_NCOL - 1) * pitch_c + ent_ink,
              B_H_TOT - B_LEG_TOP - leg_bot, leg_fw, leg_bot))

  # ---- "top 60" bracket + its caption -----------------------------
  nk      <- st[st$gene == "NKAP", ]
  br_x    <- min(nk$x) - 1.15 * cd
  br_lo   <- BODY_T[["NKAP"]]
  # Set on ONE line whenever the room left of the bracket takes it. The
  # two-line break was sized for the old narrow pocket; on one line the
  # caption spans the width the legend block above it does, and its
  # smaller height buys the bracket back the space it costs.
  cap1_txt <- B_CAP1_TXT
  cap1_x   <- br_x - 0.55 * cd
  .one     <- gsub("\n", " ", B_CAP1_TXT, fixed = TRUE)
  if (str_w(.one, G$note_pt) <= cap1_x - 0.02) cap1_txt <- .one
  cap1_h   <- lengths(strsplit(cap1_txt, "\n", fixed = TRUE)) *
              1.15 * G$note_pt / 72
  # The bracket and its caption share the upper-left with the legend, so
  # the top is clamped to clear it BY MEASUREMENT: cap1 is centred on the
  # bracket, so its ink top is (br_lo + br_hi)/2 + cap1_h/2.
  br_hi   <- max(br_lo + 0.5 * cd,
                 min(BODY_T[["NKAP"]] + max(nk$n) * cd + 0.30 * max(nk$ink),
                     2 * (leg_bot - 0.05 - cap1_h / 2) - br_lo,
                     leg_bot - 0.05))
  br_tip  <- 0.45 * cd
  bracket <- data.frame(
    x    = c(br_x + br_tip, br_x, br_x,  br_x),
    xend = c(br_x,          br_x, br_x,  br_x + br_tip),
    y    = c(br_hi,         br_hi, br_lo, br_lo),
    yend = c(br_hi,         br_lo, br_lo, br_lo))
  cap1 <- data.frame(
    x = cap1_x, y = (br_lo + br_hi) / 2,
    lab = cap1_txt,
    size = G$note_pt / .pt, stringsAsFactors = FALSE)

  # ---- below-the-body caption + its two leaders --------------------
  # The caption hangs from its own two circles: it is centred on them
  # (section 5) and its ink top is exactly one leader below them, so the
  # two leaders come out the SAME LENGTH and the SAME SLOPE, mirrored
  # about the caption's centre. Neither is parallel to a page edge, and
  # neither is long: previously one was a 0.23 in vertical riser and the
  # other a 0.61 in near-horizontal reach to a circle half an inch
  # outside the caption, and the pair read as two different ideas.
  # They cross the aa-axis row, which carries no ink at those x -- the
  # numbers sit on the GROUP EDGES, 0 / 272 / 415.
  isN    <- circ_b$y0 > 0.5 * (BODY_B[["NKAP"]] + BODY_T[["PRPF8"]])
  .o     <- order(circ_b$x0[isN])
  nbx    <- circ_b$x0[isN][.o]
  nby    <- (circ_b$y0[isN] - circ_b$r[isN])[.o]
  cap2_x <- G$capx1
  # The leaders stop at lead_end; the caption's top ink starts one
  # stub_pad BELOW that, so no stroke reaches the text.
  lead_end <- min(nby) - G$stub_dy
  cap_top <- lead_end - G$stub_pad
  cap2_y <- cap_top - G$note_h / 2
  cap2  <- data.frame(
    x = cap2_x, y = cap2_y,
    lab = B_CAP2_TXT,
    size = G$note_pt / .pt, stringsAsFactors = FALSE)
  # Sorted left-to-right, the left half of the circles lean right and
  # the right half lean left: a converging pair whose vertex would be
  # the caption. Taking the direction from the ORDER, not from the sign
  # of (x - centre), keeps the two mirrored even if both circles ended
  # up on the same side of the caption after the clamp.
  .dir  <- ifelse(seq_along(nbx) <= length(nbx) / 2, 1, -1)
  lead  <- data.frame(
    x    = nbx,
    y    = nby,
    xend = pmin(pmax(nbx + .dir * G$stub_dx, cap2_x + 0.03),
                cap2_x + G$note_w - 0.03),
    yend = lead_end)
  # Both prose blocks share the middle band; the footnote keeps the old
  # treatment (CENTRED in its own free slot rather than hung from the
  # band's top edge -- pinned to the top it left one 0.34 in blank strip
  # running the full panel width between it and PRPF8's labels). The
  # floor of a slot is the tallest PRPF8 label whose ink is under it.
  band_t  <- yN_note_t
  band_b  <- yP_body_t
  pr_top  <- function(a, b) {
    ix <- which(st$gene == "PRPF8")
    ok <- st$ink_x2[ix] > a - G$lab_pad & st$ink_x1[ix] < b + G$lab_pad
    if (!any(ok)) band_b else band_b + max(st$ink_t[ix][ok])
  }
  slot_y  <- function(a, b, h) {
    fl <- pr_top(a, b) + G$lab_pad
    max(fl + h / 2, (fl + band_t) / 2)
  }
  # The WXS/RNA-seq footnote belongs between the two cartoons -- the
  # starred circles are in both -- and it is the only furniture that can
  # fill the left of the middle band, which PRPF8's own short labels
  # cannot reach. It is therefore CENTRED in the slot, across as well as
  # down: hung hard left it left the whole column between its right end
  # and the caption empty, and that column was the panel's largest blank
  # rectangle (1.05 x 0.85 in). Centred it straddles that column, and
  # the residual splits in two.
  foot_x <- max(0.02, (cap2_x - G$lab_pad - G$foot_w) / 2)
  foot <- data.frame(
    x = foot_x,
    y = slot_y(foot_x, foot_x + G$foot_w, 1.15 * G$note_pt / 72),
    lab = B_FOOT_TXT,
    size = G$note_pt / .pt, stringsAsFactors = FALSE)

  # ================================================================
  # 8. Assemble.
  # ================================================================
  p <- ggplot() +
    # protein bodies
    geom_rect(data = boxes, aes(xmin = x1, xmax = x2, ymin = ymin, ymax = ymax),
              fill = "white", colour = "black", linewidth = 0.25) +
    geom_segment(data = conn, aes(x = x, xend = xend, y = y, yend = yend),
                 colour = "black", linewidth = 0.4) +
    geom_text(data = boxes, aes(x = xm, y = ym, label = lab, size = size),
              family = FONT_FAMILY, colour = "black", vjust = 0.5) +
    # aa axis
    geom_segment(data = ticks, aes(x = x, xend = x, y = y, yend = yend),
                 colour = "black", linewidth = 0.25) +
    geom_text(data = ticks, aes(x = x, y = ytx, label = lab, hjust = hj),
              size = G$tick_pt / .pt, family = FONT_FAMILY, vjust = 1) +
    geom_text(data = gname, aes(x = x, y = y, label = lab, size = size),
              family = FONT_FAMILY, fontface = "italic", hjust = 0, vjust = 1) +
    # lollipops
    ggforce::geom_circle(data = circ_all, aes(x0 = x0, y0 = y0, r = r, fill = fill),
                         colour = B_STROKE, linewidth = 0.22) +
    # shape 8's drawn diameter is ~0.0296 in per unit of `size`, so the
    # asterisk is sized as a fraction of the lollipop it sits on.
    geom_point(data = stars, aes(x = x, y = y), shape = 8,
               size = cd / 0.0296 * 0.42, stroke = 0.40, colour = "black") +
    # rotated stack labels
    geom_text(data = labs_df,
              aes(x = x, y = y, label = lab, angle = angle, size = size),
              family = FONT_FAMILY, hjust = 0, vjust = 0.5,
              lineheight = B_LAB_LH) +
    # annotations
    geom_segment(data = bracket, aes(x = x, xend = xend, y = y, yend = yend),
                 colour = "black", linewidth = 0.4, lineend = "square") +
    geom_text(data = cap1, aes(x = x, y = y, label = lab, size = size),
              family = FONT_FAMILY, hjust = 1, vjust = 0.5, lineheight = 1.15) +
    geom_segment(data = lead, aes(x = x, xend = xend, y = y, yend = yend),
                 colour = "black", linewidth = 0.3, lineend = "square") +
    geom_text(data = cap2, aes(x = x, y = y, label = lab, size = size),
              family = FONT_FAMILY, hjust = 0, vjust = 0.5, lineheight = 1.15) +
    # legend
    ggforce::geom_circle(data = leg, aes(x0 = x0, y0 = y0, r = r, fill = fill),
                         colour = B_STROKE, linewidth = 0.22) +
    geom_text(data = leg_tx, aes(x = x, y = y, label = lab, size = size),
              family = FONT_FAMILY, hjust = 0, vjust = 0.5) +
    geom_text(data = foot, aes(x = x, y = y, label = lab, size = size),
              family = FONT_FAMILY, hjust = 0, vjust = 0.5) +
    scale_fill_identity() +
    scale_size_identity() +
    coord_cartesian(xlim = c(0, B_W_TOT), ylim = c(0, B_H_TOT),
                    expand = FALSE, clip = "off") +
    labs(title = gene_md("Mutations cluster in specific regions of\nNKAP and PRPF8")) +
    theme_publication() +
    theme(
      legend.position = "none",
      axis.line = element_blank(), axis.ticks = element_blank(),
      axis.text = element_blank(), axis.title = element_blank(),
      panel.background = element_blank(),
      plot.title = md_title(margin = margin(0, 0, 2, 0)),   # italic genes
      plot.margin = margin(B_MAR[["t"]], B_MAR[["r"]], B_MAR[["b"]],
                           B_MAR[["l"]], "pt"))

  wrap_elements(full = p)
}))
