# Figure 2, panel c. Sourced by ../compose.R, which supplies the
# shared setup (theme tokens, FLASHIER/CACHE paths, ROW_H, panel_cached,
# missing_panel, sf3b1_info). Must leave a `panel_c` object defined.
#
# Rebuild just this panel:   Rscript compose.R c

# ---- c. GSP x covariate correlation heatmap ----------------------
if (!have_cache("panel_c_matrix.tsv", "panel_c_cols.tsv",
                "panel_c_rows.tsv")) {
  panel_c <- missing_panel("c  GSP x covariate heatmap",
                           "build_cache_panel_c.R")
} else if (!identical(REBUILD, "ALL") && !("c" %in% REBUILD) &&
           file.exists(file.path(PANEL_DIR, "c.rds"))) {
  cat("[panel] c (cached)\n"); panel_c <- readRDS(file.path(PANEL_DIR, "c.rds"))
} else {
  cat("[panel] c building ...\n")
  c_mat  <- fread(file.path(CACHE, "panel_c_matrix.tsv"), data.table = FALSE)
  c_cols <- fread(file.path(CACHE, "panel_c_cols.tsv"),   data.table = FALSE)
  c_rows <- fread(file.path(CACHE, "panel_c_rows.tsv"),   data.table = FALSE)

  # Geometry of the drawn panel, needed to convert the label text's
  # point size into row units for collision spacing.
  #
  # The drawn height of the matrix is not knowable here -- patchwork
  # only settles it when the whole figure is laid out -- so these two
  # fractions of the cell are MEASURED off the render. At 400 dpi this
  # row's cell is 3.935 in, of which the tile block is 2.397 in (0.609)
  # and the x-axis label band below it 0.535 in (0.136).
  #
  # RE-MEASURE whenever this row's chrome changes -- title lines,
  # ROW_H, the x-label ANGLE (45 deg labels are much shorter than the
  # 90 deg ones they replaced, which alone moved the split by 0.13 in)
  # or the longest x-label. The scan:
  #   pdftoppm -r 400 -png -singlefile out/figure2.pdf /tmp/f2
  #   then find, in a vertical slice through the middle of the heatmap,
  #   the first/last rows that are mostly tiles and the last dark row.
  C_CELL_H_IN  <- ROW_H[2]   # ROW_H is in inches; see ../compose.R
  # RE-MEASURED after the y-expansion and title-margin fixes below grew
  # the matrix: the tile block is now 2.645 in of this 3.937 in cell
  # (was 2.400 in / 0.609).
  C_VP_H_IN    <- C_CELL_H_IN * 0.672   # tile block (the matrix itself)
  C_BELOW_IN   <- C_CELL_H_IN * 0.136   # x-axis label band under it
  C_LAB_W      <- 16   # x-units reserved for the label text
  C_LEAD_GAP   <- 2.2  # matrix edge -> label text; room for the leaders
  c_cols <- c_cols[order(c_cols$display_order), ]
  c_rows <- c_rows[order(c_rows$display_order), ]
  n_gsp  <- nrow(c_rows)

  # The cache stores raw covariate names in the row labels
  # ("GSP 3 [avg_len]"); swap in the prettified column labels so the
  # right-hand annotations read like the x-axis ("read length").
  pretty <- setNames(c_cols$label, c_cols$column)
  c_rows$row_label <- vapply(c_rows$row_label, function(lb) {
    m <- regmatches(lb, regexec("^(GSP \\d+) \\[(.*)\\]$", lb))[[1]]
    if (length(m) != 3) return(lb)
    key <- m[3]
    sub_lab <- if (key %in% names(pretty)) unname(pretty[key]) else key
    # Match the reference's shorter wording ("[read length]", not
    # "[Avg read length]") -- the long form pushed these labels into
    # the legend column.
    sub_lab <- sub("^Avg +", "", sub_lab)
    sub_lab <- sub("^junction cov$", "junction coverage", sub_lab)
    sprintf("%s [%s]", m[2], sub_lab)
  }, character(1))

  # Long form. y is flipped so GSP 1 sits at the TOP.
  c_long <- c_mat %>%
    pivot_longer(-gsp, names_to = "column", values_to = "r") %>%
    mutate(x = match(column, c_cols$column),
           y = n_gsp - match(gsp, c_rows$gsp) + 1)

  # Block boundaries + the Group strip above the tiles.
  grp_of   <- c_cols$group[match(seq_len(nrow(c_cols)), c_cols$display_order)]
  div_at   <- which(diff(match(grp_of, c("Tissue", "Technical", "Summary"))) != 0) + 0.5
  GROUP_COL <- c(Tissue = "#F1A340", Technical = "#998EC3", Summary = "#7FBF7B")
  STRIP_H   <- 0.7
  strip_y   <- n_gsp + 0.85 + STRIP_H / 2          # 0.35 row gap above matrix
  strip_top <- strip_y + STRIP_H / 2
  c_strip  <- data.frame(x = seq_len(nrow(c_cols)),
                          y = strip_y,
                          group = factor(grp_of, levels = names(GROUP_COL)))
  # Block dividers span matrix + strip only (a full-height vline ran up
  # through the block-name row and touched "Summary").
  c_div <- data.frame(x = div_at, y0 = 0.5, y1 = strip_top)

  # Block names: centred one text-height above the strip top, so they
  # clear it by a fixed 1 pt regardless of panel size.
  # One matrix row in points. Uses the MEASURED panel viewport height
  # (C_VP_H_IN), not the whole cell: the earlier estimate treated 78%
  # of the cell as plot area when the real figure is 722/1181 = 61%,
  # which made row_pt ~28% too large and set the GSP labels at
  # effectively zero leading.
  row_pt        <- C_VP_H_IN * 72 / n_gsp
  grp_half_rows <- (SIZE_BASE * 0.68) / row_pt / 2
  grp_y         <- strip_top + 1 / row_pt + grp_half_rows
  c_top_expand  <- grp_y + grp_half_rows + 0.3 - n_gsp

  # Right-hand GSP labels: each label stays LEVEL WITH ITS OWN ROW, in
  # one of two staggered columns (near / far), so a leader is a short
  # horizontal tick instead of a long diagonal (Quinn, 2026-10-05: "I
  # don't want GSP21's label to be all the way down nearby GSP35").
  #
  # Placement is a greedy 2-D search, top row first. Candidates are
  # (column, dy) with dy stepping away from the row; a candidate must
  # not overlap any placed label box, its leader must not cross any
  # placed label, and no placed leader may cross it. Cost = |dy| plus a
  # small premium for the far column, so a label only moves vertically
  # when neither column is free at its own row.
  #
  # The far column does NOT have to clear the longest near label: the
  # 2-D test only forbids overlaps that actually happen, so a short
  # near label leaves the far slot beside it usable. Units: label width
  # (in) -> x units through the panel's drawn width, derived from the
  # cell (PAGE_W, COL_W) minus margins and the legend, all measured on
  # cairo. x range and drawn width depend on each other, so iterate.
  c_lab <- c_rows %>%
    mutate(y_row = n_gsp - display_order + 1) %>%
    filter(nzchar(row_label)) %>%
    arrange(desc(y_row))
  x_edge <- nrow(c_cols) + 0.5
  LAB_PT <- SIZE_BASE * 0.62
  c_w_in <- with_cairo_metrics(vapply(c_lab$row_label, function(l)
    grid::convertWidth(grid::grobWidth(grid::textGrob(l, gp = grid::gpar(
      fontsize = LAB_PT, fontfamily = FONT_FAMILY))), "in", TRUE), 0))
  leg_w_in <- with_cairo_metrics(max(
    grid::convertWidth(grid::grobWidth(grid::textGrob("Pearson r", gp = grid::gpar(
      fontsize = SIZE_BASE * 0.75, fontfamily = FONT_FAMILY))), "in", TRUE),
    (2.2 + 1) / 25.4 + grid::convertWidth(grid::grobWidth(grid::textGrob("-0.4",
      gp = grid::gpar(fontsize = SIZE_BASE * 0.7, fontfamily = FONT_FAMILY))), "in", TRUE)))
  C_MARGIN_LR_PT <- c(22, 2)      # keep in sync with plot.margin below
  c_panel_w_in <- PAGE_W * COL_W$bc[["heatmap"]] / sum(COL_W$bc) -
    sum(C_MARGIN_LR_PT) / 72 - leg_w_in - 1 / 25.4 - 0.04   # 0.04: patchwork gutter
  C_GAP_IN   <- 0.02                     # leader -> text
  C_LEAD_IN  <- 0.045                    # matrix edge -> near column
  text_h_rows <- LAB_PT * 1.12 / row_pt  # box height incl. a little air
  C_FLOOR <- 0.5 - (C_BELOW_IN * 72 / row_pt) * 0.82
  c_top   <- n_gsp + 0.5

  # Column assignment. Labels are ~1.5 rows tall, so labels on adjacent
  # rows must alternate columns. Split the top-down label list into runs
  # of rows closer than one label height; within a run alternate
  # near/far, choosing the phase that sends the LONG labels far -- the
  # far column's offset has to clear the widest label left in the near
  # column, so this keeps the label block narrow.
  run_id <- cumsum(c(TRUE, -diff(c_lab$y_row) >= text_h_rows + 0.1))
  c_lab$pref <- NA_character_
  for (r in unique(run_id)) {
    k <- which(run_id == r)
    if (length(k) == 1) { c_lab$pref[k] <- "near"; next }
    ph <- list(ifelse(seq_along(k) %% 2 == 1, "near", "far"),
               ifelse(seq_along(k) %% 2 == 1, "far", "near"))
    wmax <- vapply(ph, function(p) max(c_w_in[k][p == "near"]), 0)
    c_lab$pref[k] <- ph[[which.min(wmax)]]
  }
  # Far offset only has to clear the near labels that actually sit
  # BESIDE a far label (within one label height of it), not the widest
  # near label anywhere -- the label block packs tight and the heatmap
  # gets the freed width (Quinn, 2026-10-05: "labels really right up
  # against each other").
  nb_w <- vapply(which(c_lab$pref == "far"), function(i) {
    k <- which(c_lab$pref == "near" & abs(c_lab$y_row - c_lab$y_row[i]) < text_h_rows)
    if (length(k)) max(c_w_in[k]) else 0
  }, 0)
  C_FAR_IN <- max(nb_w) + C_LEAD_IN + 0.02

  seg_hits_box <- function(x0, y0, x1, y1, bx) {   # segment vs box (L R B T)
    if (max(x0, x1) < bx[1] || min(x0, x1) > bx[2] ||
        max(y0, y1) < bx[3] || min(y0, y1) > bx[4]) return(FALSE)
    ts <- seq(0, 1, length.out = 25)
    xs <- x0 + ts * (x1 - x0); ys <- y0 + ts * (y1 - y0)
    any(xs >= bx[1] & xs <= bx[2] & ys >= bx[3] & ys <= bx[4])
  }
  x_max <- x_edge + 30
  for (iter in 1:4) {
    u <- (x_max - 0.5) / c_panel_w_in               # x units per inch
    x_near <- x_edge + C_LEAD_IN * u
    x_far  <- x_near + C_FAR_IN * u
    boxes <- list(); leads <- list()
    place <- data.frame(x_lab = numeric(0), y_lab = numeric(0), col = character(0))
    for (i in seq_len(nrow(c_lab))) {
      yr <- c_lab$y_row[i]; w <- c_w_in[i] * u
      best <- NULL
      for (dy in c(0, rbind(-seq(0.25, 12, 0.25), seq(0.25, 12, 0.25)))) {
        for (cl in c("near", "far")) {
          xl <- if (cl == "near") x_near else x_far
          yl <- yr + dy
          if (yl + text_h_rows / 2 > c_top + 1.5 || yl - text_h_rows / 2 < C_FLOOR) next
          bx <- c(xl - 0.15, xl + w, yl - text_h_rows / 2, yl + text_h_rows / 2)
          ld <- c(x_edge, yr, xl - C_GAP_IN * u, yl)
          ok <- all(vapply(boxes, function(o) bx[2] < o[1] || bx[1] > o[2] ||
                             bx[4] < o[3] || bx[3] > o[4], TRUE)) &&
                !any(vapply(boxes, function(o) seg_hits_box(ld[1], ld[2], ld[3], ld[4], o), TRUE)) &&
                !any(vapply(leads, function(l) seg_hits_box(l[1], l[2], l[3], l[4], bx), TRUE))
          cost <- 3 * abs(dy) + if (cl != c_lab$pref[i]) 1 else 0
          if (ok && (is.null(best) || cost < best$cost))
            best <- list(cost = cost, xl = xl, yl = yl, cl = cl, bx = bx, ld = ld)
        }
        if (!is.null(best) && best$cost <= 3 * abs(dy)) break
      }
      if (is.null(best)) stop("panel c: no slot for ", c_lab$row_label[i])
      boxes[[i]] <- best$bx; leads[[i]] <- best$ld
      place[i, ] <- list(best$xl, best$yl, best$cl)
    }
    new_max <- max(vapply(boxes, `[`, 0, 2)) + 0.05
    if (abs(new_max - x_max) < 0.2) { x_max <- new_max; break }
    x_max <- new_max
  }
  c_lab$x_lab  <- place$x_lab
  c_lab$y_lab  <- place$y_lab
  c_lab$x_edge <- x_edge
  c_lab$x_end  <- place$x_lab - C_GAP_IN * u
  cat(sprintf("[panel c] widths %.2f..%.2f in (q90 %.2f), text_h %.2f rows, row_pt %.2f\n", min(c_w_in), max(c_w_in), quantile(c_w_in, .9), text_h_rows, row_pt))
  cat(sprintf("[panel c] %d labels: %d near / %d far, max |dy| %.2f rows, x_max %.1f, panel %.2f in\n",
              nrow(c_lab), sum(place$col == "near"), sum(place$col == "far"),
              max(abs(c_lab$y_lab - c_lab$y_row)), x_max, c_panel_w_in))

  # The label block does not fit inside the matrix: 34 labels need
  # 34 * MIN_DY ~ 66 row-heights and there are only n_gsp = 50, so the
  # sweep above spills the bottom ~16 rows' worth BELOW the matrix.
  # That spill must NOT be allowed to drive the y scale -- ggplot sizes
  # the panel from every layer's data range, so the overflowing text
  # stretched the panel by ~0.6 in, leaving dead space under the matrix
  # and pushing the x-axis labels that far down with it (which is what
  # made them look "way displaced").
  #
  # Pin the panel to the matrix + strip instead and let the spill draw
  # outside it (clip = "off"). It lands in the x-axis label band, but
  # at x > x_lab, i.e. to the RIGHT of the last column -- the rotated
  # tick labels only occupy the matrix's own x range, so the two never
  # collide.
  C_YLIM <- c(0.5 - 0.6, n_gsp + 0.5 + c_top_expand)

  # Block names sit ON the strip rather than in a legend -- a legend
  # here cost ~25% of panel c's width for three words.
  c_grp_lab <- c_strip %>%
    group_by(group) %>%
    summarise(x0 = min(x), x1 = max(x), .groups = "drop") %>%
    mutate(y     = grp_y,
           wide  = (x1 - x0) >= 6,
           # Narrow blocks: left-anchor just inside the divider (0.5 is
           # the divider itself) so the text does not touch the line.
           x     = ifelse(wide, (x0 + x1) / 2, x0 - 0.15),
           hjust = ifelse(wide, 0.5, 0))

  panel_c <- wrap_elements(full =
    ggplot() +
      geom_tile(data = c_long, aes(x = x, y = y, fill = r)) +
      geom_tile(data = c_strip, aes(x = x, y = y, fill = NULL),
                fill = GROUP_COL[as.character(c_strip$group)], height = 0.7) +
      geom_text(data = c_grp_lab,
                aes(x = x, y = y, label = group, hjust = hjust),
                size = SIZE_BASE * 0.68 / .pt, family = FONT_FAMILY,
                fontface = "bold", color = "grey20") +
      geom_segment(data = c_div, aes(x = x, xend = x, y = y0, yend = y1),
                   color = "black", linewidth = 0.3) +
      geom_segment(data = c_lab,
                   aes(x = x_edge, xend = x_end,
                       y = y_row, yend = y_lab),
                   color = "grey60", linewidth = 0.16) +
      geom_text(data = c_lab, aes(x = x_lab, y = y_lab, label = row_label),
                hjust = 0, vjust = 0.5, size = SIZE_BASE * 0.62 / .pt,
                family = FONT_FAMILY, color = "grey15") +
      scale_fill_gradientn(colours = heatmap_diverging(),
                           limits = c(-0.6, 0.6), oob = scales::squish,   # was +-0.4: saturated (Quinn, 2026-10-05)
                           name = "Pearson r",
                           breaks = c(-0.6, 0, 0.6)) +
      scale_x_continuous(breaks = seq_len(nrow(c_cols)),
                         labels = c_cols$label, expand = c(0, 0)) +
      # NO y expansion. `coord_cartesian(ylim = C_YLIM)` below already
      # carries the headroom (C_YLIM is 0.5 - 0.6 .. n_gsp + 0.5 +
      # c_top_expand), and ggplot applies the SCALE expansion on top of
      # the coord limits -- so the old expansion(add = c(0.6,
      # c_top_expand)) added both bands a SECOND time. That was ~4.3
      # rows of pure air, 7% of the drawn panel: an empty band between
      # the title and the Group strip, plus slack under the x labels.
      # Zeroing it grows the matrix into that space without moving the
      # strip, the block names or the label block, all of which are
      # positioned in row units off the matrix.
      scale_y_continuous(expand = expansion(add = c(0, 0))) +
      coord_cartesian(xlim = c(0.5, x_max), ylim = C_YLIM,
                      clip = "off") +
      labs(x = NULL, y = NULL,
           title = "GSPs represent a mix of tissue-driven\nand technical effects") +
      theme_publication() +
            # Zero bottom margin on the title: an ink scan of the
            # render showed a 0.13 in blank band between the title and
            # the Group strip that was pure default title padding.
            # Only the title's BOTTOM edge moves, so c's title top --
            # which is aligned with b's by the 12 pt plot.margin below
            # -- stays exactly where it is.
      theme(plot.title      = element_text(hjust = 0.5,
                                            margin = margin(0, 0, 0, 0)),
            # Diagonal, not vertical -- see angled_x_labels() in
            # ../../theme.R. Applied after this theme() below so it
            # wins; the size comes from here.
            axis.text.x     = element_text(size = SIZE_BASE * 0.55),
            axis.text.y     = element_blank(),
            axis.ticks      = element_blank(),
            axis.line       = element_blank(),
            legend.position = "right",
            legend.key.width  = unit(2.2, "mm"),
            legend.key.height = unit(5.5, "mm"),
            legend.text     = element_text(size = SIZE_BASE * 0.7),
            legend.title    = element_text(size = SIZE_BASE * 0.75),
            legend.margin      = margin(0, 0, 0, 0),
            legend.box.spacing = unit(1, "mm"),
            # LEFT 22 pt: panel c has no y-axis text, so nothing but
            # this margin holds the matrix off the cell's left edge --
            # at 2 pt the heatmap started LEFT of panel c's own "c" tag.
            # TOP 12 pt: panel b is a patchwork inside its
            # wrap_elements (it carries an inset_element), so it gets
            # patchwork's own plot margin ABOVE its title while panel
            # c, a plain ggplot, does not. Measured at 400 dpi, c's
            # title sat 54 px (~9.7 pt) above b's. Re-measure with
            # scratch/titletop.R if either panel's wrapper changes.
            # 12 -> 9.7 pt after COL_W$bc widened c: re-measured b 67 / c 80 px.
            plot.margin     = margin(9.7, 2, 2, 22, "pt")) +
      # 65 deg, not the default 45: at ~4.4 pt per column a 45 deg label
      # has only 3.1 pt of perpendicular pitch -- less than its own cap
      # height -- and the ~45 column names overprinted into a smear (QA).
      # 70 deg plus a wider cell (COL_W$bc 1.22/1.28 -> 1.12/1.38 in
      # compose.R) gives a clear gap, and is still diagonal, not vertical.
      angled_x_labels(angle = 70, size = SIZE_BASE * 0.55))
  saveRDS(panel_c, file.path(PANEL_DIR, "c.rds"))
}

