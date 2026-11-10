# Supplementary figure 15: NKAP extended sashimi plots (top 15 NKAP-GSP
# cryptic alt-3'ss events), 5 rows x 3 cols, no panel letters.
#   Each block: mean per-bp coverage of the top-15 NKAP-GSP samples (blue)
#   and 150 matched controls (grey) over a 150-nt window (50 nt exon +
#   100 nt intron around the annotated 3'ss), sashimi arcs from the
#   cryptic and annotated acceptors toward the off-panel donor (width ~
#   sqrt(group-mean PSI), label = PSI %), a zoom trapezoid onto the MANE
#   gene model (green star = start codon, red star = stop codon).
# Data: cache/s14/<GENE>.rds, built by cache/build_s14_cache.R from
#   unified_results nkap 20v60 df_final + recount3 bigwigs (events.tsv, recount3
#   bigwigs, leafcutter PSI tables, MANE v1.4). Redraw of that dir's
#   custom_plot.R + assemble_metaplot.R.
source("common.R")

CACHE <- file.path(FIG_ROOT, "supplementary", "cache", "s14")
GENES <- data.table::fread(file.path(CACHE, "events.tsv"))$gene   # rank order, from build_s14_cache.R

COND_COL <- "#4C7FBF"; CTRL_COL <- "grey55"; HL_COL <- "#B0C4FF"
PAGE_W <- 8.5; PAGE_H <- 11
LAB_PT <- SIZE_BASE * 0.8            # PSI labels / small annotation text
HEIGHTS <- c(1, 1, 0.28, 0.62)       # cond, ctrl, trapezoid, gene model
# Drawn track-panel size (in), measured off the render; used to convert
# label extents to data units for the collision search below.
TRACK_W_IN <- 2.25; TRACK_H_IN <- 0.50
HEAD <- 1.45                          # panel top = HEAD * ymax
ARC_END <- 1.25                       # upper arc leaves the panel edge at ARC_END * ymax
ARC_GAP <- 0.40                       # lower arc exits ARC_GAP * ymax below it

# ---- PSI label placement --------------------------------------------
# Each label sits ON its own arc. Search arc fractions nearest the middle
# of the visible arc for a pair whose boxes miss each other, the coverage
# bars beneath, and the panel top.
place_labels <- function(arcs, cov, x_rng, ymax) {
  y_top <- HEAD * ymax
  w <- (0.07 + 0.056 * 4) / TRACK_W_IN * diff(x_rng)   # ~"38.4" + padding
  h <- (LAB_PT / 72 + 0.05) / TRACK_H_IN * y_top
  box <- function(a, f) {
    x <- a$x0 + f * (a$x1 - a$x0); y <- a$y0 + (a$y1 - a$y0) * sin(pi / 2 * f)
    c(x = x, y = y, l = x - w / 2, r = x + w / 2, b = y - h / 2, t = y + h / 2)
  }
  ok1 <- function(b) {
    under <- cov$cov[cov$pos >= b["l"] & cov$pos <= b["r"]]
    b["t"] <= y_top && b["b"] >= max(c(0, under)) + 0.02 * y_top &&
      b["l"] >= x_rng[1] && b["r"] <= x_rng[2]
  }
  gx <- 0.05 / TRACK_W_IN * diff(x_rng); gy <- 0.03 / TRACK_H_IN * y_top
  overl <- function(b1, b2) !(b1["r"] + gx < b2["l"] || b2["r"] + gx < b1["l"] ||
                              b1["t"] + gy < b2["b"] || b2["t"] + gy < b1["b"])
  fr <- seq(0.15, 0.85, by = 0.05)
  cand <- expand.grid(f1 = fr, f2 = fr)
  cand$cost <- abs(cand$f1 - 0.5) + abs(cand$f2 - 0.5)
  cand <- cand[order(cand$cost), ]
  # y of arc `a` at genomic x (NA outside its span)
  arc_y <- function(a, x) {
    f <- (x - a$x0) / (a$x1 - a$x0)
    ifelse(f < 0 | f > 1, NA, a$y0 + (a$y1 - a$y0) * sin(pi / 2 * f))
  }
  # Pass 1 also keeps each label off the OTHER arc, so a label cannot
  # read as belonging to the wrong junction; pass 2 drops that.
  clear_of <- function(b, a) {
    y <- arc_y(a, seq(b["l"], b["r"], length.out = 7)); y <- y[!is.na(y)]
    !length(y) || all(y < b["b"] - gy | y > b["t"] + gy)
  }
  for (strict in c(TRUE, FALSE)) for (i in seq_len(nrow(cand))) {
    b1 <- box(arcs[[1]], cand$f1[i]); b2 <- box(arcs[[2]], cand$f2[i])
    if (ok1(b1) && ok1(b2) && !overl(b1, b2) &&
        (!strict || (clear_of(b1, arcs[[2]]) && clear_of(b2, arcs[[1]]))))
      { if (!strict) message("  relaxed label placement: ", CUR_TRACK); return(rbind(b1, b2)) }
  }
  warning("no collision-free label placement; using f = 0.35 / 0.75")
  rbind(box(arcs[[1]], 0.35), box(arcs[[2]], 0.75))
}

coverage_track <- function(d, grp, ymax, title = NULL, show_ylab = TRUE) {
  CUR_TRACK <<- paste(d$gene, grp)
  cov <- dplyr::filter(d$cov, group == grp)
  x_rng <- d$win
  edge_x <- if (d$donor > mean(x_rng)) x_rng[2] else x_rng[1]
  cov_at <- function(p) { v <- cov$cov[match(p, cov$pos)]; ifelse(is.na(v), 0, v) }
  arcs <- lapply(seq_len(nrow(d$psi)), function(i) {
    a <- d$psi$acceptor[i]
    list(x0 = a, x1 = edge_x, y0 = cov_at(a), y1 = ARC_END * ymax,
         psi = d$psi[[grp]][i])
  })
  # Arcs leave the panel edge at staggered heights (the one starting
  # higher exits higher) so they run parallel instead of converging on
  # one point, leaving room for a label on each.
  lo_arc <- which.min(sapply(arcs, `[[`, "y0"))
  arcs[[lo_arc]]$y1 <- (ARC_END - ARC_GAP) * ymax
  arc_df <- bind_rows(lapply(seq_along(arcs), function(i) {
    a <- arcs[[i]]; f <- seq(0, 1, length.out = 100)
    data.frame(x = a$x0 + f * (a$x1 - a$x0),
               y = a$y0 + (a$y1 - a$y0) * sin(pi / 2 * f),
               grp = i, lw = max(0.2, sqrt(a$psi) / 11))
  }))
  lb <- as.data.frame(place_labels(arcs, cov, x_rng, ymax))
  lb$lab <- sprintf("%.1f", sapply(arcs, `[[`, "psi"))
  ylab <- if (grp == "cond") "Top 15\nNKAP" else "150\ncontrols"

  p <- ggplot(cov, aes(pos, cov)) +
    geom_col(fill = if (grp == "cond") COND_COL else CTRL_COL,
             colour = NA, width = 1) +
    geom_path(data = arc_df, aes(x, y, group = grp, linewidth = lw),
              colour = "black", lineend = "round", inherit.aes = FALSE) +
    scale_linewidth_identity() +
    geom_label(data = lb, aes(x, y, label = lab), inherit.aes = FALSE,
               fill = "white", label.size = 0,
               label.padding = unit(0.6, "pt"), size = LAB_PT / .pt,
               family = FONT_FAMILY) +
    coord_cartesian(xlim = x_rng, ylim = c(0, HEAD * ymax),
                    expand = FALSE, clip = "off") +
    scale_y_continuous(breaks = c(0, ymax), labels = c("0", format(ymax)),
                       guide = guide_axis(cap = "both")) +
    scale_x_continuous(breaks = NULL) +
    labs(x = NULL, y = ylab, title = title) +
    theme_publication() +
    theme(axis.line.x = element_blank(),
          # Group names only on the left column. Other columns keep the
          # title but transparent, so every track has the same drawn width
          # (TRACK_W_IN, which the label search is calibrated to).
          axis.title.y = element_text(size = LAB_PT, face = "plain", lineheight = 0.9,
                                      angle = 90, margin = margin(r = 2),
                                      colour = if (show_ylab) "grey10" else "transparent"),
          axis.text.y = element_text(size = LAB_PT),
          plot.title = md_title(hjust = 0.5, margin = margin(b = 3)),
          plot.margin = margin(1, 4, 1, 2))
  p
}

gene_model <- function(d) {
  ex <- d$exons[order(d$exons$start), ]
  lo <- min(ex$start); hi <- max(ex$end); span <- hi - lo
  # Window highlight widened to >= 1.5% of the gene so it stays visible.
  mid <- mean(d$win); hw <- max(diff(d$win), 0.015 * span) / 2
  hl <- c(max(lo, mid - hw), min(hi, mid + hw))
  intr <- data.frame(s = head(ex$end, -1), e = tail(ex$start, -1)) %>%
    dplyr::filter(e - s > 0.04 * span) %>%
    mutate(m = (s + e) / 2, h = 0.012 * span,
           x = if (d$strand == "-") m + h else m - h,
           xend = if (d$strand == "-") m - h else m + h)
  cod <- d$codons
  cod$col <- ifelse(cod$type == "start_codon", "forestgreen", "red")
  brks <- pretty(c(lo, hi), n = 3); brks <- brks[brks > lo & brks < hi]
  g <- ggplot() +
    annotate("rect", xmin = hl[1], xmax = hl[2], ymin = -0.4, ymax = 1,
             fill = HL_COL, alpha = 0.55) +
    annotate("segment", x = lo, xend = hi, y = 0, yend = 0,
             colour = "grey30", linewidth = 0.35) +
    geom_segment(data = intr, aes(x = x, xend = xend, y = 0, yend = 0),
                 arrow = arrow(length = unit(0.035, "in"), type = "open"),
                 colour = "grey30", linewidth = 0.3) +
    geom_rect(data = ex, aes(xmin = start, xmax = pmax(end, start + 0.002 * span),
                             ymin = -0.25, ymax = 0.25),
              fill = "grey30", colour = NA) +
    geom_point(data = cod, aes(pos, 0.6, colour = col), shape = 8,
               size = 1.6, stroke = 0.7) +
    scale_colour_identity() +
    annotate("text", x = lo, y = -0.8, label = d$chr, hjust = 0,
             size = LAB_PT / .pt, family = FONT_FAMILY) +
    annotate("text", x = hi, y = -0.8, label = d$gene, hjust = 1,
             size = SIZE_BASE / .pt, family = FONT_FAMILY, fontface = "italic") +
    coord_cartesian(xlim = c(lo, hi), ylim = c(-1.2, 1), expand = FALSE,
                    clip = "off") +
    scale_x_continuous(breaks = brks, labels = scales::comma) +
    scale_y_continuous(breaks = NULL) +
    labs(x = NULL, y = NULL) +
    theme_publication() +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = LAB_PT),
          plot.margin = margin(0, 4, 2, 2)) +
    angled_x_labels()
  trap <- ggplot() +
    annotate("polygon", x = c(lo, hi, hl[2], hl[1]), y = c(1, 1, 0, 0),
             fill = HL_COL, alpha = 0.55) +
    coord_cartesian(xlim = c(lo, hi), ylim = c(0, 1), expand = FALSE,
                    clip = "off") +
    theme_void() + theme(plot.margin = margin(0, 4, 0, 2))
  list(trap = trap, gene = g)
}

gene_block <- function(g, left = TRUE) {
  d <- readRDS(file.path(CACHE, paste0(g, ".rds")))
  # Shared y-scale within a gene; tick at a 2-significant-figure max.
  ymax <- signif(max(d$cov$cov), 2)
  ymax <- if (ymax < max(d$cov$cov)) signif(max(d$cov$cov) * 1.05, 2) else ymax
  ttl <- gene_md(sprintf("%s | %s nt | %s", g,
                         sub("-", "−", sprintf("%+d", d$shift)), d$region),
                 genes = g)
  gm <- gene_model(d)
  coverage_track(d, "cond", ymax, ttl, left) / coverage_track(d, "ctrl", ymax, NULL, left) /
    gm$trap / gm$gene + plot_layout(heights = HEIGHTS)
}

# Key for the gene-model stars, one line under the grid.
key <- ggplot() +
  geom_point(aes(x = c(0, 0.16), y = 0), shape = 8, size = 1.6, stroke = 0.6,
             colour = c("forestgreen", "red")) +
  annotate("text", x = c(0.012, 0.172), y = 0, hjust = 0, size = LAB_PT / .pt,
           family = FONT_FAMILY, label = c("Start codon", "Stop codon")) +
  coord_cartesian(xlim = c(-0.02, 1), ylim = c(-1, 1), expand = FALSE) +
  theme_void()

fig <- with_cairo_metrics(
  wrap_plots(lapply(seq_along(GENES), function(i) gene_block(GENES[i], left = i %% 3 == 1)),
             ncol = 3, byrow = TRUE) / key + plot_layout(heights = c(1, 0.012)))
save_supp(fig, "figS15", PAGE_W, PAGE_H)
