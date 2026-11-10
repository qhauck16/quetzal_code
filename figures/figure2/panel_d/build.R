# Figure 2, panel d. Sourced by ../compose.R, which supplies the
# shared setup (theme tokens, FLASHIER/CACHE paths, ROW_H, panel_cached,
# missing_panel, sf3b1_info). Must leave a `panel_d` object defined.
#
# Rebuild just this panel:   Rscript compose.R d

# ------------------------------------------------------------------
# d. GSP #10 corresponds to SF3B1 hotspot mutations
#    Left: full CDF with the zoom region boxed. Right: the zoom, as a
#    labelled dot plot. Both share PAL$sf3b1_mut.
# ------------------------------------------------------------------

# Width split between the left column (CDF over its legend) and the
# zoom. Declared up here because the label-width maths below needs it.
D_SPLIT <- c(1, 2.1)

MUT_LEVELS <- names(PAL$sf3b1_mut)
d_cum <- sf3b1_info %>%
  mutate(mutation_label = factor(plot_color, levels = MUT_LEVELS)) %>%
  arrange(logk10) %>%
  group_by(mutation_label) %>%
  mutate(cumulative = row_number() / n() * 100) %>%
  ungroup()

D_ZOOM_FROM <- -4        # the zoom panel's left edge, boxed on the CDF

# Both halves share this y range so their 0/25/50/75/100 gridlines sit
# at the same height across the panel. The headroom above 100% is set
# by the zoom's label repel (see d_ylim below); the CDF just follows.
D_YLIM <- c(-2, 110)

# Legend text, with consistent sentence case (the data's class names mix
# "Mutation Near ... Hotspot" with "Other ... Mutation"). Display only --
# the class values themselves are untouched, so colours and filters
# keyed on them still match.
D_LEG_LAB <- c("SF3B1 hotspot mutation"      = "SF3B1 hotspot mutation",
               "Mutation Near SF3B1 Hotspot" = "Mutation near SF3B1 hotspot",
               "Other SF3B1 Mutation"        = "Other SF3B1 mutation",
               "No SF3B1 Mutation"           = "No SF3B1 mutation")

# Zoom callout (same idiom as figure 3 panel a): two short lines from the
# CDF box's right corners, splaying out to the zoom panel's corners.
# Ends in the CDF's own data units, set against the render.
# MEASURED 2026-10-05 with ../../scripts/measure_zoom_connectors.py (box corners
# -> zoom panel's left axis line, in this CDF's data units). Re-run it after ANY
# layout change to this row; hand-tuned values drifted off the corners.
D_CALL_DX  <- 1.7708
D_CALL_TOP <- 109.433
D_CALL_BOT <- -83.142

panel_d_cdf <- ggplot(d_cum, aes(x = logk10, y = cumulative,
                                  color = mutation_label)) +
  geom_step(linewidth = 0.55) +
  # The box marks the region expanded in the dot plot to its right.
  annotate("rect", xmin = D_ZOOM_FROM, xmax = max(d_cum$logk10) + 0.25,
           ymin = 0, ymax = 102, fill = NA, color = "black",
           linewidth = 0.4) +
  annotate("segment", x = max(d_cum$logk10) + 0.25,
           xend = max(d_cum$logk10) + 0.25 + D_CALL_DX,
           y = c(102, 0), yend = c(D_CALL_TOP, D_CALL_BOT),
           linewidth = 0.35, colour = "black") +
  scale_color_manual(values = PAL$sf3b1_mut, drop = FALSE, name = NULL,
                     labels = function(x) gene_pm(unname(D_LEG_LAB[x]))) +  # italic genes (plotmath; see theme.R)
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                      breaks = c(0, 25, 50, 75, 100)) +
  # xlim pinned too: the callout segments reach past the data, and a
  # layer's data would otherwise stretch the panel to include them.
  coord_cartesian(xlim = range(d_cum$logk10), ylim = D_YLIM,
                  expand = FALSE, clip = "off") +
  labs(x = "log(GSP #10)", y = "Cumulative percentage") +
  theme_publication() +
  theme(legend.position = "none",
        axis.title  = element_text(size = SIZE_BASE),
        plot.margin = margin(4, 6, 4, 22, "pt"))

d_zoom <- d_cum %>%
  filter(logk10 >= D_ZOOM_FROM) %>%
  mutate(proteinChange = ifelse(proteinChange %in% c("NA", NA), "",
                                 proteinChange))

# ------------------------------------------------------------------
# Mutation labels.
#
# The reference figure places every label immediately beside its own
# point, alternating to the LEFT and RIGHT of the point cloud, so the
# labels track the diagonal in two hugging columns with leaders only a
# few points long. Reproduce that instead of parking labels in fixed
# columns (which is what produced the long crossing leaders):
#
#   * side alternates along each class's own rank in y, so consecutive
#     points -- which sit ~2% of the axis apart -- never compete for the
#     same slot; each side then has to resolve half the density.
#   * `direction = "y"` pins a label's x to its point's x plus one small
#     gap, so a leader can only ever be `gap` wide horizontally.
#   * the red and light-blue classes share ONE repel layer: their point
#     clouds interleave near the top of the zoom, and only a single
#     layer can guarantee no cross-class overlap. Green sits alone in
#     the empty top-left corner and gets its own layer.
#
# Per-point side/gap/hjust are carried as DATA COLUMNS (the panel is
# serialised by panel_cached(), so an aes() must never close over a
# local scalar).
# ------------------------------------------------------------------

# On-page pt size of a mutation label. 52 hotspot labels alternate into
# two columns, so each column needs 26 x 1.15 x D_LAB_PT of vertical
# room; row 3 was grown (see ../compose.R) until that fits with slack,
# which is what lets these be legible rather than 4.8 pt.
D_LAB_PT <- SIZE_BASE * 0.65

# Geometry of the zoom's drawing area, used to convert label widths
# into data units: the page width x panel d's share of row 3 x the
# zoom's share of panel d, less ~0.7 in of y axis and plot margins.
D_ZOOM_W_IN <- PAGE_W * (COL_W$de[["cdf"]] / sum(COL_W$de)) *
  (D_SPLIT[2] / sum(D_SPLIT)) - 0.70

d_x_lo  <- min(d_zoom$logk10)
d_x_hi  <- max(d_zoom$logk10)
d_x_rng <- d_x_hi - d_x_lo

# One label character is ~0.6 em of Helvetica; express that in x units.
D_CHAR_X <- 0.60 * (D_LAB_PT / 72) * d_x_rng / D_ZOOM_W_IN
# Leader gap: the point marker's radius plus a visible hairline of
# white space, so a label never sits on the point ribbon it labels.
# The corridor between the left- and right-hand label columns is 2x
# this gap, so it has to clear the point ribbon AND keep the two
# columns from touching where the diagonal is steepest ("K700E|K700E"
# printed flush at 2.5 pt).
D_PT_SIZE <- 0.9                                   # geom_point size below
D_GAP_X   <- (D_PT_SIZE * .pt / 2 + 4.5) / 72 * d_x_rng / D_ZOOM_W_IN

# Side assignment. Red (hotspot) and blue (near-hotspot) share one
# repel layer because their clouds interleave, so alternating BOTH
# classes left/right let a blue label land inside the dense red column
# near the 100% row -- "H662R" was printed straight through by
# "R625H". Blue is pinned to the LEFT flank instead: its points sit to
# the left of the red cloud at every y, so the two label columns can
# no longer occupy the same x band. Red keeps alternating.
D_LEFT_ONLY <- "Mutation Near SF3B1 Hotspot"
# Above this cumulative %, red labels go RIGHT only: that is where the
# blue near-hotspot labels occupy the left flank, and a left-flank red
# label printed flush into "N626Y" (QA review).
D_RED_RIGHT_ABOVE <- 101   # i.e. OFF: routing the top red labels right (86) just
                            # moved the overprints to the right flank.
d_lab <- d_zoom %>%
  filter(nzchar(proteinChange)) %>%
  arrange(mutation_label, cumulative, logk10) %>%
  group_by(mutation_label) %>%
  mutate(side = ifelse(as.character(mutation_label) == D_LEFT_ONLY, -1,
                ifelse(cumulative > D_RED_RIGHT_ABOVE, 1,
                       ifelse(row_number() %% 2 == 1, -1, 1)))) %>%
  ungroup() %>%
  mutate(nx = side * D_GAP_X,
         hj = ifelse(side < 0, 1, 0),
         # x reached by this label's far edge, for the axis gutters
         reach = logk10 + side * (D_GAP_X + nchar(proteinChange) * D_CHAR_X))

# Axis range: exactly the point cloud plus whatever the outermost
# label on each flank needs. No dead space, no clipped labels.
d_lim_lo <- min(d_lab$reach, d_x_lo) - D_GAP_X
d_lim_hi <- max(d_lab$reach, d_x_hi) + D_GAP_X
# Headroom above the 100% row. EVERY class tops out at exactly 100%
# (it is a cumulative percentage), so the 100% row is the single most
# crowded band in the panel -- and the band where the blue and red
# clouds are closest in x (H662R at -0.18 vs R625H at 0.00). With the
# old 105/102.5 ceiling those labels had 2.5 y-units to move into and
# a solid red column below, so repel simply gave up and overprinted.
# The extra headroom is empty plot area, not extra axis: breaks still
# stop at 100%.
d_ylim   <- c(-2, 110)

# The band from D_TOP_BAND to 100 is where the unlabelled "No SF3B1
# Mutation" cloud lives: a solid row of navy points with its step line
# running the full width of the panel. ggrepel cannot see that layer
# (it only avoids its OWN data), so any label left in the band gets a
# line drawn straight through it -- this is what struck out "H662R".
# Labels are therefore partitioned into two repel layers with DISJOINT
# y ranges: everything at or above the band is lifted into the
# headroom, everything below is capped under it. Disjoint ranges also
# mean the two layers cannot collide with each other.
D_TOP_BAND   <- 97
d_repel_ylim <- c(-2, D_TOP_BAND - 1)      # below the navy band
d_top_ylim   <- c(101.5, 108)              # above the 100% line

# Vertical counterpart of D_GAP_X: one label's height plus the marker,
# in y units. The zoom's drawing area is row 3 of the figure
# (ROW_H[3] inches) less the panel title, the x axis and margins.
D_ZOOM_H_IN <- ROW_H[3] - 0.85
D_ROW_CLEAR <- (D_LAB_PT + D_PT_SIZE * .pt / 2 + 1) / 72 *
  diff(d_ylim) / D_ZOOM_H_IN

# Three label layers: above the navy band (all classes), and below it
# split by class (red+blue share one layer because their clouds
# interleave in x; green sits alone on the left).
d_lab_top  <- filter(d_lab, cumulative >= D_TOP_BAND)
d_lab_low  <- filter(d_lab, cumulative <  D_TOP_BAND)
d_lab_diag <- filter(d_lab_low, mutation_label != "Other SF3B1 Mutation")
d_lab_oth  <- filter(d_lab_low, mutation_label == "Other SF3B1 Mutation")

# ggrepel only avoids ITS OWN data, so the step lines and unlabelled
# points were invisible to it and labels landed on them ("N626Y" sat on
# the blue step; Quinn, 2026-10-05). Sample every class's step path
# (geom_step "hv": horizontal run, then vertical rise) plus every point
# into blank-label obstacle rows appended to each repel layer.
d_step_obs <- d_zoom %>%
  arrange(mutation_label, logk10, cumulative) %>%
  group_by(mutation_label) %>%
  group_modify(function(g, k) {
    if (nrow(g) < 2) return(tibble(logk10 = g$logk10, cumulative = g$cumulative))
    bind_rows(lapply(seq_len(nrow(g) - 1), function(i) {
      x0 <- g$logk10[i]; x1 <- g$logk10[i + 1]
      y0 <- g$cumulative[i]; y1 <- g$cumulative[i + 1]
      nh <- max(2, ceiling((x1 - x0) / (D_CHAR_X * 1.5)))
      nv <- max(2, ceiling((y1 - y0) / 1.5))
      tibble(logk10     = c(seq(x0, x1, length.out = nh), rep(x1, nv)),
             cumulative = c(rep(y0, nh), seq(y0, y1, length.out = nv)))
    }))
  }) %>% ungroup() %>%
  bind_rows(d_zoom %>% dplyr::select(mutation_label, logk10, cumulative)) %>%
  distinct(mutation_label, logk10 = round(logk10, 3), cumulative = round(cumulative, 2)) %>%
  mutate(proteinChange = "", nx = 0, hj = 0.5)
with_obs <- function(lab, ylim) bind_rows(lab,
  filter(d_step_obs, cumulative >= ylim[1], cumulative <= ylim[2]))
# Near-hotspot (blue) labels get their OWN layer: sharing the red
# layer's vertical-only repel squeezed the top blue cluster (N626D,
# N626Y x2) into each other and onto the blue step. The blue layer may
# move both ways into the empty band up-left of the staircase, but its
# xlim stops at the left edge of the red labels' left column over the
# blue labels' y span, so the two layers cannot meet.
d_lab_blue <- filter(d_lab_diag, mutation_label == D_LEFT_ONLY)
d_lab_diag <- filter(d_lab_diag, mutation_label != D_LEFT_ONLY)
d_red_left <- filter(d_lab_diag, side < 0,
                     cumulative >= min(d_lab_blue$cumulative) - D_ROW_CLEAR)
D_BLUE_XMAX <- min(d_red_left$reach) - D_GAP_X
# Bottom-most K666N: its default left nudge put the leader straight
# through the G740V point just above-left of it (Quinn, 2026-10-05).
# Push it further left so the leader clears G740V.
k_low <- with(d_lab_diag, which(proteinChange == "K666N")[which.min(cumulative[proteinChange == "K666N"])])
if (length(k_low)) d_lab_diag$nx[k_low] <- -4 * D_GAP_X
d_lab_diag <- with_obs(d_lab_diag, d_repel_ylim)
# Deterministic placement (repel put the two N626Y labels on one spot):
# label right edge just left of its point (capped at D_BLUE_XMAX), at the
# point's own height; working bottom-up, a label that would overlap the
# one below in x is lifted just clear of it -- the staircase rises to the
# right, so the space up-left of each point is empty. Widths measured on
# cairo; heights from the label pt size.
D_BLUE_DY <- (D_LAB_PT * 1.15) / 72 * diff(d_ylim) / D_ZOOM_H_IN
# A blue point within one label height of the navy band (N626Y at 96.6%)
# has no room below it: its label joins the headroom layer above 100%.
d_lab_top  <- bind_rows(d_lab_top, filter(d_lab_blue, cumulative >= D_TOP_BAND - D_BLUE_DY / 2))
d_lab_blue <- d_lab_blue %>% filter(cumulative < D_TOP_BAND - D_BLUE_DY / 2) %>%
  arrange(cumulative, desc(logk10))
d_lab_blue$w <- with_cairo_metrics(vapply(d_lab_blue$proteinChange, function(l)
  grid::convertWidth(grid::grobWidth(grid::textGrob(l, gp = grid::gpar(
    fontsize = D_LAB_PT, fontfamily = FONT_FAMILY))), "in", TRUE), 0)) *
  d_x_rng / D_ZOOM_W_IN
# Cap each label at the left edge of the red left-column labels near
# ITS OWN height (+-2 label heights), not the whole blue span.
d_blue_cap <- vapply(d_lab_blue$cumulative, function(y) {
  r <- d_red_left$reach[abs(d_red_left$cumulative - y) < 2 * D_BLUE_DY]
  if (length(r)) min(r) - D_GAP_X else Inf
}, 0)
d_lab_blue$x1 <- pmin(d_lab_blue$logk10 - D_GAP_X, d_blue_cap)
d_lab_blue$x0 <- d_lab_blue$x1 - d_lab_blue$w
d_lab_blue$yl <- d_lab_blue$cumulative
for (i in seq_len(nrow(d_lab_blue))[-1]) for (j in seq_len(i - 1)) {
  xo <- d_lab_blue$x0[i] < d_lab_blue$x1[j] + D_GAP_X / 2 &&
        d_lab_blue$x0[j] < d_lab_blue$x1[i] + D_GAP_X / 2
  if (xo && d_lab_blue$yl[i] < d_lab_blue$yl[j] + D_BLUE_DY)
    d_lab_blue$yl[i] <- d_lab_blue$yl[j] + D_BLUE_DY
}
stopifnot(max(d_lab_blue$yl) < D_TOP_BAND - 1)
d_lab_oth  <- with_obs(d_lab_oth,  d_repel_ylim)
cat(sprintf("[panel d] %d step/point obstacles\n", nrow(d_step_obs)))

REPEL_COMMON <- list(size = D_LAB_PT / .pt, seed = 20260922,
                     # NB max.time is a WALL-CLOCK budget: leave it
                     # effectively unbounded or the layout depends on how
                     # busy the machine is and the seed buys nothing.
                     max.overlaps = Inf, max.iter = 50000, max.time = 1e6,
                     min.segment.length = 0, segment.size = 0.12,
                     segment.color = "grey60", segment.alpha = 0.9,
                     point.padding = unit(1, "pt"), force_pull = 0.4,
                     family = FONT_FAMILY, show.legend = FALSE,
                     xlim = c(d_lim_lo, d_lim_hi))

panel_d_zoom <- ggplot(d_zoom, aes(x = logk10, y = cumulative,
                                    color = mutation_label)) +
  geom_step(linewidth = 0.35, alpha = 0.5) +
  geom_point(size = D_PT_SIZE) +
  # Hotspot (red): vertical-only repel so each label stays in its
  # point's own column. Blue points/steps are obstacles in this layer.
  do.call(geom_text_repel, c(list(
    data = d_lab_diag,
    mapping = aes(label = proteinChange, hjust = hj),
    nudge_x = d_lab_diag$nx, direction = "y",
    ylim = d_repel_ylim,
    # box.padding 0.02 -> 0.07: ggrepel measures text with narrower
    # metrics than cairo draws it, so boxes that "just touch" in its model
    # overlap on the page -- "N626Y" printed flush into "R625H" (QA).
    force = 2.5, box.padding = 0.07), REPEL_COMMON)) +
  # Near-hotspot (blue): placed deterministically (see d_lab_blue).
  geom_segment(data = d_lab_blue, inherit.aes = FALSE,
               aes(x = logk10, y = cumulative, xend = x1 + D_GAP_X * 0.3, yend = yl),
               color = "grey60", linewidth = 0.12, alpha = 0.9) +
  geom_text(data = d_lab_blue, aes(x = x1, y = yl, label = proteinChange),
            hjust = 1, size = D_LAB_PT / .pt, family = FONT_FAMILY,
            show.legend = FALSE) +
  # Other SF3B1 mutations (green) sit alone on the left, below the band.
  do.call(geom_text_repel, c(list(
    data = d_lab_oth,
    mapping = aes(label = proteinChange, hjust = hj),
    nudge_x = d_lab_oth$nx, nudge_y = -D_ROW_CLEAR,
    ylim = d_repel_ylim, direction = "both",
    force = 0.8, box.padding = 0.12), REPEL_COMMON)) +
  # Everything sitting in the navy band, lifted into the headroom above
  # the 100% line. Six labels across the full panel width, so "both"
  # lets them settle sideways into their own gaps.
  do.call(geom_text_repel, c(list(
    data = d_lab_top,
    mapping = aes(label = proteinChange, hjust = hj),
    nudge_x = d_lab_top$nx,
    nudge_y = mean(d_top_ylim) - d_lab_top$cumulative,
    ylim = d_top_ylim, direction = "both",
    force = 1.2, box.padding = 0.12), REPEL_COMMON)) +
  scale_color_manual(values = PAL$sf3b1_mut, drop = FALSE, name = NULL,
                     labels = function(x) gene_pm(x)) +  # italic genes (plotmath; see theme.R)
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                      breaks = c(0, 25, 50, 75, 100)) +
  coord_cartesian(xlim = c(d_lim_lo, d_lim_hi), ylim = d_ylim,
                   expand = FALSE, clip = "off") +
  labs(x = "log(GSP #10)", y = NULL) +
  theme_publication() +
  theme(legend.position = "none",
        axis.title  = element_text(size = SIZE_BASE),
        # No y tick text on the zoom, as in figure 3 panel a: the zoom
        # callout's lower line runs to this panel's bottom-left corner
        # and cut through "25%" / "0%" when they were printed here.
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(6, 6, 4, 4, "pt"))

# Shared legend drawn once, below both halves.
panel_d_legend <- with_cairo_metrics(
  cowplot::get_plot_component(
  panel_d_cdf +
    theme(legend.position  = "bottom",
          legend.key.size   = unit(2.6, "mm"),
          legend.text       = element_text(size = SIZE_BASE * 0.9, family = FONT_FAMILY),
          legend.box.spacing = unit(0, "mm")) +
    # One entry per row, stacked under the CDF (as in the reference) --
    # a 2x2 grid ran the long class names into each other.
    guides(color = guide_legend(ncol = 1,
                                 override.aes = list(linewidth = 1.2))),
  "guide-box-bottom", return_all = TRUE))

# Left column = CDF above its legend; right column = the zoom.
panel_d_left <- (panel_d_cdf / wrap_elements(full = panel_d_legend)) +
  plot_layout(heights = c(1, 0.42))

panel_d <- panel_cached("d", wrap_elements(full =
  ((panel_d_left | panel_d_zoom) + plot_layout(widths = D_SPLIT)) +
    # Two lines, like every other panel title in the figure: a
    # one-line title here left panel d's content sitting 39 px lower
    # than panel e's across the same row.
    plot_annotation(title = gene_md("GSP #10 corresponds to\nSF3B1 hotspot mutations"),
                    # Panel d is the only panel titled via
                    # plot_annotation() rather than labs(), so its top
                    # spacing has to be set by hand to match the rest
                    # of the figure. patchwork's default 5.5 pt plot
                    # margin PLUS a title band put it 39 px too low;
                    # zeroing both put it 51 px too high (measured at
                    # 400 dpi against panel e). 9 pt on the title lands
                    # it on panel e's title line.
                    theme = theme(plot.title = ggtext::element_markdown(halign = 0.5, lineheight = 0.9,
                      size = SIZE_TITLE, face = "bold", hjust = 0.5,
                      family = FONT_FAMILY,
                      margin = margin(9, 0, 2, 0)),
                      plot.margin = margin(0, 0, 0, 0)))))
