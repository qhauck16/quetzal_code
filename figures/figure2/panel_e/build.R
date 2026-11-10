# Figure 2, panel e. Sourced by ../compose.R, which supplies the
# shared setup (theme tokens, FLASHIER/CACHE paths, ROW_H, panel_cached,
# missing_panel, sf3b1_info). Must leave a `panel_e` object defined.
#
# Rebuild just this panel:   Rscript compose.R e

# ------------------------------------------------------------------
# e. Cryptic 3'ss distance histogram + event-type pie inset
#    (the reference pasted the pie in Illustrator; composed in R here)
# ------------------------------------------------------------------
if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
fact10 <- fread(file.path(QROOT, "gsp_analysis/fact_10_top_junctions.tsv"),
                 data.table = FALSE)

e_types <- fact10 %>%
  group_by(snaptron_id) %>% slice_head(n = 1) %>% ungroup() %>%
  count(type) %>%
  mutate(prop  = n / sum(n),
         label = percent(prop, accuracy = 0.1),
         type  = factor(type, levels = names(PAL$splice_type)))

e_dist <- fact10 %>%
  filter(type == "alt 3'ss",
         distance_to_cryptic > -100, distance_to_cryptic < 100)

# ------------------------------------------------------------------
# Where the inset may live.
#
# The histogram's ink is extremely uneven: every bin is <= 4 counts
# except the peak at -25..-10, which runs to the top of the panel.
# Measured in panel fractions the peak occupies x 0.42-0.52, so the
# band above ~0.25 of panel height is free on BOTH sides of it. That
# is why the pie and the legend are two SEPARATE inset_elements
# straddling the peak (the reference figure does the same) rather than
# one plot with a right-hand legend: a single inset wide enough to hold
# pie + legend side by side necessarily runs through the peak.
#
# Two insets also decouple the legend from the pie's coord_equal, so
# neither is letterboxed by the other's aspect ratio.
# ------------------------------------------------------------------
# 2026-11-09 (Quinn): pie shrunk, % labels enlarged, legend pulled left toward the pie.
E_PIE_BOX <- c(left = 0.025, bottom = 0.52, right = 0.290, top = 1.05)
E_LEG_BOX <- c(left = 0.300, bottom = 0.56, right = 0.640, top = 1.03)

# Pie drawn in CARTESIAN coords (geom_polygon wedges + coord_equal)
# rather than coord_polar: polar coords give no control over where an
# outside label lands, and the two narrow wedges (other, skipped exon)
# are adjacent at 12 o'clock so their labels collided.
#
# Wedge order is clockwise from 12 o'clock in REVERSE level order,
# reproducing coord_polar(theta = "y") + position_stack (first level on
# top of the stack = last clockwise) and the reference figure. `thm` is
# the angular form of the stack midpoint 1 - (cumsum(prop) - prop/2):
# ONE mid-angle per wedge, shared by that wedge's leader and its text,
# so the two can never drift apart.
e_types <- e_types %>% arrange(type)
e_ord   <- e_types[nrow(e_types):1, ]                # clockwise order
e_ord   <- e_ord %>%
  mutate(th0  = c(0, head(cumsum(prop), -1)) * 2 * pi,   # clockwise, from top
         th1  = cumsum(prop) * 2 * pi,
         thm  = (th0 + th1) / 2,
         wide = prop >= 0.15)
# math angle (radians, counter-clockwise from +x) for a clockwise angle
# from 12 o'clock:  a = pi/2 - theta
cw2xy <- function(theta, r) cbind(x = r * cos(pi / 2 - theta),
                                  y = r * sin(pi / 2 - theta))

e_wedge <- do.call(rbind, lapply(seq_len(nrow(e_ord)), function(i) {
  th <- seq(e_ord$th0[i], e_ord$th1[i], length.out = 120)
  xy <- cw2xy(th, 1)
  data.frame(type = e_ord$type[i],
             x = c(0, xy[, "x"]), y = c(0, xy[, "y"]))
}))

# ---- inside labels (wide wedges) --------------------------------
# On the wedge's mid-ray. The radius follows from the wedge's angular
# span, and the two spans here pull in opposite directions:
#   * 67.5% spans 243 deg, so its mid-ray has enormous angular
#     clearance and the only binding constraint is the RIM -- a
#     horizontal label centred at 0.62 R puts its far corner at
#     ~0.91 R, visibly kissing the edge, so it is pulled in to 0.52 R
#     (far corner ~0.82 R).
#   * 18.9% spans only 68 deg, so near the apex the chord is narrower
#     than the text; the label is pushed out until its inner corner
#     clears the upper straight edge. 0.575 R balances the two: the
#     inner corner sits ~8 deg inside the edge and the outer corner at
#     ~0.91 R.
E_R_IN_WIDE   <- 0.520  # wedge > half the circle: rim-limited
E_R_IN_NARROW <- 0.575  # wedge < half the circle: chord-limited
e_in <- e_ord %>% filter(wide)
e_in <- cbind(e_in, cw2xy(e_in$thm,
                          ifelse(e_in$prop > 0.5, E_R_IN_WIDE, E_R_IN_NARROW)))

# ---- outside labels (narrow wedges) -----------------------------
# A radial leader from just off the rim out to R_ELB; the text is then
# anchored AT the leader's tip with hjust/vjust chosen so the text box
# grows strictly AWAY from the tip. That makes "leader crosses its own
# text" geometrically impossible without knowing the string's size.
#
# The anchor side is read off the mid-ray: a ray within 30 deg of
# vertical (|cos| < 0.5) gets its label stacked above/below it, any
# other ray gets it set beside, on the side the ray points. Here
# skipped exon (mid-ray 9 deg from 12 o'clock) goes straight up and
# other (34 deg) goes out to the right, which separates the two labels
# in BOTH axes even though their wedges are adjacent.
E_R_TIP <- 1.03   # leader start, just off the rim
E_R_ELB <- 1.20   # leader end = text anchor
E_R_PAD <- 0.05   # tip-to-text gap, in pie radii

e_out <- e_ord %>% filter(!wide)
e_out <- cbind(e_out,
               setNames(as.data.frame(cw2xy(e_out$thm, E_R_TIP)), c("x0", "y0")),
               setNames(as.data.frame(cw2xy(e_out$thm, E_R_ELB)), c("x1", "y1")))
e_out$upright <- abs(cos(pi / 2 - e_out$thm)) < 0.5
e_out$xlab <- with(e_out, ifelse(upright, x1, x1 + sign(x1) * E_R_PAD))
e_out$ylab <- with(e_out, ifelse(upright, y1 + sign(y1) * E_R_PAD, y1))
e_out$hj   <- with(e_out, ifelse(upright, 0.5, ifelse(x1 >= 0, 0, 1)))
e_out$vj   <- with(e_out, ifelse(upright, ifelse(y1 >= 0, 0, 1), 0.5))

# Limits: the pie plus its outside labels, so the whole inset -- text
# included -- is inside E_PIE_BOX and therefore clear of the peak.
# The headroom past each anchor is the rendered extent of the longest
# string at E_LAB_SIZE, measured off the PDF (0.54 R wide, 0.18 R tall
# at this pie size): +0.14 R to the right of the "8.4%" anchor and
# +0.19 R above the "5.2%" anchor. clip = "off" is kept so a label can
# never be chopped if a future string is longer -- but the limits, not
# clipping, are what keep it off the histogram.
E_PIE_XLIM <- c(-1.06, 1.27)
E_PIE_YLIM <- c(-1.04, 1.50)

E_LAB_SIZE <- SIZE_BASE * 0.85          # pie % labels and legend text

# The pie now uses the SHARED helper, in the style Quinn asked every pie
# to match (Figure3-NKAP_improved.pdf panel d): big wedges plain text,
# small wedges in a rounded chip filled with their own colour, sitting on
# the pie. The hand-tuned outside-label geometry above (e_in / e_out /
# E_PIE_XLIM) is no longer drawn. Wedge order is unchanged: the helper
# draws clockwise in REVERSE of levels_order, and this panel's order was
# clockwise in level order; passing the levels as-is reproduces it (verified by eye against the old render).
# E_PIE_R_IN is the drawn radius, MEASURED off the render and iterated
# to a fixed point -- see ../scripts/pie_chart.R.
source("../scripts/pie_chart.R")
E_PIE_LAB_SIZE <- SIZE_BASE * 1.0      # % labels on the wedges (were 0.85x)
E_PIE_R_IN <- 0.477   # re-measured 2026-11-09 after the shrink
panel_e_pie <- pie_labels_on_wedges(
  e_types, "type", PAL$splice_type,
  levels_order = as.character(e_types$type),
  pie_r_in = E_PIE_R_IN, label_pt = E_PIE_LAB_SIZE, family = FONT_FAMILY)

# ---- legend, extracted as its own grob --------------------------
# ggplot draws it, so key/text metrics and spacing stay ggplot's; it is
# then placed by itself and centres at its natural size inside its box.
e_leg_src <- ggplot(e_types, aes(x = 1, y = prop, fill = type)) +
  geom_col() +
  scale_fill_manual(values = PAL$splice_type, drop = FALSE, name = NULL,
                    breaks = names(PAL$splice_type)) +
  theme_void(base_size = SIZE_BASE) +
  theme(legend.text        = element_text(size = E_LAB_SIZE,
                                          family = FONT_FAMILY),
        legend.key.size    = unit(3.0, "mm"),
        legend.key.spacing.y = unit(0.6, "mm"),
        legend.margin      = margin(0, 0, 0, 0),
        legend.box.spacing = unit(0, "mm")) +
  guides(fill = guide_legend(ncol = 1))
panel_e_leg <- cowplot::get_plot_component(e_leg_src, "guide-box-right",
                                           return_all = TRUE)

panel_e_hist <- ggplot(e_dist, aes(x = -distance_to_cryptic)) +
  geom_histogram(bins = 60, fill = PAL$fill_neutral, color = "white",
                 linewidth = 0.2) +
  labs(x = "Cryptic 3'ss distance from MANE 3'ss", y = "Count",
       title = gene_md("Unannotated junctions driving GSP #10\nmatch known SF3B1 phenotype")) +
  theme_publication() +
  theme(plot.title = md_title())

panel_e <- panel_cached("e", wrap_elements(full =
  panel_e_hist +
    inset_element(panel_e_pie,
                  left   = E_PIE_BOX[["left"]],  bottom = E_PIE_BOX[["bottom"]],
                  right  = E_PIE_BOX[["right"]], top    = E_PIE_BOX[["top"]],
                  align_to = "panel") +
    inset_element(panel_e_leg,
                  left   = E_LEG_BOX[["left"]],  bottom = E_LEG_BOX[["bottom"]],
                  right  = E_LEG_BOX[["right"]], top    = E_LEG_BOX[["top"]],
                  align_to = "panel")))
