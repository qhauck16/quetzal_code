# Supplementary figure: GSP #10 comparison to Alsafadi et al. (2016).
#   a  exact-coordinate overlap of the GSP #10 alt 3'ss driving junctions with
#      Alsafadi's SF3B1-MUT vs WT uveal-melanoma junctions, up vs down; exact
#      binomial test of up vs down among the hits (p0 = 0.5).
#   b  Venn: Alsafadi significant upregulated alt 3'ss vs the same set
# Inputs from cache/build_overlap_cache.R.
source("common.R")
cells <- fread(file.path("cache", "s6_alsafadi_cells.tsv"), data.table = FALSE) %>%
  filter(set == "gsp")
n_up <- cells$n_overlap[cells$direction == "upregulated"]
n_dn <- cells$n_overlap[cells$direction == "downregulated"]
binom_p <- binom.test(n_up, n_up + n_dn, p = 0.5)$p.value
cells$fisher_p <- binom_p   # read only by overlap_heatmap()'s caption, nulled below
vn    <- fread(file.path("cache", "s6_venn.tsv"), data.table = FALSE)

p_a <- overlap_heatmap(cells, "GSP #10 alt 3'ss in top 250\ndriving unannotated junctions",
                       xlab = "Alsafadi et al. (2016) direction", pct_digits = 1) +
  # p-value as a second line of the x title, centred under the tiles. No
  # coord_fixed: square tiles are too narrow for "upregulated" /
  # "downregulated" at this height.
  labs(caption = NULL,
       x = sprintf("Alsafadi et al. (2016) direction\nBinomial test, up vs down: p = %.1e",
                   binom_p)) +
  theme(legend.key.height = unit(4, "mm"),
        axis.title.x = element_text(lineheight = 1.1))

r <- 1; sep <- 1.15
TOP_DX <- 0.62   # titles sit over the OUTER part of each circle so the two do not collide
circ <- data.frame(x0 = c(-sep / 2, sep / 2), y0 = 0, r = r, set = c("a", "b"))
nums <- data.frame(x = c(-sep / 2 - 0.55 * r, 0, sep / 2 + 0.55 * r), y = 0,
                   lab = c(vn$only_alsafadi, vn$shared, vn$only_gsp))
# Three-line titles anchored at the OUTER edges of the pair (hjust 0 / 1):
# two-line versions are each wider than half the pair and collide.
tops <- data.frame(x = c(-sep / 2 - r, sep / 2 + r), y = r + 0.12,
                   lab = c("Alsafadi et al.\nupregulated\nalt 3'ss",
                           "GSP #10 alt 3'ss\nin top 250 driving\nunannotated junctions"))
ns <- data.frame(x = c(-sep / 2, sep / 2), y = -r - 0.12,
                 lab = sprintf("(n = %d)", c(vn$only_alsafadi + vn$shared, vn$only_gsp + vn$shared)))
p_b <- ggplot() +
  ggforce::geom_circle(data = circ, aes(x0 = x0, y0 = y0, r = r, fill = set),
                       alpha = 0.45, colour = "grey20", linewidth = 0.4) +
  geom_text(data = nums, aes(x, y, label = lab), size = SIZE_BASE * 1.2 / .pt,
            family = FONT_FAMILY) +
  geom_text(data = tops, aes(x, y, label = lab, hjust = c(0, 1)), vjust = 0, lineheight = 0.95,
            size = SIZE_BASE / .pt, family = FONT_FAMILY) +
  geom_text(data = ns, aes(x, y, label = lab), vjust = 1,
            size = SIZE_BASE / .pt, family = FONT_FAMILY) +
  scale_fill_manual(values = c(a = "#4C72B0", b = "#D62728"), guide = "none") +
  coord_equal(clip = "off") + theme_void() +
  theme(plot.margin = margin(34, 10, 4, 10))

fig <- (p_a | p_b) + plot_layout(widths = c(1.1, 1)) + plot_annotation_publication() &
  # theme.R's annotation theme does not reach per-plot tags
  theme(plot.tag = element_text(size = SIZE_PANEL_LABEL, face = "bold", family = FONT_FAMILY))
save_supp(fig, "figS6", 8.5, 3.0)
