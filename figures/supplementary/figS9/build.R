# Supplementary figure: overlap of splice-factor GSPs and LeafCutter
# results. Rows GSP #10 / #13 / #30, columns LeafCutter top-N vs 3N matched-control
# runs 10v30 / 20v60 / 30v90. One tile row per panel: the GSP driving junctions'
# up/down LeafCutter overlap; p = binomial test of up vs down (p0 = 0.5).
# Inputs from cache/build_overlap_cache.R.
source("common.R")
cells <- fread(file.path("cache", "s7_leafcutter_cells.tsv"), data.table = FALSE)
lab <- c(`10` = "GSP #10 alt 3'ss in top 250\ndriving unannotated junctions",
         `13` = "GSP #13 alt 3'ss in top 250\ndriving unannotated junctions",
         `30` = "GSP #30 alt 3'/5'ss in top 250\ndriving unannotated junctions")
ttl <- function(p) gene_md(sub("^(\\w+) ", "\\1 ", sub("^(SF3B1|NKAP|YLPM1)", "\\1", p)))
pn_all <- unique(cells$panel)
# p as markdown with a real superscript exponent, replacing
# overlap_heatmap()'s plain caption.
p_md <- function(p) {
  e <- floor(log10(p)); m <- p / 10^e
  if (sprintf("%.1f", m) == "10.0") { m <- 1; e <- e + 1 }
  sprintf("<i>p</i> = %.1f &times; 10<sup>&minus;%d</sup>", m, -e)
}
panels <- lapply(seq_along(pn_all), function(k) {
  pn <- pn_all[k]
  d <- cells[cells$panel == pn, ]
  first_col <- (k - 1) %% 3 == 0
  d$direction <- d$direction   # levels up/down relabelled below
  # p goes on the title line (plain weight) rather than in the caption:
  # a caption sits under the x title, so the bottom row's p would hang
  # lower than the other rows'.
  # binom_p: binom.test(up, up + down, p = 0.5) among the driving junctions that
  # are LeafCutter hits. overlap_heatmap() reads `fisher_p` only for its own
  # caption, which is nulled below.
  overlap_heatmap(d %>% mutate(fisher_p = binom_p), lab[[as.character(d$fact[1])]],
                  title = paste0("<b>", gene_md(pn), "</b>\u2003",
                                 p_md(d$binom_p[1])),
                  xlab = "LeafCutter direction") +
    labs(caption = NULL) +
    scale_x_discrete(expand = c(0, 0),
                     labels = c(upregulated = "Up", downregulated = "Down")) +
    # Row labels only on the left column: repeated in every column they
    # take most of each panel's width.
    theme(axis.text.y = if (first_col) element_text() else element_blank()) +
    theme(plot.title = ggtext::element_markdown(size = SIZE_BASE, face = "plain",
                                                family = FONT_FAMILY, hjust = 0))
})
# axis_titles = "collect": one "LeafCutter direction" per column under
# the bottom row instead of one per panel.
fig <- wrap_plots(panels, ncol = 3) +
  plot_layout(guides = "collect", axis_titles = "collect") &
  theme(legend.position = "right")
save_supp(fig, "figS9", 8.5, 4.6)
