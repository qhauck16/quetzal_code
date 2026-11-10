# Supplementary Note 1, Figure 3: the IDH1 / SPOP LeafCutter results are
# representative and overlap their GSPs.
#   a: signature score in HELD-OUT LGG / PRAD tumours, mutant vs wild-type
#      (one-sided Wilcoxon), and PC1 of the signature junctions' scaled logit
#      PSI in four groups (mutant / wild-type x in LeafCutter run / held out;
#      two-sided Wilcoxon on the held-out pair) -- built by signature_plots.R.
#   b: S9-style overlap tiles: GSP driving alt 3'/5'ss junctions (from the top 250
#      unannotated, mutation_ds/gsp_top) against the LeafCutter up / down junctions
#      (cluster adj. P < 0.05, no effect-size floor); exact binomial test of up vs
#      down among the hits (p0 = 0.5) in each title. The PSI/RPM-matched background
#      row and its Fisher test were dropped 2026-11-10, as in S9. IDH1 R132 vs GSP #1 and
#      each SPOP residue vs GSP #12 (mutation_ds/out/overlap_cells.tsv).
# Run from supplementary/:  Rscript --vanilla note1/fig3.R
source("common.R")
N1_COMPOSE <- TRUE
source("note1/signature_plots.R")
source("note1/n1_common.R")

cells <- fread(file.path(N1_MDS, "out", "overlap_cells.tsv"), data.table = FALSE) %>%
  filter(min_dpsi == 0, min_fc == 1, types == "alt 3'/5'ss",
         (panel == "IDH1_R132" & gsp == 1) | (grepl("^SPOP_", panel) & gsp == 12))
COMP_LAB <- c(IDH1_R132 = "IDH1 R132 (20 v 60)", SPOP_F133 = "SPOP F133 (20 v 60)",
              SPOP_W131 = "SPOP W131 (14 v 42)", SPOP_F102 = "SPOP F102 (8 v 24)",
              SPOP_Y87 = "SPOP Y87 (7 v 21)")
p_md <- function(p) {
  e <- floor(log10(p)); m <- p / 10^e
  if (sprintf("%.1f", m) == "10.0") { m <- 1; e <- e + 1 }
  sprintf("<i>p</i> = %.1f &times; 10<sup>&minus;%d</sup>", m, -e)
}
heat <- lapply(seq_along(COMP_LAB), function(k) {
  d <- cells[cells$panel == names(COMP_LAB)[k] & cells$set == "gsp", ]
  up <- d$n_overlap[d$direction == "upregulated"]; dn <- d$n_overlap[d$direction == "downregulated"]
  d$fisher_p <- binom.test(up, up + dn, p = 0.5)$p.value   # binomial p (overlap_heatmap reads fisher_p)
  # short row labels on EVERY panel: IDH1 (GSP #1) and SPOP (GSP #12) panels
  # share rows of the layout, so one shared label would mislabel the others
  overlap_heatmap(d, sprintf("GSP #%d driving\njunctions", d$gsp[1]),
                  title = paste0("<b>", gene_md(COMP_LAB[[k]]), "</b><br>", p_md(d$fisher_p[1])),
                  xlab = "LeafCutter direction") +
    labs(caption = NULL) +
    scale_x_discrete(expand = c(0, 0), labels = c(upregulated = "Up", downregulated = "Down")) +
    theme(plot.title = ggtext::element_markdown(size = SIZE_BASE, face = "plain", family = FONT_FAMILY,
                                                hjust = 0, lineheight = 1.1))
})
p_b <- wrap_plots(c(heat, list(guide_area())), ncol = 3) + plot_layout(guides = "collect")

p_a <- (sig_panels$box_idh1 | sig_panels$box_spop) / (sig_panels$pc1_idh1 | sig_panels$pc1_spop)
fig <- wrap_elements(full = p_a) / wrap_elements(full = p_b) + plot_layout(heights = c(2.1, 1.3)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = SIZE_PANEL_LABEL, face = "bold", family = FONT_FAMILY))
ggsave_pub(file.path(N1_OUT, "SN1_Fig3.pdf"), fig, width = 10, height = 13)
cat("[fig3] wrote", file.path(N1_OUT, "SN1_Fig3.pdf"), "\n")
