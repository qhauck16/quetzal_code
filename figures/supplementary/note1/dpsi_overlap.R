# Supplementary Note 1 (draft): IDH1 / SPOP mutant-vs-wild-type LeafCutter
# comparisons recover their GSP's driving junctions, but as small PSI shifts.
# Data: mutation_ds/out/{junction_dpsi,overlap_cells}.tsv (mutation_ds/overlap.R).
#   a: ECDF of LeafCutter dPSI (mutant - control) for the GSP's top driving
#      alt 3'/5'ss junctions vs the PSI/RPM-matched background, per comparison;
#      SF3B1 20v60 (GSP #10, S7) as reference. Dashed: |dPSI| = 0.05.
#      Label: one-sided Wilcoxon rank-sum, GSP > matched dPSI.
#   b: % of each set called LeafCutter-up (cluster p.adjust < 0.05, dPSI > 0)
#      under four filters: none (as S7), |dPSI| > 0.05, PSI fold change
#      (mutant / control) >= 2.5, and both; Fisher p on top. FC >= 2.5 asks
#      whether the GSP junctions are newly used (low-PSI junctions switched on)
#      or already-used junctions nudged a little.
if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
source("common.R")
MDS <- file.path(dirname(dirname(FIG_ROOT)), "quetzal_code", "mutation_ds", "out")
if (!dir.exists(MDS)) MDS <- file.path(QROOT, "mutation_ds/out")

COMPS <- tribble(~comparison, ~gsp, ~lab,
  "IDH1_R132",   1L,  "IDH1 R132 (20 v 60)\nGSP #1",
  "SPOP_F133",   12L, "SPOP F133 (20 v 60)\nGSP #12",
  "SPOP_W131",   12L, "SPOP W131 (14 v 42)\nGSP #12",
  "SPOP_F102",   12L, "SPOP F102 (8 v 24)\nGSP #12",
  "SPOP_Y87",    12L, "SPOP Y87 (7 v 21)\nGSP #12",
  "SF3B1 20v60", 10L, "SF3B1 GSP-high (20 v 60)\nGSP #10 (reference)")
COMPS$lab_md <- gene_md(COMPS$lab)
COMPS$short_md <- gene_md(c("IDH1 R132", "SPOP F133", "SPOP W131", "SPOP F102", "SPOP Y87", "SF3B1 (ref.)"))
SET_COL <- c(gsp = "#1f77b4", matched = "grey60")
SET_LAB <- c(gsp = "GSP driving junctions", matched = "Matched background")

jd <- fread(file.path(MDS, "junction_dpsi.tsv"), data.table = FALSE) %>%
  inner_join(COMPS, by = c("comparison", "gsp")) %>%
  mutate(lab_md = factor(lab_md, levels = COMPS$lab_md))
wt <- jd %>% group_by(lab_md) %>%
  summarise(p = wilcox.test(deltapsi[set == "gsp"], deltapsi[set == "matched"], alternative = "greater")$p.value,
            n_g = sum(set == "gsp"), n_m = sum(set == "matched"),
            med_g = median(deltapsi[set == "gsp"]), .groups = "drop") %>%
  mutate(lab = sprintf("Wilcoxon P = %s\nmedian GSP dPSI = %+.3f\nn = %d vs %d",
                       formatC(p, format = "e", digits = 1), med_g, n_g, n_m))
print(as.data.frame(wt))

p_a <- ggplot(jd, aes(deltapsi, colour = set)) +
  geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.3) +
  geom_vline(xintercept = c(-0.05, 0.05), linetype = "dashed", colour = "grey75", linewidth = 0.3) +
  stat_ecdf(linewidth = 0.6) +
  geom_text(data = wt, aes(x = 0.155, y = 0.04, label = lab), inherit.aes = FALSE,
            hjust = 1, vjust = 0, size = SIZE_BASE * 0.8 / .pt, family = FONT_FAMILY, lineheight = 0.9) +
  facet_wrap(~ lab_md, nrow = 2) +
  coord_cartesian(xlim = c(-0.1, 0.16)) +          # a view: keeps every junction in the ECDF
  scale_colour_manual(values = SET_COL, labels = SET_LAB, name = NULL) +
  labs(x = "LeafCutter dPSI (mutant or GSP-high minus control)", y = "Cumulative fraction of junctions") +
  theme_publication() +
  theme(strip.background = element_blank(),
        strip.text = ggtext::element_markdown(size = SIZE_STRIP, family = FONT_FAMILY, lineheight = 0.95),
        legend.position = "top")

oc <- fread(file.path(MDS, "overlap_cells.tsv"), data.table = FALSE) %>%
  filter(direction == "upregulated") %>%
  inner_join(COMPS, by = c("panel" = "comparison", "gsp")) %>%
  filter(types == ifelse(gsp == 10, "alt 3'ss", "alt 3'/5'ss")) %>%
  mutate(pct = 100 * n_overlap / n_set,
         thr = factor(case_when(min_dpsi == 0 & min_fc == 1 ~ "Any (as S7)",
                                min_dpsi > 0 & min_fc == 1  ~ "|dPSI| > 0.05",
                                min_dpsi == 0               ~ "FC >= 2.5",
                                TRUE                        ~ "|dPSI| > 0.05\n& FC >= 2.5"),
                      levels = c("Any (as S7)", "|dPSI| > 0.05", "FC >= 2.5", "|dPSI| > 0.05\n& FC >= 2.5")),
         short_md = factor(short_md, levels = COMPS$short_md))
pl <- oc %>% group_by(short_md, thr) %>%
  summarise(y = max(pct), p = dplyr::first(fisher_p), .groups = "drop") %>%
  mutate(lab = ifelse(p < 0.05, sprintf("P = %s", formatC(p, format = "e", digits = 1)), "n.s."))
p_b <- ggplot(oc, aes(thr, pct, fill = set)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_text(aes(label = n_overlap, group = set), position = position_dodge(width = 0.75),
            vjust = -0.3, size = SIZE_BASE * 0.75 / .pt, family = FONT_FAMILY) +
  # alternate heights so neighbouring P labels can't collide
  geom_text(data = pl, aes(thr, y + ifelse(as.integer(thr) %% 2 == 0, 15, 7), label = lab), inherit.aes = FALSE,
            size = SIZE_BASE * 0.8 / .pt, family = FONT_FAMILY) +
  facet_wrap(~ short_md, nrow = 1) +
  scale_fill_manual(values = SET_COL, labels = SET_LAB, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "% of set LeafCutter-up\n(cluster adj. P < 0.05)") +
  theme_publication() + angled_x_labels() +
  theme(strip.background = element_blank(),
        strip.text = ggtext::element_markdown(size = SIZE_STRIP, family = FONT_FAMILY, lineheight = 0.95),
        legend.position = "none")

fig <- p_a / p_b + plot_layout(heights = c(1.6, 1)) + plot_annotation_publication() &
  theme(plot.tag = element_text(size = SIZE_PANEL_LABEL, face = "bold", family = FONT_FAMILY))
save_supp(fig, "note1/note1_dpsi_overlap", 11, 8.5)
