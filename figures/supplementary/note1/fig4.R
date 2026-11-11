# Supplementary Note 1, Figure 4: how many junctions each LeafCutter
# comparison changes at increasing effect-size thresholds.
# Count = junctions in significant clusters (adj. P <= 0.05) that go UP in the
# test group (mutant / GSP-high) with dPSI > {0.01, 0.05, 0.10} and PSI fold
# change (test / control) > {1, 2.5, 5}, counted exactly as Fig. S8
# (cache/leafcutter_counts_lib.R: junctions in genes with a MANE transcript). Bars per comparison, coloured by the
# GSP the comparison belongs to; the wild-type-vs-wild-type null runs (LGG,
# PRAD), when present, are grey. dPSI is test - control.
# Run from supplementary/:  Rscript --vanilla note1/fig4.R
source("common.R")
source("note1/n1_common.R")
DS <- file.path(QROOT, "leafcutter_analyses/diff_splicing/results")

COMPS <- bind_rows(
  expand_grid(g = c("sf3b1", "nkap", "ylpm1"), cmp = c("10v30", "20v60", "30v90")) %>%
    transmute(comp = paste(toupper(g), cmp), base = file.path(DS, paste0(g, "_", cmp)),
              gsp = c(sf3b1 = "GSP #10 (SF3B1)", nkap = "GSP #13 (NKAP)", ylpm1 = "GSP #30 (YLPM1)")[g]),
  tibble(comp = c("IDH1 R132 20v60", "SPOP F133 20v60", "SPOP W131 14v42", "SPOP F102 8v24", "SPOP Y87 7v21"),
         base = file.path(N1_MDS, "results", c("IDH1_R132", "SPOP_F133", "SPOP_W131", "SPOP_F102", "SPOP_Y87")),
         gsp = c("GSP #1 (IDH1)", rep("GSP #12 (SPOP)", 4))),
  tibble(comp = c("Null LGG 20v60", "Null PRAD 20v60"),
         base = file.path(N1_MDS, "results", c("NULL_LGG", "NULL_PRAD")), gsp = "Null (wild-type vs wild-type)")) %>%
  filter(file.exists(file.path(base, "leafcutter_ds_cluster_significance.txt")))

source("cache/leafcutter_counts_lib.R")   # same junction counting as Fig. S8
TH <- expand_grid(dpsi = c(0.01, 0.05, 0.10), fc = c(1, 2.5, 5))
counts <- bind_rows(lapply(seq_len(nrow(COMPS)), function(i) {
  x <- load_full_leaf_sig(COMPS$base[i])
  TH %>% mutate(n = mapply(function(a, b) nrow(assign_splice_types(x, a, b)), dpsi, fc),
                comp = COMPS$comp[i], gsp = COMPS$gsp[i])
}))
fwrite(counts, file.path(N1_OUT, "SN1_Fig4_counts.tsv"), sep = "\t")
print(counts %>% filter(dpsi == 0.05) %>% pivot_wider(names_from = fc, values_from = n, names_prefix = "fc>") %>%
        as.data.frame())

GSP_COL <- c("GSP #10 (SF3B1)" = "#C0392B", "GSP #13 (NKAP)" = "#8E44AD", "GSP #30 (YLPM1)" = "#E67E22",
             "GSP #1 (IDH1)" = "#1F77B4", "GSP #12 (SPOP)" = "#2CA02C",
             "Null (wild-type vs wild-type)" = "grey60")
counts <- counts %>%
  mutate(comp = factor(comp, levels = COMPS$comp), gsp = factor(gsp, levels = names(GSP_COL)),
         dlab = factor(sprintf("dPSI > %g", dpsi), levels = sprintf("dPSI > %g", c(0.01, 0.05, 0.10))),
         flab = factor(sprintf("FC > %g", fc), levels = sprintf("FC > %g", c(1, 2.5, 5))))
fig <- ggplot(counts, aes(comp, n, fill = gsp)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = scales::comma(n)), vjust = -0.3, size = SIZE_BASE * 0.6 / .pt, family = FONT_FAMILY) +
  facet_wrap(~ dlab + flab, ncol = 3, scales = "free_y",   # own y per panel: FC > 5 counts are small
             labeller = labeller(.multi_line = FALSE)) +
  scale_fill_manual(values = GSP_COL, name = NULL, labels = gene_md(names(GSP_COL)), drop = TRUE) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "Junctions up in test group\n(cluster adj. P < 0.05)") +
  theme_publication() + angled_x_labels() +
  theme(strip.background = element_blank(), strip.text = element_text(size = SIZE_STRIP),
        axis.text.x = ggtext::element_markdown(size = SIZE_BASE * 0.8, angle = 45, hjust = 1, family = FONT_FAMILY),
        legend.text = ggtext::element_markdown(size = SIZE_LEGEND, family = FONT_FAMILY),
        legend.position = "bottom")
fig <- fig + scale_x_discrete(labels = function(x) gene_md(x))
ggsave_pub(file.path(N1_OUT, "SN1_Fig4.pdf"), fig, width = 12, height = 8)
cat("[fig4] wrote", file.path(N1_OUT, "SN1_Fig4.pdf"), "\n")
