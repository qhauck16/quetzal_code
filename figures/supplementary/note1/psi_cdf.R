# Supplementary Note 1: is GSP #1 / #12 driven by tissue (a covariate) or by
# the mutation's splicing? PSI of each GSP's top driving junctions
# (cache/psi_group_means.tsv, from note1/build_psi_cache.R).
#   row 1, GSP #1:  brain vs non-brain | IDH1 R132 vs other brain | per-junction dPSI of both
#   row 2, GSP #12: PRAD vs non-PRAD   | SPOP mutant vs other PRAD | per-junction dPSI of both
# Columns 1-2: CDF over junctions of the group's mean PSI (sqrt axis -- most
# driving junctions are low-PSI). Column 3: CDF of the per-junction difference
# for the tissue contrast (wild-type tissue vs other tissues) vs the mutation
# contrast (mutant vs wild-type tissue); Wilcoxon signed-rank
# of each difference vs 0, and the median.
# Run from supplementary/:  Rscript --vanilla note1/psi_cdf.R
source("common.R")
d <- fread("note1/cache/psi_group_means.tsv", data.table = FALSE) %>%
  filter(is.finite(mean_psi)) %>% mutate(junc = paste(gene_name, start, end))

CONTRASTS <- tribble(~gsp, ~kind, ~a, ~b,
  1,  "Tissue",   "Brain (LGG + GBM)",  "Non-brain",
  1,  "Mutation", "IDH1 R132 (brain)",  "Other brain",
  12, "Tissue",   "PRAD",               "Non-PRAD",
  12, "Mutation", "SPOP mutant (PRAD)", "Other PRAD")
COL <- c(a = "#d62728", b = "grey50")
KIND_COL <- c(Tissue = "#7570b3", Mutation = "#d62728")

ns <- d %>% group_by(gsp, group) %>% summarise(n = max(n), .groups = "drop")
psi_panel <- function(g, k) {
  cc <- filter(CONTRASTS, gsp == g, kind == k)
  x <- d %>% filter(gsp == g, group %in% c(cc$a, cc$b)) %>%
    group_by(junc) %>% filter(n() == 2) %>% ungroup() %>%
    mutate(group = factor(group, levels = c(cc$a, cc$b)))
  lab <- ns %>% filter(gsp == g, group %in% levels(x$group)) %>%
    mutate(group = factor(group, levels = levels(x$group))) %>% arrange(group)
  ggplot(x, aes(mean_psi, colour = group)) +
    stat_ecdf(linewidth = 0.6) +
    scale_x_continuous(trans = "sqrt", breaks = c(0, 0.01, 0.05, 0.1, 0.25, 0.5, 1)) +
    scale_y_continuous(labels = function(v) paste0(v * 100, "%")) +
    scale_colour_manual(values = setNames(unname(COL), levels(x$group)),
                        labels = gene_md(sprintf("%s (%d tumours)", lab$group, lab$n)), name = NULL) +
    labs(x = "Mean PSI per junction (sqrt scale)", y = "Cumulative % of junctions",
         title = sprintf("GSP #%d, %s: %d junctions", g, tolower(k), n_distinct(x$junc))) +
    theme_publication() + angled_x_labels() +
    theme(plot.title = element_text(size = SIZE_TITLE, face = "bold", hjust = 0.5),
          legend.text = ggtext::element_markdown(size = SIZE_LEGEND * 0.9, family = FONT_FAMILY),
          legend.position = "inside", legend.position.inside = c(0.98, 0.02),
          legend.justification = c(1, 0), legend.background = element_blank())
}
# Difference panels decompose additively: tissue = WILD-TYPE tissue - other
# tissues, mutation = mutant - wild-type tissue (sum = mutant - other
# tissues). Using all brain / all PRAD for the tissue term would count the
# mutants twice (412 of the 660 brain tumours are R132).
DCON <- tribble(~gsp, ~kind, ~a, ~b,
  1,  "Tissue",   "Other brain",        "Non-brain",
  1,  "Mutation", "IDH1 R132 (brain)",  "Other brain",
  12, "Tissue",   "Other PRAD",         "Non-PRAD",
  12, "Mutation", "SPOP mutant (PRAD)", "Other PRAD")
delta <- d %>% inner_join(DCON, by = "gsp", relationship = "many-to-many") %>%
  filter(group %in% c(a, b)) %>%
  group_by(gsp, kind, junc) %>%
  summarise(dpsi = mean_psi[group == a[1]] - mean_psi[group == b[1]], .groups = "drop") %>%
  filter(is.finite(dpsi))
dst <- delta %>% group_by(gsp, kind) %>%
  summarise(n = n(), median = median(dpsi), mean_abs = mean(abs(dpsi)),
            frac_up = mean(dpsi > 0), p = wilcox.test(dpsi)$p.value, .groups = "drop")
print(as.data.frame(dst), digits = 3)
fwrite(dst, file.path(OUT_DIR, "note1", "note1_psi_contrast_stats.tsv"), sep = "\t")

delta_panel <- function(g) {
  x <- delta %>% filter(gsp == g) %>% mutate(kind = factor(kind, levels = c("Tissue", "Mutation")))
  cc <- filter(DCON, gsp == g)
  s <- dst %>% filter(gsp == g) %>% mutate(kind = factor(kind, levels = levels(x$kind))) %>% arrange(kind)
  keys <- sprintf("%s: %s - %s<br>median %+.3f, P = %s", s$kind, cc$a[match(s$kind, cc$kind)],
                  cc$b[match(s$kind, cc$kind)], s$median, formatC(s$p, format = "e", digits = 1))
  ggplot(x, aes(dpsi, colour = kind)) +
    geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
    stat_ecdf(linewidth = 0.6) +
    scale_colour_manual(values = KIND_COL, labels = gene_md(keys), name = NULL) +
    scale_y_continuous(labels = function(v) paste0(v * 100, "%")) +
    coord_cartesian(xlim = quantile(x$dpsi, c(0.01, 0.99))) +   # a view; every junction stays in the ECDF
    labs(x = "Per-junction difference in mean PSI", y = "Cumulative % of junctions",
         title = sprintf("GSP #%d: tissue vs mutation shift", g)) +
    theme_publication() + angled_x_labels() +
    theme(plot.title = element_text(size = SIZE_TITLE, face = "bold", hjust = 0.5),
          legend.text = ggtext::element_markdown(size = SIZE_LEGEND * 0.85, family = FONT_FAMILY, lineheight = 1.1),
          legend.position = "bottom", legend.direction = "vertical",
          legend.justification = c(0, 1), legend.margin = margin(0, 0, 0, 0),
          legend.key.height = unit(18, "pt"))
}
fig <- (psi_panel(1, "Tissue") | psi_panel(1, "Mutation") | delta_panel(1)) /
       (psi_panel(12, "Tissue") | psi_panel(12, "Mutation") | delta_panel(12)) +
  plot_annotation_publication() &
  theme(plot.tag = element_text(size = SIZE_PANEL_LABEL, face = "bold", family = FONT_FAMILY))
save_supp(fig, "note1/note1_psi_cdf", 13, 9)
