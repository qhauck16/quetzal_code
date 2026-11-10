# Supplementary Figure S8. Summary of LeafCutter results: splice-event
# types among differentially spliced junctions, per GSP gene and its
# cohort-proportion-matched control run, across top-vs-bottom sample
# splits (columns) and dPSI / FC thresholds (rows).
# Source: diff_splicing/unified_results/leafcutter_meta_analysis.Rmd;
# counts from cache/build_leafcutter_counts.R (that Rmd's code verbatim).
source("common.R")
cn <- fread(file.path("cache", "s8_counts.tsv"), data.table = FALSE)
G <- c(sf3b1 = "SF3B1", nkap = "NKAP", ylpm1 = "YLPM1")
# House splice-type colours (PAL$splice_type, as in figure 2e) so a type
# looks the same across the paper; "canonical" (absent from PAL) gets
# Set2 yellow (Set2's 5th, lime, read too close to alt 3'ss teal). The
# published S8's plain Set2 order painted canonical in figure 2's
# "other" lavender.
ST <- c("alt 3'ss", "alt 5'ss", "canonical", "skipped exon", "other")
STC <- c(PAL$splice_type, canonical = RColorBrewer::brewer.pal(6, "Set2")[6])[ST]
LV <- c("SF3B1", "SF3B1 ctrl", "NKAP", "NKAP ctrl", "YLPM1", "YLPM1 ctrl")
TH <- expand_grid(dpsi = c(0.01, 0.05), fc = c(1, 2.5, 5)) %>%
  mutate(l = sprintf("ΔPSI > %g\nFC > %g", dpsi, fc)) %>% pull(l)
cn <- cn %>% mutate(
  cond = factor(ifelse(condition == "control", paste0(G[gene], " ctrl"), G[gene]), levels = LV),
  type = factor(type, levels = ST),
  comparison = factor(comparison, levels = c("10v30", "20v60", "30v90")),
  th = factor(sprintf("ΔPSI > %g\nFC > %g", dpsi, fc), levels = TH))
tot <- cn %>% group_by(th, comparison, cond) %>% summarise(t = sum(n), .groups = "drop") %>%
  mutate(lab = format(t, big.mark = ",", trim = TRUE))   # total per bar group

p <- ggplot(cn, aes(cond, n, fill = type)) +
  geom_col(position = position_dodge2(width = 0.85, preserve = "single"), width = 0.8) +
  # Totals sit in the band that the lower y expansion opens between
  # the zero baseline and the axis line (they used to print ON the
  # axis line). 6 pt = the paper's minimum.
  geom_text(data = tot, aes(cond, 0, label = lab), inherit.aes = FALSE, vjust = 1.4,
            size = 6 / .pt, family = FONT_FAMILY, colour = "grey25") +
  facet_grid(th ~ comparison, scales = "free_y") +
  scale_fill_manual(values = STC, name = "Splice type", drop = FALSE) +
  scale_x_discrete(labels = function(x) paste0("<i>", sub(" ctrl", "</i> ctrl", ifelse(grepl("ctrl", x), x, paste0(x, "</i>"))))) +
  scale_y_continuous(labels = scales::label_comma(), expand = expansion(mult = c(0.11, 0.06))) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "Differentially spliced junctions (n)") +
  theme_publication() +
  theme(axis.text.x = ggtext::element_markdown(angle = 45, hjust = 1, vjust = 1,
                                               size = SIZE_BASE * 0.85,
                                               margin = margin(t = 2)),
        legend.position = "bottom",
        plot.caption = element_text(size = SIZE_BASE * 0.85, hjust = 0, colour = "grey25"),
        plot.caption.position = "plot",
        strip.background = element_rect(fill = "grey93", colour = NA),
        strip.text = element_text(size = SIZE_BASE),
        # Threshold strips horizontal, not rotated to vertical.
        strip.text.y = element_text(size = SIZE_BASE, angle = 0, hjust = 0),
        panel.spacing.y = unit(1.1, "lines"), panel.spacing.x = unit(0.6, "lines"))
save_supp(p, "figS8", 8.5, 11)
