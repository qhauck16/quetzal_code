# Supplementary Figure S7. SF3B1 p.Y765C mutants.
# Same CDF as figure 2 panel d's left half (GSP #10, per-class
# cumulative percentage), with the two p.Y765C carriers labelled.
# Data: SF3B1_full_plus_variants (vt_path()).
source("common.R")

d <- fread(vt_path("SF3B1_full_plus_variants"), data.table = FALSE)
has <- function(x) str_detect(coalesce(x, ""), "Y765C")
y765c <- d$rail_id[has(d$gene_mutation) | has(d$wxs_mutation) | has(d$rnaseq_mutation)]
stopifnot(length(y765c) == 2)          # the caption says "the two samples"

LEV <- names(PAL$sf3b1_mut)
LEG <- c("SF3B1 hotspot mutation"      = "SF3B1 hotspot mutation",
         "Mutation Near SF3B1 Hotspot" = "Mutation near SF3B1 hotspot",
         "Other SF3B1 Mutation"        = "Other SF3B1 mutation",
         "No SF3B1 Mutation"           = "No SF3B1 mutation")

cum <- d %>%
  filter(!is.na(logk10)) %>%
  mutate(cls = factor(plot_color, levels = LEV)) %>%
  arrange(logk10) %>%
  group_by(cls) %>% mutate(cumulative = row_number() / n() * 100) %>% ungroup()
lab <- cum %>% filter(rail_id %in% y765c) %>%
  mutate(txt = paste0("p.Y765C (", project_id, ")"),
         # Label anchors placed by hand in empty space, clear of
         # the curves: BRCA right of the steep rise, below the blue
         # plateau; LUSC above its point on the blue plateau.
         lx = logk10 + ifelse(project_id == "TCGA-BRCA", 1.4, 0),
         ly = cumulative + ifelse(project_id == "TCGA-BRCA", 9, 11))

p <- ggplot(cum, aes(logk10, cumulative, colour = cls)) +
  geom_step(linewidth = 0.6) +
  geom_segment(data = lab, aes(x = logk10, y = cumulative, xend = lx, yend = ly),
               colour = "grey45", linewidth = 0.3) +
  geom_point(data = lab, size = 2.2, shape = 21, fill = "white", stroke = 0.9) +
  geom_label(data = lab, aes(lx, ly, label = txt), colour = "grey10",
             size = SIZE_BASE / .pt, family = FONT_FAMILY, hjust = 0, vjust = 0,
             label.size = 0, label.padding = unit(1, "pt"), fill = "white") +
  scale_colour_manual(values = PAL$sf3b1_mut, limits = LEV, name = NULL,
                      labels = function(x) gene_pm(unname(LEG[x]))) +
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     breaks = c(0, 25, 50, 75, 100)) +
  coord_cartesian(ylim = c(0, 102), clip = "off") +
  labs(x = "log(GSP #10)", y = "Cumulative percentage") +
  theme_publication() +
  theme(legend.position = "inside",
        # Above the red curve's long bottom run (2%), clear of the
        # steep hotspot rise at the far right.
        legend.position.inside = c(0.84, 0.10),
        legend.justification = c(1, 0),
        legend.background = element_blank(),
        legend.key = element_blank()) +
  guides(colour = guide_legend(override.aes = list(linewidth = 1.2, shape = NA)))

save_supp(p, "figS7", 5.5, 4.2)
