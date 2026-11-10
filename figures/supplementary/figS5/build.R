# Supplementary Figure S5. Extended GSP x covariate heatmap: every one of
# the 150 GSPs (rows, in GSP order, top = GSP 1) against the 33 TCGA
# cohorts (alphabetical), the technical covariates and the four summary
# columns. Same statistic, palette and group strip as main Figure 2c;
# only N_TOP (150 instead of 60) and the ordering differ. Matrix built by
# cache/build_s5_heatmap_cache.R (a copy of
# ../figure2/scripts/build_cache_panel_c.R with N_TOP = 150).
source("common.R")

m    <- fread(file.path("cache", "s5_matrix.tsv"), data.table = FALSE)
cols <- fread(file.path("cache", "s5_cols.tsv"),   data.table = FALSE)

GROUPS <- c("Tissue", "Technical", "Summary")
cols <- cols %>%
  mutate(g = match(group, GROUPS),
         key = ifelse(group == "Tissue", column, sprintf("%02d", display_order))) %>%
  arrange(g, key) %>%
  mutate(x = row_number(),
         label = ifelse(group == "Tissue", column, label))
n_gsp <- nrow(m)

long <- m %>%
  pivot_longer(-gsp, names_to = "column", values_to = "r") %>%
  mutate(x = cols$x[match(column, cols$column)],
         y = n_gsp - gsp + 1)

GROUP_COL <- c(Tissue = "#F1A340", Technical = "#998EC3", Summary = "#7FBF7B")
STRIP_H <- 3.2                       # in rows (150 rows -> ~0.02 per row)
strip_y <- n_gsp + 1.2 + STRIP_H / 2
strip   <- cols %>% transmute(x, y = strip_y, group)
grp_lab <- strip %>% group_by(group) %>%
  summarise(x0 = min(x), x1 = max(x), .groups = "drop") %>%
  mutate(x = (x0 + x1) / 2, y = strip_y + STRIP_H / 2 + 1.2)
div_at  <- which(diff(match(cols$group, GROUPS)) != 0) + 0.5

# Row ticks every 10 GSPs on the left so a reader can find a GSP.
y_brk <- c(1, seq(10, 150, 10))

p <- ggplot() +
  geom_tile(data = long, aes(x, y, fill = r)) +
  geom_tile(data = strip, aes(x, y), height = STRIP_H,
            fill = GROUP_COL[strip$group]) +
  geom_text(data = grp_lab, aes(x, y, label = group), vjust = 0,
            size = SIZE_BASE * 0.8 / .pt, family = FONT_FAMILY,
            fontface = "bold", colour = "grey20") +
  annotate("segment", x = div_at, xend = div_at, y = 0.5,
           yend = strip_y + STRIP_H / 2, colour = "black", linewidth = 0.3) +
  scale_fill_gradientn(colours = heatmap_diverging(), limits = c(-0.6, 0.6),
                       oob = scales::squish, breaks = c(-0.6, 0, 0.6),
                       name = "Pearson r") +
  scale_x_continuous(breaks = cols$x, labels = cols$label, expand = c(0, 0)) +
  scale_y_continuous(breaks = n_gsp - y_brk + 1, labels = y_brk,
                     expand = c(0, 0), name = "GSP") +
  coord_cartesian(ylim = c(0.5, strip_y + STRIP_H / 2 + 6), clip = "off") +
  labs(x = NULL) +
  theme_publication() +
  theme(axis.line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(linewidth = 0.25),
        axis.text.y = element_text(size = SIZE_BASE * 0.8),   # >= 6 pt
        panel.grid = element_blank(),
        legend.position = "right",
        legend.key.width  = unit(2.5, "mm"),
        legend.key.height = unit(9, "mm"),
        plot.margin = margin(4, 4, 4, 4)) +
  angled_x_labels(size = SIZE_BASE * 0.75)

save_supp(p, "figS5", 8.5, 9)
