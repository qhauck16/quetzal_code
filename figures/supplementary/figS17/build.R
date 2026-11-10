# Supplementary Figure S17. Intron coverage across the top 100 YLPM1-axis
# (GSP #30) disrupted introns: one row per gene (its top alt 3'/5'ss
# intron by delta-PSI, non-coding/UTR junctions excluded), left = mean
# coverage in the top 15 GSP #30 samples, right = project-matched
# controls. Each row is [exon]--intron--[exon] drawn 5' -> 3'; the maroon
# fill is intron coverage / mean coverage of the 20 flanking exonic bp
# (1 = exon-box height, bins < 0.1 zeroed, clipped at 2 x box). Rows
# sorted by the centre of mass of the top-15 trace (3'-most on top).
# Traces from cache/build_s16_cache.py (--n_genes 100).
source("common.R")

tr <- fread(file.path("cache", "s16_traces.tsv"), data.table = FALSE)
N_ROWS <- 100
N_CTRL <- as.integer(sub(".*n_ctrl ", "",
  grep("^wrote", readLines(file.path("cache", "build_s16_cache.log")), value = TRUE)))

genes <- tr %>% distinct(dpsi_rank, gene, com) %>% arrange(dpsi_rank) %>%
  slice_head(n = N_ROWS) %>%
  arrange(desc(ifelse(is.na(com), -Inf, com))) %>%
  mutate(row = row_number())

# Gene labels at 6 pt (the text floor); the 11 in page gives a ~7.7 pt
# row pitch, enough for 6 pt italic without rows touching.
LABEL_PT <- 6
EF <- 0.05; Y0 <- 0.2; EH <- 0.6; YMAX <- 1.4; FILL <- "#800000"
N_BIN <- max(tr$bin)
d <- tr %>% inner_join(genes %>% select(gene, row), by = "gene") %>%
  mutate(base = (N_ROWS - row) * YMAX,
         x    = EF + (bin - 1) / (N_BIN - 1) * (1 - 2 * EF),
         ymin = base + Y0,
         ymax = pmin(base + Y0 + EH * ratio, base + YMAX)) %>%
  # No baseline: runs of zeroed bins become NA so the ribbon breaks there
  # instead of drawing a hairline under every intron (a baseline makes
  # the enrichment look weaker). Zero bins touching coverage stay,
  # so each island tapers to the floor.
  group_by(gene, column) %>% arrange(bin, .by_group = TRUE) %>%
  mutate(nz = ratio > 0,
         ymax = ifelse(nz | lag(nz, default = FALSE) | lead(nz, default = FALSE),
                       ymax, NA)) %>%
  ungroup() %>%
  mutate(
         column = factor(column, levels = c("top", "ctrl"),
           labels = c("Top 15 samples on GSP #30",
                      sprintf("%d matched controls", N_CTRL))))
ex <- expand_grid(genes %>% select(row), side = c(0, 1 - EF)) %>%
  mutate(base = (N_ROWS - row) * YMAX)
labs_df <- genes %>% mutate(y = (N_ROWS - row) * YMAX + Y0 + EH / 2,
                            column = factor(levels(d$column)[1], levels(d$column)))

p <- ggplot() +
  geom_ribbon(data = d, aes(x = x, ymin = ymin, ymax = ymax,
              group = interaction(gene, column)), fill = FILL, colour = NA) +
  geom_rect(data = expand_grid(ex, column = factor(levels(d$column), levels(d$column))),
            aes(xmin = side, xmax = side + EF, ymin = base + Y0, ymax = base + Y0 + EH),
            fill = "black") +
  geom_text(data = labs_df, aes(x = -0.015, y = y, label = gene), hjust = 1,
            fontface = "italic", size = LABEL_PT / .pt, family = FONT_FAMILY) +
  facet_wrap(~column, nrow = 1) +
  scale_x_continuous(breaks = c(EF, 1 - EF), labels = c("5′ss", "3′ss"),
                     expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, N_ROWS * YMAX), clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_publication() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = SIZE_BASE * 0.8),
        strip.background = element_blank(),
        strip.text = element_text(size = SIZE_TITLE, face = "plain"),
        panel.spacing = unit(6, "mm"),
        plot.margin = margin(4, 6, 4, 46))

save_supp(p, "figS17", 8.5, 11)
