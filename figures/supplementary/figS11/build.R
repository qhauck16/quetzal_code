# Supplementary Figure S11. Other NKAP mutants.
# Empirical CDF of GSP #13 as in figure 3 panel a: any NKAP variant vs
# every other sample, with every NKAP mutant affected at amino acid 329
# or earlier labelled.
# Source: tcga/genome_wide/flashier/cdf_supplements.Rmd (chunk
# "NKAP aa <= 329"); data NKAP_full_plus_variants.tsv. The mutation
# parsing is that chunk's, verbatim: it reads all three mutation
# columns (TCGA call, WXS re-call, RNA-seq re-call).
source("common.R")

d <- fread(vt_path("NKAP_full_plus_variants"), data.table = FALSE)
MCOL <- c("gene_mutation", "wxs_mutation", "rnaseq_mutation")
nkap_changes <- function(row) {
  out <- character(0)
  for (col in MCOL) {
    s <- row[[col]]
    if (is.na(s) || s %in% c("", "NA")) next
    for (b in regmatches(s, gregexpr("NKAP:[^;]+", s))[[1]])
      out <- c(out, strsplit(sub("^NKAP:", "", b), ",")[[1]])
  }
  unique(trimws(out[out != ""]))
}
ch  <- lapply(seq_len(nrow(d)), function(i) nkap_changes(d[i, ]))
pos <- lapply(ch, function(x) suppressWarnings(as.numeric(regmatches(x, regexpr("\\d+", x)))))
d$cls <- ifelse(lengths(ch) > 0, "NKAP variant", "Other")
d$lab <- vapply(seq_along(ch), function(i) {
  k <- which(!is.na(pos[[i]]) & pos[[i]] <= 329)
  # Samples that also carry an R402 hotspot change (e.g. R402Q + L313M)
  # sit high because of R402, not the early variant -- don't label them.
  if (any(grepl("^R402", ch[[i]]))) k <- integer(0)
  if (length(k)) paste(ch[[i]][k], collapse = ", ") else ""
}, character(1))

COL <- c("NKAP variant" = unname(VARIANT_COL_NKAP <- "#8E44AD"), "Other" = "grey70")
cum <- d %>% filter(!is.na(logk13)) %>%
  mutate(cls = factor(cls, levels = names(COL))) %>%
  arrange(logk13) %>%
  group_by(cls) %>% mutate(cumulative = row_number() / n() * 100) %>% ungroup()
lab <- cum %>% filter(nzchar(lab))
cat("[S11]", nrow(lab), "labelled NKAP mutants (aa <= 329)\n")
# The low-loading mutants (log GSP #13 < -9) pile up on the steep part
# of the curve, where repelled labels landed on both curves. They go in
# a stacked column in the empty lower-right quadrant instead, in the
# same vertical order as their points so leader lines do not cross.
COL_X <- -9.3
left <- lab %>% filter(logk13 < -9) %>% arrange(cumulative, logk13) %>%
  mutate(ly = seq(3, 44, length.out = n()))
right <- lab %>% filter(logk13 >= -9)

p <- ggplot(cum, aes(logk13, cumulative, colour = cls)) +
  geom_step(linewidth = 0.6) +
  geom_segment(data = left, aes(x = logk13, y = cumulative, xend = COL_X - 0.15, yend = ly),
               colour = "grey55", linewidth = 0.25) +
  geom_point(data = lab, size = 1.6) +
  geom_text(data = left, aes(COL_X, ly, label = lab), colour = "grey10", hjust = 0,
            size = SIZE_BASE * 0.8 / .pt, family = FONT_FAMILY) +
  geom_text_repel(data = right, aes(label = lab), colour = "grey10",
                  size = SIZE_BASE * 0.8 / .pt, family = FONT_FAMILY,
                  direction = "both", force = 2, max.overlaps = Inf,
                  nudge_x = -0.6, nudge_y = 3,
                  min.segment.length = 0, segment.size = 0.25,
                  segment.colour = "grey55", box.padding = 0.25, seed = 1) +
  scale_colour_manual(values = COL, name = NULL,
                      labels = function(x) gene_pm(x)) +
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     breaks = c(0, 25, 50, 75, 100)) +
  coord_cartesian(ylim = c(0, 102), clip = "off") +
  labs(x = "log(GSP #13)", y = "Cumulative percentage") +
  theme_publication() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.98, 0.04),
        legend.justification = c(1, 0),
        legend.background = element_blank(),
        legend.key = element_blank()) +
  guides(colour = guide_legend(override.aes = list(linewidth = 1.2, shape = NA)))

save_supp(p, "figS11", 5.5, 4.2)
