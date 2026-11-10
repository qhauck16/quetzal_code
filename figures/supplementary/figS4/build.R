# Supplementary Figure S4. Junctions driving GSP #3 (read-length GSP).
# Driving junctions: within a factor of two of the top-scoring junction in
# their single-gene group, over the top 200 groups on GSP #3. Left: those
# with mean PSI <= 0.1% in <= 50 bp samples; right: the same for >= 75 bp.
# Per junction, mean PSI over each cohort's samples in each read-length
# bucket; cohorts with >= 50 samples in both buckets.
# Inputs: cache/build_readlength_cache.R.
source("common.R")
COL <- c(`48 bp` = "#1F77B4", `76 bp` = "#FF7F0E")
# One row of three cohort CDFs per read-length condition. Row titles
# carry the condition; one shared legend at the bottom.
mk <- function(f, title) {
  d <- fread(file.path("cache", f), data.table = FALSE) %>%
    mutate(rl_group = factor(rl_group, levels = names(COL)))
  ggplot(d, aes(mean_psi, colour = rl_group)) +
    stat_ecdf(geom = "step", linewidth = 0.5) +
    facet_wrap(~ project_id, nrow = 1) +
    coord_cartesian(xlim = c(0, 0.5), ylim = c(0.75, 1)) +
    scale_x_continuous(breaks = c(0, 0.25, 0.5)) +
    scale_colour_manual(values = COL, name = "Read length") +
    labs(x = "Mean PSI", y = "Empirical CDF", title = title) +
    theme_publication() +
    theme(plot.title = element_text(size = SIZE_TITLE, face = "bold", hjust = 0),
          strip.background = element_blank(), strip.text = element_text(size = SIZE_BASE),
          aspect.ratio = 1, panel.spacing = unit(1, "lines"))
}
fig <- (mk("s4_low48_per_cohort.tsv", "Junctions with low PSI in samples with 48 bp reads") /
        mk("s4_low76_per_cohort.tsv", "Junctions with low PSI in samples with 76 bp reads")) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
save_supp(fig, "figS4", 6.5, 5.6)
