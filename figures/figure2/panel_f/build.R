# Figure 2, panel f. Sourced by ../compose.R, which supplies the
# shared setup (theme tokens, FLASHIER/CACHE paths, ROW_H, panel_cached,
# missing_panel, sf3b1_info). Must leave a `panel_f` object defined.
#
# Rebuild just this panel:   Rscript compose.R f

# ------------------------------------------------------------------
# f. PSI fold-change vs binned GSP #10
# ------------------------------------------------------------------
f_psi <- fread(vt_path("sf3b1_top500_alt3_fc5_psi"),
               data.table = FALSE)
f_psi$mean_psi_alt3 <- rowMeans(as.matrix(f_psi[, -1]), na.rm = TRUE)

f_breaks <- c(-Inf, seq(-8, 0, by = 0.5))
f_labels <- c("(-Inf,-8]",
              paste0("(", head(seq(-8, 0, by = 0.5), -1), ",",
                     tail(seq(-8, 0, by = 0.5), -1), "]"))

f_df <- f_psi %>%
  select(rail_id, mean_psi_alt3) %>%
  inner_join(sf3b1_info, by = "rail_id") %>%
  mutate(fc_vs_global = mean_psi_alt3 / mean(mean_psi_alt3, na.rm = TRUE),
         logk10_bin   = cut(logk10, breaks = f_breaks, labels = f_labels,
                            right = TRUE))

panel_f <- panel_cached("f", wrap_elements(full =
  ggplot(f_df, aes(x = logk10_bin, y = fc_vs_global)) +
    # colour = grey25 (not the default black): the left-hand bins are so
    # short that top/bottom/median strokes touch, and at full black they
    # collapsed into solid bars that hid the steelblue fill. The same
    # colour softens the dense outlier stack in the "(-Inf,-8]" bin so it
    # reads as a cloud of points rather than a black spike.
    geom_boxplot(fill = PAL$fill_neutral, colour = "grey25", alpha = 0.6,
                 linewidth = 0.25, outlier.size = 0.2, outlier.alpha = 0.3) +
    geom_hline(yintercept = 1, linetype = "dotted", color = "grey40",
               linewidth = 0.4) +
    # Two lines: single-line, the y title was taller than the y axis and
    # overran into the header band beside the "f" tag. Wording unchanged.
    labs(x = "binned log(GSP #10)",
         y = "(sample mean PSI /\nglobal mean PSI)",
         title = paste0("Quetzal results strongly correlate with\n",
                        "LeafCutter results from top-loaded samples")) +
    theme_publication() +
    theme(plot.title  = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.margin = margin(4, 6, 4, 16, "pt"))))

