# Figure 2, panel g. Sourced by ../compose.R, which supplies the
# shared setup (theme tokens, FLASHIER/CACHE paths, ROW_H, panel_cached,
# missing_panel, sf3b1_info). Must leave a `panel_g` object defined.
#
# Rebuild just this panel:   Rscript compose.R g

# ---- g. Mutation enrichment: CMH vs linear ------------------------
if (!have_cache("panel_g_both_tests.tsv", "panel_g_thresholds.tsv")) {
  panel_g <- missing_panel("g  mutation enrichment", "build_cache_panel_g.R")
} else if (!identical(REBUILD, "ALL") && !("g" %in% REBUILD) &&
           file.exists(file.path(PANEL_DIR, "g.rds"))) {
  cat("[panel] g (cached)\n"); panel_g <- readRDS(file.path(PANEL_DIR, "g.rds"))
} else {
  cat("[panel] g building ...\n")
  g_df  <- fread(file.path(CACHE, "panel_g_both_tests.tsv"), data.table = FALSE)
  g_thr <- fread(file.path(CACHE, "panel_g_thresholds.tsv"), data.table = FALSE)
  g_cut <- g_thr$value[g_thr$quantity == "neg_log10_bonf_09"]

  g_df <- g_df %>%
    filter(n_mut_top25 > 2) %>%
    mutate(x = -log10(linear_pval),
           y = -log10(cmh_pval),
           # "13-NKAP" -> "GSP #13-NKAP" (the reference formatted this
           # in Illustrator); en-dash to match the other panels.
           lab = paste0("GSP #", gsp_num, "–", gene_name),
           # plotmath twin of `lab` with the gene italic: ggrepel has
           # no rich-text support, but it does parse plotmath.
           lab_pm = sprintf('paste("GSP #%s\u2013", italic("%s"))', gsp_num, gene_name))

  # Label the top 5 by LINEAR significance. That rule reproduces the
  # reference's label set exactly (SF3B1, NKAP, IDH1, SPOP, YLPM1;
  # #6 KIT at x=68 is unlabelled there too). Note it is not "passes
  # both thresholds" -- IDH1 is labelled despite a near-zero CMH
  # score, and HRAS/RBM10/PRKDC clear both cutoffs but are not shown.
  N_LABEL <- 5
  g_lab <- g_df %>% arrange(desc(x)) %>% head(N_LABEL)
  # A label whose point sits within ~a text height of the horizontal
  # Bonferroni line would be drawn straddling that dotted line (IDH1,
  # y = 0.9 vs cut 7.05). Start those labels well above it so the
  # repel keeps them clear; everything else starts at its point.
  # REVISED (QA review): lifting IDH1's label above the line made its
  # leader CROSS the dotted line. Instead a label near the line stays at
  # its point's height and moves sideways (right), so it sits entirely
  # below the line. NKAP's label drops below its point, clear of the
  # colourbar in the top-left.
  near <- abs(g_lab$y - g_cut) < 8
  g_lab$nudge_y <- ifelse(near, 0, ifelse(grepl("NKAP", g_lab$lab), -14, 0))
  g_lab$nudge_x <- ifelse(near, 28, 0)

  panel_g <- wrap_elements(full =
    ggplot(g_df, aes(x = x, y = y, color = n_mut_top25)) +
      geom_vline(xintercept = g_cut, linetype = "dotted",
                 color = "grey45", linewidth = 0.35) +
      geom_hline(yintercept = g_cut, linetype = "dotted",
                 color = "grey45", linewidth = 0.35) +
      geom_point(size = 1.4) +
      # Fixed dark text, NOT the viridis point colour: the top of the
      # scale is near-white yellow (SF3B1) / yellow-green (IDH1), which
      # is illegible as text on white. The reference sets these labels
      # in black for the same reason; the leader line carries the
      # label-to-point association instead.
      geom_text_repel(data = g_lab, aes(label = lab_pm), parse = TRUE,
                      size = SIZE_BASE * 0.85 / .pt,
                      family = FONT_FAMILY, color = "grey15",
                      seed = 1, box.padding = 0.45,
                      nudge_y = g_lab$nudge_y, nudge_x = g_lab$nudge_x,
                      min.segment.length = 0.2, segment.size = 0.25,
                      segment.color = "grey60", show.legend = FALSE) +
      scale_color_viridis_c(name = "Number of mutants in\ntop 25 samples on GSP",
                            breaks = c(5, 10, 15, 20)) +
      # `*` not `~`: plotmath's `~` inserts a full word-space, which
      # left "log10  (linear p-value)" oddly gapped. `*` juxtaposes.
      labs(x = expression(-log[10] * "(linear p-value)"),
           y = expression(-log[10] * "(CMH p-value)"),
           title = paste0("Mutation enrichment discovers additional\n",
                          "mutation-GSP relationships")) +
      theme_publication() +
      theme(plot.title        = element_text(hjust = 0.5),
            plot.margin       = margin(4, 6, 4, 14, "pt"),
            # x = 0.07, not 0.02: the vertical Bonferroni line sits at
            # ~0.045 npc (x = 7.05 on a 7-308 range), and at 0.02 the
            # dotted line struck straight through the legend title.
            legend.position.inside = c(0.07, 0.98),
            legend.position   = "inside",
            legend.justification = c(0, 1),
            legend.direction  = "horizontal",
            legend.key.height = unit(2.2, "mm"),
            legend.key.width  = unit(5, "mm"),
            legend.title      = element_text(size = SIZE_BASE * 0.8),
            legend.text       = element_text(size = SIZE_BASE * 0.75),
            legend.background = element_blank()))
  saveRDS(panel_g, file.path(PANEL_DIR, "g.rds"))
}

