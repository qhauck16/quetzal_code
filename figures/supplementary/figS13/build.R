# Supplementary figure: GSP #13 and #30 characterization.
#   a, b  splice types of the driving unannotated junctions of the top 250
#         single-gene junction groups on GSP #13 / #30
#         (fact_{13,30}_top_junctions.tsv, one row per snaptron_id).
#   c, d  per-sample mean PSI of the top 500 LeafCutter junctions,
#         as fold change over the global mean, binned by log loading
if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
source("common.R")

pie_df <- function(fact) {
  fread(file.path(QROOT, "gsp_analysis", sprintf("fact_%d_top_junctions.tsv", fact)),
        data.table = FALSE) %>%
    group_by(snaptron_id) %>% slice_head(n = 1) %>% ungroup() %>%
    dplyr::count(type) %>%
    mutate(prop = n / sum(n), label = sprintf("%.1f%%", prop * 100))
}
# Drawn pie radius (in), measured off the render and iterated to a fixed
# point, as for the main-figure pies (../scripts/pie_chart.R).
PIE_R_IN <- 1.05
mk_pie <- function(fact) {
  d <- pie_df(fact)
  lv <- intersect(names(PAL$splice_type), d$type)
  pie_labels_on_wedges(d, "type", PAL$splice_type, levels_order = lv,
                       pie_r_in = PIE_R_IN, label_pt = SIZE_BASE,
                       family = FONT_FAMILY) +
    labs(title = paste0("GSP #", fact)) +
    theme(plot.title = element_text(size = SIZE_TITLE, face = "bold",
                                    hjust = 0.5, family = FONT_FAMILY))
}

BR  <- c(-Inf, seq(-8, 0, by = 0.5))
BLB <- c("(-Inf,-8]", paste0("(", head(seq(-8, 0, 0.5), -1), ",",
                             tail(seq(-8, 0, 0.5), -1), "]"))
mk_box <- function(psi_file, info_file, kcol, title) {
  psi <- fread(psi_file, data.table = FALSE)
  psi$mean_psi <- rowMeans(as.matrix(psi[, -1]), na.rm = TRUE)
  info <- fread(vt_path(info_file), data.table = FALSE)
  d <- psi %>% dplyr::select(rail_id, mean_psi) %>% inner_join(info, by = "rail_id")
  d$fc  <- d$mean_psi / mean(d$mean_psi, na.rm = TRUE)
  d$bin <- cut(d[[paste0("log", kcol)]], breaks = BR, labels = BLB, right = TRUE)
  ggplot(d, aes(bin, fc)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey55",
               linewidth = 0.3) +
    geom_boxplot(outlier.size = 0.3, fill = PAL$fill_neutral, alpha = 0.7,
                 linewidth = 0.3) +
    scale_x_discrete(drop = FALSE) +
    labs(x = paste0("log(GSP #", sub("k", "", kcol), ") bin"),
         y = "Fold change\n(sample mean PSI / global mean PSI)",
         title = gene_md(title)) +
    theme_publication() +
    theme(plot.title = md_title(hjust = 0)) +
    angled_x_labels(size = SIZE_BASE * 0.8)
}
DS <- file.path(QROOT, "leafcutter_analyses/diff_splicing")
p_c <- mk_box(vt_path("nkap_top500_alt3_fc5_psi"),
              "NKAP_full_plus_variants", "k13",
              "NKAP 20v60 top 500 alt 3'ss")
p_d <- mk_box(vt_path("ylpm1_top500_alt3_5_fc5_psi"),
              "YLPM1_full_plus_variants", "k30",
              "YLPM1 6v18 top 500 alt 3'ss + alt 5'ss")

leg <- wrap_elements(full = pie_legend_grob(
  pie_df(13), "type", PAL$splice_type,
  levels_order = intersect(names(PAL$splice_type), pie_df(13)$type),
  label_pt = SIZE_BASE, family = FONT_FAMILY), ignore_tag = TRUE)

fig <- (wrap_elements(full = mk_pie(13)) | wrap_elements(full = mk_pie(30)) |
          leg) + plot_layout(widths = c(1, 1, 0.45))
fig <- (fig / (p_c | p_d)) + plot_layout(heights = c(1, 1.15))
# plot_annotation(theme=) does not style tags (patchwork takes each tag
# from its own panel's theme), so set the bold tag on every panel too.
fig <- fig + plot_annotation_publication() &
  theme(plot.tag = element_text(size = SIZE_PANEL_LABEL, face = "bold",
                                family = FONT_FAMILY))
save_supp(fig, "figS13", 8.5, 7.2)
