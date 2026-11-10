# Supplementary figure: GSP #30 vs YLPM1 expression.
#   a  top 100 samples on GSP #30: log loading vs YLPM1 expression
#      (TCGA cohort-scaled RPM), YLPM1 mutants labelled
#      (ylpm1_lof_investigation.Rmd)
#   b  YLPM1 LoF samples split by SVM phenotype status
#   c  YLPM1 expression: phenotype-positive, the 25 phenotype-negative
#      intermediates nearest the boundary, and everyone else
#      (diff_splicing/ylpm1/top500_alt3_fc_vs_loading.Rmd, last two
#       chunks). Wilcoxon rank-sum p-values.
# NB LoF = protein change matching "*", "fs" or "splice", the rule panel
# b's source uses AND figure 4 panel a uses. The original panel a used
# the table's plot_color column, which leaves X1783_splice as "other";
# it is LoF here so the three panels and figure 4 agree.
source("common.R")

v <- fread(vt_path("YLPM1_full_plus_variants"), data.table = FALSE)
er <- fread(file.path("cache", "ylpm1_expression_row.tsv"), data.table = FALSE, header = TRUE)
ex <- tibble(rail_id = as.integer(colnames(er)[-1]), expr = as.numeric(er[1, -1]))
phen <- fread(vt_path("ylpm1_phenotype_axis_scores"),
              data.table = FALSE)
LOF <- "\\*|fs|[Ss]plice"
pc <- function(x) ifelse(is.na(x) | x %in% c("", "NA"), "", gsub("YLPM1:", "", x))

## ---- a -----------------------------------------------------------------
CLS <- c("No YLPM1 mutation", "Other YLPM1 mutation", "YLPM1 frameshift/stop gained")
COL <- setNames(c("grey65", "#1F77B4", "#D62728"), CLS)
a_df <- v %>% filter(!is.na(k30)) %>% slice_max(k30, n = 100, with_ties = FALSE) %>%
  inner_join(ex, by = "rail_id") %>%
  mutate(pc = pc(proteinChange),
         cls = factor(case_when(grepl(LOF, pc) ~ CLS[3], nzchar(pc) ~ CLS[2], TRUE ~ CLS[1]),
                      levels = CLS))
p_a <- ggplot(a_df, aes(log(k30), expr, colour = cls)) +
  geom_point(size = 1.4, alpha = 0.9) +
  # Repel over ALL points (unlabelled ones get "") -- ggrepel only avoids
  # the points in its own layer's data, so labels otherwise sit on greys.
  geom_text_repel(data = a_df, aes(label = pc),
                  size = SIZE_BASE * 0.75 / .pt, family = FONT_FAMILY,
                  max.overlaps = Inf, min.segment.length = 0, segment.size = 0.25,
                  segment.colour = "grey50", box.padding = 0.35, point.padding = 0.3, max.time = 5, max.iter = 1e5,
                  point.size = 1.4, force = 4, force_pull = 0.5, seed = 7,
                  # the solver leaves M1180Cfs*27 sitting on its own point
                  # (boxed in by greys); start it below the cloud instead.
                  nudge_y = ifelse(a_df$pc == "M1180Cfs*27", -0.5, 0),
                  nudge_x = ifelse(a_df$pc == "M1180Cfs*27", 0.25, 0),
                  show.legend = FALSE) +
  scale_colour_manual(values = COL, name = NULL, labels = function(x) gene_pm(x)) +
  labs(x = "log(GSP #30)", y = "YLPM1 expression (scaled RPM)",
       title = gene_md("YLPM1 expression vs GSP #30 loading")) +
  theme_publication() +
  theme(plot.title = md_title(hjust = 0), axis.title.y = md_text(),
        legend.position = "bottom") +
  guides(colour = guide_legend(override.aes = list(size = 2.2)))
p_a$labels$y <- gene_md("YLPM1 expression (scaled RPM)")

## ---- b, c -------------------------------------------------------------
fmt_p <- function(p) if (p < 1e-4) sprintf("p = %.1e", p) else sprintf("p = %.3g", p)
bracket <- function(x1, x2, y, tick, lab) list(
  annotate("segment", x = x1, xend = x2, y = y, yend = y, linewidth = 0.3),
  annotate("segment", x = c(x1, x2), xend = c(x1, x2), y = y, yend = y - tick, linewidth = 0.3),
  annotate("text", x = (x1 + x2) / 2, y = y + tick, label = lab, vjust = 0,
           size = SIZE_BASE * 0.85 / .pt, family = FONT_FAMILY))
BOX <- function(d, x, fillv, ylab, title) ggplot(d, aes({{ x }}, expr, fill = {{ x }})) +
  geom_boxplot(outlier.size = 0.4, alpha = 0.8, linewidth = 0.3) +
  scale_fill_manual(values = fillv, guide = "none") +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = gene_md(ylab), title = gene_md(title)) +
  theme_publication() + theme(plot.title = md_title(hjust = 0), axis.title.y = md_text()) +
  angled_x_labels()

pj <- phen %>% select(-any_of("k30")) %>% inner_join(ex, by = "rail_id") %>% filter(!is.na(expr))
lof <- pj %>% inner_join(v, by = "rail_id") %>% filter(grepl(LOF, proteinChange)) %>%
  mutate(g = factor(ifelse(score > 0, "pos", "neg"), levels = c("pos", "neg")))
wb <- wilcox.test(expr ~ g, data = lof)$p.value
nb <- table(lof$g)
yb <- max(lof$expr); rb <- diff(range(lof$expr))
p_b <- BOX(lof, g, c(pos = "#D62728", neg = "#4C72B0"), "YLPM1 expression (scaled RPM)",
           "YLPM1 LoF samples by\nSVM phenotype status") +
  scale_x_discrete(labels = c(pos = sprintf("Phen+ (n = %d)", nb[["pos"]]),
                              neg = sprintf("Phen− (n = %d)", nb[["neg"]]))) +
  bracket(1, 2, yb + 0.08 * rb, 0.03 * rb, fmt_p(wb))

n25 <- pj %>% filter(group == "intermediate", score <= 0) %>% arrange(desc(score)) %>%
  slice_head(n = 25) %>% pull(rail_id)
cd <- pj %>% mutate(g = factor(case_when(score > 0 ~ "pos", rail_id %in% n25 ~ "n25", TRUE ~ "rest"),
                               levels = c("pos", "n25", "rest")))
rest <- cd$expr[cd$g == "rest"]
w1 <- wilcox.test(cd$expr[cd$g == "pos"], rest)$p.value
w2 <- wilcox.test(cd$expr[cd$g == "n25"], rest)$p.value
yc <- max(cd$expr); rc <- diff(range(cd$expr))
p_c <- BOX(cd, g, c(pos = "#D62728", n25 = "#F0A400", rest = "#4C72B0"),
           "YLPM1 expression (scaled RPM)", "YLPM1 expression\nby SVM group") +
  scale_x_discrete(labels = c(pos = "Phen+", n25 = "Next 25 intermediates", rest = "Other samples")) +
  bracket(1, 3, yc + 0.06 * rc, 0.02 * rc, fmt_p(w1)) +
  bracket(2, 3, yc + 0.17 * rc, 0.02 * rc, fmt_p(w2))
cat(sprintf("[S17] b p = %.2g (n %d/%d); c p = %.2g, %.2g\n", wb, nb[1], nb[2], w1, w2))

# plot_annotation(theme=) does not style tags (patchwork takes each tag
# from its own panel's theme), so set the bold tag on every panel too.
fig <- p_a / (p_b | p_c) + plot_layout(heights = c(1.25, 1)) + plot_annotation_publication() &
  theme(plot.tag = element_text(size = SIZE_PANEL_LABEL, face = "bold", family = FONT_FAMILY))
save_supp(fig, "figS18", 8.5, 9)
