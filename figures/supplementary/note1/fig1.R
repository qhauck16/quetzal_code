# Supplementary Note 1, Figure 1: within-cohort GSP loading CDFs, mutant vs
# the rest of the same cohort.
#   a: TCGA-UVM, GSP #10 -- SF3B1 hotspot (Fig. 2d class) vs other UVM
#   b: TCGA-UCEC, GSP #13 -- NKAP C-term (300-413) / PRPF8 RNase H (1570-1900)
#      variant (Fig. 3a classes) vs other UCEC
#   c: TCGA-LGG, GSP #1 -- IDH1 mutant vs IDH1 wild-type (tumours with calls)
#   d: TCGA-PRAD, GSP #12 -- SPOP mutant vs SPOP wild-type (tumours with calls)
# x = log(GSP loading) (natural log, as Fig. 3a / 4a), y = cumulative %. A
# dotted segment at 50% joins the two medians and is labelled with the fold
# change exp(median log mutant - median log other). One-sided Wilcoxon
# rank-sum, mutant > other, in the legend.
# Run from supplementary/:  Rscript --vanilla note1/fig1.R
source("common.R")
source("note1/n1_common.R")
COL <- c("#d62728", "grey45")

two_cdf <- function(x, mut_lab, oth_lab, gsp, title) {
  x <- x %>% filter(is.finite(logk)) %>%
    mutate(grp = factor(ifelse(mut, mut_lab, oth_lab), levels = c(mut_lab, oth_lab)))
  p  <- wilcox.test(x$logk[x$mut], x$logk[!x$mut], alternative = "greater")$p.value
  md <- c(median(x$logk[x$mut]), median(x$logk[!x$mut]))
  fold <- exp(md[1] - md[2])
  n  <- table(x$grp)
  cat(sprintf("[fig1] %s: %d vs %d, P = %s, medians %.2f / %.2f, %s\n", title, n[[1]], n[[2]],
              fmt_p(p), md[1], md[2], fmt_fold(fold)))
  d <- x %>% arrange(logk) %>% group_by(grp) %>% mutate(cumulative = row_number() / n() * 100) %>% ungroup()
  keys <- gene_md(c(sprintf("%s (n = %d), P = %s", mut_lab, n[[1]], fmt_p(p)),
                    sprintf("%s (n = %d)", oth_lab, n[[2]])))
  ggplot(d, aes(logk, cumulative, colour = grp)) +
    geom_step(linewidth = 0.6) +
    annotate("segment", x = md[2], xend = md[1], y = 50, yend = 50, linetype = "dotted", linewidth = 0.5) +
    annotate("point", x = md, y = 50, size = 1.2) +
    annotate("text", x = mean(md), y = 53, label = fmt_fold(fold), vjust = 0,
             size = SIZE_BASE * 0.9 / .pt, family = FONT_FAMILY) +
    scale_colour_manual(values = setNames(COL, levels(d$grp)), breaks = levels(d$grp),
                        labels = setNames(keys, levels(d$grp)), name = NULL) +
    scale_y_continuous(labels = function(v) paste0(v, "%"), breaks = c(0, 25, 50, 75, 100)) +
    labs(x = sprintf("log(GSP #%d)", gsp), y = "Cumulative percentage", title = gene_md(title)) +
    theme_publication() +
    theme(plot.title = md_title(),
          legend.text = ggtext::element_markdown(size = SIZE_LEGEND * 0.9, family = FONT_FAMILY),
          legend.position = "inside", legend.position.inside = c(0.02, 0.98),
          legend.justification = c(0, 1), legend.background = element_blank(),
          legend.key.width = unit(16, "pt"), legend.key.height = unit(9, "pt"))
}

sf <- n1_sf3b1() %>% filter(cohort == "TCGA-UVM") %>%
  transmute(logk = logk10, mut = plot_color == "SF3B1 hotspot mutation")
nk <- n1_nkap() %>% filter(cohort == "TCGA-UCEC") %>%
  transmute(logk = logk13, mut = variant_class %in% c("NKAP variant", "PRPF8 variant"))
tm <- n1_tumours()
ih <- tm %>% filter(cohort == "TCGA-LGG")  %>% transmute(logk = logk1,  mut = idh1_any)
sp <- tm %>% filter(cohort == "TCGA-PRAD") %>% transmute(logk = logk12, mut = spop_any)

p_a <- two_cdf(sf, "SF3B1 hotspot", "Other UVM", 10, "TCGA-UVM: SF3B1 and GSP #10")
p_b <- two_cdf(nk, "NKAP C-term / PRPF8 RNase H variant", "Other UCEC", 13, "TCGA-UCEC: NKAP / PRPF8 and GSP #13")
p_c <- two_cdf(ih, "IDH1 mutant", "IDH1 wild-type", 1,  "TCGA-LGG: IDH1 and GSP #1")
p_d <- two_cdf(sp, "SPOP mutant", "SPOP wild-type", 12, "TCGA-PRAD: SPOP and GSP #12")

fig <- (p_a | p_b) / (p_c | p_d) + plot_annotation_publication() &
  theme(plot.tag = element_text(size = SIZE_PANEL_LABEL, face = "bold", family = FONT_FAMILY))
ggsave_pub(file.path(N1_OUT, "SN1_Fig1.pdf"), fig, width = 10, height = 8)
cat("[fig1] wrote", file.path(N1_OUT, "SN1_Fig1.pdf"), "\n")
