# Supplementary Note 1: the IDH1 / SPOP LeafCutter splicing signatures
# (mutation_ds/signature_check.R).
#   a/b: signature score (mean PSI of mutant-up junctions - mean PSI of
#        mutant-down junctions; cluster adj. P < 0.05, |dPSI| > 0.05) in the
#        HELD-OUT tumours of the home cohort, i.e. those not in the LeafCutter
#        comparison; one-sided Wilcoxon, mutant > wild-type.
#   c/d: PC1 of a PCA of the same junctions' scaled logit(PSI) over every
#        tumour of the cohort, as four violins: mutant / wild-type x in the
#        LeafCutter run / held out; two-sided Wilcoxon (PC1's sign is
#        arbitrary) of held-out mutant vs held-out wild-type.
# Run from supplementary/:  Rscript --vanilla note1/signature_plots.R
if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
source("common.R")
MDS <- file.path(QROOT, "mutation_ds")
MUT_COL <- c(mutant = "#d62728", `wild-type` = "grey45")

sc <- fread(file.path(MDS, "out", "signature_scores.tsv"), data.table = FALSE,
            colClasses = c(rail_id = "character"))
in_run <- function(comp) as.character(fread(file.path(MDS, "samples", paste0(comp, ".tsv")))$rail_id)
fmt_p <- function(p) formatC(p, format = "e", digits = 1)

SETS <- list(
  idh1 = list(comp = "IDH1_R132", cohort = "TCGA-LGG", score = "idh1_score", gene = "IDH1",
              mut = function(x) x$r132, wt = function(x) !x$idh1_any, mut_lab = "IDH1 R132"),
  spop = list(comp = "SPOP_F133", cohort = "TCGA-PRAD", score = "spop_score", gene = "SPOP",
              mut = function(x) x$spop_any, wt = function(x) !x$spop_any, mut_lab = "SPOP mutant"))

# ---- a/b: held-out boxplots -------------------------------------------------
box_panel <- function(s) {
  x <- sc %>% filter(cohort == s$cohort, !rail_id %in% in_run(s$comp)) %>%
    mutate(status = case_when(s$mut(.) ~ "mutant", s$wt(.) ~ "wild-type")) %>%
    filter(!is.na(status), is.finite(.data[[s$score]])) %>%
    mutate(score = .data[[s$score]], status = factor(status, levels = names(MUT_COL)))
  n <- table(x$status)
  p <- wilcox.test(x$score[x$status == "mutant"], x$score[x$status == "wild-type"],
                   alternative = "greater")$p.value
  labs_x <- gene_md(c(sprintf("%s<br>(n = %d)", s$mut_lab, n[["mutant"]]),
                      sprintf("%s wild-type<br>(n = %d)", s$gene, n[["wild-type"]])))
  top <- max(x$score); rng <- diff(range(x$score))
  ggplot(x, aes(status, score, colour = status)) +
    geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.4) +
    geom_jitter(width = 0.15, height = 0, size = 0.6, alpha = 0.55) +
    annotate("segment", x = 1, xend = 2, y = top + 0.06 * rng, yend = top + 0.06 * rng, linewidth = 0.3) +
    annotate("text", x = 1.5, y = top + 0.09 * rng, label = sprintf("P = %s", fmt_p(p)), vjust = 0,
             size = SIZE_BASE * 0.9 / .pt, family = FONT_FAMILY) +
    scale_colour_manual(values = MUT_COL, guide = "none") +
    scale_x_discrete(labels = labs_x) +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.12))) +
    labs(x = NULL, y = "Signature score\n(mean PSI up - mean PSI down)",
         title = gene_md(sprintf("%s signature, held-out %s tumours", s$gene, sub("TCGA-", "", s$cohort)))) +
    theme_publication() +
    theme(plot.title = md_title(), axis.text.x = ggtext::element_markdown(size = SIZE_BASE, family = FONT_FAMILY))
}

# ---- c/d: PCA ----------------------------------------------------------------
pca_panel <- function(s) {
  r <- readRDS(file.path(MDS, "out", sprintf("signature_psi_%s.rds", s$comp)))
  M <- r$M[r$cohort == s$cohort, , drop = FALSE]
  keep_j <- colMeans(is.na(M)) <= 0.2                       # junctions measured in >= 80% of tumours
  M <- M[, keep_j, drop = FALSE]
  M <- apply(M, 2, function(v) { v[is.na(v)] <- median(v, na.rm = TRUE); v })
  L <- qlogis(pmin(pmax(M, 0.001), 0.999))                  # logit, clipped away from 0 / 1
  L <- L[, apply(L, 2, sd) > 0, drop = FALSE]
  pc <- prcomp(L, center = TRUE, scale. = TRUE)
  ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
  meta <- sc[match(rownames(L), sc$rail_id), ]
  d <- tibble(PC1 = pc$x[, 1], PC2 = pc$x[, 2],
              status = case_when(s$mut(meta) ~ "mutant", s$wt(meta) ~ "wild-type"),
              run = ifelse(rownames(L) %in% in_run(s$comp), "In LeafCutter run", "Held out")) %>%
    filter(!is.na(status)) %>% mutate(status = factor(status, levels = names(MUT_COL)))
  cat(sprintf("[pca] %s: %d tumours x %d junctions; PC1 %.1f%%, PC2 %.1f%%\n",
              s$comp, nrow(L), ncol(L), ve[1], ve[2]))
  # four groups: status x (in LeafCutter run | held out); the test is the
  # held-out pair (PC1's sign is arbitrary, so two-sided)
  lv <- c("mutant.In LeafCutter run", "wild-type.In LeafCutter run", "mutant.Held out", "wild-type.Held out")
  d <- d %>% mutate(grp = factor(paste(status, run, sep = "."), levels = lv))
  n <- table(d$grp)
  ho <- d %>% filter(run == "Held out")
  p <- wilcox.test(ho$PC1[ho$status == "mutant"], ho$PC1[ho$status == "wild-type"])$p.value
  top <- max(d$PC1); rng <- diff(range(d$PC1))
  xl <- gene_md(sprintf("%s<br>%s<br>(n = %d)",
                        c(s$mut_lab, paste(s$gene, "wild-type"), s$mut_lab, paste(s$gene, "wild-type")),
                        c("in LeafCutter run", "in LeafCutter run", "held out", "held out"), as.integer(n)))
  ggplot(d, aes(grp, PC1, colour = status)) +
    geom_violin(fill = NA, linewidth = 0.4, scale = "width", width = 0.75) +
    geom_boxplot(width = 0.12, outlier.shape = NA, linewidth = 0.35) +
    geom_jitter(width = 0.12, height = 0, size = 0.6, alpha = 0.55) +
    annotate("segment", x = 3, xend = 4, y = top + 0.06 * rng, yend = top + 0.06 * rng, linewidth = 0.3) +
    annotate("text", x = 3.5, y = top + 0.09 * rng, label = sprintf("held out: P = %s", fmt_p(p)), vjust = 0,
             size = SIZE_BASE * 0.9 / .pt, family = FONT_FAMILY) +
    scale_colour_manual(values = MUT_COL, guide = "none") +
    scale_x_discrete(labels = xl) +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.12))) +
    labs(x = NULL, y = sprintf("PC1 (%.1f%% of variance)", ve[1]),
         title = gene_md(sprintf("%s: PC1 of %d signature junctions (scaled logit PSI)",
                                 sub("TCGA-", "", s$cohort), ncol(L)))) +
    theme_publication() +
    theme(plot.title = md_title(),
          axis.text.x = ggtext::element_markdown(size = SIZE_BASE * 0.9, family = FONT_FAMILY, lineheight = 1.1))
}

sig_panels <- list(box_idh1 = box_panel(SETS$idh1), box_spop = box_panel(SETS$spop),
                   pc1_idh1 = pca_panel(SETS$idh1), pc1_spop = pca_panel(SETS$spop))
# fig3.R sets N1_COMPOSE and lays these out itself
if (!exists("N1_COMPOSE")) {
  fig <- (sig_panels$box_idh1 | sig_panels$box_spop) / (sig_panels$pc1_idh1 | sig_panels$pc1_spop) +
    plot_layout(heights = c(1, 1.1)) + plot_annotation_publication() &
    theme(plot.tag = element_text(size = SIZE_PANEL_LABEL, face = "bold", family = FONT_FAMILY))
  save_supp(fig, "note1/note1_signature", 10, 9)
}
