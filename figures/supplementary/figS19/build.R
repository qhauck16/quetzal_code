# Supplementary Figure 19: loss of heterozygosity around YLPM1. Per 100 kb bin
# across YLPM1 +/- 50 Mb (hg19, clipped to chr14), the % of samples whose
# PanCanAtlas ABSOLUTE allelic segment keeps both homologues (minor HSCN > 0), for:
#   - the 11 single-LoF samples among the top 13 on GSP #30 (the other 2 carry two LoF mutations);
#   - the 25 intermediates closest to the decision boundary (score <= 0, highest scores);
#   - background: every other TCGA tumour with ABSOLUTE calls (incl. the other phenotype-positive samples;
#     5/16 LOH, not enriched).
# P = one-sided Fisher's exact test (greater) of LOH exactly at the YLPM1 locus (segment covering the YLPM1
# midpoint), each group vs background. Inputs: PanCanAtlas TCGA_mastercalls.abs_segtabs.fixed.txt,
# figure4/cache/ylpm1_phenotype_axis_scores.tsv, tcga_v2_samples.tsv (rail_id -> barcode).
if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
source("common.R")
PA  <- file.path(QROOT, "inputs/pancanatlas")
YL  <- c(75230069, 75326186); MID <- mean(YL)                                  # YLPM1, hg19
WIN <- c(max(0, YL[1] - 5e7), min(107349540, YL[2] + 5e7))
TOP13  <- c(882990, 869984, 441042, 461535, 870567, 939375, 440985, 865440, 866724, 434773, 222488, 434710, 443515)
DOUBLE <- c(461535, 866724)
sc   <- fread("../figure4/cache/ylpm1_phenotype_axis_scores.tsv")
meta <- fread(file.path(QROOT, "inputs/tcga_v2_samples.tsv"), select = c("rail_id", "tcga_barcode"))
arr  <- function(r) substr(meta$tcga_barcode[match(r, meta$rail_id)], 1, 15)
near <- head(sc[group == "intermediate" & score <= 0][order(-score)], 25)$rail_id
G <- list(single = arr(setdiff(TOP13, DOUBLE)), near = arr(near))
seg <- fread(file.path(PA, "TCGA_mastercalls.abs_segtabs.fixed.txt"))[
  Chromosome == 14 & End >= WIN[1] & Start <= WIN[2], .(Sample, Start, End, het = pmin(Modal_HSCN_1, Modal_HSCN_2) > 0)]
G$bg <- setdiff(unique(seg$Sample), c(arr(TOP13), G$near))
G <- lapply(G, function(s) intersect(s, seg$Sample))                           # samples with ABSOLUTE calls
at <- seg[Start <= MID & End >= MID]
loh <- function(s) sum(!at[Sample %in% s]$het)
fp  <- sapply(c("single", "near"), function(g) {
  a <- loh(G[[g]]); b <- length(G[[g]]) - a; c <- loh(G$bg); d <- length(G$bg) - c
  fisher.test(matrix(c(a, b, c, d), 2), alternative = "greater")$p.value })
fmtp <- function(p) if (p < 1e-3) sprintf("p = %.1e", p) else sprintf("p = %.2g", p)
LAB <- c(single = sprintf("Top GSP #30, single LoF (n = %d)", length(G$single)),
         near   = sprintf("25 nearest the decision boundary (n = %d)", length(G$near)),
         bg     = sprintf("TCGA background (n = %d)", length(G$bg)))
cat(sprintf("LOH at YLPM1: single %d/%d (%s); near %d/%d (%s); background %.1f%%\n", loh(G$single), length(G$single),
            fmtp(fp["single"]), loh(G$near), length(G$near), fmtp(fp["near"]), 100 * loh(G$bg) / length(G$bg)))
bins <- data.table(Start = seq(WIN[1], WIN[2], 1e5)); bins[, End := Start]
setkey(seg, Start, End)
prof <- rbindlist(lapply(names(G), function(g) {
  o <- foverlaps(bins, seg[Sample %in% G[[g]]], type = "within", nomatch = NULL)
  o[, .(het = mean(het)), by = .(pos = i.Start)][, grp := LAB[[g]]]
}))
prof[, grp := factor(grp, LAB)]
# p-values written just above each tested line. x per line: the single-LoF
# label sits over the flat stretch past 87 Mb, clear of the step at ~84 Mb.
PX <- c(single = 88e6, near = MID + 3e6)
pl <- data.table(grp = factor(LAB[c("single", "near")], LAB), x = PX,
                 y = mapply(function(l, x) prof[grp == l][which.min(abs(pos - x))]$het, LAB[c("single", "near")], PX) + 0.035,
                 lab = sapply(fp[c("single", "near")], fmtp))
COL <- setNames(c("#D62728", "#9467BD", "#1F77B4"), LAB)
p <- ggplot(prof, aes(pos, het, colour = grp)) +
  geom_vline(xintercept = MID, colour = "#E6A800", linewidth = 1.6, alpha = 0.55) +
  annotate("text", x = MID, y = 1.05, label = "YLPM1", fontface = "italic", size = SIZE_BASE / .pt, family = FONT_FAMILY) +
  geom_step(linewidth = 0.8) +
  geom_text(data = pl, aes(x, y, label = lab, colour = grp), hjust = 0, vjust = 0, size = SIZE_BASE / .pt,
            family = FONT_FAMILY, show.legend = FALSE) +
  scale_colour_manual(values = COL, name = NULL) +
  scale_y_continuous(limits = c(0, 1.07), breaks = seq(0, 1, 0.25), labels = scales::percent) +
  scale_x_continuous(expand = c(0, 0), labels = function(v) sprintf("%.0f Mb", v / 1e6)) +
  labs(x = "chr14 (hg19)", y = "Samples heterozygous") +
  theme_publication() + theme(legend.position = "bottom", legend.direction = "vertical")
save_supp(p, "figS19", 8, 4.3)
