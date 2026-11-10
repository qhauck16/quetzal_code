# Supplementary figure: gene-expression effects of mis-spliced genes
# driving GSP #13, across PSI and FC cutoffs. Same as figure 3 panel e
# (frameshifted vs in-frame alt 3'ss, cryptic site -20..0 nt from the
# annotated 3'ss, mean scaled RPM in the top 20 samples on GSP #13) for
# every (dPSI, FC) cell. Source: diff_splicing/combined_expression_
# analysis.Rmd, plot_faceted_cdf("NKAP", nkap_df, "frameshift",
# dist_window = c(-20, 0), require_alt3 = TRUE). One-sided Wilcoxon
# (frameshifted < in-frame).
# Expression: figure 3 panel e's cache (the 13.4 GB matrix sliced to the
# same 20 top samples), so the full matrix is never read here.
# NB the LeafCutter table is the 20v60 dpsi0.01_fc2.5 instance, so every
# junction already has FC > 2.5 -- the FC > 1 row necessarily equals the
# FC > 2.5 row, exactly as in the original figure.
if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
source("common.R")

# header = TRUE: fread autodetect treats this header as data (its first
# cell pattern-matches the gene symbols in column 1).
em <- fread(file.path(FIG_ROOT, "figure3", "cache", "panel_e_scaled.tsv"),
            data.table = FALSE, header = TRUE)
rownames(em) <- em$gene_name; em <- as.matrix(em[, -1])
df <- fread(file.path(QROOT, "leafcutter_analyses/diff_splicing/unified_results/nkap/20v60/dpsi0.01_fc2.5/df_final.tsv"), data.table = FALSE)

PSI <- c(0.01, 0.02, 0.05, 0.10); FC <- c(1, 2.5, 5)
cells <- expand.grid(psi = PSI, fc = FC)
curves <- list(); ann <- list()
for (i in seq_len(nrow(cells))) {
  ps <- cells$psi[i]; fc <- cells$fc[i]
  b <- df %>% filter(deltapsi > ps, test / control > fc, type == "alt 3'ss",
                     !is.na(distance_to_cryptic),
                     distance_to_cryptic > -20, distance_to_cryptic < 0)
  tg <- intersect(unique(b$gene_name[b$distance_to_cryptic %% 3 != 0]), rownames(em))
  cp <- setdiff(intersect(unique(b$gene_name[b$distance_to_cryptic %% 3 == 0]), rownames(em)), tg)
  tx <- rowMeans(em[tg, , drop = FALSE], na.rm = TRUE)
  cx <- rowMeans(em[cp, , drop = FALSE], na.rm = TRUE)
  curves[[i]] <- bind_rows(tibble(x = tx, g = "Frameshifted"), tibble(x = cx, g = "In-frame")) %>%
    mutate(psi = ps, fc = fc)
  pv <- wilcox.test(tx, cx, alternative = "less")$p.value
  ann[[i]] <- tibble(psi = ps, fc = fc,
                     lab = sprintf("Frameshifted: %d\nIn-frame: %d\np = %.2g", length(tg), length(cp), pv))
}
cv <- bind_rows(curves); an <- bind_rows(ann)
lab_psi <- function(x) paste0("ΔPSI > ", x); lab_fc <- function(x) paste0("FC > ", x)

p <- ggplot(cv, aes(x, colour = g)) +
  stat_ecdf(linewidth = 0.55) +
  geom_text(data = an, aes(x = Inf, y = 0.03, label = lab), inherit.aes = FALSE,
            hjust = 1.03, vjust = 0, size = SIZE_BASE * 0.8 / .pt,
            family = FONT_FAMILY, lineheight = 0.95) +
  facet_grid(fc ~ psi, labeller = labeller(psi = lab_psi, fc = lab_fc)) +
  scale_colour_manual(values = c(Frameshifted = "#B22222", `In-frame` = "#4682B4"), name = NULL) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  labs(x = "Mean scaled RPM in top-20 samples on GSP #13", y = "Cumulative fraction") +
  theme_publication() +
  theme(legend.position = "bottom", panel.spacing = unit(0.6, "lines"),
        strip.background = element_rect(fill = "grey93", colour = NA),
        strip.text = element_text(size = SIZE_BASE, face = "plain")) +
  guides(colour = guide_legend(override.aes = list(linewidth = 1.2)))
save_supp(p, "figS14", 8.5, 5.8)
print(an)
