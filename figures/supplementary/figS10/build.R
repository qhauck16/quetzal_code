# Supplementary Figure S10. Enrichment of SUGP1 mutants on GSP #10.
# Samples ranked 40-100 on GSP #10, coloured by SF3B1 / SUGP1 status;
# SUGP1 carriers labelled with their protein change. p_adj is the
# cohort-stratified CMH for SUGP1 among the top-100 non-SF3B1 samples,
# Bonferroni over every gene (flash_fact_char.Rmd, "SUGP1 along GSP 10").
# Inputs are built by cache/build_mutation_cache.R.
source("common.R")

d   <- fread(file.path("cache", "s9_top100_gsp10.tsv"), data.table = FALSE)
cmh <- fread(file.path("cache", "s9_cmh_sf3b1_excluded.tsv"), data.table = FALSE)
padj <- cmh$padj[cmh$gene == "SUGP1"]

COL <- c(PAL$sf3b1_mut[c("SF3B1 hotspot mutation", "Mutation Near SF3B1 Hotspot",
                         "Other SF3B1 Mutation")],
         "SUGP1 carrier" = "#8E44AD",
         "No SF3B1/SUGP1 Mutation" = unname(PAL$sf3b1_mut["No SF3B1 Mutation"]))
LEG <- c("SF3B1 hotspot mutation" = "SF3B1 hotspot mutation",
         "Mutation Near SF3B1 Hotspot" = "Mutation near SF3B1 hotspot",
         "Other SF3B1 Mutation" = "Other SF3B1 mutation",
         "SUGP1 carrier" = "SUGP1 carrier",
         "No SF3B1/SUGP1 Mutation" = "No SF3B1/SUGP1 mutation")

pd <- d %>% filter(rank >= 40) %>%
  mutate(category = factor(category, levels = names(COL)),
         lab = ifelse(category == "SUGP1 carrier", sugp1_label, NA))

# p_adj as plotmath so the exponent is a real superscript.
e  <- floor(log10(padj)); mant <- padj / 10^e
p_lab <- sprintf('italic("SUGP1")~~italic(p)[adj] == %.1f %%*%% 10^%d', mant, e)

p <- ggplot(pd, aes(rank, k10, fill = category)) +
  geom_col(width = 0.9) +
  geom_text_repel(aes(label = lab), colour = COL[["SUGP1 carrier"]],
                  size = SIZE_BASE * 0.85 / .pt, family = FONT_FAMILY,
                  # "both" (was "y"): adjacent carriers at ranks 94/95
                  # stacked vertically and 94's leader ran through 95's label.
                  nudge_y = 0.05, direction = "both", min.segment.length = 0,
                  segment.size = 0.3, segment.colour = COL[["SUGP1 carrier"]],
                  box.padding = 0.2, na.rm = TRUE, seed = 1) +
  annotate("text", x = 72, y = 0.30, label = p_lab, parse = TRUE,
           size = SIZE_BASE * 1.1 / .pt, family = FONT_FAMILY, hjust = 0) +
  scale_fill_manual(values = COL, name = NULL, breaks = names(COL),
                    labels = function(x) gene_pm(unname(LEG[x]))) +
  scale_x_continuous(breaks = seq(40, 100, 10), expand = expansion(add = 0.8)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Rank along GSP #10 (ranks 40-100)", y = "GSP #10 loading") +
  theme_publication() +
  theme(legend.position = "inside", legend.position.inside = c(0.99, 0.99),
        legend.justification = c(1, 1), legend.background = element_blank(),
        legend.key.size = unit(3.5, "mm"))
save_supp(p, "figS10", 6.5, 4)
