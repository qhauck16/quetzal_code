# Supplementary figure: enrichment of mutations in non-NKAP mutants.
# Co-occurrence of AlphaMissense-filtered mutations among the top 50
# samples on GSP #13 that carry no NKAP mutation, for the genes that
# pass the cohort-stratified CMH (Bonferroni p_adj < 0.05). Diagonal =
# samples mutated in that gene (outlined); off-diagonal = samples
# mutated in both. Source: flash_fact_char.Rmd, "NKAP sig-gene pairwise
# co-occurrence heatmap". Inputs from cache/build_mutation_cache.R.
source("common.R")

co  <- fread(file.path("cache", "s12_cooccurrence.tsv"), data.table = FALSE)
cmh <- fread(file.path("cache", "s12_cmh_nkap_excluded.tsv"), data.table = FALSE)
g   <- unique(co$row)                         # already in p_adj order
padj <- setNames(cmh$padj, cmh$gene)[g]
# p_adj with a real superscript exponent (as S9), not "2.23e-15".
p_md <- function(p) {
  e <- floor(log10(p)); m <- p / 10^e
  if (sprintf("%.1f", m) == "10.0") { m <- 1; e <- e + 1 }
  sprintf("%.1f &times; 10<sup>&minus;%d</sup>", m, -e)
}
ylab <- sprintf("<i>%s</i><br><i>p</i><sub>adj</sub> = %s", g,
                vapply(padj, p_md, character(1)))

h <- co %>% mutate(i = match(row, g), j = match(col, g)) %>% filter(i <= j) %>%
  mutate(diag = i == j, x = factor(col, levels = g), y = factor(row, levels = rev(g)))

p <- ggplot(h, aes(x, y, fill = count)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_tile(data = filter(h, diag), fill = NA, colour = "black", linewidth = 0.7) +
  geom_text(aes(label = count), size = SIZE_BASE * 1.1 / .pt,
            family = FONT_FAMILY, fontface = "bold") +
  scale_fill_gradient(low = "grey94", high = "steelblue", name = "Samples (n)",
                      breaks = c(0, 4, 8)) +
  scale_x_discrete(position = "top", labels = function(x) paste0("<i>", x, "</i>")) +
  scale_y_discrete(labels = setNames(ylab, g)) +
  coord_equal() +
  labs(x = NULL, y = NULL) +
  theme_publication() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x.top = md_text(), axis.text.y = md_text(hjust = 1),
        # Legend in the empty lower-left triangle instead of a side
        # column of whitespace; page narrowed to match.
        legend.position = "inside", legend.position.inside = c(0.03, 0.05),
        legend.justification = c(0, 0), legend.direction = "horizontal",
        legend.title.position = "top",
        legend.key.width = unit(6, "mm"), legend.key.height = unit(2.5, "mm"))
save_supp(p, "figS12", 4.6, 3.75)
