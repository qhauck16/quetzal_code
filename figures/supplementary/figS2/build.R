# Supplementary Figure S2. Quetzal's recovery of differentially spliced
# junctions, as in figure 1 panel b, over a range of fold-change
# requirements (rows) for the three cohort comparisons (columns).
# Recovery logic and cache are figure 1's (figure1/compose.R
# build_panel_b()), looped over FC thresholds. FC >= 1 is the plain
# direction test.
if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
source("common.R")
suppressPackageStartupMessages({ library(yaml); library(GenomicRanges) })
CACHE <- file.path(FIG_ROOT, "figure1", "panel_b", "recovery_cache")

summed_props <- function(fls_obj) {
  res <- fls_obj$summary; jids <- str_split_i(colnames(res), "-", 1); uniq <- unique(jids)
  out <- matrix(0, 2, length(uniq), dimnames = list(NULL, uniq))
  for (j in uniq) out[, j] <- rowSums(res[, jids == j, drop = FALSE])
  out
}
hits_ref <- function(mdf) mdf %>% arrange(desc(abs(deltapsi))) %>%
  mutate(junc_id = paste0(start, "-", end)) %>% distinct(junc_id, .keep_all = TRUE)
run_recovery_fc <- function(sp, hr, sk, psi_vals, fc) {
  sapply(psi_vals, function(psi) {
    cn <- sum(abs(hr$deltapsi) > psi) + 1
    ranks <- as.numeric(str_split_i(colnames(sp), ":", 1)); ok <- ranks < cn
    cand <- sp[, ok, drop = FALSE]; cs <- sign(hr$deltapsi[ranks[ok]])
    c1 <- cand[1, ]; c2 <- cand[2, ]; hi <- pmax(c1, c2); lo <- pmin(c1, c2)
    ratio <- ifelse(lo > 0, hi / lo, ifelse(hi > 0, Inf, 1))
    succ <- sum(ratio >= fc & ((cs > 0 & c2 > c1) | (cs < 0 & c1 > c2)))
    sr <- as.numeric(str_split_i(sk$no_factor_match, ":", 1))
    succ / (succ + ncol(cand) - succ + sum(sr < cn))
  })
}
gencode <- readRDS(file.path(QROOT, "inputs/gencode_v46_granges.RDS"))
g2c <- data.frame(gene_name = gencode$gene_name[gencode$type == "gene"],
                  chrom = as.character(seqnames(gencode[gencode$type == "gene"]))) %>%
  distinct(gene_name, .keep_all = TRUE)

COMPS <- c(KICH_vs_LUAD = "KICH vs LUAD", KICH_vs_TGCT = "KICH vs TGCT",
           TGCT_vs_LUAD = "TGCT vs LUAD")
FCS <- c(1, 1.05, 1.1, 1.25, 1.5, 2)
psi_vals <- c(50, 47.5, 45, 42.5, 40, 39:10) / 100
rows <- list()
for (comp in names(COMPS)) {
  ft <- readRDS(file.path(CACHE, paste0(comp, "_fls_test.RDS")))
  fm <- readRDS(file.path(CACHE, paste0(comp, "_fls_matched.RDS")))
  mt <- fread(file.path(CACHE, paste0(comp, "_matched_test.tsv")), data.table = FALSE) %>%
    left_join(g2c, by = "gene_name") %>% filter(!is.na(chrom))
  mb <- fread(file.path(CACHE, paste0(comp, "_matched_background.tsv")), data.table = FALSE) %>%
    dplyr::select(-deltapsi) %>% left_join(mt %>% dplyr::select(subclass, deltapsi), by = "subclass") %>%
    left_join(g2c, by = "gene_name") %>% filter(!is.na(chrom)) %>% arrange(desc(abs(deltapsi)))
  tsp <- summed_props(ft); thr <- hits_ref(mt); msp <- summed_props(fm); mhr <- hits_ref(mb)
  for (fc in FCS) rows[[length(rows) + 1]] <- tibble(
    comp = COMPS[[comp]], fc = fc, psi = psi_vals,
    `LeafCutter DS` = run_recovery_fc(tsp, thr, ft$skipped, psi_vals, fc),
    `Matched control` = run_recovery_fc(msp, mhr, fm$skipped, psi_vals, fc))
}
d <- bind_rows(rows) %>% pivot_longer(c(`LeafCutter DS`, `Matched control`),
                                      names_to = "set", values_to = "rec") %>%
  mutate(fc = factor(paste0("FC ≥ ", fc), levels = paste0("FC ≥ ", FCS)),
         comp = factor(comp, levels = COMPS))

p <- ggplot(d, aes(psi, rec, colour = set)) +
  geom_line(linewidth = 0.4) + geom_point(size = 0.8) +
  facet_grid(fc ~ comp) +
  scale_x_reverse(breaks = c(0.5, 0.4, 0.3, 0.2, 0.1)) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1),
                     labels = scales::percent_format(accuracy = 1)) +
  scale_colour_manual(values = c(`LeafCutter DS` = PAL$leafcutter,
                                 `Matched control` = PAL$matched), name = NULL) +
  labs(x = "|ΔPSI| cutoff", y = "Percent recovery") +
  theme_publication() +
  theme(legend.position = "bottom", panel.spacing = unit(0.8, "lines"),
        # gridlines on in the supplement
        panel.grid.major = element_line(colour = "grey90", linewidth = 0.25),
        strip.background = element_rect(fill = "grey93", colour = NA),
        strip.text = element_text(size = SIZE_BASE),
        # FC strips read horizontally instead of rotated 90 deg
        strip.text.y = element_text(size = SIZE_BASE, angle = 0)) +
  guides(colour = guide_legend(override.aes = list(size = 2, linewidth = 0.8)))
save_supp(p, "figS2", 6.5, 8)
