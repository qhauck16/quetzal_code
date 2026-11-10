# Supplementary Note 1, Figure 2: GSP loading CDFs outside the home cohort, drawn
# like the NKAP (Fig 3a) and YLPM1 (Fig 4a) CDFs: x = log(GSP loading)
# (natural log), y = cumulative percentage, step curves.
#   a: GSP #1 -- LGG IDH1 mutants vs LGG IDH1 wild-type, plus IDH1 R132
#      carriers vs IDH1 wild-type pooled over every NON-BRAIN cohort (GBM out).
#   b: GSP #12 -- PRAD SPOP mutants vs PRAD wild-type, plus every tumour of
#      the other cohorts ("Other cohorts", SPOP mutants included); the 3
#      hotspot (F133/W131/F102/Y87) carriers outside PRAD are labelled points
#      ON that curve (cohort + change), not a curve.
# Samples: tumours in the flashier fit with mutation calls (AlphaMissense-
# merged PanCanAtlas + GDC table); home-cohort mutant = any call in the gene;
# wild-type = no call in the gene. All curves solid (Quinn 2026-10-09).
# Per-cohort stats for every cohort are still written to note1_cdf_stats.tsv.
# Loadings: shared_cache/gsp_loadings.tsv (the figures' normalised L).
# Test: one-sided Wilcoxon rank-sum, mutant > wild-type.
# Run from supplementary/:  Rscript --vanilla note1/fig2.R
if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
source("common.R")

L <- fread(file.path(FIG_ROOT, "shared_cache", "gsp_loadings.tsv"), select = c("rail_id", "k1", "k12"),
           data.table = FALSE)
meta <- fread(file.path(QROOT, "inputs/tcga_v2_samples.tsv"), select = c("rail_id", "tcga_barcode", "gdc_cases.project.project_id"),
              data.table = FALSE) %>%
  transmute(rail_id, cohort = gdc_cases.project.project_id,
            sampleId = vapply(strsplit(tcga_barcode, "-"), function(x) {
              x[4] <- sub("[A-Za-z]+$", "", x[4]); paste(x[1:4], collapse = "-") }, character(1))) %>%
  filter(rail_id %in% L$rail_id)
am <- fread(file.path(QROOT, "mutation_processing/merged_grch38_with_am.tsv"), select = c("sampleId", "hugoGeneSymbol", "proteinChange"),
            data.table = FALSE)
meta <- meta %>% filter(sampleId %in% am$sampleId) %>% left_join(L, by = "rail_id")

status <- function(gene, hot_re, kcol) {
  g <- am %>% filter(hugoGeneSymbol == gene)
  meta %>% mutate(group = ifelse(sampleId %in% g$sampleId, "mutant", "wild-type"),
                  hotspot = sampleId %in% g$sampleId[grepl(hot_re, g$proteinChange)],
                  hot_change = vapply(sampleId, function(x) paste(unique(
                    g$proteinChange[g$sampleId == x & grepl(hot_re, g$proteinChange)]), collapse = ", "),
                    character(1)),
                  logk = log(.data[[kcol]]), gene = gene)
}
idh  <- status("IDH1", "R132", "k1")
spop <- status("SPOP", "(F133|W131|F102|Y87)", "k12")

stats_for <- function(d) d %>% group_by(gene, cohort) %>%
  filter(sum(group == "mutant") >= 1) %>%
  summarise(n_mut = sum(group == "mutant"), n_hot = sum(hotspot), n_wt = sum(group == "wild-type"),
            p = if (n_mut >= 2) wilcox.test(logk[group == "mutant"], logk[group == "wild-type"],
                                            alternative = "greater")$p.value else NA_real_,
            med_mut = median(logk[group == "mutant"]), med_wt = median(logk[group == "wild-type"]),
            .groups = "drop")
st <- bind_rows(stats_for(idh), stats_for(spop)) %>% arrange(gene, desc(n_mut))
print(as.data.frame(st %>% filter(n_mut >= 5)), digits = 3)
fwrite(st, file.path(OUT_DIR, "note1", "note1_cdf_stats.tsv"), sep = "\t")

fmt_p <- function(p) ifelse(p < 1e-3, formatC(p, format = "e", digits = 1), formatC(p, format = "f", digits = 3))

# a/b: four curves per gene -- home-cohort mutants and wild-type, plus the
# hotspot carriers (IDH1 R132; SPOP F133/W131/F102/Y87) and the gene-wild-type
# tumours of every OTHER cohort, pooled. Mutants of the gene outside the home
# cohort that are not hotspots are left out of both "other" curves.
# `drop`: cohorts left out of the "other cohorts" curves entirely (IDH1: the
# other brain cohort, GBM -- Quinn 2026-10-09, "for IDH1, just do non-brain",
# since GSP #1 is brain-associated and GBM confounds the pooled comparison).
# all_other: the 4th curve is EVERY tumour outside the home cohort, mutants
# included, labelled just "Other cohorts"; the hotspot carriers are then
# returned as `pts` to be labelled on that curve (SPOP, Quinn 2026-10-09).
four_way <- function(d, gene, home, hot, drop = character(0), other_lab = "other cohorts",
                     all_other = FALSE) {
  lv <- c(sprintf("%s mutant (%s)", gene, sub("TCGA-", "", home)),
          sprintf("%s wild-type (%s)", gene, sub("TCGA-", "", home)),
          sprintf("%s %s, %s", gene, hot, other_lab),
          if (all_other) "Other cohorts" else sprintf("%s wild-type, %s", gene, other_lab))
  d <- d %>% filter(!cohort %in% drop)
  if (all_other) {
    pts <- d %>% filter(cohort != home, hotspot)
    d <- d %>% mutate(curve = case_when(cohort == home & group == "mutant" ~ lv[1],
                                        cohort == home ~ lv[2], TRUE ~ lv[4])) %>%
      mutate(curve = factor(curve, levels = lv))
    wt <- function(x, y) wilcox.test(x, y, alternative = "greater")$p.value
    s <- d %>% count(curve, .drop = FALSE) %>%
      mutate(p = c(wt(d$logk[d$curve == lv[1]], d$logk[d$curve == lv[2]]), NA, NA, NA))
    return(list(d = d, s = s, lv = lv, pts = pts))
  }
  d <- d %>% mutate(curve = case_when(
      cohort == home & group == "mutant" ~ lv[1],
      cohort == home                    ~ lv[2],
      hotspot                           ~ lv[3],
      group == "wild-type"              ~ lv[4])) %>%
    filter(!is.na(curve)) %>% mutate(curve = factor(curve, levels = lv))
  wt <- function(x, y) wilcox.test(x, y, alternative = "greater")$p.value
  s <- d %>% count(curve, .drop = FALSE) %>%
    mutate(p = c(wt(d$logk[d$curve == lv[1]], d$logk[d$curve == lv[2]]), NA,
                 if (sum(d$curve == lv[3]) >= 2) wt(d$logk[d$curve == lv[3]], d$logk[d$curve == lv[4]]) else NA, NA))
  list(d = d, s = s, lv = lv)
}
# label_hot: draw the hotspot carriers as labelled points ON the other-cohort
# wild-type curve instead of as a curve of their own (SPOP: only 3 carriers
# outside PRAD -- too few for a CDF, Quinn 2026-10-09).
four_plot <- function(fw, gsp, title, label_hot = FALSE) {
  d <- fw$d %>% arrange(logk) %>% group_by(curve) %>% mutate(cumulative = row_number() / n() * 100) %>% ungroup()
  key <- fw$s %>% mutate(k = ifelse(is.na(p), sprintf("%s (n = %d)", curve, n),
                                    sprintf("%s (n = %d), P = %s", curve, n, fmt_p(p))))
  show <- if (label_hot) fw$lv[-3] else fw$lv
  cols <- setNames(c("#d62728", "grey45", "#ff7f0e", "grey70"), fw$lv)[show]
  key  <- key %>% filter(curve %in% show)
  p <- ggplot(filter(d, curve %in% show), aes(logk, cumulative, colour = curve)) +
    geom_step(linewidth = 0.6) +
    scale_colour_manual(values = cols, breaks = show, labels = gene_md(key$k), name = NULL) +
    scale_y_continuous(labels = function(x) paste0(x, "%"), breaks = c(0, 25, 50, 75, 100)) +
    labs(x = sprintf("log(GSP #%d)", gsp), y = "Cumulative percentage", title = gene_md(title)) +
    theme_publication() +
    theme(plot.title = md_title(),
          legend.text = ggtext::element_markdown(size = SIZE_LEGEND * 0.9, family = FONT_FAMILY),
          legend.position = "inside", legend.position.inside = c(0.02, 0.98),
          legend.justification = c(0, 1), legend.background = element_blank(),
          legend.key.width = unit(16, "pt"), legend.key.height = unit(9, "pt"))
  if (label_hot) {
    ref <- d %>% filter(curve == fw$lv[4])
    pts <- fw$pts %>% mutate(y = 100 * ecdf(ref$logk)(logk))
    print(pts %>% select(cohort, hot_change, logk, y))
    p <- p +
      geom_point(data = pts, aes(logk, y), inherit.aes = FALSE, colour = "#ff7f0e", size = 2) +
      ggrepel::geom_text_repel(data = pts, aes(logk, y, label = sprintf("%s: %s %s",
                                 sub("TCGA-", "", cohort), fw$gene, hot_change)),
                               inherit.aes = FALSE, colour = "black", size = SIZE_BASE * 0.85 / .pt,
                               family = FONT_FAMILY, nudge_x = 1.8, nudge_y = -12, direction = "y",
                               hjust = 0, segment.size = 0.3, segment.colour = "grey40",
                               min.segment.length = 0, seed = 1)
  }
  p
}
fw_i <- four_way(idh,  "IDH1", "TCGA-LGG",  "R132", drop = "TCGA-GBM", other_lab = "non-brain cohorts")
fw_s <- four_way(spop, "SPOP", "TCGA-PRAD", "hotspot", all_other = TRUE)
fw_i$gene <- "IDH1"; fw_s$gene <- "SPOP"
print(fw_i$s); print(fw_s$s)
p_a <- four_plot(fw_i, 1,  "IDH1 and GSP #1")
p_b <- four_plot(fw_s, 12, "SPOP and GSP #12", label_hot = TRUE)

fig <- (p_a | p_b) + plot_annotation_publication() &
  theme(plot.tag = element_text(size = SIZE_PANEL_LABEL, face = "bold", family = FONT_FAMILY))
ggsave_pub(file.path(OUT_DIR, "note1", "SN1_Fig2.pdf"), fig, width = 11, height = 4.6)
cat("[fig2] wrote", file.path(OUT_DIR, "note1", "SN1_Fig2.pdf"), "\n")
