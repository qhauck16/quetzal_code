# Supplementary Note 1: the five GSP CDFs on ONE shared x scale, to show how
# far each gene's mutants sit from the non-mutant loadings.
#   SF3B1 / GSP #10, NKAP / GSP #13, YLPM1 / GSP #30: pan-cancer, with the
#     exact classes of their main-figure CDFs (Fig 2d, 3a, 4a).
#   IDH1 / GSP #1 (TCGA-LGG only) and SPOP / GSP #12 (TCGA-PRAD only):
#     mutant (any call in the gene) vs wild-type, as in note1/cdf_loadings.R.
# Every panel: x = log(GSP loading) over the same limits, y = cumulative %.
# Run from supplementary/:  Rscript --vanilla note1/cdf_scale_compare.R
if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
source("common.R")
source(file.path(FIG_ROOT, "scripts", "vt_paths.R"))
source(file.path(FIG_ROOT, "figure3", "scripts", "nkap_variants.R"))   # classify_variants, VARIANT_*

# ---- main-figure classes ---------------------------------------------
sf <- fread(vt_path("SF3B1_full_info"), data.table = FALSE) %>%           # Fig 2d
  transmute(logk = logk10, class = plot_color)
nk <- classify_variants(fread(vt_path("NKAP_full_plus_variants"), data.table = FALSE)) %>%  # Fig 3a
  transmute(logk = logk13, class = variant_class)
Y_LOF <- "YLPM1 frameshift/early stop"; Y_NOT <- "No YLPM1 frameshift/early stop"   # Fig 4a
yl <- fread(vt_path("YLPM1_full_plus_variants"), data.table = FALSE) %>%
  transmute(logk = logk30,
            class = ifelse(!is.na(proteinChange) & grepl("fs|\\*|splice", proteinChange), Y_LOF, Y_NOT))

# ---- IDH1 (LGG) / SPOP (PRAD) -----------------------------------------
L <- fread(file.path(FIG_ROOT, "shared_cache", "gsp_loadings.tsv"), select = c("rail_id", "k1", "k12"),
           data.table = FALSE)
meta <- fread(file.path(QROOT, "inputs/tcga_v2_samples.tsv"), select = c("rail_id", "tcga_barcode", "gdc_cases.project.project_id"),
              data.table = FALSE) %>%
  transmute(rail_id, cohort = gdc_cases.project.project_id,
            sampleId = vapply(strsplit(tcga_barcode, "-"), function(x) {
              x[4] <- sub("[A-Za-z]+$", "", x[4]); paste(x[1:4], collapse = "-") }, character(1))) %>%
  filter(rail_id %in% L$rail_id)
am <- fread(file.path(QROOT, "mutation_processing/merged_grch38_with_am.tsv"), select = c("sampleId", "hugoGeneSymbol"), data.table = FALSE)
meta <- meta %>% filter(sampleId %in% am$sampleId) %>% left_join(L, by = "rail_id")
home <- function(gene, cohort, kcol) meta %>% filter(cohort == !!cohort) %>%
  transmute(logk = log(.data[[kcol]]),
            class = ifelse(sampleId %in% am$sampleId[am$hugoGeneSymbol == gene],
                           paste(gene, "mutant"), paste(gene, "wild-type")))
ih <- home("IDH1", "TCGA-LGG", "k1")
sp <- home("SPOP", "TCGA-PRAD", "k12")

PANELS <- list(
  list(d = sf, gsp = 10, title = "SF3B1 (pan-cancer)", col = PAL$sf3b1_mut),
  list(d = nk, gsp = 13, title = "NKAP (pan-cancer)",  col = VARIANT_COL),
  list(d = yl, gsp = 30, title = "YLPM1 (pan-cancer)", col = setNames(c("#E74C3C", "#34495E"), c(Y_LOF, Y_NOT))),
  list(d = ih, gsp = 1,  title = "IDH1 (TCGA-LGG)",   col = c("IDH1 mutant" = "#E74C3C", "IDH1 wild-type" = "#34495E")),
  list(d = sp, gsp = 12, title = "SPOP (TCGA-PRAD)",  col = c("SPOP mutant" = "#E74C3C", "SPOP wild-type" = "#34495E")))
XLIM <- range(unlist(lapply(PANELS, function(p) p$d$logk)), na.rm = TRUE)
cat(sprintf("shared x range: %.2f to %.2f\n", XLIM[1], XLIM[2]))

panel <- function(p) {
  d <- p$d %>% filter(!is.na(logk), class %in% names(p$col)) %>%
    mutate(class = factor(class, levels = names(p$col))) %>%
    arrange(logk) %>% group_by(class) %>% mutate(cumulative = row_number() / n() * 100) %>% ungroup()
  n <- d %>% count(class, .drop = FALSE)
  print(as.data.frame(n))
  ggplot(d, aes(logk, cumulative, colour = class)) +
    geom_step(linewidth = 0.55) +
    scale_colour_manual(values = p$col, drop = TRUE, name = NULL,
                        labels = gene_md(sprintf("%s (n = %d)", n$class, n$n)) %>% setNames(n$class)) +
    scale_x_continuous(breaks = seq(-15, 0, 5)) +
    scale_y_continuous(labels = function(x) paste0(x, "%"), breaks = c(0, 25, 50, 75, 100)) +
    coord_cartesian(xlim = XLIM) +
    labs(x = sprintf("log(GSP #%d)", p$gsp), y = "Cumulative percentage", title = gene_md(p$title)) +
    theme_publication() +
    theme(plot.title = md_title(),
          legend.text = ggtext::element_markdown(size = SIZE_LEGEND * 0.85, family = FONT_FAMILY),
          # Below the panel: the pan-cancer wild-type curves rise steeply at the
          # left, so there is no empty corner on a shared scale.
          legend.position = "bottom", legend.direction = "vertical",
          legend.justification = c(0, 1), legend.key.height = unit(8, "pt"),
          legend.margin = margin(0, 0, 0, 0)) +
    guides(colour = guide_legend(ncol = 1))
}
ps <- lapply(PANELS, panel)
fig <- wrap_plots(ps, nrow = 2) + plot_annotation_publication() &
  theme(plot.tag = element_text(size = SIZE_PANEL_LABEL, face = "bold", family = FONT_FAMILY))
save_supp(fig, "note1/note1_cdf_scale_compare", 11, 8.2)
