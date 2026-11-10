# Supplementary Figure S1. Representative single-gene Quetzal plots.
# Left: structure plots of a gene's TCGA factorisation over a subset of
# cohorts. Right: PSI of the factor's unique (DE) junctions -- CDFs over
# the cohort(s) of interest vs matched GTEx tissue for the cohort-specific
# genes, and a read-length boxplot of the top DE junction for the
# read-length genes.
# Source: tcga/gene_level_poisson/showcase_genes/supplementary_plots
# (Snakefile + scripts/make_{structure_plot,cdf_plot,psi_box}.R +
# config.yaml); the logic of each script is ported here as a function
# returning a ggplot, restyled to the paper theme.
if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
source("common.R")
suppressPackageStartupMessages({ library(yaml); library(fastTopics) })
SP  <- file.path(QROOT, "gene_level_analysis/showcase_plots")
cfg <- read_yaml(file.path(SP, "config.yaml"))
P   <- function(p) normalizePath(file.path(SP, p), mustWork = FALSE)
TH  <- cfg$thresholds
d3  <- ggsci::pal_d3("category10")(10)

tmeta <- fread(P(cfg$paths$tcga_metadata), data.table = FALSE) %>%
  transmute(rail_id = as.integer(rail_id), project_id = gdc_cases.project.project_id,
            avg_len = as.numeric(avg_len))
gmeta <- fread(P(cfg$paths$gtex_metadata), data.table = FALSE) %>%
  transmute(rail_id = as.integer(rail_id), study)
entry <- function(g) if (g %in% names(cfg$read_length))
  list(grp = "read_length", e = cfg$read_length[[g]]) else
  list(grp = "tcga_cohort", e = cfg$tcga_cohort[[g]])
gdir  <- function(e, g) P(file.path(cfg$paths$tcga_chr_dir, e$chr, g))
read_psi <- function(dir, e, g) {
  d <- fread(P(file.path(dir, e$chr, g, "psi_table.tsv")), sep = "\t", data.table = FALSE)
  rownames(d) <- d$V1; d$V1 <- NULL; d
}
de_juncs <- function(res, de, fact) {
  k <- as.integer(sub("^k", "", fact))
  pass <- de$svalue[, k] < TH$sval_thresh & de$postmean[, k] > TH$lfc_thresh &
          res$F[, k] > TH$abs_thresh
  list(names = rownames(res$F)[pass],
       top = rownames(res$F)[order(de$postmean[, k] * pass, decreasing = TRUE)[1]])
}

## ---- structure plot ---------------------------------------------------
struct_plot <- function(g, title) {
  le <- entry(g); e <- le$e; coh <- unlist(e$cohorts)
  res <- readRDS(file.path(gdir(e, g), "res.RDS"))
  md <- tmeta %>% filter(rail_id %in% as.integer(rownames(res$L)), project_id %in% coh)
  short <- sub("^TCGA-", "", md$project_id)
  md$lab <- if (le$grp == "read_length") paste0(short, "-", round(md$avg_len), "bp") else short
  keep <- rownames(res$L) %in% as.character(md$rail_id)
  r <- res; r$L <- r$L[keep, , drop = FALSE]
  if (!is.null(r$Ln)) r$Ln <- r$Ln[keep, , drop = FALSE]
  if (!is.null(r$s) && length(r$s) == length(keep)) r$s <- r$s[keep]
  grp <- md$lab[match(rownames(r$L), as.character(md$rail_id))]
  lv <- unique(grp)
  if (le$grp == "read_length") {
    co <- sub("^TCGA-", "", coh)
    lv <- lv[order(match(sub("-\\d+bp$", "", lv), co), as.integer(sub(".*-(\\d+)bp$", "\\1", lv)))]
  } else lv <- sub("^TCGA-", "", coh)
  nf <- ncol(r$F)
  pal <- if (nf <= 10) ggsci::pal_d3("category10")(nf) else ggsci::pal_d3("category20")(nf)
  set.seed(1)   # structure_plot orders samples stochastically -- seed it
  # read-length rows have 12 groups, some tiny (PRAD-50bp): a wider gap
  # keeps their 45-deg labels from touching
  fastTopics::structure_plot(r, grouping = factor(grp, levels = lv),
                             gap = if (le$grp == "read_length") 220 else 100,
                             colors = pal, n = 5000, verbose = FALSE) +
    labs(title = gene_md(title), y = "Topic proportion", fill = "Topic", colour = "Topic") +
    theme_publication() +
    theme(plot.title = md_title(hjust = 0), axis.line.x = element_blank(),
          axis.ticks.x = element_blank(), axis.title.x = element_blank(),
          legend.key.size = unit(2.4, "mm"),
          legend.text = element_text(size = SIZE_BASE * 0.75),
          legend.title = element_text(size = SIZE_BASE * 0.8)) +
    angled_x_labels(size = SIZE_BASE * 0.8)
}

## ---- DE-junction CDFs (cohort genes) ------------------------------------
cdf_plot <- function(g) {
  e <- entry(g)$e
  coi <- trimws(strsplit(e$cohort_of_interest, ",")[[1]])
  tissue <- unique(unlist(lapply(coi, function(x) cfg$tcga_to_gtex[[x]])))
  res <- readRDS(file.path(gdir(e, g), "res.RDS")); de <- readRDS(file.path(gdir(e, g), "de_res.RDS"))
  j <- de_juncs(res, de, e$factor)$names
  tp <- read_psi(cfg$paths$tcga_psi_dir, e, g); gp <- read_psi(cfg$paths$gtex_psi_dir, e, g)
  long <- function(ps, rids, lab) bind_rows(lapply(j, function(jj) {
    col <- names(ps)[endsWith(names(ps), paste0(":", jj))][1]
    if (is.na(col)) return(NULL)
    v <- suppressWarnings(as.numeric(ps[as.character(rids), col])); v <- v[is.finite(v)]
    if (!length(v)) return(NULL)
    data.frame(psi = v, source = lab, junction = jj)
  }))
  tl <- paste("TCGA", sub("^TCGA-", "", coi)); gl <- paste("GTEx", tissue)
  d <- bind_rows(
    bind_rows(lapply(seq_along(coi), function(i)
      long(tp, tmeta$rail_id[tmeta$project_id == coi[i] & tmeta$rail_id %in% as.integer(rownames(tp))], tl[i]))),
    long(gp, gmeta$rail_id[gmeta$study == tissue & gmeta$rail_id %in% as.integer(rownames(gp))], gl)) %>%
    mutate(source = factor(source, levels = c(tl, gl)),
           junction = factor(junction, levels = sort(unique(junction))))
  cols <- setNames(c(d3[c(1, 2)][seq_along(tl)], d3[3]), c(tl, gl))
  ggplot(d, aes(psi, colour = source)) +
    stat_ecdf(geom = "step", linewidth = 0.45) +
    facet_wrap(~ junction, ncol = min(2, nlevels(d$junction))) +
    scale_colour_manual(values = cols, name = NULL) +
    guides(colour = guide_legend(nrow = if (length(cols) > 2) 2 else 1)) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    scale_x_continuous(breaks = scales::breaks_extended(3)) +
    labs(x = "PSI", y = "Cumulative fraction",
         title = paste0(e$factor, "-unique junctions")) +
    theme_publication() +
    theme(plot.title = element_text(size = SIZE_BASE, face = "plain", hjust = 0.5),
          strip.background = element_blank(),
          strip.text = element_text(size = SIZE_BASE * 0.8, face = "plain"),
          legend.position = "bottom", legend.text = element_text(size = SIZE_BASE * 0.8),
          legend.key.height = unit(3, "mm"),
          panel.spacing = unit(1, "lines"))
}

## ---- read-length boxplot (read-length genes) -----------------------------
box_plot <- function(g) {
  e <- entry(g)$e; coh <- unlist(e$cohorts)
  res <- readRDS(file.path(gdir(e, g), "res.RDS")); de <- readRDS(file.path(gdir(e, g), "de_res.RDS"))
  tj <- de_juncs(res, de, e$factor)$top
  ps <- read_psi(cfg$paths$tcga_psi_dir, e, g)
  col <- names(ps)[endsWith(names(ps), paste0(":", tj))][1]
  d <- data.frame(rail_id = as.integer(rownames(ps)), psi = suppressWarnings(as.numeric(ps[[col]]))) %>%
    filter(is.finite(psi)) %>%
    inner_join(tmeta %>% filter(project_id %in% coh), by = "rail_id") %>%
    mutate(rl = case_when(round(avg_len) <= 50 ~ "≤ 50 bp", round(avg_len) >= 75 ~ "≥ 75 bp")) %>%
    filter(!is.na(rl))
  ggplot(d, aes(rl, psi)) +
    geom_boxplot(width = 0.5, fill = "white", outlier.size = 0.4, linewidth = 0.3) +
    labs(x = "Sample read length", y = "PSI",
         title = gene_md(sprintf("%s (%s) top DE junction<br>%s", g, e$factor, tj))) +
    theme_publication() + theme(plot.title = md_title(hjust = 0.5, size = SIZE_BASE))
}

# wrap_elements(): the structure plot's tall angled tick labels would otherwise push
# the right plot's x title and legend far below its axis.
rowp <- function(left, right) (left | wrap_elements(full = right)) + plot_layout(widths = c(2.3, 1.35))
hdr  <- function(t) wrap_elements(full = grid::textGrob(t, gp = grid::gpar(
  fontsize = SIZE_TITLE + 2, fontface = "bold", fontfamily = FONT_FAMILY)))
R <- list(
  hdr("Examples of cohort/tissue-specific splicing captured by Quetzal"),
  rowp(struct_plot("C2",    "TCGA factorisation of <i>C2</i> (k10 is prostate-specific)"),    cdf_plot("C2")),
  rowp(struct_plot("MVK",   "TCGA factorisation of MVK (k6 is kidney-specific)"),      cdf_plot("MVK")),
  rowp(struct_plot("NHEJ1", "TCGA factorisation of NHEJ1 (k2 is colon/rectum-specific)"), cdf_plot("NHEJ1")),
  hdr("Examples of read-length associated variation captured by Quetzal"),
  rowp(struct_plot("MAP2K2", "TCGA factorisation of MAP2K2 (k5 only appears in 75/76 bp samples)"), box_plot("MAP2K2")),
  rowp(struct_plot("KLHDC4", "TCGA factorisation of KLHDC4 (k6 only appears in 75/76 bp samples)"), box_plot("KLHDC4")))
fig <- wrap_plots(R, ncol = 1, heights = c(0.25, 1.4, 1.4, 2.2, 0.25, 1.4, 1.4))
save_supp(fig, "figS1", 8.5, 11)
