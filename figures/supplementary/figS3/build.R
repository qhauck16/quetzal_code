# Supplementary Figure S3. Cancer-gene mis-splicing found by single-gene
# Quetzal fits: (a) TCF12 GSP k10 in TCGA-LGG and TCGA-GBM, (b) SRC k8 in
# TCGA-COAD, (c) SETD2 k2 in TCGA-LAML. Each panel = loading violin
# (cohort vs every other cohort, Wilcoxon), one PSI CDF per DE junction
# (cohort, its matched GTEx tissue, pooled other TCGA, pooled other
# GTEx) and the gene model with the DE junctions drawn as arcs, numbered
# to match the CDF titles; red asterisk = non-coding junction.
# Same drawing code as main Figure 1d (../figure1/compose.R, EZH2).
# gene_data.RDS per gene from cache/s2/run.sh (make_gene_outputs.R).
if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
source("common.R")
suppressPackageStartupMessages(library(yaml))

LDE <- new.env()
sys.source(file.path(FIG_ROOT, "figure1/panel_d_e/scripts/lib.R"), envir = LDE)
`%||%` <- function(a, b) if (is.null(a)) b else a
load_gd <- function(dir) {
  gd <- readRDS(file.path("cache/s2", dir, "gene_data.RDS"))
  gd$paths$tcga_metadata <- file.path(QROOT, "inputs/tcga_v2_samples.tsv")
  gd$paths$gtex_metadata <- file.path(QROOT, "inputs/gtexv2_metadata.tsv") %||% gd$paths$gtex_metadata
  gd$paths$psi_dir       <- file.path(QROOT, "leafcutter_analyses/psi_files")      %||% gd$paths$psi_dir
  gd$jc_tsv <- file.path("cache/s2", dir, "junctions_classified.tsv")
  gd
}

# One colour per TCGA cohort, as in the original figure (LGG blue, GBM teal,
# COAD orange, LAML pink); the cohort's violin and its CDF curve share it.
COHORT_COL <- c("TCGA-LGG" = "#2E6FB7", "TCGA-GBM" = "#6FB3C3",
                "TCGA-COAD" = "#F28E2B", "TCGA-LAML" = "#E377C2")
cohort_col <- function(gd) unname(COHORT_COL[gd$cancer_type] %||% PAL$tcga_test)
CDF_COL <- function(gd) setNames(
  c(PAL$gtex, "#a8d5a2", cohort_col(gd), PAL$tcga_other),
  c(paste0("GTEx-", gd$gtex_study), "Other GTEx", gd$cancer_type, "Other TCGA"))

P_SIZE <- SIZE_BASE * 0.85 / .pt   # every p-value label on the page
fmt_p <- function(p) {
  if (is.na(p)) return("")
  if (p == 0) return("p < 1e-300")
  if (p >= 0.01) return(sprintf("p = %.2f", p))
  e <- floor(log10(p)); sprintf("p = %.1fe%d", p / 10^e, e)
}

# -- violin: same as fig 1d ------------------------------------------------
violin <- function(gd) {
  mx <- max(gd$res$L[, gd$fact_idx])
  v <- LDE$plot_loading_violins(gd)
  # Drop ggpubr's bracket (its "p < 2.22e-16" is a floor and was drawn at
  # ~10 pt, larger than every other label): redraw it with the exact
  # Wilcoxon p in the same format/size as the CDF p-values.
  v$layers <- Filter(function(l) !inherits(l$geom, "GeomSignif"), v$layers)
  lv <- levels(v$data$group)
  pw <- wilcox.test(v$data$loading[v$data$group == lv[1]],
                    v$data$loading[v$data$group == lv[2]])$p.value
  yb <- mx * 1.06; tk <- mx * 0.03
  v +
    annotate("path", x = c(1, 1, 2, 2), y = c(yb - tk, yb, yb, yb - tk), linewidth = 0.3) +
    annotate("text", x = 1.5, y = yb + tk, label = fmt_p(pw), vjust = 0,
             size = P_SIZE, family = FONT_FAMILY) +
    scale_fill_manual(values = setNames(c(cohort_col(gd), PAL$tcga_other),
                                        c(gd$cancer_type, "Other\nCohorts"))) +
    coord_cartesian(ylim = c(-0.24 * mx, mx * 1.20), clip = "off") +
    scale_y_continuous(breaks = pretty(c(0, mx), 3)) +
    labs(title = NULL) +
    theme_publication() +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          plot.margin = margin(12, 4, 2, 2),
          axis.text.x = element_text(angle = 0, hjust = 0.5, size = SIZE_BASE * 0.9))
}

# -- CDFs with the four cohort pools (lib's plot_psi_cdfs_full_static) ----
cdf_data <- function(gd) {
  tm <- fread(gd$paths$tcga_metadata, data.table = FALSE)
  gm <- fread(gd$paths$gtex_metadata, data.table = FALSE)
  o_tcga <- tm$rail_id[tm$gdc_cases.project.project_id != gd$cancer_type &
                         tm$rail_id %in% rownames(gd$dat)]
  o_gtex <- gm$rail_id[gm$study != gd$gtex_study & gm$rail_id %in% rownames(gd$gtex_dat)]
  rd <- function(src, ids) {
    x <- fread(file.path(gd$paths$psi_dir, src, gd$chr, gd$gene_name, "psi_table.tsv"),
               data.table = FALSE)
    x[x$V1 %in% ids, , drop = FALSE]
  }
  ot <- rd("tcga", o_tcga); og <- rd("gtex", o_gtex)
  lab <- names(CDF_COL(gd))
  lapply(gd$junctions_to_compare, function(j) {
    col <- function(m) which(grepl(j, colnames(m)))[1]
    v <- list(gd$gtex_norm_dat[, col(gd$gtex_norm_dat)], og[, col(og)],
              gd$tcga_norm_dat[, col(gd$tcga_norm_dat)], ot[, col(ot)])
    data.frame(psi = unlist(v), cohort = factor(rep(lab, lengths(v)), levels = lab))
  })
}

cdf_plots <- function(gd, legend = TRUE) {
  dl <- cdf_data(gd)
  lapply(seq_along(dl), function(i) {
    j <- gd$junctions_to_compare[i]
    ggplot(dl[[i]], aes(psi, colour = cohort)) +
      stat_ecdf(geom = "step", linewidth = 0.6) +
      scale_colour_manual(values = CDF_COL(gd), name = NULL, drop = FALSE) +
      scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
      scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
      annotate("text", x = 1, y = 0.45, label = if (isTRUE(gd$junction_coefs[[j]] > 0))
                 fmt_p(gd$junction_pvals[[j]]) else "ns",
               hjust = 1, size = P_SIZE, family = FONT_FAMILY) +
      labs(title = paste("Junction", i), x = "Junction PSI", y = "Empirical CDF") +
      theme_publication() +
      theme(legend.position = if (legend && i == length(dl)) "right" else "none",
            legend.key.size = unit(3, "mm"),
            legend.text = element_text(size = SIZE_BASE * 0.8),
            plot.title = element_text(size = SIZE_STRIP, hjust = 0.5, face = "bold"),
            plot.margin = margin(2, 6, 6, 4))
  })
}

# -- gene model + numbered arcs + non-coding stars (fig 1d EZH2 code) ----
gene_arcs <- function(gd, arc_col = "#1F77B4FF") {
  js <- gd$junctions_to_compare
  arcs <- do.call(rbind, lapply(seq_along(js), function(i) {
    pp <- as.numeric(strsplit(js[i], "-")[[1]])
    data.frame(x_start = pp[1], x_end = pp[2],
               width = gd$res$F[js[i], gd$fact_idx], number = i)
  }))
  jc <- read.table(gd$jc_tsv, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  jc$j1 <- paste0(as.integer(sub("-.*", "", jc$junction)) + 1, "-", sub(".*-", "", jc$junction))
  nc <- setNames(!jc$Coding, jc$j1)[paste0(arcs$x_start, "-", arcs$x_end)]
  arcs$noncoding <- !is.na(nc) & nc

  span <- gd$max_coord - gd$min_coord
  p <- gene_body_gg(gd$t_models, gd$chr, gd$min_coord, gd$max_coord,
                    gene_label = gd$gene_name, exon_h = 0.35, n_arrows = 30,
                    label_size = 3.4, label_family = FONT_FAMILY) +
    scale_y_continuous(limits = c(0.3, 2.55), expand = c(0, 0))
  amp <- 0.4; mw <- max(arcs$width)
  for (k in seq_len(nrow(arcs))) {
    xs <- seq(arcs$x_start[k], arcs$x_end[k], length.out = 200)
    p <- p + geom_path(data = data.frame(x = xs,
                         y = 1.15 + sin((xs - xs[1]) / (xs[200] - xs[1]) * pi) * amp),
                       aes(x, y), inherit.aes = FALSE, colour = arc_col,
                       linewidth = 0.15 + 0.75 * arcs$width[k] / mw)
  }
  num <- data.frame(x = (arcs$x_start + arcs$x_end) / 2, y = 1.15 + amp + 0.40,
                    label = as.character(arcs$number))
  o <- order(num$x); x <- num$x[o]; sep <- 0.022 * span
  for (it in 1:20) for (k in seq_along(x)[-1]) if (x[k] - x[k - 1] < sep) {
    d <- (sep - (x[k] - x[k - 1])) / 2; x[k - 1] <- x[k - 1] - d; x[k] <- x[k] + d
  }
  num$x[o] <- x
  star <- num[arcs$noncoding, , drop = FALSE]
  p + geom_text(data = num, aes(x, y, label = label), inherit.aes = FALSE,
                size = 2.8, family = FONT_FAMILY, fontface = "bold", colour = "grey15") +
    { if (nrow(star)) geom_text(data = star, aes(x, y + 0.35), label = "*",
                                inherit.aes = FALSE, size = 5, family = FONT_FAMILY,
                                fontface = "bold", colour = "red") } +
    theme(plot.margin = margin(2, 30, 4, 6))
}

# Fixed panel widths (inches) between two flexible spacers, so each
# violin + CDF row is centred over its full-width gene track (Quinn: no
# dead white band on the right). CDFs use aspect.ratio = 1.
VIOLIN_W <- 1.35; CDF_W <- 1.45
# Panel c's CDFs are a 2 x 2 block (~2.9 in tall); its violin is made
# roughly square (VIOLIN_C_W wide, vertically centred) instead of a tall
# 1.35 in sliver stretched to the block height.
VIOLIN_C_W <- 2.1
row_of <- function(gd, cdf_ncol = NULL, legend = TRUE) {
  cp <- cdf_plots(gd, legend)
  if (is.null(cdf_ncol)) {
    n <- length(cp)
    wrap_plots(c(list(plot_spacer(), violin(gd)), cp, list(plot_spacer())), nrow = 1,
               widths = unit(c(1, VIOLIN_W, rep(CDF_W, n), 1),
                             c("null", rep("in", n + 1), "null")))
  } else {
    cdfs <- wrap_plots(cp, ncol = cdf_ncol)
    vio  <- plot_spacer() / violin(gd) / plot_spacer() + plot_layout(heights = c(1, 6, 1))
    (plot_spacer() | vio | cdfs | plot_spacer()) +
      plot_layout(widths = unit(c(1, VIOLIN_C_W, cdf_ncol * CDF_W + 0.5, 1),
                                c("null", "in", "in", "null")))
  }
}

tcf_lgg <- load_gd("TCF12_LGG/TCF12"); tcf_gbm <- load_gd("TCF12_GBM/TCF12")
src     <- load_gd("SRC_COAD/SRC");   setd2   <- load_gd("SETD2_LAML/SETD2")

pa <- wrap_elements(full = row_of(tcf_lgg) / row_of(tcf_gbm) / gene_arcs(tcf_lgg) +
                      plot_layout(heights = c(1, 1, 0.45)))
pb <- wrap_elements(full = row_of(src) / gene_arcs(src) + plot_layout(heights = c(1, 0.45)))
pc <- wrap_elements(full = row_of(setd2, cdf_ncol = 2) / gene_arcs(setd2) +
                      plot_layout(heights = c(2, 0.45)))

key <- ggplot() +
  annotate("text", x = 0, y = 0, label = "*", colour = "red", size = 5,
           fontface = "bold", family = FONT_FAMILY, vjust = 0.7) +
  annotate("text", x = 0.02, y = 0, label = "Non-coding junction", hjust = 0,
           size = SIZE_BASE * 0.85 / .pt, family = FONT_FAMILY) +
  coord_cartesian(xlim = c(0, 1), ylim = c(-1, 1)) + theme_void()

p <- pa / pb / pc / wrap_elements(key, ignore_tag = TRUE) +
  plot_layout(heights = c(2.4, 1.4, 2.6, 0.12)) +
  plot_annotation_publication() &
  # theme.R's annotation theme does not reach the tags of nested plots
  # (they rendered regular weight); set the bold tag on every plot
  theme(plot.tag = element_text(size = SIZE_PANEL_LABEL, face = "bold", family = FONT_FAMILY))
save_supp(p, "figS3", 8.5, 11)
