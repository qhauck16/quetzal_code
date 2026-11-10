.libPaths('/scratch/midway2/qhauck/conda_env/rstudio-server/lib/R/library/')

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(bit64)
  library(yaml)
  library(ggbio)
  library(patchwork)
  library(GenomicRanges)
  library(rtracklayer)
  library(ggsci)
  library(ggpubr)
})

# ---- Argument parsing -----------------------------------------------------

parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  flags <- args[grepl("^--", args)]
  out <- list()
  for (f in flags) {
    idx <- which(args == f)
    if (length(idx) == 0) next
    out[[sub("^--", "", f)]] <- args[idx + 1]
  }
  out
}

require_arg <- function(args, name) {
  if (is.null(args[[name]]) || !nzchar(args[[name]]))
    stop(sprintf("Missing required --%s", name))
  args[[name]]
}

# ---- TCGA -> matched GTEx tissue lookup -----------------------------------

tcga_to_gtex_lookup <- c(
  "TCGA-ACC"   = "ADRENAL_GLAND",
  "TCGA-BLCA"  = "BLADDER",
  "TCGA-BRCA"  = "BREAST",
  "TCGA-CESC"  = "CERVIX_UTERI",
  "TCGA-COAD"  = "COLON",
  "TCGA-DLBC"  = "SPLEEN",
  "TCGA-ESCA"  = "ESOPHAGUS",
  "TCGA-GBM"   = "BRAIN",
  "TCGA-HNSC"  = "SALIVARY_GLAND",
  "TCGA-KICH"  = "KIDNEY",
  "TCGA-KIRC"  = "KIDNEY",
  "TCGA-KIRP"  = "KIDNEY",
  "TCGA-LAML"  = "BLOOD",
  "TCGA-LGG"   = "BRAIN",
  "TCGA-LIHC"  = "LIVER",
  "TCGA-LUAD"  = "LUNG",
  "TCGA-LUSC"  = "LUNG",
  "TCGA-MESO"  = "LUNG",
  "TCGA-OV"    = "OVARY",
  "TCGA-PAAD"  = "PANCREAS",
  "TCGA-PCPG"  = "ADRENAL_GLAND",
  "TCGA-PRAD"  = "PROSTATE",
  "TCGA-READ"  = "COLON",
  "TCGA-SARC"  = "MUSCLE",
  "TCGA-SKCM"  = "SKIN",
  "TCGA-STAD"  = "STOMACH",
  "TCGA-TGCT"  = "TESTIS",
  "TCGA-THCA"  = "THYROID",
  "TCGA-UCEC"  = "UTERUS",
  "TCGA-UCS"   = "UTERUS"
)

palet <- c(
  'red', 'blue', 'green', 'purple', 'orange', 'cyan2', 'goldenrod4', 'pink',
  'magenta', 'grey', 'brown', 'darkolivegreen',
  'skyblue2', 'yellow', 'chartreuse4', 'black'
)

round_down_to_thousand <- function(value) floor(value   / 1000) * 1000
round_up_to_thousand   <- function(value) ceiling(value / 1000) * 1000

# ---- Core data loader -----------------------------------------------------
# Loads everything load_gene_data() in updated_single_gene_figures.Rmd loads,
# but driven by an explicit `paths` list (read from config.yaml).

load_gene_data <- function(gene_name, chr, fact, cancer_type,
                           paths, thresholds) {

  tcga_metadata    <- fread(paths$tcga_metadata,    data.table = FALSE)
  gtex_metadata    <- fread(paths$gtex_metadata)
  mane_transcripts <- fread(paths$mane_transcripts, data.table = FALSE)
  gencode          <- readRDS(paths$gencode_granges)

  tcga_dir <- file.path(paths$tcga_chr_dir, chr, gene_name)
  res     <- readRDS(file.path(tcga_dir, "res.RDS"))
  de_res  <- readRDS(file.path(tcga_dir, "de_res.RDS"))
  corr_info <- fread(file.path(tcga_dir, "results.tsv"))

  dat <- as.matrix(read.table(file.path(tcga_dir, "dat.tsv")))
  colnames(dat) <- rownames(res$F)

  gene_data <- gencode[gencode$gene_name == gene_name & gencode$type %in% c("exon")]
  mane_tx <- mane_transcripts$Ensembl_nuc[mane_transcripts$symbol == gene_name]
  if (length(mane_tx) > 0) {
    mane_exons <- gene_data[gene_data$transcript_id %in% mane_tx]
    if (length(mane_exons) > 0) gene_data <- mane_exons
  }
  t_models  <- rtracklayer::split(gene_data, gene_data$transcript_id)
  min_coord <- round_down_to_thousand(min(c(start(gene_data@ranges), end(gene_data@ranges))))
  max_coord <- round_up_to_thousand(max(c(start(gene_data@ranges), end(gene_data@ranges))))
  gene_strand <- as.character(t_models@unlistData@strand[1])

  gtex_dat_path <- file.path(paths$gtex_dat_dir, chr, "FastTopics_output", gene_name, "dat.tsv")
  gtex_dat <- fread(gtex_dat_path, data.table = FALSE) %>% column_to_rownames("rownames")
  gtex_dat <- gtex_dat[, colnames(gtex_dat) %in% colnames(dat)]
  for (col in setdiff(colnames(dat), colnames(gtex_dat))) gtex_dat[[col]] <- 0
  gtex_dat <- as.matrix(gtex_dat[, c(colnames(dat), setdiff(colnames(gtex_dat), colnames(dat)))])

  gtex_study   <- tcga_to_gtex_lookup[cancer_type]
  dat_cancer   <- dat[rownames(dat) %in% filter(tcga_metadata, gdc_cases.project.project_id == cancer_type)$rail_id, ]
  dat_gtex_tissue <- gtex_dat[rownames(gtex_dat) %in% filter(gtex_metadata, study == gtex_study)$rail_id, ]

  fact_idx <- as.integer(sub("^k", "", fact))
  junctions_to_compare <- rownames(res$F)[
    de_res$svalue[, fact_idx] < thresholds$sval_thresh &
    de_res$postmean[, fact_idx] > thresholds$lfc_thresh &
    res$F[, fact_idx]          > thresholds$abs_thresh
  ]

  gtex_norm_dat <- fread(file.path(paths$psi_dir, "gtex", chr, gene_name, "psi_table.tsv"),
                         data.table = FALSE) %>%
    filter(V1 %in% rownames(dat_gtex_tissue)) %>% column_to_rownames("V1")
  tcga_norm_dat <- fread(file.path(paths$psi_dir, "tcga", chr, gene_name, "psi_table.tsv"),
                         data.table = FALSE) %>%
    filter(V1 %in% rownames(dat_cancer)) %>% column_to_rownames("V1")

  junction_pvals    <- c()
  junction_coefs    <- c()
  junction_z        <- c()   # Wald z; log10 P is recomputed from it when P underflows to 0
  junction_psi_difs <- c()
  junction_psi_fcs  <- c()
  for (junc in junctions_to_compare) {
    clust <- str_split_i(names(tcga_norm_dat)[grepl(junc, names(tcga_norm_dat))], ":", 1)[1]
    tcga_juncs <- str_split_i(names(tcga_norm_dat)[grepl(clust, names(tcga_norm_dat))], ":", 2)
    gtex_juncs <- str_split_i(names(gtex_norm_dat)[grepl(clust, names(gtex_norm_dat))], ":", 2)

    tcga_vals <- dat[rownames(dat) %in% rownames(tcga_norm_dat), colnames(dat) == junc, drop = FALSE]
    tcga_sums <- rowSums(dat[rownames(dat) %in% rownames(tcga_norm_dat),
                             colnames(dat) %in% tcga_juncs, drop = FALSE])
    gtex_vals <- gtex_dat[rownames(gtex_dat) %in% rownames(gtex_norm_dat),
                          colnames(gtex_dat) == junc, drop = FALSE]
    gtex_sums <- rowSums(gtex_dat[rownames(gtex_dat) %in% rownames(gtex_norm_dat),
                                  colnames(gtex_dat) %in% gtex_juncs, drop = FALSE])

    cluster_test <- data.frame(
      junc_count  = c(tcga_vals, gtex_vals),
      cluster_sum = c(tcga_sums, gtex_sums),
      study = c(rep("TCGA", length(tcga_vals)), rep("GTEx", length(gtex_vals)))
    ) %>% filter(cluster_sum > 0)

    if (length(unique(cluster_test$study)) < 2) {
      junction_z        <- c(junction_z, NA)
      junction_coefs    <- c(junction_coefs, NA)
      junction_pvals    <- c(junction_pvals, NA)
      junction_psi_difs <- c(junction_psi_difs, NA)
      junction_psi_fcs  <- c(junction_psi_fcs, NA)
      next
    }

    # Beta-binomial GLM, the SAME model as the genome-wide search
    # (tcga/gene_level_poisson/cancer_specific_factors.R fit_bb): reads for
    # the junction vs total reads in its cluster ~ study, one shared
    # overdispersion (aod::betabin, random = ~ 1). Replaced a quasibinomial
    # GLM on 2026-10-08 so the paper uses a single junction-level test.
    cluster_test$study <- factor(cluster_test$study, levels = c("GTEx", "TCGA"))
    co <- tryCatch({
      fit <- suppressWarnings(aod::betabin(
        cbind(junc_count, cluster_sum - junc_count) ~ study,
        random = ~ 1, data = cluster_test))
      summary(fit)@Coef
    }, error = function(e) NULL)
    if (is.null(co) || !"studyTCGA" %in% rownames(co)) {
      glm_res <- c(NA_real_, NA_real_); zval <- NA_real_
    } else {
      p_col   <- grep("^Pr|p.value|pvalue", colnames(co), ignore.case = TRUE, value = TRUE)[1]
      glm_res <- c(co["studyTCGA", "Estimate"], co["studyTCGA", p_col])
      zval    <- co["studyTCGA", "Estimate"] / co["studyTCGA", "Std. Error"]
    }
    junction_z <- c(junction_z, zval)

    junction_coefs <- c(junction_coefs, glm_res[1])
    junction_pvals <- c(junction_pvals, glm_res[2])

    psi_junction_id <- paste0(clust, ":", junc)
    junction_psi_difs <- c(junction_psi_difs,
      (median(tcga_norm_dat[, psi_junction_id], na.rm = TRUE) + 1e-4) -
      (median(gtex_norm_dat[, psi_junction_id], na.rm = TRUE) + 1e-4))
    junction_psi_fcs <- c(junction_psi_fcs,
      (median(tcga_norm_dat[, psi_junction_id], na.rm = TRUE) + 1e-4) /
      (median(gtex_norm_dat[, psi_junction_id], na.rm = TRUE) + 1e-4))
  }
  names(junction_pvals) <- junctions_to_compare
  names(junction_coefs) <- junctions_to_compare
  names(junction_z)     <- junctions_to_compare
  # log10 of the two-sided Wald P, exact even where junction_pvals underflows to 0
  junction_log10p <- (pnorm(-abs(junction_z), log.p = TRUE) + log(2)) / log(10)

  list(
    gene_name = gene_name, chr = chr, fact = fact, fact_idx = fact_idx,
    cancer_type = cancer_type, gtex_study = gtex_study, gene_strand = gene_strand,
    dat = dat, gtex_dat = gtex_dat,
    res = res, de_res = de_res, corr_info = corr_info,
    dat_cancer = dat_cancer, dat_gtex_tissue = dat_gtex_tissue,
    tcga_norm_dat = tcga_norm_dat, gtex_norm_dat = gtex_norm_dat,
    junctions_to_compare = junctions_to_compare,
    junction_pvals = junction_pvals, junction_coefs = junction_coefs,
    junction_z = junction_z, junction_log10p = junction_log10p,
    junction_psi_difs = junction_psi_difs, junction_psi_fcs = junction_psi_fcs,
    t_models = t_models, min_coord = min_coord, max_coord = max_coord,
    paths = paths, thresholds = thresholds
  )
}

# ---- Junction classifications ---------------------------------------------

junction_classifications <- function(gd) {
  leaf2_class <- fread(gd$paths$leaf2_class, data.table = FALSE)
  leaf2_class$coord <- sub('^.*:', '', leaf2_class$Intron_coord)

  data.frame(junction = gd$junctions_to_compare) %>%
    mutate(start = str_split_i(junction, '-', 1),
           end   = str_split_i(junction, '-', 2)) %>%
    mutate(start = as.numeric(start) - 1) %>%
    mutate(junction = paste0(start, '-', end)) %>%
    left_join(leaf2_class, by = c("junction" = "coord")) %>%
    select(junction, Gene_name, Intron_coord, Strand, Annot, Coding, UTR, GencodePC) %>%
    mutate(
      junction_pval = gd$junction_pvals[
        match(junction, paste0(as.numeric(str_split_i(names(gd$junction_pvals), '-', 1)) - 1, '-',
                               str_split_i(names(gd$junction_pvals), '-', 2)))
      ],
      junction_coef = gd$junction_coefs[
        match(junction, paste0(as.numeric(str_split_i(names(gd$junction_coefs), '-', 1)) - 1, '-',
                               str_split_i(names(gd$junction_coefs), '-', 2)))
      ]
    )
}

# ---- Static (PDF) factor-junction arc plot --------------------------------
# Returns a ggplot (no plotly), suitable for ggsave to PDF.

plot_factor_junctions_static <- function(gd, height_scaling = 0.7,
                                         arc_thickness = 0.15,
                                         exon_rect_h = 0.4,
                                         # When NULL (default), the arc
                                         # y-baseline is
                                         # `max_y + factor_idx * 4`.
                                         # Setting to a small integer
                                         # (e.g. 2) pulls arcs down so
                                         # they sit just above the gene
                                         # body -- used at compose time
                                         # for tight panel layouts.
                                         arc_y_offset = NULL) {
  fact_col <- gd$fact_idx
  num_factors <- ncol(gd$res$L)

  factor_info <- data.frame(gd$res$F, rownames(gd$res$F))
  colnames(factor_info) <- c(seq_len(num_factors), "junction")

  factor_plot <- factor_info %>%
    mutate(start = as.numeric(str_split_i(junction, '-', 1)),
           end   = as.numeric(str_split_i(junction, '-', 2)))

  whole_line_df <- data.frame()
  for (j in seq_len(num_factors)) {
    this_df <- data.frame(
      x = factor_plot$start, y = 5,
      xend = factor_plot$end,  yend = 5,
      width = factor_plot[, j], factor = as.character(j)
    ) %>% filter(width > gd$thresholds$abs_thresh / 10)
    whole_line_df <- rbind(whole_line_df, this_df)
  }
  whole_line_df$label <- paste0(whole_line_df$x, "-", whole_line_df$xend)
  plot_df <- whole_line_df %>%
    filter(label %in% gd$junctions_to_compare,
           factor == as.character(fact_col))

  p <- ggplot() +
    ggbio::geom_alignment(gd$t_models, gap.geom = getOption("figure1.gap_geom", "arrow"), label = FALSE,
                          exon.rect.h = exon_rect_h) +
    theme_minimal() +
    scale_x_continuous(breaks = seq(gd$min_coord, gd$max_coord, by = 10000))

  if (nrow(plot_df) > 0) {
    max_y     <- 5
    max_width <- max(plot_df$width)
    for (i in seq_len(nrow(plot_df))) {
      y_base <- if (is.null(arc_y_offset)) {
        max_y + as.numeric(plot_df$factor[i]) * 4
      } else {
        max_y + arc_y_offset
      }
      x_vals <- seq(plot_df$x[i], plot_df$xend[i], length.out = 100)
      arc <- sin((x_vals - plot_df$x[i]) /
                 (plot_df$xend[i] - plot_df$x[i]) * pi)
      scaled_thickness <- arc_thickness * (plot_df$width[i] / max_width)
      y_top    <- arc * height_scaling + y_base + scaled_thickness
      y_bottom <- arc * height_scaling + y_base
      polygon_df <- data.frame(x = c(x_vals, rev(x_vals)),
                               y = c(y_top, rev(y_bottom)))
      p <- p + geom_polygon(
        data = polygon_df, aes(x = x, y = y),
        fill = palet[as.numeric(plot_df$factor[i])], alpha = 1
      )
    }
  }

  p + scale_y_continuous(breaks = c(0, 10), labels = c("", "")) +
    labs(title = paste0(gd$gene_name, " ", gd$fact, " (", gd$cancer_type, ")"))
}

# ---- TCGA vs GTEx PSI ECDFs (static PDF) ----------------------------------

plot_psi_cdfs_full_static <- function(gd, ncol = 2) {
  other_tcga_ids <- {
    tcga_metadata <- fread(gd$paths$tcga_metadata, data.table = FALSE)
    tcga_metadata$rail_id[
      tcga_metadata$gdc_cases.project.project_id != gd$cancer_type &
      tcga_metadata$rail_id %in% rownames(gd$dat)
    ]
  }
  other_gtex_ids <- {
    gtex_metadata <- fread(gd$paths$gtex_metadata)
    gtex_metadata$rail_id[
      gtex_metadata$study != gd$gtex_study &
      gtex_metadata$rail_id %in% rownames(gd$gtex_dat)
    ]
  }

  tcga_other_psi <- fread(file.path(gd$paths$psi_dir, "tcga", gd$chr,
                                    gd$gene_name, "psi_table.tsv"),
                          data.table = FALSE) %>%
    filter(V1 %in% other_tcga_ids) %>% column_to_rownames("V1")
  gtex_other_psi <- fread(file.path(gd$paths$psi_dir, "gtex", gd$chr,
                                    gd$gene_name, "psi_table.tsv"),
                          data.table = FALSE) %>%
    filter(V1 %in% other_gtex_ids) %>% column_to_rownames("V1")

  plots <- list(); k <- 1
  for (junc in gd$junctions_to_compare) {
    gtex_col <- which(grepl(junc, colnames(gd$gtex_norm_dat)))[1]
    tcga_col <- which(grepl(junc, colnames(gd$tcga_norm_dat)))[1]
    other_tcga_col <- which(grepl(junc, colnames(tcga_other_psi)))[1]
    other_gtex_col <- which(grepl(junc, colnames(gtex_other_psi)))[1]
    if (is.na(gtex_col) || is.na(tcga_col)) next

    to_plot <- data.frame(
      junction = c(
        gd$gtex_norm_dat[, gtex_col],
        gd$tcga_norm_dat[, tcga_col],
        if (!is.na(other_tcga_col)) tcga_other_psi[, other_tcga_col] else numeric(0),
        if (!is.na(other_gtex_col)) gtex_other_psi[, other_gtex_col] else numeric(0)
      ),
      cohort_id = c(
        rep(gd$gtex_study,  nrow(gd$gtex_norm_dat)),
        rep(gd$cancer_type, nrow(gd$tcga_norm_dat)),
        if (!is.na(other_tcga_col)) rep("Other TCGA", nrow(tcga_other_psi)) else character(0),
        if (!is.na(other_gtex_col)) rep("Other GTEx", nrow(gtex_other_psi)) else character(0)
      )
    )

    pval <- gd$junction_pvals[[junc]]
    coef <- gd$junction_coefs[[junc]]
    show_pval <- !is.na(pval) && !is.na(coef) && coef > 0
    pval_label <- if (show_pval) sprintf("p = %.3g", pval) else ""

    p <- ggplot(to_plot, aes(x = junction, color = cohort_id)) +
      stat_ecdf(geom = "step", linewidth = 1.2) +
      scale_color_manual(values = setNames(
        c("#1b9e77", "#d95f02", "#7570b3", "#e7298a"),
        c(gd$gtex_study, gd$cancer_type, "Other TCGA", "Other GTEx")
      )) +
      labs(title = paste0(gd$gene_name, ":", junc),
           x = "Leafcutter PSI", y = "Empirical CDF") +
      annotate("text", x = 1, y = 0, label = pval_label,
               hjust = 1, vjust = -0.3, size = 4) +
      coord_fixed(ratio = 1) +
      theme_minimal(base_size = 14) +
      theme(legend.position = "right",
            legend.title    = element_blank(),
            panel.grid      = element_blank(),
            axis.line       = element_line(color = "black", linewidth = 0.8),
            plot.title      = element_text(face = "bold", hjust = 0.5, size = 15),
            plot.margin     = margin(10, 20, 10, 10))

    plots[[k]] <- p; k <- k + 1
  }

  if (length(plots) == 0) return(ggplot() + theme_void() +
                                 labs(title = "no DE junctions"))
  wrap_plots(plots, ncol = ncol)
}

# ---- Factor-loading violin plot (static PDF) ------------------------------
# Highlighted TCGA cohorts vs an "Other Cohorts" pool. Default highlight is
# the sample's own cancer_type. Highlighted colors come from the d3
# category20 palette (ggsci::pal_d3("category20")); Other is gray80.

plot_loading_violins <- function(gd, cancer_types = NULL) {
  if (is.null(cancer_types) || length(cancer_types) == 0) {
    cancer_types <- gd$cancer_type
  }

  tcga_metadata <- fread(gd$paths$tcga_metadata, data.table = FALSE)

  fact_col <- gd$fact_idx
  loadings <- data.frame(
    rail_id = as.numeric(rownames(gd$res$L)),
    loading = gd$res$L[, fact_col]
  ) %>%
    left_join(tcga_metadata %>% dplyr::select(rail_id, gdc_cases.project.project_id),
              by = "rail_id") %>%
    mutate(group = ifelse(gdc_cases.project.project_id %in% cancer_types,
                          gdc_cases.project.project_id, "Other\nCohorts"))

  group_order <- c(sort(cancer_types), "Other\nCohorts")
  loadings$group <- factor(loadings$group, levels = group_order)

  n_high <- length(cancer_types)
  d3_20  <- ggsci::pal_d3("category20")(20)
  group_colors <- c(d3_20[seq_len(n_high)], "gray80")
  names(group_colors) <- group_order

  n_labels <- loadings %>%
    group_by(group) %>%
    summarize(N = n(), .groups = "drop")
  # Place the "N=" annotations a fixed FRACTION below zero rather than
  # a hard -0.08: genes whose loadings span a much larger range (e.g.
  # NF1, max ~1.0 vs EZH2 ~0.4) otherwise get the labels crammed onto
  # the axis line.
  n_label_y <- -0.15 * max(loadings$loading, na.rm = TRUE)

  comparisons <- lapply(cancer_types, function(ct) c("Other\nCohorts", ct))

  theme_violin <- theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.4),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 10, margin = margin(r = 4)),
      axis.text = element_text(size = 9, color = "black"),
      axis.text.x = element_text(angle = 30, hjust = 1),
      plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
      legend.position = "none",
      plot.margin = margin(4, 4, 4, 4)
    )

  ggplot(loadings, aes(x = group, y = loading, fill = group)) +
    geom_violin(width = 0.8, color = "black", linewidth = 0.3,
                scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.12, fill = "white", color = "black",
                 outlier.size = 0.4, linewidth = 0.3) +
    ggpubr::stat_compare_means(
      method       = "wilcox.test",
      label        = "p.format",
      comparisons  = comparisons,
      hide.ns      = FALSE,
      size         = 3.4,
      tip.length   = 0.01,
      # Extra vertical padding between the p-value text and the
      # comparison bracket line below it.
      vjust        = -1.5
    ) +
    geom_text(data = n_labels,
              aes(x = group, y = n_label_y, label = paste0("N=", N)),
              size = 2.8, color = "grey40") +
    scale_fill_manual(values = group_colors) +
    theme_violin +
    labs(y = "Loading", title = paste0(gd$gene_name, " ", gd$fact))
}

# ---- MANE gene body plot (static PDF) -------------------------------------
# Mirrors the NKAP MANE gene model block in pygenome_from_GSPs.Rmd: a thin
# ggbio::geom_alignment of the MANE transcript exons, optionally with a
# translucent red rect highlighting `region` (a "chr:start-end" string).

plot_mane_gene_model_static <- function(gd, region = NULL) {
  p <- ggplot() +
    ggbio::geom_alignment(gd$t_models, gap.geom = getOption("figure1.gap_geom", "arrow"), label = FALSE,
                          exon.rect.h = 0.05) +
    theme_minimal() +
    ggtitle(paste0(gd$gene_name, " (MANE transcript)"))

  if (!is.null(region) && nzchar(region)) {
    coords <- gsub(",", "", region)
    region_start <- suppressWarnings(as.numeric(sub(".*:(\\d+)-.*", "\\1", coords)))
    region_end   <- suppressWarnings(as.numeric(sub(".*-(\\d+)$",  "\\1", coords)))
    if (!is.na(region_start) && !is.na(region_end)) {
      p <- p + annotate("rect", xmin = region_start, xmax = region_end,
                        ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "red")
    }
  }
  p
}

# ---- pyGenomeTracks input prep --------------------------------------------
# Writes junctions.bed, sample_bigwigs.tsv, groups_file.txt, region.txt
# into out_dir for one region.

prep_pygenome_inputs <- function(gd, out_dir, region = NULL) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  all_norm_dat <- bind_rows(
    gd$tcga_norm_dat %>% as.data.frame() %>% rownames_to_column("sample_id"),
    gd$gtex_norm_dat %>% as.data.frame() %>% rownames_to_column("sample_id")
  ) %>% column_to_rownames("sample_id")
  all_norm_dat[is.na(all_norm_dat)] <- 0

  df_vals <- as.data.frame(t(all_norm_dat * 100))
  df_vals$junc_id <- rownames(df_vals)
  df_meta <- df_vals %>%
    separate(junc_id, into = c("clu", "range"), sep = ":", remove = FALSE) %>%
    separate(range,   into = c("start", "end"), sep = "-", convert = TRUE)
  df_meta$start    <- df_meta$start - 1
  df_meta$`#Chrom` <- gd$chr
  df_meta$strand   <- gd$gene_strand
  df_meta$junc <- paste0(df_meta$`#Chrom`, ":", df_meta$start, ":", df_meta$end,
                         ":clu_", df_meta$clu, "_", df_meta$strand)
  df_meta$gid  <- paste0(df_meta$`#Chrom`, "_clu_", df_meta$clu, "_", df_meta$strand)

  final_bed <- df_meta %>%
    relocate(`#Chrom`, start, end, junc, gid, strand) %>%
    select(-clu, -junc_id, -junc) %>%
    group_by(start, end) %>% slice_head(n = 1) %>% ungroup()

  num_cols <- sapply(final_bed, is.numeric)
  final_bed[, num_cols] <- lapply(final_bed[, num_cols, drop = FALSE],
                                  function(x) round(as.double(x), 4))

  write.table(final_bed, file.path(out_dir, "junctions.bed"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  tcga_samples <- rownames(gd$tcga_norm_dat)
  gtex_samples <- rownames(gd$gtex_norm_dat)

  sample_df <- bind_rows(
    tibble(
      sample = tcga_samples,
      bigwigPath = file.path(gd$paths$tcga_bw_dir,
                             str_split_i(gd$cancer_type, '-', 2),
                             paste0(tcga_samples, ".bw")),
      Group_label = gd$cancer_type, strand = "."
    ),
    tibble(
      sample = gtex_samples,
      bigwigPath = file.path(gd$paths$gtex_bw_dir, gd$gtex_study,
                             paste0(gtex_samples, ".bw")),
      Group_label = paste0("GTEx-", gd$gtex_study), strand = "."
    )
  )
  # AggregateBigwigsForPlotting.py crashes with UnboundLocalError when one of
  # the per-group bigwigs is missing (the loop's last iteration leaves `bw`
  # unset). Drop any rows whose bigwigPath doesn't exist on disk.
  exists_mask <- file.exists(sample_df$bigwigPath)
  if (any(!exists_mask)) {
    cat(sprintf("[prep_pygenome_inputs] dropping %d sample(s) with missing bigwigs\n",
                sum(!exists_mask)))
    sample_df <- sample_df[exists_mask, , drop = FALSE]
  }
  write.table(sample_df, file.path(out_dir, "sample_bigwigs.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  group_df <- data.frame(
    Group_label = c(gd$cancer_type, paste0("GTEx-", gd$gtex_study)),
    color       = c("#969696", "#7fc97f"),
    Bedgz = "", Supergroup = "", stringsAsFactors = FALSE
  )
  write.table(group_df, file.path(out_dir, "groups_file.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  if (is.null(region) || !nzchar(region)) {
    junc_starts <- suppressWarnings(as.numeric(str_split_i(gd$junctions_to_compare, '-', 1)))
    junc_ends   <- suppressWarnings(as.numeric(str_split_i(gd$junctions_to_compare, '-', 2)))
    gene_starts <- start(gd$t_models@unlistData@ranges)
    gene_ends   <- end(gd$t_models@unlistData@ranges)
    coords_lo <- c(junc_starts, gene_starts)
    coords_hi <- c(junc_ends,   gene_ends)
    coords_lo <- coords_lo[is.finite(coords_lo)]
    coords_hi <- coords_hi[is.finite(coords_hi)]
    if (!length(coords_lo) || !length(coords_hi)) {
      stop(sprintf(
        "Cannot auto-compute region for %s: no finite junction or gene coords (n_juncs=%d, n_gene_ranges=%d)",
        gd$gene_name, length(gd$junctions_to_compare), length(gene_starts)))
    }
    region_start <- min(coords_lo) - 100
    region_end   <- max(coords_hi) + 100
    region <- paste0(gd$chr, ":",
                     format(region_start, big.mark = ",", scientific = FALSE), "-",
                     format(region_end,   big.mark = ",", scientific = FALSE))
  }

  writeLines(region, file.path(out_dir, "region.txt"))
  invisible(list(out_dir = out_dir, region = region))
}
