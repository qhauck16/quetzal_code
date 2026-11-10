
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(bit64)
  library(yaml)
  library(ggbio)
  library(patchwork)
  library(GenomicRanges)
  library(rtracklayer)
  library(fastTopics)
  library(ggsci)
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

round_down_to_thousand <- function(value) floor(value   / 1000) * 1000
round_up_to_thousand   <- function(value) ceiling(value / 1000) * 1000

# ---- Palette helpers ------------------------------------------------------
# scale_color_d3 / pal_d3() default is the "category10" palette (10 colors).
# Use "category20" to allow more factors. Both come from ggsci.

factor_palette <- function(num_factors) {
  if (num_factors <= 10) {
    ggsci::pal_d3("category10")(num_factors)
  } else if (num_factors <= 20) {
    ggsci::pal_d3("category20")(num_factors)
  } else {
    # Fall back: tile category20 to cover everything (>20 factors is unusual).
    rep_len(ggsci::pal_d3("category20")(20), num_factors)
  }
}

# ---- Core data loader (TCGA-only, multi-cohort) ---------------------------

load_gene_data <- function(gene_name, chr, fact, cohorts, paths, thresholds) {
  tcga_metadata    <- fread(paths$tcga_metadata,    data.table = FALSE)
  mane_transcripts <- fread(paths$mane_transcripts, data.table = FALSE)
  gencode          <- readRDS(paths$gencode_granges)

  tcga_dir <- file.path(paths$tcga_chr_dir, chr, gene_name)
  res    <- readRDS(file.path(tcga_dir, "res.RDS"))
  de_res <- readRDS(file.path(tcga_dir, "de_res.RDS"))

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

  fact_idx <- as.integer(sub("^k", "", fact))
  junctions_to_compare <- rownames(res$F)[
    de_res$svalue[, fact_idx]   < thresholds$sval_thresh &
    de_res$postmean[, fact_idx] > thresholds$lfc_thresh &
    res$F[, fact_idx]           > thresholds$abs_thresh
  ]

  # rail_id -> cancer_type map for samples in this gene's dat.
  cohort_lookup <- tcga_metadata %>%
    dplyr::filter(rail_id %in% rownames(dat),
                  gdc_cases.project.project_id %in% cohorts) %>%
    dplyr::select(rail_id, cohort = gdc_cases.project.project_id)

  list(
    gene_name = gene_name, chr = chr, fact = fact, fact_idx = fact_idx,
    cohorts = cohorts, gene_strand = gene_strand,
    dat = dat, res = res, de_res = de_res,
    cohort_lookup = cohort_lookup,
    junctions_to_compare = junctions_to_compare,
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
    dplyr::select(junction, Gene_name, Intron_coord, Strand, Annot, Coding, UTR, GencodePC)
}

# ---- Structure plot (PDF) -------------------------------------------------
# Subset to rail_ids in the specified cohorts, strip "TCGA-" from labels,
# color factors with ggsci::pal_d3.

plot_structure_static <- function(gd) {
  rail_ids <- as.character(gd$cohort_lookup$rail_id)
  if (length(rail_ids) == 0) {
    stop(sprintf("No samples found in cohorts: %s",
                 paste(gd$cohorts, collapse = ", ")))
  }

  # fastTopics::structure_plot accepts a (multinom) fit + a grouping vector,
  # but selects by rownames(L). Subset every per-sample component of the fit
  # (L, Ln, s, ...) so verify.fit() stays consistent.
  res_sub <- gd$res
  keep <- rownames(res_sub$L) %in% rail_ids
  res_sub$L <- res_sub$L[keep, , drop = FALSE]
  if (!is.null(res_sub$Ln)) res_sub$Ln <- res_sub$Ln[keep, , drop = FALSE]
  if (!is.null(res_sub$s)  && length(res_sub$s) == length(keep))
    res_sub$s <- res_sub$s[keep]

  grp <- gd$cohort_lookup$cohort[
    match(rownames(res_sub$L), as.character(gd$cohort_lookup$rail_id))
  ]
  grp_short <- sub("^TCGA-", "", grp)
  grp_short <- factor(grp_short, levels = sub("^TCGA-", "", gd$cohorts))

  num_factors <- ncol(res_sub$F)
  pal <- factor_palette(num_factors)

  # structure_plot's default loadings_order = "embed" orders samples by
  # a t-SNE embedding, and n = 5000 subsamples -- both unseeded. Without
  # this the panel reshuffles on every run (a rerun moved 6034 px at
  # 110 dpi while every other panel was byte-identical), which makes any
  # figure diff useless for spotting real regressions.
  set.seed(20260922)
  fastTopics::structure_plot(
    res_sub,
    grouping = grp_short,
    gap      = 300,
    colors   = pal,
    n        = 5000
  ) +
    theme(aspect.ratio = 1 / 5) +
    labs(title = paste0("TCGA Factorization of ", gd$gene_name,
                        " (", gd$fact, ")"))
}

# ---- Static (PDF) factor-junction arc plot --------------------------------

plot_factor_junctions_static <- function(gd, height_scaling = 0.7,
                                         arc_thickness = 0.15,
                                         exon_rect_h = 0.4,
                                         arc_y_offset = NULL) {
  fact_col    <- gd$fact_idx
  num_factors <- ncol(gd$res$L)
  pal         <- factor_palette(num_factors)

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
        fill = pal[as.numeric(plot_df$factor[i])], alpha = 1
      )
    }
  }

  p + scale_y_continuous(breaks = c(0, 10), labels = c("", "")) +
    labs(title = paste0(gd$gene_name, " ", gd$fact))
}

# ---- MANE gene body plot (static PDF) -------------------------------------

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

# ---- pyGenomeTracks input prep (multi-cohort, no GTEx) --------------------
# Writes junctions.bed, sample_bigwigs.tsv, groups_file.txt, region.txt into
# out_dir for one region.

prep_pygenome_inputs <- function(gd, out_dir, region = NULL) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # PSI tables, restricted to samples in the requested cohorts.
  rail_ids <- as.character(gd$cohort_lookup$rail_id)
  norm_dat <- fread(file.path(gd$paths$tcga_chr_dir, "..",
                              "leafcutter/leaf2_try/new_psi_files/tcga",
                              gd$chr, gd$gene_name, "psi_table.tsv"),
                    data.table = FALSE) %>%
    dplyr::filter(V1 %in% rail_ids) %>%
    tibble::column_to_rownames("V1")
  norm_dat[is.na(norm_dat)] <- 0

  df_vals <- as.data.frame(t(norm_dat * 100))
  df_vals$junc_id <- rownames(df_vals)
  df_meta <- df_vals %>%
    tidyr::separate(junc_id, into = c("clu", "range"),
                    sep = ":", remove = FALSE) %>%
    tidyr::separate(range,   into = c("start", "end"),
                    sep = "-", convert = TRUE)
  df_meta$start    <- df_meta$start - 1
  df_meta$`#Chrom` <- gd$chr
  df_meta$strand   <- gd$gene_strand
  df_meta$junc <- paste0(df_meta$`#Chrom`, ":", df_meta$start, ":", df_meta$end,
                         ":clu_", df_meta$clu, "_", df_meta$strand)
  df_meta$gid  <- paste0(df_meta$`#Chrom`, "_clu_", df_meta$clu, "_", df_meta$strand)

  final_bed <- df_meta %>%
    relocate(`#Chrom`, start, end, junc, gid, strand) %>%
    dplyr::select(-clu, -junc_id, -junc) %>%
    group_by(start, end) %>% slice_head(n = 1) %>% ungroup()

  num_cols <- sapply(final_bed, is.numeric)
  final_bed[, num_cols] <- lapply(final_bed[, num_cols, drop = FALSE],
                                  function(x) round(as.double(x), 4))

  write.table(final_bed, file.path(out_dir, "junctions.bed"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  # One row per (rail_id) with its cohort assignment.
  sample_df <- gd$cohort_lookup %>%
    mutate(
      cohort_short = sub("^TCGA-", "", cohort),
      bigwigPath   = file.path(gd$paths$tcga_bw_dir, cohort_short,
                               paste0(rail_id, ".bw")),
      Group_label  = cohort_short,
      strand       = "."
    ) %>%
    dplyr::rename(sample = rail_id) %>%
    dplyr::select(sample, bigwigPath, Group_label, strand)

  exists_mask <- file.exists(sample_df$bigwigPath)
  if (any(!exists_mask)) {
    cat(sprintf("[prep_pygenome_inputs] dropping %d sample(s) with missing bigwigs\n",
                sum(!exists_mask)))
    sample_df <- sample_df[exists_mask, , drop = FALSE]
  }
  write.table(sample_df, file.path(out_dir, "sample_bigwigs.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  # Color each cohort with a distinct d3 swatch (preserves config order).
  cohort_short_order <- sub("^TCGA-", "", gd$cohorts)
  cohort_colors      <- factor_palette(length(cohort_short_order))
  group_df <- data.frame(
    Group_label = cohort_short_order,
    color       = cohort_colors,
    Bedgz       = "",
    Supergroup  = "",
    stringsAsFactors = FALSE
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
