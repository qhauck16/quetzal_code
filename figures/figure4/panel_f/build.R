# ------------------------------------------------------------------
# Figure 4, panel f -- GIST PSI heatmap + row dendrogram.
#
# Drawn natively so it can sit in the patchwork: geom_tile for the
# matrix and the three annotation bars, dendextend::as.ggdend()
# segments for the dendrogram. The clustering follows heatmaply's
# defaults (euclidean / complete / OLO seriation), reproduced with
# stats + dendextend so no plotly is needed at render time.
# ------------------------------------------------------------------

if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
panel_f <- panel_cached("f", local({

  suppressPackageStartupMessages({
    library(dendextend)
    library(grid)
  })

  F_CACHE_COL  <- file.path(CACHE, "panel_f_col_order.rds")
  F_CACHE_JUNC <- file.path(CACHE, "panel_f_tcga_cluster_ids.rds")
  F_CACHE_DATA <- file.path(CACHE, "panel_f_data.rds")
  dir.create(CACHE, showWarnings = FALSE)

  GIST_DIR <- file.path(QROOT, "gist")
  PTN_DIR  <- file.path(GIST_DIR, "paired_tumor_normal")

  F_RAMP <- colorRampPalette(c("navy", "white", "firebrick3"))(50)

  # ----------------------------------------------------------------
  # 1. Column order: the TCGA heatmap's own column clustering.
  #    This is the ONE place heatmaply/plotly is used.  It is
  #    deterministic, so it runs once and is cached; after that the
  #    figure has no plotly dependency at all.
  # ----------------------------------------------------------------
  if (!file.exists(F_CACHE_COL)) {
    message("[panel f] building ", F_CACHE_COL, " (needs heatmaply + plotly)")
    f_tcga_mat <- as.matrix(data.table::fread(
      file.path(CACHE, "panel_f_tcga_mat.tsv"),   # scripts/build_cache_panel_f_tcga.R
      data.table = FALSE, check.names = FALSE))
    f_hm <- heatmaply::heatmaply(f_tcga_mat, scale = "none", colors = F_RAMP)
    saveRDS(plotly::plotly_build(f_hm)$x$layout$xaxis$ticktext, F_CACHE_COL)
  }
  tcga_ordered_juncs <- readRDS(F_CACHE_COL)
  if (length(tcga_ordered_juncs) != 512 || tcga_ordered_juncs[1] != "PTPRO:15595641-15596089")
    message("[panel f] TCGA column order differs from the published one: ",
            length(tcga_ordered_juncs), " junctions, first ", tcga_ordered_juncs[1])

  # ----------------------------------------------------------------
  # 2. Matrix + phenotype bars.
  # ----------------------------------------------------------------
  if (!file.exists(F_CACHE_DATA)) {
    message("[panel f] building ", F_CACHE_DATA, " (reads the 6 MB counts.gz)")

    gist_junc_counts <- data.table::fread(
      file.path(PTN_DIR, "ylpm1_perind_numers.counts.gz"), data.table = FALSE)
    gist_junc_counts$strand <- str_split_i(gist_junc_counts$V1, "_", 3)
    gist_junc_counts <- gist_junc_counts %>%
      separate(V1, into = c("chrom", "start", "end", "cluster_id"), sep = ":") %>%
      mutate(end   = as.numeric(end)   - 1,
             start = as.numeric(start) + 1) %>%
      dplyr::select(-cluster_id)

    high_psi_hits <- readr::read_tsv(file.path(CACHE, "panel_f_high_psi_hits.tsv"),
                                     show_col_types = FALSE)

    # The LeafCutter cluster the PSI denominator is taken over is the
    # TCGA cluster, not the GIST one (GIST's cluster_id is dropped and
    # the cluster is taken from the TCGA psi_table headers). Using GIST's
    # own clusters changes the matrix enough to permute the dendrogram.
    if (!file.exists(F_CACHE_JUNC)) {
      message("[panel f] building ", F_CACHE_JUNC, " (psi_table.tsv headers)")
      f_pairs <- unique(high_psi_hits[, c("chrom", "gene_name")])
      f_hdr <- function(chr, gene) {
        # TCGA PSI tables.
        f <- file.path(file.path(QROOT, "leafcutter_analyses/psi_files"),
                       "tcga", chr, gene, "psi_table.tsv")
        if (!file.exists(f)) return(NULL)
        data.table::data.table(chr = chr, gene_name = gene,
                               column_name = names(data.table::fread(f, nrows = 0)))
      }
      suppressWarnings(juncs_for_ylpm1 <- data.table::rbindlist(
        Map(f_hdr, f_pairs$chrom, f_pairs$gene_name), use.names = TRUE))
      saveRDS(juncs_for_ylpm1, F_CACHE_JUNC)
    }
    tcga_base <- readRDS(F_CACHE_JUNC) %>%
      filter(column_name != "V1") %>%
      separate(column_name, into = c("cluster_id", "junction_coord"), sep = ":") %>%
      separate(junction_coord, into = c("start", "end"), sep = "-") %>%
      mutate(start = as.numeric(start), end = as.numeric(end))

    gist_psi <- tcga_base %>%
      right_join(gist_junc_counts, by = c("chr" = "chrom", "start", "end")) %>%
      group_by(cluster_id) %>%
      mutate(across(starts_with("HRS"), ~ .x / sum(.x), .names = "{.col}_PSI")) %>%
      ungroup() %>%
      dplyr::select(-gene_name)

    to_use <- high_psi_hits %>%
      dplyr::select(gene_name, chrom, start, end, strand) %>%
      distinct() %>%
      left_join(gist_psi, by = c("chrom" = "chr", "start", "end", "strand")) %>%
      group_by(chrom, start, end, strand) %>%
      slice_head(n = 1) %>%
      ungroup() %>%
      mutate(junc_id = paste0(gene_name, ":", start, "-", end),
             ord     = match(junc_id, tcga_ordered_juncs)) %>%
      filter(!is.na(ord)) %>%
      arrange(ord) %>%
      dplyr::select(-ord)

    f_mat <- to_use[, grepl("PSI", colnames(to_use))] %>% as.matrix() %>% t()
    f_all_na <- apply(f_mat, 2, function(x) all(is.na(x)))
    f_mat[, !f_all_na] <- scale(qlogis(pmin(pmax(f_mat[, !f_all_na], 1e-6),
                                            1 - 1e-6)))
    colnames(f_mat) <- paste0(to_use$gene_name, ": ", to_use$start, "-", to_use$end)
    f_mat[is.na(f_mat)] <- 0

    # ---- phenotype bars ------------------------------------------
    gist_mutant_status <- readr::read_tsv(
      file.path(GIST_DIR, "YLPM1_mutant_status.txt"), show_col_types = FALSE)
    gist_experiment_info <- readr::read_tsv(
      file.path(GIST_DIR, "accession_to_experiment.txt"), show_col_types = FALSE) %>%
      filter(`BioSample accession` %in% list.files(PTN_DIR)) %>%
      filter(grepl("Whole-transcriptomic-sequencing", `Experiment title`)) %>%
      left_join(gist_mutant_status, by = c("BioSample name" = "Case")) %>%
      mutate(
        `YLPM1 genomic`  = ifelse(is.na(`YLPM1 genomic`), "NO", `YLPM1 genomic`),
        individual       = str_split_i(`BioSample name`, "[A-Z]", 1),
        tumor_or_normal  = ifelse(grepl("T", `BioSample name`), "Tumor", "Normal"),
        row_key          = paste0(`BioSample accession`, "_PSI"))

    f_pheno <- gist_experiment_info %>%
      dplyr::select(row_key, sample = `BioSample accession`,
                    biosample = `BioSample name`,
                    mut = `YLPM1 genomic`, individual, tn = tumor_or_normal) %>%
      as.data.frame()
    rownames(f_pheno) <- f_pheno$row_key
    f_pheno <- f_pheno[rownames(f_mat), ]

    saveRDS(list(mat = f_mat, pheno = f_pheno), F_CACHE_DATA)
  }
  f_dat   <- readRDS(F_CACHE_DATA)
  f_mat   <- f_dat$mat
  f_pheno <- f_dat$pheno
  if (!identical(dim(f_mat), c(17L, 512L)))
    message("[panel f] matrix is ", paste(dim(f_mat), collapse = " x "), " (published: 17 x 512)")

  # ----------------------------------------------------------------
  # 3. Clustering -- heatmaply's defaults, reproduced natively.
  #    heatmapr() does: dist(x) -> hclust(complete) -> as.dendrogram ->
  #    dendextend::seriate_dendrogram(method = "OLO").  Its displayed
  #    top-to-bottom row order is order.dendrogram() of that dendrogram
  #    (heatmapr$matrix$data comes back REVERSED, because plotly draws
  #    bottom-up).  Deterministic, but seeded anyway.
  # ----------------------------------------------------------------
  set.seed(1)
  f_d    <- stats::dist(f_mat)                      # euclidean
  f_hc   <- stats::hclust(f_d, method = "complete")
  f_dend <- dendextend::seriate_dendrogram(stats::as.dendrogram(f_hc),
                                           f_d, method = "OLO")
  f_ord  <- stats::order.dendrogram(f_dend)         # display order, top -> bottom

  F_EXPECT <- c("HRS989008", "HRS989023", "HRS989059", "HRS988935", "HRS988930",
                "HRS988939", "HRS988918", "HRS988965", "HRS989010", "HRS989065",
                "HRS989085", "HRS989084", "HRS989064", "HRS988968", "HRS988995",
                "HRS989006", "HRS989058")
  f_got <- f_pheno$sample[f_ord]
  if (!identical(f_got, F_EXPECT)) {
    warning("[panel f] row order does NOT match the published dendrogram:\n  ",
            paste(f_got, collapse = " "))
  }
  if (!identical(f_pheno$mut[f_ord], c(rep("NO", 11), rep("YES", 6))))
    warning("[panel f] mutants no longer form one bottom block: ",
            paste(f_pheno$mut[f_ord], collapse = " "))

  N_ROW <- nrow(f_mat); N_COL <- ncol(f_mat)

  # ----------------------------------------------------------------
  # 4. Geometry.  ONE ggplot, one coordinate system: heatmap, strips,
  #    dendrogram, leaders and labels all live in the same x/y units,
  #    which is what guarantees the dendrogram lines up row-for-row.
  #    Row i (1 = top) is drawn at y = N_ROW + 1 - i.
  # ----------------------------------------------------------------
  F_HM_X0   <- 0.5                       # heatmap left  edge
  F_HM_X1   <- N_COL + 0.5               # heatmap right edge
  F_GAP1    <- 5                         # heatmap -> strips
  F_STRIP_W <- 21                        # one annotation bar
  F_GAP2    <- 13                        # strips -> dendrogram
  F_DEND_W  <- 150
  F_STRIP_X0 <- F_HM_X1 + F_GAP1
  F_DEND_X0  <- F_STRIP_X0 + 3 * F_STRIP_W + F_GAP2
  F_XMAX     <- F_DEND_X0 + F_DEND_W

  F_HM_Y0 <- 0.5; F_HM_Y1 <- N_ROW + 0.5
  F_BAND  <- 7.0                          # label + legend band below
  F_YMIN  <- F_HM_Y0 - F_BAND

  # Drawn panel size, measured from the composed figure. Used ONLY to
  # turn point sizes into data units, so text metrics drive the leader
  # geometry. Update if ROW_H[4] or COL_W$fg change.
  F_PANEL_W_IN <- 4.07
  F_PANEL_H_IN <- 2.52
  f_x_per_in <- (F_XMAX - F_HM_X0 + 1) / F_PANEL_W_IN
  f_y_per_in <- (F_HM_Y1 - F_YMIN)     / F_PANEL_H_IN

  # ---- heatmap -----------------------------------------------------
  # as.vector() of a samples x junctions matrix walks DOWN each column,
  # i.e. the sample index varies fastest -- so y must be the `times`
  # sequence and x the `each` one (a transpose would still render a
  # plausible-looking texture).
  f_hm_df <- data.frame(
    x = rep(seq_len(N_COL), each = N_ROW),
    y = rep(N_ROW + 1 - seq_len(N_ROW), times = N_COL),
    z = as.vector(f_mat[f_ord, ]))
  stopifnot(
    identical(f_hm_df$z[f_hm_df$x == 1     & f_hm_df$y == N_ROW],
              unname(f_mat[f_ord[1],     1])),
    identical(f_hm_df$z[f_hm_df$x == N_COL & f_hm_df$y == 1],
              unname(f_mat[f_ord[N_ROW], N_COL])),
    identical(f_hm_df$z[f_hm_df$x == 7     & f_hm_df$y == N_ROW - 4],
              unname(f_mat[f_ord[5],     7])))

  # ---- annotation strips ------------------------------------------
  # heatmaply melts all three bars into ONE factor and spreads a single
  # palette across the union of their 18 levels, so the bars SHARE a
  # palette, hardcoded here.
  F_ANN_COL <- c(
    "104" = "#A6CEE3", "14" = "#1F78B4", "17" = "#B2DF8A", "19" = "#33A02C",
    "35"  = "#FB9A99", "37" = "#E31A1C", "56" = "#FDBF6F", "63" = "#FF7F00",
    "64"  = "#CAB2D6", "65" = "#6A3D9A", "7"  = "#FFFF99", "73" = "#B15928",
    "87"  = "#66C2A5", "91" = "#FC8D62",
    "NO"  = "#8DA0CB", "Normal" = "#E78AC3", "Tumor" = "#A6D854",
    "YES" = "#FFD92F")

  F_BAR_LAB <- c("YLPM1 mutant status", "Individual", "Tumor or Normal")
  f_bar_val <- list(f_pheno$mut[f_ord], f_pheno$individual[f_ord], f_pheno$tn[f_ord])

  f_ann <- do.call(rbind, lapply(seq_along(f_bar_val), function(j) {
    data.frame(bar   = j,
               x     = F_STRIP_X0 + F_STRIP_W * (j - 0.5),
               y     = N_ROW + 1 - seq_len(N_ROW),
               color = unname(F_ANN_COL[as.character(f_bar_val[[j]])]),
               stringsAsFactors = FALSE)
  }))
  stopifnot(!any(is.na(f_ann$color)))

  # ---- dendrogram --------------------------------------------------
  # as.ggdend()'s x is the LEAF POSITION along order.dendrogram(), i.e.
  # exactly our display row index; y is merge height.  Swap the axes and
  # the leaves land dead centre on their heatmap rows by construction.
  f_seg <- dendextend::as.ggdend(f_dend)$segments
  f_hmax <- max(f_seg$y, f_seg$yend)
  f_seg <- data.frame(
    x    = F_DEND_X0 + (f_seg$y    / f_hmax) * F_DEND_W,
    xend = F_DEND_X0 + (f_seg$yend / f_hmax) * F_DEND_W,
    y    = N_ROW + 1 - f_seg$x,
    yend = N_ROW + 1 - f_seg$xend)

  # ---- leader lines to the three bars ------------------------------
  # Label widths come from the font metrics, so the label block's left
  # edge is derived, not placed: the longest label's right edge sits
  # one pad to the left of the first bar.
  f_tmpdev <- tempfile(fileext = ".pdf")
  grDevices::cairo_pdf(f_tmpdev, width = 4, height = 3); on.exit(unlink(f_tmpdev))
  f_lab_w_in <- vapply(F_BAR_LAB, function(s)
    grid::convertWidth(grid::grobWidth(grid::textGrob(
      s, gp = grid::gpar(fontfamily = FONT_FAMILY, fontsize = SIZE_BASE))),
      "in", valueOnly = TRUE), numeric(1))
  grDevices::dev.off()

  # The labels are LEFT-aligned at F_LAB_L, so the longest one has the
  # shortest leader.  Widening the pad alone cannot help it -- it shifts
  # the whole block -- so the block is set back by a minimum horizontal
  # run as well, which every leader then inherits, so the top label's
  # leader is not a bare riser.
  F_PAD_X   <- 0.055 * f_x_per_in                   # text -> leader gap
  F_MIN_RUN <- 0.10  * f_x_per_in                   # shortest leader run
  f_lab_w_x <- f_lab_w_in * f_x_per_in
  F_LAB_L   <- (F_STRIP_X0 + F_STRIP_W * 0.5) - F_PAD_X - F_MIN_RUN -
               max(f_lab_w_x)

  F_LINE_Y  <- (SIZE_BASE / 72) * f_y_per_in        # one text line, in y units
  # First label 0.85 lines below the heatmap to clear its bottom edge.
  f_lab_y   <- F_HM_Y0 - 0.85 * F_LINE_Y - (seq_along(F_BAR_LAB) - 1) * 1.55 * F_LINE_Y

  # ---- colour key for the two categorical bars ---------------------
  # "Individual" needs no key -- its colours only mark which rows come
  # from the same patient. Placed in
  # the empty band under the dendrogram, one entry per line.
  F_KEY <- data.frame(
    # The bar records mutant status (YES/NO in the GIST table), so the key
    # says Mutant / Wild-type rather than claiming loss of function.
    lab = c("Mutant", "Wild-type", "Tumor", "Normal"),
    col = c("#FFD92F", "#8DA0CB", "#A6D854", "#E78AC3"))
  F_KEY$y <- F_HM_Y0 - 0.85 * F_LINE_Y - (seq_len(nrow(F_KEY)) - 1) * 1.30 * F_LINE_Y
  F_KEY$x <- F_DEND_X0 + 0.02 * f_x_per_in
  f_key_sq <- 0.07 * f_x_per_in                     # swatch width, x units

  f_lab_df <- data.frame(bar = seq_along(F_BAR_LAB), label = F_BAR_LAB,
                         x = F_LAB_L, y = f_lab_y,
                         x_right = F_LAB_L + f_lab_w_x + F_PAD_X,
                         x_bar   = F_STRIP_X0 + F_STRIP_W * (seq_along(F_BAR_LAB) - 0.5),
                         stringsAsFactors = FALSE)
  # leader = horizontal run from the label to the bar's x, then a riser
  # up to the bottom edge of the strips.  Both carried as COLUMNS --
  # an aes() closing over a local scalar renders blank after readRDS.
  f_lead_df <- rbind(
    data.frame(x = f_lab_df$x_right, xend = f_lab_df$x_bar,
               y = f_lab_df$y,       yend = f_lab_df$y),
    data.frame(x = f_lab_df$x_bar,   xend = f_lab_df$x_bar,
               y = f_lab_df$y,       yend = F_HM_Y0))

  # ---- borders -----------------------------------------------------
  f_box_df <- data.frame(
    xmin = c(F_HM_X0, F_STRIP_X0),
    xmax = c(F_HM_X1, F_STRIP_X0 + 3 * F_STRIP_W),
    ymin = F_HM_Y0, ymax = F_HM_Y1)

  F_LIM <- range(f_mat)

  p <- ggplot() +
    # geom_tile, not geom_raster: a raster is one bitmap that PDF
    # viewers resample with smoothing (blurry); tiles stay crisp vector
    # rectangles at any zoom.
    # The 1.03 overlap hides the hairline seams viewers anti-alias
    # between exactly-abutting vector rectangles.
    geom_tile(data = f_hm_df, aes(x = x, y = y, fill = z),
              width = 1.03, height = 1.03) +
    geom_tile(data = f_ann, aes(x = x, y = y),
              fill = f_ann$color, width = F_STRIP_W, height = 1) +
    geom_rect(data = f_box_df,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = NA, color = "grey25", linewidth = 0.30) +
    geom_segment(data = f_seg, aes(x = x, y = y, xend = xend, yend = yend),
                 color = "grey15", linewidth = 0.34, lineend = "round") +
    geom_segment(data = f_lead_df, aes(x = x, y = y, xend = xend, yend = yend),
                 color = "grey25", linewidth = 0.30) +
    # Rich text so "YLPM1" in the first label can be italic.
    ggtext::geom_richtext(data = transform(f_lab_df, label = gene_md(label)),
              aes(x = x, y = y, label = label),
              hjust = 0, vjust = 0.5, family = FONT_FAMILY,
              size = SIZE_BASE / .pt, color = "grey10",
              fill = NA, label.colour = NA,
              label.padding = grid::unit(0, "pt"),
              label.r = grid::unit(0, "pt")) +
    geom_tile(data = F_KEY, aes(x = x + f_key_sq / 2, y = y),
              fill = F_KEY$col, width = f_key_sq, height = 0.8 * F_LINE_Y,
              inherit.aes = FALSE) +
    ggtext::geom_richtext(data = transform(F_KEY, lab = gene_md(lab)),
              aes(x = x + f_key_sq * 1.5, y = y, label = lab),
              hjust = 0, vjust = 0.5, family = FONT_FAMILY,
              size = SIZE_BASE * 0.85 / .pt, color = "grey10",
              fill = NA, label.colour = NA,
              label.padding = grid::unit(0, "pt"),
              label.r = grid::unit(0, "pt"), inherit.aes = FALSE) +
    scale_fill_gradientn(colours = F_RAMP, limits = F_LIM,
                         breaks = c(-2, 0, 2), labels = c("-2", "0", "2"),
                         name = "Scaled logit(PSI)") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(F_HM_X0 - 1, F_XMAX + 10), ylim = c(F_YMIN, F_HM_Y1),
                    clip = "off") +
    guides(fill = guide_colourbar(
      title.position = "top", title.hjust = 0,
      frame.colour = "grey25", frame.linewidth = 0.25)) +
    labs(x = NULL, y = NULL,
         title = gene_md(paste0("GIST tumors with YLPM1 LoF recapitulate\n",
                        "differential splicing results identified in TCGA"))) +
    theme_publication() +
    theme(
      plot.title             = md_title(),
      axis.text              = element_blank(),
      axis.ticks             = element_blank(),
      axis.line              = element_blank(),
      panel.grid             = element_blank(),
      legend.position        = "inside",
      legend.position.inside = c(0.005, (F_HM_Y0 - 0.9 * F_LINE_Y - F_YMIN) /
                                         (F_HM_Y1 - F_YMIN)),
      legend.justification   = c(0, 1),
      legend.direction       = "vertical",
      legend.background      = element_blank(),
      legend.title           = element_text(size = SIZE_BASE, hjust = 0,
                                            margin = margin(b = 2, unit = "pt")),
      legend.text            = element_text(size = SIZE_BASE),
      # guide_colourbar multiplies legend.key.height by 5 for a
      # vertical bar (ggplot2 >= 3.5), so this is the bar's drawn
      # height divided by 5 -- NOT a typo.
      legend.key.width       = unit(2.6, "mm"),
      legend.key.height      = unit(12.0 / 5, "mm"),
      legend.margin          = margin(0, 0, 0, 0),
      plot.margin            = margin(4, 4, 2, 4, unit = "pt")
    )

  wrap_elements(full = p)
}))
