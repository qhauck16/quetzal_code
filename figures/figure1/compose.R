# Figure 1 assembly. Sources theme.R, regenerates every R-native panel
# from the cached gene_data.RDS + recovery_cache TSVs/RDSes, rasterizes
# the sashimi PDFs (and the panel-a schematic) via pdftocairo + autocrop,
# lays everything out with patchwork, writes out/figure1.pdf.
#
# Edit theme.R (font, palette, sizes, arrowhead style, panel-label size)
# and re-run this script -- no other file needs touching.
#
#   Rscript compose.R
#
# ~30-60 s cold, <10 s on repeat (external PDF rasters cached to raster_cache/).


# ------------------------------------------------------------------
# 0. Setup: paths + theme.
# ------------------------------------------------------------------
if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
HERE <- normalizePath(dirname(sub("--file=", "",
    commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))])),
    mustWork = FALSE)
if (!nzchar(HERE) || !dir.exists(HERE)) HERE <- getwd()
setwd(HERE)

source("../theme.R")
options(figure1.gap_geom = GAP_GEOM)

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(patchwork)
  library(cowplot)
  library(grid)
  library(png)
  library(ggplot2)
  library(ggbio)
  library(yaml)
  library(GenomicRanges)
})
select <- dplyr::select
rename <- dplyr::rename

dir.create("out",           showWarnings = FALSE)
dir.create("raster_cache",  showWarnings = FALSE)

# Custom coverage / gene body / sashimi helpers built from raw bigwigs
# via rtracklayer + hand-rolled ggplot geoms. See scripts/coverage_plot.R.
source("../scripts/coverage_plot.R")

# ------------------------------------------------------------------
# 1. External PDF -> autocropped raster grob helpers.
#    Shared with the other figures so the rasterize/crop behaviour
#    stays identical across the paper.
# ------------------------------------------------------------------
source("../scripts/raster_utils.R")
source("../scripts/pdf_font_subst.R")

# ------------------------------------------------------------------
# Red-star overlay for non-coding junctions on plot_factor_junctions_static.
# Reads junctions_classified.tsv (the DE-junction annotations output by
# make_gene_outputs.R), keeps rows with Coding == FALSE and junction in
# gd$junctions_to_compare, and drops a red asterisk above each arc.
# ------------------------------------------------------------------
# plot_factor_junctions_static() places arcs at y = 5 + fact_idx*4 and
# the gene body at y = 5 -- most of the plot is empty air above the arcs
# and below them. Crop the y-range so both fit tightly.
compress_arcs_ylim <- function(p, gd, top_pad = 3.5, bot = 3.5,
                               arc_y_offset = NULL) {
  base_offset <- if (is.null(arc_y_offset)) as.numeric(gd$fact_idx) * 4
                 else                        arc_y_offset
  y_top <- 5 + base_offset + top_pad
  p + coord_cartesian(ylim = c(bot, y_top), expand = FALSE, clip = "off")
}

# Label each of gd$junctions_to_compare with an integer 1..N above its
# arc, so panel-d/e CDFs can be titled with the matching number.
# Returns list(plot, mapping) where mapping is a data.frame with
# junction id -> integer label.
add_junction_numbers <- function(p, gd, y_pad = 2.6, size = 3,
                                  arc_y_offset = NULL) {
  juncs <- gd$junctions_to_compare
  if (length(juncs) == 0) return(list(plot = p, mapping = NULL))
  parts <- strsplit(juncs, "-")
  x_mid <- vapply(parts,
                  function(pp) (as.numeric(pp[1]) + as.numeric(pp[2])) / 2,
                  numeric(1))
  base_offset <- if (is.null(arc_y_offset)) as.numeric(gd$fact_idx) * 4
                 else                        arc_y_offset
  y_lab <- 5 + base_offset + y_pad
  df <- data.frame(x = x_mid, y = y_lab,
                   label = as.character(seq_along(juncs)))
  # ggrepel auto-nudges labels apart when consecutive junctions cluster
  # within a few hundred bp (as EZH2 1/2/3 do).
  p2 <- p + ggrepel::geom_text_repel(
    data = df, aes(x = x, y = y, label = label),
    inherit.aes    = FALSE,
    size           = size,
    fontface       = "bold",
    color          = "black",
    direction      = "x",
    box.padding    = 0.35,
    point.size     = NA,
    segment.color  = NA,
    max.overlaps   = Inf,
    seed           = 1L
  )
  list(plot = p2, mapping = data.frame(junction = juncs,
                                        number   = seq_along(juncs)))
}

add_noncoding_stars <- function(p, gd, jc_tsv, y_pad = 1.6, star_size = 3.2,
                                 arc_y_offset = NULL) {
  if (!file.exists(jc_tsv)) return(p)
  jc <- data.table::fread(jc_tsv, data.table = FALSE)
  if (!"Coding" %in% names(jc)) return(p)
  nc <- jc[!jc$Coding & jc$junction %in% gd$junctions_to_compare, , drop = FALSE]
  if (nrow(nc) == 0) return(p)
  parts <- strsplit(nc$junction, "-")
  x_star <- vapply(parts,
                   function(pp) (as.numeric(pp[1]) + as.numeric(pp[2])) / 2,
                   numeric(1))
  base_offset <- if (is.null(arc_y_offset)) as.numeric(gd$fact_idx) * 4
                 else                        arc_y_offset
  y_star <- 5 + base_offset + y_pad
  df <- data.frame(x = x_star, y = y_star)
  p + geom_point(data = df, aes(x = x, y = y),
                 inherit.aes = FALSE, shape = 8,
                 color = "red", size = star_size, stroke = 0.8)
}

# ------------------------------------------------------------------
# 2. Panel 1a -- Illustrator schematic.
#    .ai files are PDF-compatible; pdftocairo reads them.
# ------------------------------------------------------------------
# preserve_aspect = FALSE fills the cell exactly so the schematic is
# flush with panel b's top; the row-1 cell aspect (~1.18) is within ~1%
# of the image's (~1.19), so the stretch is imperceptible. Re-check if
# the row-1 geometry changes.
# wrap_elements() gives panel a the same tag treatment as panel b.
# pdf_subst_font() re-sets the schematic's Arial text in FONT_FAMILY
# (see ../scripts/pdf_font_subst.R).
panel_a <- wrap_elements(full =
  gg_from_pdf(pdf_subst_font("panel_a/gene_level_workflow.pdf"),
              dpi = 600,
              preserve_aspect = FALSE,
              # Drop the light-gray backing boxes (fill level ~0.906
              # + its antialiased edge ~0.95). The exon rectangles
              # sit at ~0.796 and are left untouched.
              whiten_gray = c(0.88, 0.96)) +
    theme(plot.margin = margin(0, 0, 0, 0, "pt")))

# ------------------------------------------------------------------
# 3. Panel 1b -- recovery curves at FC = 1.5. Rebuild ggplots from the
#    recovery_cache/ files (fast; no SLURM needed).
# ------------------------------------------------------------------
build_panel_b <- function() {
  cache_dir <- "panel_b/recovery_cache"

  # Helpers copied from panel_b_recovery_fc1.5.R (self-contained here so
  # compose.R stays a single entry point).
  summed_props <- function(fls_obj) {
    res  <- fls_obj$summary
    jids <- str_split_i(colnames(res), "-", 1)
    uniq <- unique(jids)
    out  <- matrix(0, nrow = 2, ncol = length(uniq),
                   dimnames = list(NULL, uniq))
    for (j in uniq) out[, j] <- rowSums(res[, jids == j, drop = FALSE])
    out
  }
  hits_ref <- function(mdf) mdf %>%
    arrange(desc(abs(deltapsi))) %>%
    mutate(junc_id = paste0(start, "-", end)) %>%
    distinct(junc_id, .keep_all = TRUE)

  run_recovery <- function(sp, hr, sk, psi_vals) {
    storage <- matrix(0, length(psi_vals), 2); storage[, 1] <- psi_vals
    for (i in seq_along(psi_vals)) {
      psi <- psi_vals[i]
      cn <- sum(abs(hr$deltapsi) > psi) + 1
      col_ranks <- as.numeric(str_split_i(colnames(sp), ":", 1))
      in_cut <- col_ranks < cn
      cand <- sp[, in_cut, drop = FALSE]
      cs <- sign(hr$deltapsi[col_ranks[in_cut]])
      succ <- sum((cs > 0 & cand[2,] > cand[1,]) | (cs < 0 & cand[1,] > cand[2,]))
      sr <- as.numeric(str_split_i(sk$no_factor_match, ":", 1))
      fail <- ncol(cand) - succ + sum(sr < cn)
      storage[i, 2] <- succ / (succ + fail)
    }
    storage
  }

  run_recovery_fc <- function(sp, hr, sk, psi_vals, fc) {
    storage <- matrix(0, length(psi_vals), 3); storage[, 1] <- psi_vals; storage[, 3] <- fc
    for (i in seq_along(psi_vals)) {
      psi <- psi_vals[i]
      cn <- sum(abs(hr$deltapsi) > psi) + 1
      col_ranks <- as.numeric(str_split_i(colnames(sp), ":", 1))
      in_cut <- col_ranks < cn
      cand <- sp[, in_cut, drop = FALSE]
      cs <- sign(hr$deltapsi[col_ranks[in_cut]])
      c1 <- cand[1,]; c2 <- cand[2,]
      hi <- pmax(c1, c2); lo <- pmin(c1, c2)
      ratio <- ifelse(lo > 0, hi/lo, ifelse(hi > 0, Inf, 1))
      above_fc <- ratio >= fc
      dir_ok <- (cs > 0 & c2 > c1) | (cs < 0 & c1 > c2)
      succ <- sum(above_fc & dir_ok)
      sr <- as.numeric(str_split_i(sk$no_factor_match, ":", 1))
      fail <- ncol(cand) - succ + sum(sr < cn)
      storage[i, 2] <- succ / (succ + fail)
    }
    storage
  }

  # Gene-to-chrom map for the matched-background frame construction.
  gencode <- readRDS(file.path(QROOT, "inputs/gencode_v46_granges.RDS"))
  gene_to_chrom <- as.data.frame(mcols(gencode)) %>%
    filter(type == "gene") %>%
    dplyr::select(gene_name) %>%
    mutate(chrom = as.character(seqnames(gencode[gencode$type == "gene"]))) %>%
    distinct(gene_name, .keep_all = TRUE)

  comps <- list(
    KICH_vs_LUAD = "KICH vs LUAD",
    KICH_vs_TGCT = "KICH vs TGCT",
    TGCT_vs_LUAD = "TGCT vs LUAD"
  )

  plots <- list()
  for (comp in names(comps)) {
    fls_test    <- readRDS(file.path(cache_dir, paste0(comp, "_fls_test.RDS")))
    fls_matched <- readRDS(file.path(cache_dir, paste0(comp, "_fls_matched.RDS")))
    m_test <- fread(file.path(cache_dir, paste0(comp, "_matched_test.tsv")),
                     data.table = FALSE) %>%
      left_join(gene_to_chrom, by = "gene_name") %>% filter(!is.na(chrom))
    m_bg   <- fread(file.path(cache_dir, paste0(comp, "_matched_background.tsv")),
                     data.table = FALSE) %>%
      dplyr::select(-deltapsi) %>%
      left_join(m_test %>% dplyr::select(subclass, deltapsi), by = "subclass") %>%
      left_join(gene_to_chrom, by = "gene_name") %>%
      filter(!is.na(chrom)) %>%
      arrange(desc(abs(deltapsi)))

    tsp <- summed_props(fls_test);   thr <- hits_ref(m_test)
    msp <- summed_props(fls_matched); mhr <- hits_ref(m_bg)
    psi_vals <- c(50, 47.5, 45, 42.5, 40, 39:10) / 100

    # base storage (for x-axis range parity with the paper)
    base_stor <- run_recovery(tsp, thr, fls_test$skipped, psi_vals)
    ts <- run_recovery_fc(tsp, thr, fls_test$skipped, psi_vals, 1.5)
    ms <- run_recovery_fc(msp, mhr, fls_matched$skipped, psi_vals, 1.5)

    df <- data.frame(psi_cutoff = psi_vals,
                     leafcutter = ts[, 2],
                     matched    = ms[, 2]) %>%
      pivot_longer(c(leafcutter, matched), names_to = "set", values_to = "recovery")

    p <- ggplot(df, aes(x = psi_cutoff, y = recovery, color = set)) +
      geom_line(linewidth = 0.7) +
      geom_point(size = 1.3) +
      scale_x_reverse(limits = c(max(psi_vals), min(psi_vals))) +
      # Percent labels to match the "Percent recovery" axis title.
      scale_y_continuous(limits = c(0, 1),
                         labels = scales::percent_format(accuracy = 1)) +
      scale_color_manual(values = c(leafcutter = PAL$leafcutter,
                                     matched    = PAL$matched),
                          labels = c(leafcutter = "LeafCutter DS",
                                     matched    = "Matched control"),
                          name   = NULL) +
      labs(x = expression(paste("|", Delta, "PSI| cutoff")),
           y = "Percent recovery",
           subtitle = comps[[comp]]) +
      theme_publication() +
      theme(legend.position = "none",
            aspect.ratio    = 1,
            axis.text.x     = element_text(angle = 45, hjust = 1, vjust = 1,
                                            size = SIZE_BASE),
            plot.subtitle   = element_text(size = SIZE_STRIP, face = "plain",
                                            family = FONT_FAMILY,
                                            hjust = 0.5,
                                            margin = margin(0, 0, 1, 0, "pt")),
            # Tight margins: subplot content sits flush at the top of
            # its cell, and the two 2x2 rows are pulled closer
            # together (less whitespace between them).
            plot.margin     = margin(0, 4, 0, 4, "pt"))

    plots[[comp]] <- p
  }

  # Legend goes below the row of three (shared).
  legend_p <- ggplot(data.frame(x = 1, y = 1, set = c("leafcutter", "matched")),
                     aes(x, y, color = set)) +
    geom_point(size = 2) +
    scale_color_manual(values = c(leafcutter = PAL$leafcutter,
                                    matched    = PAL$matched),
                        labels = c(leafcutter = "LeafCutter DS",
                                    matched    = "Matched control"),
                        name   = NULL) +
    theme_publication() +
    theme(legend.position = "bottom")
  plots$legend <- cowplot::get_legend(legend_p)

  plots
}

panel_b_plots <- build_panel_b()

# 2x2 grid: top row = KICH_vs_LUAD | KICH_vs_TGCT
#            bottom row = TGCT_vs_LUAD | legend
# The legend is drawn by hand in normalized [0,1] coords rather than as
# a ggplot guide, which sizes itself from the text and would overflow
# this narrow cell.
panel_b_legend_keys <- data.frame(
  y     = c(0.60, 0.40),
  col   = c(PAL$leafcutter, PAL$matched),
  label = c("LeafCutter DS", "Matched control"),
  stringsAsFactors = FALSE)
panel_b_legend_only <- ggplot() +
  geom_point(data = panel_b_legend_keys,
             aes(x = 0.05, y = y), colour = panel_b_legend_keys$col,
             size = 2.6) +
  geom_text(data = panel_b_legend_keys,
            aes(x = 0.13, y = y, label = label),
            hjust = 0, vjust = 0.5, size = 2.5,
            family = FONT_FAMILY, colour = "grey15") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1),
                   expand = FALSE, clip = "off") +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0, "pt"))
panel_b <- wrap_elements(full = cowplot::plot_grid(
  panel_b_plots$KICH_vs_LUAD, panel_b_plots$KICH_vs_TGCT,
  panel_b_plots$TGCT_vs_LUAD, panel_b_legend_only,
  ncol = 2, align = "hv", rel_widths = c(1, 1), rel_heights = c(1, 1)) +
  # Left gutter so the "b" panel tag doesn't collide with the
  # top-left subplot's title.
  theme(plot.margin = margin(0, 0, 0, 16, "pt")))

# ------------------------------------------------------------------
# 4. Panel 1c -- BDH1 structure_plot + junction_arcs + sashimi.
#    Load gene_data cache, source panel_c/lib.R in its own env,
#    re-render ggplots + theme, rasterize sashimi PDF.
# ------------------------------------------------------------------
LC <- new.env()
sys.source("panel_c/scripts/lib.R", envir = LC)
gd_repath <- function(gd) {
  gd$paths$gencode_granges  <- file.path(QROOT, "inputs/gencode_v46_granges.RDS")
  gd$paths$tcga_metadata    <- file.path(QROOT, "inputs/tcga_v2_samples.tsv")
  gd$paths$mane_transcripts <- file.path(QROOT, "inputs/MANE.GRCh38.v1.4.summary.txt")
  gd$paths$tcga_chr_dir     <- file.path(QROOT, "quetzal_runs/gene_level")
  gd$paths$leaf2_class      <- file.path(QROOT, "leafcutter_analyses/leaf2_classification/clustering/Leaf2_junction_classifications.txt")
  gd$paths$tcga_bw_dir      <- file.path(QROOT, "inputs/bigwigs/tcga")
  # panel_d_e-only paths (safe to add for panel_c too; unused there)
  gd$paths$gtex_metadata    <- file.path(QROOT, "inputs/gtexv2_metadata.tsv")     %||% gd$paths$gtex_metadata
  gd$paths$gtex_dat_dir     <- file.path(QROOT, "leafcutter_analyses/junction_tables_gtex") %||% gd$paths$gtex_dat_dir
  gd$paths$psi_dir          <- file.path(QROOT, "leafcutter_analyses/psi_files")          %||% gd$paths$psi_dir
  gd$paths$gtex_bw_dir      <- file.path(QROOT, "inputs/bigwigs/gtex")   %||% gd$paths$gtex_bw_dir
  gd
}
`%||%` <- function(a, b) if (is.null(a)) b else a
bdh1_gd <- gd_repath(readRDS("panel_c/output/BDH1_k4/gene_data.RDS"))

# Inject the shared factor palette into the sourced env before plotting,
# so plot_factor_junctions_static picks up our colors instead of its
# lib-internal `palet`.
LC$palet <- unname(sapply(as.character(seq_len(ncol(bdh1_gd$res$L))), factor_color))

# Structure plot: no gridlines or panel border, cohort labels along the
# bottom, small typography.
panel_c_struct <- LC$plot_structure_static(bdh1_gd) +
  theme_publication() +
  theme(legend.position   = "none",
        panel.grid        = element_blank(),
        panel.border      = element_blank(),
        axis.line.x       = element_line(color = "grey40", linewidth = 0.3),
        axis.ticks.y      = element_line(color = "grey40", linewidth = 0.3),
        axis.text.x       = element_text(angle = 45, hjust = 1, vjust = 1,
                                          size = SIZE_BASE),
        axis.title.x      = element_blank(),
        plot.title        = element_blank(),
        # Asymmetric L/R margin nudges the structure plot left within
        # its cell, tightening the gap to the coverage stack.
        plot.margin       = margin(2, 20, 2, 0, "pt"),
        # No fixed aspect ratio, so the height matches the coverage
        # overlay on the right.
        aspect.ratio      = NULL)

# Panel c coverage: built from raw bigwigs + a hand-rolled gene body.
# Every layer is a native ggplot / patchwork element.
bdh1_region_txt <- readLines(
  "panel_c/output/BDH1_k4/regions/default/region.txt", warn = FALSE)[1]
bdh1_region  <- gsub(",", "", bdh1_region_txt)
bdh1_chrom   <- sub(":.*", "", bdh1_region)
bdh1_x_start <- as.numeric(sub(".*:(\\d+)-.*", "\\1", bdh1_region))
bdh1_x_end   <- as.numeric(sub(".*-(\\d+)",   "\\1", bdh1_region))

COHORT_LABELS <- c("LIHC", "SKCM", "LGG", "COAD", "UVM")
# Each cohort gets its own hue family, chosen from colours the structure
# plot on the left uses little or not at all, so the coverage tracks do
# not read as related to its topics. UVM is the figures' Set2 tan.
COHORT_COLORS <- c(LIHC = "#17BECFFF", SKCM = "#C49C94FF", LGG  = "#7F7F7FFF",
                    COAD = "#98DF8AFF", UVM  = "#E5C494FF")

# The brown from the same tab20 family as COHORT_COLORS above. Used
# for TCGA-LUSC in BOTH places it appears in panel e (violin fill and
# sashimi coverage) and for TCGA-LAML in panel d (violin + CDFs).
LUSC_BROWN <- "#C49C94FF"

# k4-unique junctions -> arcs data frame (position + color per arc)
arcs_df <- do.call(rbind, lapply(bdh1_gd$junctions_to_compare, function(j) {
  parts <- as.numeric(strsplit(j, "-")[[1]])
  data.frame(x_start = parts[1], x_end = parts[2],
              width = bdh1_gd$res$F[j, bdh1_gd$fact_idx],
              color = factor_color(bdh1_gd$fact_idx))
}))

# Vertical stack of 5 coverage tracks (normalized to a shared y_max)
# with the k4 arcs drawn immediately above the top (LIHC) track.
bw_paths <- file.path(
  "panel_c/output/BDH1_k4/regions/default",
  paste0(COHORT_LABELS, "-3-..bw"))

# Each track is normalized to its OWN max so every cohort reaches the
# same height in the panel (peaks are directly readable as
# within-cohort intensity, not cross-cohort absolute coverage).
coverage_track_ggs <- lapply(seq_along(bw_paths), function(i) {
  coverage_track_gg(bw_paths[i], bdh1_chrom, bdh1_x_start, bdh1_x_end,
                    track_color = unname(COHORT_COLORS[COHORT_LABELS[i]]),
                    track_label = COHORT_LABELS[i],
                    y_max = NULL,  # per-track auto-scaled
                    # SIZE_BASE is in pt; geom_text/annotate size is in
                    # mm, so divide by .pt to land on the same rendered
                    # size as the structure plot's axis labels.
                    label_size = SIZE_BASE * 1.5 / .pt,
                    label_family = FONT_FAMILY,
                    y_axis = TRUE, y_axis_pt = SIZE_BASE * 1.2)
})

# Arcs as curved LINES (thin lines, not filled ovals) — each drawn as
# a geom_path along a half-sine curve; linewidth encodes junction
# support.
arc_row_gg <- ggplot() +
  scale_x_continuous(limits = c(bdh1_x_start, bdh1_x_end), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 2), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(plot.margin = margin(0, 44, 0, 0, "pt"))
max_w <- max(arcs_df$width, na.rm = TRUE)
# Arcs in the lower half of the panel (apex at y=0.85), so the
# caption at y=1.7 sits clearly above their tops.
for (j in seq_len(nrow(arcs_df))) {
  x1 <- arcs_df$x_start[j]; x2 <- arcs_df$x_end[j]
  xs <- seq(x1, x2, length.out = 200)
  arc_y <- sin((xs - x1) / (x2 - x1) * pi) * 0.85
  lw <- 0.6 + 0.9 * (arcs_df$width[j] / max_w)
  arc_row_gg <- arc_row_gg + geom_path(
    data = data.frame(x = xs, y = arc_y),
    aes(x = x, y = y), inherit.aes = FALSE,
    color = as.character(arcs_df$color[j]), linewidth = lw)
}
arc_row_gg <- arc_row_gg + annotate(
  "text",
  x = (min(arcs_df$x_start) + max(arcs_df$x_end)) / 2,
  y = 1.85, hjust = 0.5, vjust = 1,
  # Only the GSP name is italic.
  label = 'italic("k4")*"-unique junctions"', parse = TRUE,
  size = 3.5,
  family = FONT_FAMILY, color = "grey15")

panel_c_gene_body <- gene_body_gg(
  bdh1_gd$t_models, bdh1_chrom, bdh1_x_start, bdh1_x_end,
  gene_label = "BDH1",
  exon_h     = 0.35,
  n_arrows   = 20,
  label_size = 3.4,
  label_family = FONT_FAMILY)

all_rows    <- c(list(arc_row_gg), coverage_track_ggs, list(panel_c_gene_body))
# Arc row taller so the "k4-unique junctions" caption is clearly
# readable above the arcs.
all_heights <- c(0.9, rep(1, length(coverage_track_ggs)), 0.9)
panel_c_overlay_gg <- Reduce(`/`, all_rows) + plot_layout(heights = all_heights)

# Rasterize the whole panel-c overlay stack to a raster; cowplot
# treats it as a single opaque image (no reshaping of the inner
# grid).
panel_c_overlay <- {
  tf <- tempfile(fileext = ".pdf")
  # cairo_pdf with a white background, not pdf(): pdf() substitutes its
  # own Helvetica and rasterises its transparent page with a grey backing.
  grDevices::cairo_pdf(tf, width = 8, height = 3.5, bg = "white")
  print(panel_c_overlay_gg); dev.off()
  pf <- tempfile()
  system2("pdftocairo", c("-png", "-r", 300, "-singlefile", tf, pf))
  img_c <- autocrop_png(paste0(pf, ".png"))
  ncol_img <- ncol(img_c); nrow_img <- nrow(img_c)
  ggplot() +
    annotation_raster(img_c, 0, ncol_img, 0, nrow_img) +
    coord_cartesian(xlim = c(0, ncol_img), ylim = c(0, nrow_img),
                     expand = FALSE, clip = "off") +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0, "pt"))
}

# Side-by-side: structure_plot LEFT | coverage-with-arc+gene overlay RIGHT.
# cowplot::plot_grid rather than patchwork's `|`, whose padding around
# wrap_elements clips the bottom coverage rows.
# Both halves occupy the same vertical band (PANEL_C_H of the row) so
# the coverage stack and the structure plot read as equal height.
PANEL_C_H <- 0.82
# Raised off the row's bottom edge so the 45-deg cohort labels'
# descenders are not clipped by panel d's row.
PANEL_C_STRUCT_LIFT <- 0.04
panel_c_struct_shrunk <- cowplot::ggdraw() +
  cowplot::draw_plot(panel_c_struct, x = 0, y = PANEL_C_STRUCT_LIFT,
                     width = 1, height = PANEL_C_H - PANEL_C_STRUCT_LIFT)
# The coverage raster's content starts lower inside its region than the
# structure plot's; lift it so the two halves' content tops line up.
PANEL_C_COV_LIFT <- 0.13
panel_c_overlay <- cowplot::ggdraw() +
  cowplot::draw_plot(panel_c_overlay, x = 0, y = PANEL_C_COV_LIFT,
                     width = 1, height = PANEL_C_H)
# wrap_elements() to match panels a/b/d/e: patchwork styles tags on a
# bare cowplot gtable differently from a wrapped element.
panel_c <- wrap_elements(full = cowplot::plot_grid(
  panel_c_struct_shrunk, panel_c_overlay,
  ncol = 2, rel_widths = c(0.95, 1.45)))

# ------------------------------------------------------------------
# 5. Panel 1d -- EZH2 loading violin + gene arcs + 5 CDFs + sashimi.
# ------------------------------------------------------------------
LDE <- new.env()
sys.source("panel_d_e/scripts/lib.R", envir = LDE)
# same palette injection for the arc colors
LDE$palet <- unname(sapply(as.character(seq_len(30)), factor_color))

ezh2_gd <- gd_repath(readRDS("panel_d_e/output/EZH2_LAML/EZH2/gene_data.RDS"))

# Violin: TCGA-LAML uses LUSC_BROWN, the same as the CDFs and as TCGA-LUSC in panel e,
# and expand the y-axis upward so ggpubr's p-value bracket doesn't clip.
ezh2_violin_raw <- LDE$plot_loading_violins(ezh2_gd)
ezh2_loading_max <- max(ezh2_gd$res$L[, ezh2_gd$fact_idx])
ezh2_violin <- ezh2_violin_raw +
  scale_fill_manual(values = setNames(
    c(LUSC_BROWN, PAL$tcga_other),
    c("TCGA-LAML", "Other\nCohorts"))) +
  # Keep the visible axis range close to the data; clip="off" lets ggpubr's stat_compare_means draw its bracket + p-value ABOVE
  # the top axis without expanding the y-scale.
  coord_cartesian(ylim = c(-0.24 * ezh2_loading_max,
                            ezh2_loading_max * 1.20), clip = "off") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  labs(y = bquote(italic(.(ezh2_gd$fact)) * "-loading")) +
  theme_publication() +
  theme(legend.position = "none",
        plot.title      = element_blank(),
        axis.text.x     = element_text(angle = 0, hjust = 0.5,
                                        size = SIZE_BASE))

# EZH2 arcs over the gene body, same style as panel c, with red stars
# for non-coding junctions and numeric labels 1..N.
# Gene body sits at y=1 (per gene_body_gg row layout for a single tx);
# arcs are drawn as curved geom_paths above the gene body with
# linewidth scaled by their Poisson-factor weight, and number labels
# 1..N sit above each arc.
ezh2_chrom   <- ezh2_gd$chr
ezh2_x_start <- ezh2_gd$min_coord
ezh2_x_end   <- ezh2_gd$max_coord

ezh2_arcs_df <- do.call(rbind, lapply(seq_along(ezh2_gd$junctions_to_compare),
                                       function(i) {
  j <- ezh2_gd$junctions_to_compare[i]
  parts <- as.numeric(strsplit(j, "-")[[1]])
  data.frame(x_start = parts[1], x_end = parts[2],
             width   = ezh2_gd$res$F[j, ezh2_gd$fact_idx],
             # tab20 medium blue rather than the factor palette's
             # near-navy -- reads as blue instead of just dark.
             color   = "#1F77B4FF",
             number  = i)
}))
ezh2_max_w <- max(ezh2_arcs_df$width, na.rm = TRUE)

# Build gene body then extend its y-range upward so arcs fit above.
ezh2_gene_body <- gene_body_gg(ezh2_gd$t_models, ezh2_chrom,
                                ezh2_x_start, ezh2_x_end,
                                gene_label = "EZH2",
                                exon_h = 0.35, n_arrows = 20,
                                label_size = 3.4,
                                label_family = FONT_FAMILY) +
  # Top of the range must clear the junction numbers (y = 1.75) AND
  # the non-coding asterisks above them (y = 2.05).
  scale_y_continuous(limits = c(0.3, 2.55), expand = c(0, 0))

# Arc apex sits at y = 1.2 (just above the gene body row) + arc_amp;
# small arc_amp so arcs hug the gene body.
arc_amp <- 0.4
for (k in seq_len(nrow(ezh2_arcs_df))) {
  x1 <- ezh2_arcs_df$x_start[k]; x2 <- ezh2_arcs_df$x_end[k]
  xs <- seq(x1, x2, length.out = 200)
  # Arc starts at y=1.15 (just above gene body top at y ~ 1.175 with
  # exon_h/2 = 0.175) and rises by arc_amp.
  arc_y <- 1.15 + sin((xs - x1) / (x2 - x1) * pi) * arc_amp
  # linewidth = width / max_w, scaled to a max of ~0.9 pt. Baseline
  # 0.15 pt so the smallest arc is still visible; the max is set by
  # max_lw. Guaranteed proportional to the NMF factor F-weight
  # (ezh2_gd$res$F[j, ezh2_gd$fact_idx]).
  min_lw <- 0.15; max_lw <- 0.9
  lw <- min_lw + (max_lw - min_lw) * (ezh2_arcs_df$width[k] / ezh2_max_w)
  ezh2_gene_body <- ezh2_gene_body +
    geom_path(data = data.frame(x = xs, y = arc_y),
              aes(x = x, y = y), inherit.aes = FALSE,
              color = as.character(ezh2_arcs_df$color[k]), linewidth = lw)
}
# Number labels above each arc.
ezh2_num_df <- data.frame(
  x = (ezh2_arcs_df$x_start + ezh2_arcs_df$x_end) / 2,
  y = 1.15 + arc_amp + 0.40,
  label = as.character(ezh2_arcs_df$number))
# Spread labels of adjacent short arcs (junctions 2 and 3) to a minimum
# spacing of 2.2% of the gene's span, pushing crowded pairs apart
# symmetrically; the star data frame below reuses these x's.
local({
  sep <- 0.022 * (ezh2_x_end - ezh2_x_start)
  o <- order(ezh2_num_df$x); x <- ezh2_num_df$x[o]
  for (it in 1:20) for (k in seq_along(x)[-1]) if (x[k] - x[k - 1] < sep) {
    d <- (sep - (x[k] - x[k - 1])) / 2; x[k - 1] <- x[k - 1] - d; x[k] <- x[k] + d
  }
  ezh2_num_df$x[o] <<- x
})
# Non-coding junctions get a red asterisk above the number.
ezh2_jclass <- read.table(
  "panel_d_e/output/EZH2_LAML/EZH2/junctions_classified.tsv",
  sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# Coordinate mismatch: junctions_classified starts are 0-based
# (148815954), gd junctions are 1-based (148815955). Match by adding 1.
ezh2_jclass$junction_1b <- paste0(
  as.integer(sub("-.*", "", ezh2_jclass$junction)) + 1,
  "-", sub(".*-", "", ezh2_jclass$junction))
noncoding <- setNames(!ezh2_jclass$Coding, ezh2_jclass$junction_1b)
ezh2_arcs_df$noncoding <- unname(noncoding[
  paste0(ezh2_arcs_df$x_start, "-", ezh2_arcs_df$x_end)])
ezh2_arcs_df$noncoding[is.na(ezh2_arcs_df$noncoding)] <- FALSE
ezh2_star_df <- data.frame(
  x = ezh2_num_df$x,
  y = 1.15 + arc_amp + 0.75,
  keep = ezh2_arcs_df$noncoding)
ezh2_star_df <- ezh2_star_df[ezh2_star_df$keep, , drop = FALSE]

ezh2_arcs <- ezh2_gene_body +
  geom_text(data = ezh2_num_df, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, size = 2.8,
            family = FONT_FAMILY, fontface = "bold", color = "grey15") +
  {
    if (nrow(ezh2_star_df) > 0)
      geom_text(data = ezh2_star_df, aes(x = x, y = y),
                label = "*", inherit.aes = FALSE, size = 5,
                family = FONT_FAMILY, fontface = "bold", color = "red")
    else NULL
  } +
  theme(plot.margin = margin(2, 30, 10, 6, "pt"))

# Build the 5 CDFs directly from gd data, showing only TCGA-LAML vs
# GTEx-BLOOD (no "Other TCGA" / "Other GTEx" pools).
# Title each panel with its junction NUMBER (1..N) to match the labels
# painted on ezh2_arcs.
build_two_cohort_cdfs <- function(gd, title_prefix = "Junction ") {
  plots <- list()
  n_juncs <- length(gd$junctions_to_compare)
  for (i in seq_len(n_juncs)) {
    junc <- gd$junctions_to_compare[i]
    gtex_col <- which(grepl(junc, colnames(gd$gtex_norm_dat)))[1]
    tcga_col <- which(grepl(junc, colnames(gd$tcga_norm_dat)))[1]
    if (is.na(gtex_col) || is.na(tcga_col)) next

    df <- data.frame(
      psi       = c(gd$gtex_norm_dat[, gtex_col],
                    gd$tcga_norm_dat[, tcga_col]),
      cohort_id = c(rep(gd$gtex_study,  nrow(gd$gtex_norm_dat)),
                    rep(gd$cancer_type, nrow(gd$tcga_norm_dat)))
    )
    pval <- gd$junction_pvals[[junc]]
    # Beta-binomial Wald P (panel_d_e/scripts/lib.R). Where it underflows to
    # 0, print it from the exact log10 P stored alongside.
    l10  <- if (!is.null(gd$junction_log10p)) gd$junction_log10p[[junc]] else NA
    pval_lab <- if (is.na(pval)) "" else if (pval > 0) sprintf("p = %.3g", pval) else {
      e <- floor(l10); sprintf("p = %.2ge%d", 10^(l10 - e), e) }

    p <- ggplot(df, aes(x = psi, color = cohort_id)) +
      stat_ecdf(geom = "step", linewidth = 0.9) +
      scale_color_manual(values = setNames(
        c(PAL$gtex, LUSC_BROWN),
        c(gd$gtex_study, gd$cancer_type)),
        labels = setNames(c(paste0("GTEx-", gd$gtex_study), gd$cancer_type),
                          c(gd$gtex_study, gd$cancer_type)),
        name = NULL) +
      scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
      scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
      labs(title = paste0(title_prefix, i),
           x = "Junction PSI", y = "Empirical CDF") +
      # Middle-right: every curve is above 0.9 by PSI ~ 0.45, so that
      # region is empty. x must stay <= 1 or the scale limits drop it.
      # Junction 1 also holds the colour key, so its p-value sits higher.
      annotate("text", x = 1, y = if (i == 1) 0.62 else 0.45, label = pval_lab,
               hjust = 1, vjust = 0.5, size = 2.0, family = FONT_FAMILY) +
      theme_publication() +
      theme(legend.position = "none",
            plot.title      = element_text(size = SIZE_STRIP, hjust = 0.5,
                                           face = "bold",
                                           family = FONT_FAMILY),
            plot.margin     = margin(2, 6, 2, 4, "pt"),
            aspect.ratio    = 1)
    if (i == 1)
      p <- p + guides(colour = guide_legend(ncol = 1, reverse = TRUE)) +
        theme(legend.position          = "inside",
              legend.position.inside   = c(1, 0.04),
              legend.justification     = c(1, 0),
              legend.background        = element_blank(),
              legend.key               = element_blank(),
              legend.text              = element_text(size = 5.5, family = FONT_FAMILY),
              legend.key.height        = unit(6, "pt"),
              legend.key.width         = unit(8, "pt"),
              legend.spacing.y         = unit(0, "pt"),
              legend.margin            = margin(0, 0, 0, 0))
    plots[[length(plots) + 1L]] <- p
  }
  list(plots = plots,
       # get_plot_component(..., "guide-box-bottom"), not get_legend(),
       # which returns the first (empty) of several guide boxes.
       shared_legend = cowplot::get_plot_component(
         plots[[1]] + theme(legend.position = "bottom",
                            legend.title = element_blank()),
         "guide-box-bottom", return_all = TRUE))
}

ezh2_cdf_pack <- build_two_cohort_cdfs(ezh2_gd)
ezh2_cdfs     <- wrap_plots(ezh2_cdf_pack$plots, nrow = 1)
cdf_legend    <- ezh2_cdf_pack$shared_legend

# Layout:
#   [ gene body + arcs               ] full width
#   [ violin | CDFs (row of 5) ] widths 1 : 4  (colour key inside Junction 1)
# The colour key lives inside the Junction 1 CDF (see build_two_cohort_cdfs).
panel_d <- wrap_elements(full =
  ezh2_arcs /
    ((ezh2_violin | ezh2_cdfs) + plot_layout(widths = c(1, 4))) +
    plot_layout(heights = c(0.5, 1)))

# ------------------------------------------------------------------
# 6. Panel 1e -- NF1 loading violin + arcs + narrow sashimi + tx model.
# ------------------------------------------------------------------
nf1_gd <- gd_repath(readRDS("panel_d_e/output/NF1_LUSC/NF1/gene_data.RDS"))

nf1_violin_raw <- LDE$plot_loading_violins(nf1_gd)
nf1_loading_max <- max(nf1_gd$res$L[, nf1_gd$fact_idx])
# Lift the comparison bracket clear of the top outlier (ggpubr puts it
# just above the data maximum).
lift_bracket <- function(p, y) {
  for (i in seq_along(p$layers))
    if (inherits(p$layers[[i]]$stat, "StatSignif"))
      p$layers[[i]]$stat_params$y_position <- y
  p
}
nf1_violin_raw <- lift_bracket(nf1_violin_raw, nf1_loading_max * 1.10)
nf1_violin <- nf1_violin_raw +
  scale_fill_manual(values = setNames(
    c(LUSC_BROWN, PAL$tcga_other),
    c("TCGA-LUSC", "Other\nCohorts"))) +
  coord_cartesian(ylim = c(-0.24 * nf1_loading_max,
                            nf1_loading_max * 1.24), clip = "off") +
  labs(y = bquote(italic(.(nf1_gd$fact)) * "-loading")) +
  theme_publication() +
  theme(legend.position = "none",
        plot.title      = element_blank(),
        axis.text.x     = element_text(angle = 0, hjust = 0.5,
                                        size = SIZE_BASE))

NF1_ARC_OFFSET <- 0.8
nf1_arcs_base    <- LDE$plot_factor_junctions_static(
  nf1_gd, height_scaling = 0.8, arc_thickness = 0.35,
  exon_rect_h = 0.7, arc_y_offset = NF1_ARC_OFFSET)
nf1_arcs_starred <- add_noncoding_stars(
  nf1_arcs_base, nf1_gd,
  "panel_d_e/output/NF1_LUSC/NF1/junctions_classified.tsv",
  y_pad = 1.1, star_size = 3.5, arc_y_offset = NF1_ARC_OFFSET)
nf1_arcs <- compress_arcs_ylim(nf1_arcs_starred, nf1_gd,
                                top_pad = 2.4, bot = 4.5,
                                arc_y_offset = NF1_ARC_OFFSET) +
  theme_track() +
  theme(plot.title    = element_blank(),
        axis.text.x   = element_blank(),
        axis.ticks.x  = element_blank(),
        axis.line.x   = element_blank()) +
  annotate("text", x = nf1_gd$max_coord, y = 5, hjust = 1.05, vjust = 0.5,
           label = "NF1", fontface = "italic",
           family = FONT_FAMILY, color = "grey15",
           size = SIZE_STRIP * 0.36)

# Bespoke sashimi (no pygenometracks): TCGA-LUSC and GTEx-LUNG bigwig
# coverage tracks with per-cohort splice arcs drawn above the coverage,
# thickness ~ sqrt(PSI), staggered apexes, PSI label at each arc apex.
#
# Row heights for panel e's right column, and the resulting DRAWN
# height of one sashimi track. sashimi_track_gg needs the latter to
# convert its PSI label's point size into data units for collision
# clearance, so both are derived from one place here instead of being
# restated as magic numbers in two files.
NF1_ROW_H  <- c(2.0, 2.0, 0.35, 0.55, 0.55)   # tcga, gtex, trap, mane, 216
FIG_H_IN   <- 11
PANEL_E_H  <- 1.15                            # panel e's share of FIG_ROWS
FIG_ROWS   <- 1.95 + 1.085 + 1.30 + 1.15
NF1_TRACK_H_IN <- (NF1_ROW_H[1] / sum(NF1_ROW_H)) *
                  (PANEL_E_H / FIG_ROWS) * FIG_H_IN

nf1_narrow_dir <- "panel_d_e/output/NF1_LUSC/NF1/regions/narrow"
nf1_region_txt <- readLines(file.path(nf1_narrow_dir, "region.txt"),
                             warn = FALSE)[1]
nf1_reg  <- gsub(",", "", nf1_region_txt)
nf1_chrom   <- sub(":.*", "", nf1_reg)
# Outer splice sites of the displayed junctions (PSI >= 5) plus a small pad.
nf1_x_start <- 31248900
nf1_x_end   <- 31258700

# y max = a multiple of each track's own coverage maximum, leaving
# headroom for the arcs.
NF1_YMAX_F <- 1.8   # five arcs leave one exon, so NF1 needs more headroom than
                    # the ~1.35x used in figures 3 and 4.
NF1_YMAX <- NF1_YMAX_F * c(
  tcga = max(bigwig_binned(file.path(nf1_narrow_dir, "TCGA-LUSC-3-..bw"),
                           nf1_chrom, nf1_x_start, nf1_x_end)$coverage),
  gtex = max(bigwig_binned(file.path(nf1_narrow_dir, "GTEx-LUNG-3-..bw"),
                           nf1_chrom, nf1_x_start, nf1_x_end)$coverage))
# The junction group's unique junction(s), drawn red in both tracks.
NF1_UNIQUE <- lapply(strsplit(nf1_gd$junctions_to_compare, "-"), as.numeric)
nf1_sashimi_tcga <- sashimi_track_gg(
  bw_path      = file.path(nf1_narrow_dir, "TCGA-LUSC-3-..bw"),
  chrom = nf1_chrom, start = nf1_x_start, end = nf1_x_end,
  # tab20 brown -- same family as panel c's cohort colors, and
  # matches the TCGA-LUSC violin fill below.
  track_color  = LUSC_BROWN,
  track_label  = "TCGA-LUSC",
  links_path   = file.path(nf1_narrow_dir, "MeanPSI.TCGA-LUSC-3-..links"),
  min_psi      = 1,
  # y axis scaled from this track's own coverage maximum; arcs are drawn
  # over the coverage, anchored at its height, with a guaranteed rise.
  y_max        = NF1_YMAX[["tcga"]],
  # Containment order, outermost highest: 17% (exon 1 -> last exon)
  # contains 31% and 30%, which contain the short cryptic 21%. The
  # helper's default (shortest span highest) would force crossings.
  apex_override = c("21%" = 0.66, "31%" = 0.76, "30%" = 0.76, "17%" = 0.96),
  highlight    = NF1_UNIQUE,
  min_arc_rise = 0.10,
  apex_range   = c(0.94, 0.62),
  label_size   = SIZE_BASE * 1.5 / .pt,
  psi_size     = SIZE_BASE / .pt,
  track_height_in = NF1_TRACK_H_IN,
  label_family = FONT_FAMILY)

nf1_sashimi_gtex <- sashimi_track_gg(
  bw_path      = file.path(nf1_narrow_dir, "GTEx-LUNG-3-..bw"),
  chrom = nf1_chrom, start = nf1_x_start, end = nf1_x_end,
  track_color  = "#8FBC8F",  # muted green
  track_label  = "GTEx-Lung",
  links_path   = file.path(nf1_narrow_dir, "MeanPSI.GTEx-LUNG-3-..links"),
  min_psi      = 1,
  # As for the TCGA track.
  y_max        = NF1_YMAX[["gtex"]],
  apex_override = c("3.9%" = 0.60, "40%" = 0.70, "37%" = 0.72, "16%" = 0.96),
  highlight    = NF1_UNIQUE,
  min_arc_rise = 0.10,
  apex_range   = c(0.94, 0.62),
  label_size   = SIZE_BASE * 1.5 / .pt,
  psi_size     = SIZE_BASE / .pt,
  track_height_in = NF1_TRACK_H_IN,
  label_family = FONT_FAMILY)

# Key for the red arc, in the GTEx track's empty low band (the intron
# right of the middle exon carries almost no coverage and no arc dips
# below ~0.55 of y_max there).
NF1_KEY_X <- nf1_x_start + 0.60 * (nf1_x_end - nf1_x_start)
NF1_KEY_Y <- 0.30 * NF1_YMAX[["gtex"]]
nf1_sashimi_gtex <- nf1_sashimi_gtex +
  annotate("segment", x = NF1_KEY_X, xend = NF1_KEY_X + 0.045 * (nf1_x_end - nf1_x_start),
           y = NF1_KEY_Y, yend = NF1_KEY_Y, colour = "#D62728", linewidth = 0.9) +
  annotate("text", x = NF1_KEY_X + 0.06 * (nf1_x_end - nf1_x_start), y = NF1_KEY_Y,
           hjust = 0, vjust = 0.5, colour = "#D62728", family = FONT_FAMILY,
           size = SIZE_BASE / .pt, parse = TRUE,
           label = sprintf('%d~"unique junction ("*italic("%s")*")"',
                           length(NF1_UNIQUE), nf1_gd$fact))

# -- Own NF1 transcript models (MANE + NF1-216 NMD) with a translucent
# zoom-connector trapezoid tapering from the sashimi width (top) down
# to the highlighted narrow window inside each transcript (bottom).
# ------------------------------------------------------------------
nf1_txmodels_rds <- "panel_d_e/NF1_transcript_models/nf1_txmodels.RDS"
nf1_txmodels <- readRDS(nf1_txmodels_rds)
nf1_tx_lo <- min(sapply(nf1_txmodels, function(x) min(start(x))))
nf1_tx_hi <- max(sapply(nf1_txmodels, function(x) max(end(x))))

nf1_trap <- trapezoid_gg(nf1_tx_lo, nf1_tx_hi,
                          nf1_x_start, nf1_x_end,
                          fill = "#B0C4FF", alpha = 0.45)
nf1_tx_mane <- transcript_zoom_gg(
  nf1_txmodels[["NF1 (MANE)"]], nf1_tx_lo, nf1_tx_hi,
  highlight = c(nf1_x_start, nf1_x_end),
  tx_label  = "NF1 (MANE)",
  label_size = 2.6, label_family = FONT_FAMILY)
nf1_tx_216 <- transcript_zoom_gg(
  nf1_txmodels[["NF1-216 (NMD)"]], nf1_tx_lo, nf1_tx_hi,
  highlight = c(nf1_x_start, nf1_x_end),
  tx_label  = "NF1-216 (NMD)",
  # Left end: this isoform ends inside the highlighted window.
  label_side = "left",
  label_size = 2.6, label_family = FONT_FAMILY)

# Right column, all native ggplots. Patchwork aligns panel edges rather
# than scales, so the sashimi's narrow x-range and the transcripts' wide
# one coexist in one flat stack (flat so plot_layout addresses every row).
nf1_right_col <- (nf1_sashimi_tcga / nf1_sashimi_gtex / nf1_trap /
                   nf1_tx_mane / nf1_tx_216) +
  plot_layout(heights = NF1_ROW_H)
# Sandwich the violin between two blank spacers so its actual drawn
# height matches panel_d's violin (panel_d violin fills only ~62% of
# its row -- arcs + CDF-legend take the rest -- so panel_e's violin
# needs the same fraction here).
nf1_violin_sized <- (plot_spacer() / nf1_violin / plot_spacer()) +
  plot_layout(heights = c(0.19, 0.62, 0.19))
# Widths 1:4 match panel d's violin:CDF-row width ratio.
panel_e <- wrap_elements(full =
  (nf1_violin_sized | nf1_right_col) +
    plot_layout(widths = c(1, 4)))

# ------------------------------------------------------------------
# 7. Compose the full figure.
# ------------------------------------------------------------------
# Row 1: 1a (schematic, wider) | 1b (three recovery curves)
# Row 2: 1c (full width, taller -- coverage tracks are the star)
# Row 3: 1d (full width, arcs + violin + CDFs)
# Row 4: 1e (full width, arcs + violin + sashimi + tx model)
figure1 <- ((panel_a | panel_b) + plot_layout(widths = c(1.63, 1.37))) /
           panel_c /
           panel_d /
           panel_e +
  # Row 1 height + the (1.63, 1.37) widths above are tuned together so
  # panel b's 2x2 of square plots fills its cell.
  plot_layout(heights = c(1.95, 1.085, 1.30, 1.15)) +
  plot_annotation_publication()

# ------------------------------------------------------------------
# 8. Save.
# ------------------------------------------------------------------
out_pdf <- file.path("out", "figure1.pdf")
cat("[compose] writing", out_pdf, "\n")
ggsave_pub(out_pdf, figure1,
           width  = 8.5,   # ~ nature/paper column set
           height = 11,
           limitsize = FALSE)

cat("[compose] done.\n")
