# Figure 4, panel c. Sourced by ../compose.R, which supplies the shared
# setup (theme tokens, FLASHIER/CACHE paths, ROW_H/COL_W, panel_cached,
# missing_panel). Must leave a `panel_c` object defined.
#
# Rebuild just this panel:   Rscript compose.R c
#
# Intron-feature logistic-regression betas for the three GSP intron sets
# (GSP #10/SF3B1, #13/NKAP-PRP8, #30/YLPM1). The model is transcribed
# from diff_splicing/unified_results/scripts/fit_betas.R ("with_bp").

if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
panel_c <- panel_cached("c", local({

  C_UNIFIED <- file.path(QROOT, "leafcutter_analyses/diff_splicing/unified_results")
  C_SUB     <- file.path("20v60", "dpsi0.01_fc2.5",
                         "junctions_plus_characteristics.tsv")

  # Column key -> directory + header. Order is the reference's.
  c_cols <- tibble::tibble(
    key   = c("sf3b1", "nkap", "ylpm1"),
    # PRPF8, not PRP8 (the yeast name) -- the published header had the typo.
    label = c("SF3B1 introns", "NKAP/PRPF8 introns", "YLPM1 introns"),
    x     = 1:3
  )

  # Row order is FIXED here on purpose. The upstream script ends with
  # arrange(variable), whose result depends on locale collation -- under
  # C collation "Branchpoint_score" sorts before the lower-case names
  # and the published row order flips. Pin it instead.
  c_rows <- tibble::tibble(
    variable = c("upstream_AG_distance", "score_5ss", "score_3ss",
                 "length", "gc_content", "downstream_AG_distance",
                 "Branchpoint_score", "bp_distance_to_3ss"),
    # Prose labels: these are NOT in the pipeline code, they were typed
    # by hand in Illustrator.
    # "log(intron length)" is Quinn's confirmed wording. NB the fitted
    # model is scale(length) on raw bp (compute_features.R:41), so the
    # label and the predictor disagree. Kept as instructed and because
    # the betas here must match the published figure; if the model is
    # ever refitted on log(length) the betas in this column change.
    label = c("distance to 1st upstream AG (nt)",
              "5'ss MaxEnt score",
              "3'ss MaxEnt score",
              "log(intron length)",
              "GC%",
              "distance to 1st downstream AG (nt)",
              "branchpoint score",
              "branchpoint distance 3'ss")
  )
  c_rows$y <- rev(seq_len(nrow(c_rows)))   # first row at the TOP

  # ---- the model, verbatim from scripts/fit_betas.R (with_bp) -------
  fit_one <- function(key) {
    feats <- data.table::fread(file.path(C_UNIFIED, key, C_SUB),
                               data.table = FALSE)
    d1 <- feats %>%
      filter(!is.na(gc_content), !is.na(length),
             !is.na(score_3ss), !is.na(score_5ss),
             !is.na(downstream_AG_distance), !is.na(upstream_AG_distance))
    d2 <- d1 %>%
      filter(!is.na(Branchpoint_score), !is.na(bp_distance_to_3ss),
             bp_distance_to_3ss < 100) %>%
      group_by(subclass) %>% filter(n() == 2) %>% ungroup()
    d2$condition <- factor(d2$condition, levels = c("background", "test"))

    fit <- glm(condition ~ scale(gc_content) + scale(length) +
                 scale(score_3ss) + scale(score_5ss) +
                 scale(downstream_AG_distance) + scale(upstream_AG_distance) +
                 scale(Branchpoint_score) + scale(bp_distance_to_3ss),
               family = binomial, data = d2)

    tab <- as.data.frame(summary(fit)$coefficients) %>%
      tibble::rownames_to_column("variable") %>%
      dplyr::rename(estimate = Estimate, p_value = `Pr(>|z|)`) %>%
      filter(variable != "(Intercept)") %>%
      mutate(variable = gsub("scale\\((.*)\\)", "\\1", variable),
             sig = case_when(p_value < 0.001 ~ "***",
                             p_value < 0.01  ~ "**",
                             p_value < 0.05  ~ "*",
                             TRUE            ~ "ns"),
             # round(x, 3) prints -0.19 and 0.005 exactly as shown in
             # the reference -- do not format() to a fixed width.
             # NB `cell_label`, not `label`: BOTH c_cols and c_rows
             # carry a `label` column (the column header and the prose
             # row label), so an inner_join on a `label` here is
             # silently renamed to label.x/label.y and c_tiles$label
             # comes back NULL. The stopifnot on row count still
             # passes, so the failure only shows up later as
             # "no non-missing arguments to max".
             cell_label = paste0(round(estimate, 3), " ", sig)) %>%
      select(variable, estimate, p_value, sig, cell_label)
    list(tab = tab, n_test = sum(d2$condition == "test"))
  }

  c_fits <- lapply(c_cols$key, fit_one)
  names(c_fits) <- c_cols$key

  c_tiles <- purrr::imap_dfr(c_fits, function(f, key) {
    f$tab %>% mutate(key = key)
  }) %>%
    inner_join(c_cols, by = "key") %>%
    inner_join(c_rows %>% select(variable, y), by = "variable")

  stopifnot(nrow(c_tiles) == nrow(c_rows) * nrow(c_cols))

  # Cell text (Quinn, 2026-10-08): significant cells show beta to 2
  # decimal places plus stars; non-significant cells are blank (colour only).
  c_tiles <- c_tiles %>%
    mutate(cell_label = ifelse(sig == "ns", "", paste(sprintf("%.2f", estimate), sig)))

  c_n <- c_cols %>%
    mutate(n_test = vapply(c_fits, function(f) f$n_test, numeric(1)),
           label_n = paste0("N = ", n_test))

  cat("[panel c] N test:", paste(c_n$label_n, collapse = " | "), "\n")
  cat("[panel c] betas:",
      paste(c_tiles$cell_label[order(c_tiles$x, -c_tiles$y)], collapse = " | "),
      "\n")

  # ---- geometry -----------------------------------------------------
  # The cell is COL_W$cd[["c"]]/sum(COL_W$cd) of the page, ROW_H[2] tall.
  # The tile block's drawn width is what sizes the cell text: the
  # longest string is "-1.229 ***". Everything else (row-label gutter,
  # legend, margins) is subtracted from the cell width, and the font is
  # then chosen so that string occupies at most C_FILL of a tile.
  C_CELL_W_IN <- PAGE_W * COL_W$cd[["c"]] / sum(COL_W$cd)
  C_MARGIN_IN <- (8 + 2) / 72          # plot.margin left + right, below
  C_LEGEND_IN <- 0.62                  # key + labels + "Beta" + spacing
  C_ROWLAB_PT <- SIZE_BASE * 0.72      # row-label point size
  # Helvetica averages ~0.50 em over mixed lower-case; the longest row
  # label is 34 characters.
  C_ROWLAB_IN <- max(nchar(c_rows$label)) * 0.50 * C_ROWLAB_PT / 72 + 4 / 72
  C_TILES_IN  <- C_CELL_W_IN - C_MARGIN_IN - C_LEGEND_IN - C_ROWLAB_IN
  C_TILE_W_IN <- C_TILES_IN / nrow(c_cols)

  # Digits, spaces, "-" and "*" are all 0.556 em in Helvetica, so the
  # 10-character "-1.229 ***" is 5.56 em wide.
  C_FILL      <- 0.80                  # fraction of a tile the text may use
  C_LONGEST   <- max(nchar(c_tiles$cell_label))
  C_CELL_PT   <- C_TILE_W_IN * 72 * C_FILL / (C_LONGEST * 0.556)
  C_CELL_PT   <- min(C_CELL_PT, SIZE_BASE)   # never larger than base
  cat(sprintf("[panel c] tile %.3f in -> cell text %.2f pt (longest %d chars)\n",
              C_TILE_W_IN, C_CELL_PT, C_LONGEST))

  # Header size is capped the same way against the widest header.
  C_HDR_PT <- min(SIZE_BASE * 0.85,
                  C_TILE_W_IN * 72 * 0.98 /
                    (max(nchar(c_cols$label)) * 0.52))

  # Carry the text size as a COLUMN, not a captured scalar -- an aes()
  # (or even a bare `size =`) resolved from this local() environment is
  # gone after saveRDS/readRDS. See ../../figure2/compose.md.
  # /.pt because scale_size_identity() feeds the value straight to the
  # `size` aesthetic, which ggplot reads in MILLIMETRES -- passing
  # points there renders the text ~2.8x too large.
  c_tiles$txt_pt <- C_CELL_PT / .pt

  # The "N = ..." line is drawn as a layer inside the panel rather than
  # below it: off-panel text silently stretches the panel and drags the
  # axis with it. One row-height of bottom expansion holds it.
  C_N_PT <- SIZE_BASE * 0.78
  c_n$y  <- 0.18
  c_n$n_pt <- C_N_PT / .pt

  p <- ggplot() +
    geom_tile(data = c_tiles,
              aes(x = x, y = y, fill = estimate),
              color = "white", linewidth = 0.35) +
    geom_text(data = c_tiles,
              # `cell_label`, NOT `label`: after the join c_tiles$label
              # is c_cols' COLUMN HEADER, so binding `label` here drew
              # "SF3B1 introns" once per tile.
              aes(x = x, y = y, label = cell_label, size = txt_pt),
              family = FONT_FAMILY, color = "grey10") +
    geom_text(data = c_n,
              aes(x = x, y = y, label = label_n, size = n_pt),
              family = FONT_FAMILY, color = "grey10") +
    scale_size_identity() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, limits = c(-2, 2), name = "Beta",
                         breaks = c(-2, -1, 0, 1, 2)) +
    scale_x_continuous(position = "top", breaks = c_cols$x,
                       labels = gene_md(c_cols$label), expand = c(0, 0)) +
    scale_y_continuous(breaks = c_rows$y, labels = c_rows$label,
                       expand = expansion(add = c(0.82, 0.10))) +
    coord_cartesian(xlim = c(0.5, nrow(c_cols) + 0.5), clip = "off") +
    labs(x = NULL, y = NULL,
         # Two lines: panel c's cell is only 4.55 in wide.
         title = paste0("Introns liable to be mis-spliced show\n",
                        "distinct features across GSPs")) +
    theme_publication() +
    # Centre the title on the WHOLE panel (row labels + table + legend),
    # not just the tile area, which is ggplot's default ("panel").
    theme(plot.title    = element_text(hjust = 0.5),
          plot.title.position = "plot",
          axis.text.x.top = ggtext::element_markdown(size = C_HDR_PT, vjust = 0,
                                         margin = margin(b = 2.5, unit = "pt")),
          axis.text.y   = element_text(size = C_ROWLAB_PT, hjust = 1,
                                       margin = margin(r = 2.5, unit = "pt")),
          axis.ticks    = element_blank(),
          axis.line     = element_blank(),
          legend.position   = "right",
          legend.key.width  = unit(2.2, "mm"),
          legend.key.height = unit(5.0, "mm"),
          legend.text       = element_text(size = SIZE_BASE * 0.75),
          legend.title      = element_text(size = SIZE_BASE * 0.85),
          legend.margin      = margin(0, 0, 0, 0),
          legend.box.spacing = unit(1.2, "mm"),
          # 9.5 pt top, not 12: measured at 400 dpi, panel c's title sat
          # 14 px below panel d's across row 2. Re-measure with
          # ../scripts/title_align.R if either panel's title changes.
          plot.margin   = margin(9.5, 2, 2, 8, "pt"))

  wrap_elements(full = p)
}))
