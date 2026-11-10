# Figure 3, panel f. Sourced by ../compose.R. Must leave `panel_f`.
#
# Rebuild just this panel:   Rscript compose.R f
#
# f. Mutations spatially cluster in the C* spliceosome -- full-width row.
# Artwork panel_f/Improved_structure.pdf: three structure renders
# (overview + orange / blue zooms) with a small vector legend in Arial.
#
# Rasterised with GHOSTSCRIPT at 600 dpi (login nodes may lack poppler's
# pdftoppm/pdftocairo), then cropped to the ink bounding box so the art
# fills the cell. The legend's Arial is EMBEDDED, so it cannot be re-set
# in Nimbus Sans; at legend size the two are nearly indistinguishable.
# Title re-set in R like every panel.

F_PDF <- file.path("panel_f", "Improved_structure.pdf")
F_ART <- file.path("raster_cache", "panel_f_improved.png")

if (!file.exists(F_ART) || file.mtime(F_ART) < file.mtime(F_PDF)) {
  dir.create("raster_cache", showWarnings = FALSE)
  cat("[panel f] rasterising", F_PDF, "with ghostscript ...\n")
  tmp <- tempfile(fileext = ".png")
  st <- system2("gs", c("-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=png16m",
                        "-r600", "-dFirstPage=1", "-dLastPage=1",
                        "-dTextAlphaBits=4", "-dGraphicsAlphaBits=4",
                        paste0("-sOutputFile=", tmp), shQuote(F_PDF)))
  stopifnot(st == 0, file.exists(tmp))
  png::writePNG(autocrop_png(tmp, tol = 0.04, pad = 8), F_ART)
}

panel_f <- panel_cached("f", local({
  img <- png::readPNG(F_ART)
  h <- dim(img)[1]; w <- dim(img)[2]
  # Aspect preserved (coord_fixed on the image's own pixel ratio).
  wrap_elements(full =
    ggplot() +
      annotation_raster(img, xmin = 0, xmax = w, ymin = 0, ymax = h,
                        interpolate = TRUE) +
      coord_fixed(ratio = 1, xlim = c(0, w), ylim = c(0, h),
                  expand = FALSE) +
      labs(title = "Mutations spatially cluster in the C* spliceosome", x = NULL, y = NULL) +
      theme_void() +
      theme(plot.title  = element_text(size = SIZE_TITLE, face = "bold",
                                       family = FONT_FAMILY, hjust = 0.5,
                                       margin = margin(0, 0, 4, 0)),
            plot.margin = margin(4, 6, 2, 6)))
}))
