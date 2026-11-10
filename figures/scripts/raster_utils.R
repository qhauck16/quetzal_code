# Shared raster helpers for paper figures: rasterize an external vector
# PDF (Illustrator schematics, pyGenomeTracks output) and bring it into
# a ggplot/patchwork layout.
#
# Sourced by each figure's compose.R as `source("../scripts/raster_utils.R")`.
# The raster cache is resolved relative to the CALLING figure's working
# directory, so each figure keeps its own raster_cache/.

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
  library(png)
  library(cowplot)
})

# ------------------------------------------------------------------
# 1. External PDF -> autocropped raster grob (used for panels 1a
#    schematic, sashimi tracks, NF1 transcript model).
# ------------------------------------------------------------------
rasterize_pdf <- function(pdf_path, dpi = 300) {
  stopifnot(file.exists(pdf_path))
  # Derive a cache key from the FULL relative path (many sashimi PDFs
  # share the basename `sashimi_simple_coverage.pdf`, so basename alone
  # collides across panels).
  rel <- normalizePath(pdf_path, mustWork = TRUE)
  key <- gsub("[^A-Za-z0-9._-]+", "_", rel)
  png_path <- file.path("raster_cache", paste0(key, ".png"))
  if (!file.exists(png_path) ||
      file.info(png_path)$mtime < file.info(pdf_path)$mtime) {
    cat("[raster] ", pdf_path, " -> ", png_path, "\n", sep = "")
    system2("pdftocairo",
            c("-png", "-r", dpi, "-singlefile", pdf_path,
              tools::file_path_sans_ext(png_path)),
            stdout = TRUE, stderr = TRUE)
  }
  png_path
}

# Trim rows/columns that are ~all near-white so the composed layout
# is dense.
autocrop_png <- function(png_path, tol = 0.02, pad = 4) {
  img <- png::readPNG(png_path)
  gray <- if (length(dim(img)) == 3) 0.2989*img[,,1] + 0.587*img[,,2] + 0.114*img[,,3] else img
  mask <- gray < (1 - tol)
  if (!any(mask)) return(img)
  rr <- range(which(rowSums(mask) > 0))
  cc <- range(which(colSums(mask) > 0))
  rr <- c(max(1, rr[1] - pad), min(nrow(img), rr[2] + pad))
  cc <- c(max(1, cc[1] - pad), min(ncol(img), cc[2] + pad))
  if (length(dim(img)) == 3)
    img[rr[1]:rr[2], cc[1]:cc[2], , drop = FALSE]
  else
    img[rr[1]:rr[2], cc[1]:cc[2], drop = FALSE]
}

grob_from_pdf <- function(pdf_path, dpi = 300, crop = TRUE) {
  png_path <- rasterize_pdf(pdf_path, dpi)
  img <- if (crop) autocrop_png(png_path) else png::readPNG(png_path)
  rasterGrob(img, interpolate = TRUE)
}

# Rasterized-PDF -> ggplot with preserved aspect ratio. Uses
# annotation_raster + coord_equal so patchwork can size the outer cell
# without the image distorting.
gg_from_pdf <- function(pdf_path, dpi = 300, crop = TRUE,
                        drop_top_frac = 0, drop_bottom_frac = 0,
                        drop_left_frac = 0, drop_right_frac = 0,
                        preserve_aspect = TRUE,
                        whiten_gray = NULL) {
  png_path <- rasterize_pdf(pdf_path, dpi)
  img <- png::readPNG(png_path)
  # Knock out a narrow band of NEUTRAL gray levels (set to white).
  # Used to drop the panel-a schematic's light-gray backing boxes
  # without touching its darker gray exon rectangles / arcs. Band is
  # applied before autocrop so the crop reflects what remains.
  if (!is.null(whiten_gray) && length(dim(img)) == 3) {
    r <- img[,,1]; g <- img[,,2]; b <- img[,,3]
    neutral <- abs(r - g) < 0.03 & abs(g - b) < 0.03
    hit <- neutral & r >= whiten_gray[1] & r <= whiten_gray[2]
    img[,,1][hit] <- 1; img[,,2][hit] <- 1; img[,,3][hit] <- 1
  }
  if (crop) {
    tmp <- tempfile(fileext = ".png"); png::writePNG(img, tmp)
    img <- autocrop_png(tmp)
  }
  # Optional slice from the top/bottom/left/right in fractional units.
  # Used to trim pygenome outputs down to just the plot data area
  # (no y-labels, no track names, no auto gene body, no coord axis).
  if (drop_top_frac > 0 || drop_bottom_frac > 0) {
    h  <- nrow(img)
    r1 <- floor(h * drop_top_frac) + 1
    r2 <- h - floor(h * drop_bottom_frac)
    img <- if (length(dim(img)) == 3) img[r1:r2, , , drop = FALSE]
           else                        img[r1:r2, ,  drop = FALSE]
  }
  if (drop_left_frac > 0 || drop_right_frac > 0) {
    w  <- ncol(img)
    c1 <- floor(w * drop_left_frac) + 1
    c2 <- w - floor(w * drop_right_frac)
    img <- if (length(dim(img)) == 3) img[, c1:c2, , drop = FALSE]
           else                        img[, c1:c2,  drop = FALSE]
  }
  ncol_img <- ncol(img); nrow_img <- nrow(img)
  # For panels that need to fill their patchwork cell (e.g. cropped
  # coverage panels), coord_equal preserves aspect ratio; when stacking
  # a raster + a matching-width ggplot we usually want fill (no
  # preservation) so the two align by pixel-x.
  base <- ggplot() +
    annotation_raster(img, 0, ncol_img, 0, nrow_img) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0, "pt"))
  if (preserve_aspect) {
    base + coord_equal(xlim = c(0, ncol_img), ylim = c(0, nrow_img),
                        expand = FALSE, clip = "off")
  } else {
    base + coord_cartesian(xlim = c(0, ncol_img),
                            ylim = c(0, nrow_img),
                            expand = FALSE, clip = "off")
  }
}

# Grob-in-ggplot for legends etc. (aspect ratio doesn't matter there).
gg_wrap_grob <- function(g) {
  ggdraw() + draw_grob(g)
}
