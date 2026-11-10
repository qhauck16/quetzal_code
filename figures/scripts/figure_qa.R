# Layout QA for composed figures.
#
# Eyeballing a downscaled render misses collisions and misalignment;
# this measures them instead. Rasterises a figure PDF, builds an ink
# mask, and reports:
#
#   * per-cell ink extents, so adjacent panels' CONTENT tops can be
#     compared (bounding boxes align trivially - content often does not)
#   * contiguous all-blank row bands inside a cell, i.e. dead whitespace
#   * whether any ink reaches a cell's edge (clipping risk)
#
# Usage:
#   Rscript figure_qa.R <figure.pdf> <row_heights> <col_widths_per_row>
#
#   row_heights          comma-separated, as passed to plot_layout()
#   col_widths_per_row   semicolon-separated rows, each comma-separated
#                        e.g. "1;0.95,1.55;1.15,1;1,1"
#
# Example (figure 2):
#   Rscript ../scripts/figure_qa.R out/figure2.pdf \
#       1.00,2.33,1.26,1.75 "1;0.95,1.55;1.15,1;1,1"

suppressPackageStartupMessages(library(png))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("usage: figure_qa.R <pdf> <row_heights> <col_widths_per_row>")
pdf_path <- args[1]
row_h    <- as.numeric(strsplit(args[2], ",")[[1]])
col_spec <- lapply(strsplit(args[3], ";")[[1]],
                   function(r) as.numeric(strsplit(r, ",")[[1]]))
stopifnot(length(row_h) == length(col_spec))

DPI <- 200
tmp <- tempfile()
system2("pdftocairo", c("-png", "-r", DPI, "-singlefile", pdf_path, tmp))
img  <- readPNG(paste0(tmp, ".png"))
gray <- 0.2989*img[,,1] + 0.587*img[,,2] + 0.114*img[,,3]
ink  <- gray < 0.9
H <- nrow(ink); W <- ncol(ink)
cat(sprintf("%s  %d x %d px @ %d dpi\n\n", basename(pdf_path), H, W, DPI))

# Panel tags sit in the first ~55 px of a cell; skip them so a tag
# doesn't masquerade as panel content.
TAG_PAD <- 55

row_edges <- round(c(0, cumsum(row_h) / sum(row_h)) * H)
for (r in seq_along(row_h)) {
  r0 <- row_edges[r] + 1; r1 <- row_edges[r + 1]
  cw <- col_spec[[r]]
  col_edges <- round(c(0, cumsum(cw) / sum(cw)) * W)
  cat(sprintf("--- row %d  (px %d-%d) ---\n", r, r0, r1))
  for (k in seq_along(cw)) {
    c0 <- col_edges[k] + 1; c1 <- col_edges[k + 1]
    cell <- ink[r0:r1, min(c0 + TAG_PAD, c1):c1, drop = FALSE]
    rs <- rowSums(cell)
    if (!any(rs > 0)) { cat(sprintf("  cell %d: EMPTY\n", k)); next }
    w  <- which(rs > 0)
    top <- min(w); bot <- max(w)
    # dead horizontal bands inside the cell
    blank <- which(rs == 0)
    blank <- blank[blank > top & blank < bot]
    band <- ""
    if (length(blank)) {
      d <- c(1, diff(blank)); g <- cumsum(d != 1)
      sz <- tapply(blank, g, length)
      if (max(sz) > 15) {
        bb <- blank[g == names(which.max(sz))]
        band <- sprintf("  | largest dead band %d px (rows %d-%d)",
                        max(sz), min(bb), max(bb))
      }
    }
    edge <- if (top <= 2 || bot >= (r1 - r0 - 1)) "  | TOUCHES CELL EDGE" else ""
    cat(sprintf("  cell %d: content rows %d-%d of %d%s%s\n",
                k, top, bot, r1 - r0 + 1, band, edge))
  }
  if (length(cw) > 1) {
    tops <- vapply(seq_along(cw), function(k) {
      c0 <- col_edges[k] + 1; c1 <- col_edges[k + 1]
      cell <- ink[r0:r1, min(c0 + TAG_PAD, c1):c1, drop = FALSE]
      rs <- rowSums(cell); if (any(rs > 0)) min(which(rs > 0)) else NA_integer_
    }, integer(1))
    if (!anyNA(tops) && diff(range(tops)) > 8) {
      cat(sprintf("  ** content tops differ by %d px across this row **\n",
                  diff(range(tops))))
    }
  }
  cat("\n")
}
