.libPaths("/scratch/midway2/qhauck/conda_env/rstudio-server/lib/R/library/")
suppressPackageStartupMessages(library(png))
a <- commandArgs(trailingOnly = TRUE)
f <- a[1]; xs <- as.numeric(strsplit(a[2], ",")[[1]])   # cell x-splits in px
im <- readPNG(f); g <- if (length(dim(im)) == 3) apply(im[,,1:3], c(1,2), mean) else im
ink <- g < 0.55           # BOLD dark title text only; skips grey chrome
W <- ncol(ink)
b <- c(1, xs, W)
for (i in seq_len(length(b) - 1)) {
  sub <- ink[, b[i]:b[i+1], drop = FALSE]
  # the tag is in the leftmost 0.35 in (140 px at 400 dpi); drop it
  sub <- sub[, min(150, ncol(sub)):ncol(sub), drop = FALSE]
  rr <- which(apply(sub, 1, sum) > 120)
  cat(sprintf("cell %d: title ink rows %d-%d\n", i, rr[1], rr[length(rr)]))
}
