
suppressPackageStartupMessages({
  library(data.table)
})

args    <- commandArgs(trailingOnly = TRUE)
out     <- args[length(args)]
inputs  <- args[seq_len(length(args) - 1)]

dat <- rbindlist(lapply(inputs, fread), fill = TRUE)
fwrite(dat, out, sep = "\t")
cat("Wrote", out, "(", nrow(dat), "rows from", length(inputs), "files)\n")
