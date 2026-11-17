# Baseline for the GSP-vs-LeafCutter overlap panels (Fig. S9 and Supplementary
# Note 1 Fig. 3b): the top unannotated junction of each of the 250 gene-level
# junction groups whose weight on the GSP is closest to zero (smallest |F|),
# processed exactly like the driving set (mutation_ds/build_gsp_top.R --null ->
# mutation_ds/gsp_top/fact_<k>_null_junctions.tsv). Same splice types as the
# driving set, and never one of the GSP's top-250 junctions.
#
#   source(file.path(QROOT, "figures/supplementary/cache/null_junction_set.R"))
#   nul <- null_junction_set(fact, types, exclude_df)
null_junction_set <- function(fact, types, exclude_df) {
  f <- file.path(QROOT, "mutation_ds/gsp_top", sprintf("fact_%d_null_junctions.tsv", fact))
  d <- data.table::fread(f, data.table = FALSE)
  d <- d[d$type %in% types, ]
  ex <- paste(exclude_df$chromosome, exclude_df$start, exclude_df$end)
  d <- d[!paste(d$chromosome, d$start, d$end) %in% ex, ]
  unique(data.frame(chr = d$chromosome, strand = d$strand,
                    start = as.integer(d$start), end = as.integer(d$end)))
}
