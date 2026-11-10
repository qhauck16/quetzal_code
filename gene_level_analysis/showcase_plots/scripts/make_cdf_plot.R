#!/usr/bin/env Rscript
# DE-junction PSI CDF: TCGA samples in cohort_of_interest vs the matched
# GTEx tissue, over the gene's DE junctions for the chosen factor.
# Only runs for the `tcga_cohort` group and only when both `factor` and
# `cohort_of_interest` are set in config.

script_dir <- dirname(sub("--file=", "",
  commandArgs(trailingOnly = FALSE)[
    grep("--file=", commandArgs(trailingOnly = FALSE))]))
source(file.path(script_dir, "lib_supp.R"))

args        <- commandArgs(trailingOnly = TRUE)
gene        <- args[1]
config_path <- args[2]
out_pdf     <- args[3]

dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)

cfg <- yaml::read_yaml(config_path)
le  <- lookup_entry(cfg, gene)
if (!identical(le$group, "tcga_cohort")) {
  stop("PSI CDF is only defined for tcga_cohort genes; ", gene,
       " is in the ", le$group, " group.")
}
ent     <- le$entry
chr     <- ent$chr
fact    <- if (is.null(ent$factor))             "" else as.character(ent$factor)
coi_raw <- if (is.null(ent$cohort_of_interest)) "" else as.character(ent$cohort_of_interest)
if (!nzchar(fact) || !nzchar(coi_raw)) {
  stop("config[tcga_cohort][", gene, "]: need both `factor` and ",
       "`cohort_of_interest` set before the CDF can be built.")
}

# cohort_of_interest can be a single cohort ("TCGA-OV") or a comma-
# separated set ("TCGA-KICH, TCGA-KIRP"). Split + trim.
coi <- trimws(strsplit(coi_raw, ",", fixed = TRUE)[[1]])
coi <- coi[nzchar(coi)]
if (length(coi) == 0L) {
  stop("Empty cohort_of_interest after parsing: '", coi_raw, "'")
}

# Resolve each cohort to a GTEx tissue. Multi-cohort input must collapse
# to a single GTEx tissue (e.g. TCGA-KICH + TCGA-KIRP -> KIDNEY), so the
# two sides of the CDF are still comparable.
gtex_tissues <- unique(unlist(lapply(coi, function(c) cfg$tcga_to_gtex[[c]])))
gtex_tissues <- gtex_tissues[!is.null(gtex_tissues) & nzchar(gtex_tissues)]
missing <- coi[vapply(coi,
                       function(c) is.null(cfg$tcga_to_gtex[[c]]),
                       logical(1))]
if (length(missing) > 0L) {
  stop("No GTEx tissue mapping for: ",
       paste(missing, collapse = ", "), " in config$tcga_to_gtex.")
}
if (length(gtex_tissues) != 1L) {
  stop("cohort_of_interest cohorts map to multiple GTEx tissues (",
       paste(gtex_tissues, collapse = ", "),
       "); ambiguous control set. Use a single TCGA cohort, or cohorts ",
       "that share one GTEx tissue.")
}
gtex_tissue <- gtex_tissues[[1L]]

coi_label <- paste(sub("^TCGA-", "", coi), collapse = "+")

cat(sprintf("[psi_cdf] %s | factor=%s | TCGA={%s} vs GTEx=%s\n",
            gene, fact, paste(coi, collapse = ","), gtex_tissue))

# ---- DE junctions (coord pairs like "1018502-1049235", 1-based start) ----

res_path    <- file.path(cfg$paths$tcga_chr_dir, chr, gene, "res.RDS")
de_res_path <- file.path(cfg$paths$tcga_chr_dir, chr, gene, "de_res.RDS")
res    <- readRDS(res_path)
de_res <- readRDS(de_res_path)

juncs <- de_junction_names(res, de_res, fact, cfg$thresholds)
cat(sprintf("[psi_cdf] %d DE junctions\n", length(juncs)))
if (length(juncs) == 0L) {
  stop("No DE junctions for ", gene, " ", fact,
       " under current thresholds.")
}

# ---- PSI tables (cols are "<numeric_cluster_id>:<start>-<end>") ----

read_psi_table <- function(p) {
  if (!file.exists(p)) stop("PSI table missing: ", p)
  df <- fread(p, sep = "\t", data.table = FALSE)
  rownames(df) <- df$V1
  df$V1 <- NULL
  df
}

# Match a DE junction (start-end) to a PSI column by suffix.
psi_col_for_junc <- function(junc, cols) {
  suf <- paste0(":", junc)
  hits <- cols[endsWith(cols, suf)]
  if (length(hits) == 0L) NA_character_ else hits[1]
}

tcga_psi_path <- file.path(cfg$paths$tcga_psi_dir, chr, gene, "psi_table.tsv")
gtex_psi_path <- file.path(cfg$paths$gtex_psi_dir, chr, gene, "psi_table.tsv")
tcga_psi <- read_psi_table(tcga_psi_path)
gtex_psi <- read_psi_table(gtex_psi_path)

tcga_md <- read_tcga_meta(cfg$paths$tcga_metadata)
gtex_md <- read_gtex_meta(cfg$paths$gtex_metadata)

# Per-cohort TCGA rail_ids -- each cohort in `coi` becomes its OWN CDF
# line, not a merged TCGA pool. GTEx control stays as a single line.
tcga_rids_per_cohort <- lapply(coi, function(c) {
  tcga_md %>%
    dplyr::filter(project_id == c,
                  rail_id %in% as.integer(rownames(tcga_psi))) %>%
    dplyr::pull(rail_id)
})
names(tcga_rids_per_cohort) <- coi
gtex_rids <- gtex_md %>%
  dplyr::filter(study == gtex_tissue,
                rail_id %in% as.integer(rownames(gtex_psi))) %>%
  dplyr::pull(rail_id)
for (c in coi) {
  cat(sprintf("[psi_cdf] cohort samples: %s=%d\n",
              c, length(tcga_rids_per_cohort[[c]])))
}
cat(sprintf("[psi_cdf] cohort samples: GTEx %s=%d\n",
            gtex_tissue, length(gtex_rids)))

if (all(vapply(tcga_rids_per_cohort, length, integer(1)) == 0L)
    || length(gtex_rids) == 0L) {
  stop("No samples for cohort_of_interest / GTEx tissue in PSI table.")
}

# Long-form PSI for each set, keeping the junction identity so we can
# facet the CDF panel by junction.
psi_long <- function(psi_df, rail_ids, source_label) {
  cols    <- names(psi_df)
  matched <- vapply(juncs, psi_col_for_junc, character(1), cols = cols)
  ok      <- !is.na(matched)
  if (!any(ok)) return(NULL)
  kept_juncs   <- juncs[ok]
  matched_cols <- matched[ok]
  sub          <- psi_df[as.character(rail_ids), matched_cols, drop = FALSE]
  out <- lapply(seq_along(matched_cols), function(i) {
    v <- suppressWarnings(as.numeric(sub[[i]]))
    v <- v[is.finite(v)]
    if (length(v) == 0L) return(NULL)
    data.frame(psi      = v,
               source   = source_label,
               junction = kept_juncs[i],
               stringsAsFactors = FALSE)
  })
  out <- out[!vapply(out, is.null, logical(1))]
  if (length(out) == 0L) return(NULL)
  do.call(rbind, out)
}

tcga_labels <- vapply(coi,
                       function(c) sprintf("TCGA %s", sub("^TCGA-", "", c)),
                       character(1))
gtex_label  <- sprintf("GTEx %s", gtex_tissue)

tcga_longs <- lapply(seq_along(coi), function(i) {
  psi_long(tcga_psi, tcga_rids_per_cohort[[i]], tcga_labels[i])
})
tcga_longs <- tcga_longs[!vapply(tcga_longs, is.null, logical(1))]
gtex_long  <- psi_long(gtex_psi, gtex_rids, gtex_label)

if (length(tcga_longs) == 0L) stop("No PSI values matched in TCGA table.")
if (is.null(gtex_long))        stop("No PSI values matched in GTEx table.")

cdf_df <- do.call(rbind, c(tcga_longs, list(gtex_long)))
cdf_df$source   <- factor(cdf_df$source,
                           levels = c(tcga_labels, gtex_label))
cdf_df$junction <- factor(cdf_df$junction, levels = sort(unique(cdf_df$junction)))

n_junc     <- length(levels(cdf_df$junction))
ncol_facet <- max(1L, min(3L, n_junc))
nrow_facet <- ceiling(n_junc / ncol_facet)

p <- ggplot(cdf_df, aes(x = psi, colour = source)) +
  stat_ecdf(geom = "step", size = 0.7) +
  facet_wrap(~ junction, ncol = ncol_facet) +
  ggsci::scale_colour_d3() +
  labs(title = sprintf("%s %s DE-junction PSI -- TCGA %s vs GTEx %s",
                        gene, fact, coi_label, gtex_tissue),
       subtitle = sprintf("%d junctions, one CDF per panel", n_junc),
       x = "PSI",
       y = "Cumulative fraction",
       colour = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 8))

fig_w <- ncol_facet * 3.5
fig_h <- nrow_facet * 2.7 + 0.8
ggplot2::ggsave(out_pdf, p, width = fig_w, height = fig_h, device = cairo_pdf)
cat("[psi_cdf] wrote ", out_pdf,
    sprintf(" (%d junction panels, %dx%d)\n", n_junc, ncol_facet, nrow_facet))
