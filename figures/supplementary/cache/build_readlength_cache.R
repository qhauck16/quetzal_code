# Cache for supplementary figure S4 (GSP #3 read-length junctions).
# covariate_gsp_figures.Rmd's paths / top5 / meta / build_junction_set /
# per_cohort_cdf / build_junction_set_76 / per_cohort_cdf_76 chunks,
# VERBATIM, run from that Rmd's directory (its paths are relative). Its
# ggsave calls are disabled; only the per-cohort mean-PSI tables are
# written here. Needs the 3.7 GB flashier RDS -> build partition.
.libPaths('/scratch/midway2/qhauck/conda_env/rstudio-server/lib/R/library/')
if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
suppressPackageStartupMessages({ library(tidyverse); library(data.table); library(ggsci) })
setwd(file.path(QROOT, "inputs/productive_unproductive"))
## ---- chunk: paths
FLASH_RDS  <- "softimpute_flash_300_qc_filtered.RDS"
TCGA_META  <- file.path(QROOT, "inputs/tcga_v2_samples.tsv")
CHR_DIR    <- file.path(QROOT, "quetzal_runs/genome_wide")
PSI_DIR    <- file.path(QROOT, "leafcutter_analyses/psi_files/tcga")

GSP          <- 3L      # read-length GSP
N_TOP        <- 5L
W_THRESH     <- 0.01

# CDF plot tuning (chunks near the bottom of the doc).
N_TOP_CDF    <- 200L    # # of (gene, local_factor) groups for the CDF
LOAD_RATIO   <- 2       # within a factor of LOAD_RATIO of the top-loaded
                        # junction in that group, i.e. F[, fi] >= max(F) / 2
MAX_PSI_48   <- 0.001   # drop junctions whose MEAN PSI in 48 bp samples
                        # exceeds this -- focus on junctions that are
                        # essentially unspliced in short-read libraries
MIN_SAMPLES  <- 50L     # per-cohort filter (CDF panels): cohort must have
                        # >= this many samples in BOTH read-length buckets

## ---- chunk: top5
fit   <- readRDS(FLASH_RDS)
F_mat <- fit$F_pm
loading_3 <- F_mat[, GSP]
ord       <- order(loading_3, decreasing = TRUE)
top_rows  <- rownames(F_mat)[head(ord, N_TOP)]

parts <- str_split_fixed(top_rows, "\\.", 2)
top_df <- data.frame(
  rank          = seq_len(N_TOP),
  feature       = top_rows,
  gene          = parts[, 1],
  local_factor  = parts[, 2],
  flash_loading = loading_3[head(ord, N_TOP)],
  stringsAsFactors = FALSE
)

chrs <- paste0("chr", c(1:22, "X", "Y"))
gene_res_path <- function(chr, gene) {
  file.path(CHR_DIR, chr, "FastTopics_output", gene, "res.RDS")
}
find_chr <- function(g) {
  for (cc in chrs)
    if (file.exists(gene_res_path(cc, g))) return(cc)
  NA_character_
}
top_df$chr <- vapply(top_df$gene, find_chr, character(1))
top_df

## ---- chunk: meta
meta <- fread(TCGA_META)
meta <- meta[, .(rail_id = as.integer(rail_id),
                  avg_len = as.numeric(avg_len))]
meta <- meta[!is.na(rail_id) & !is.na(avg_len)]
meta[, rl_group := dplyr::case_when(
  round(avg_len) <= 50L ~ "48 bp",
  round(avg_len) >= 75L ~ "76 bp",
  TRUE                  ~ NA_character_
)]
meta <- meta[!is.na(rl_group)]
table(meta$rl_group)

## ---- chunk: build_junction_set
top_cdf_rows <- rownames(F_mat)[head(ord, N_TOP_CDF)]
parts_cdf    <- str_split_fixed(top_cdf_rows, "\\.", 2)
cdf_genes    <- parts_cdf[, 1]
cdf_lfacts   <- parts_cdf[, 2]
cdf_chrs     <- vapply(cdf_genes, find_chr, character(1))

# meta already has rl_group from the breakdown chunk.
meta_lookup <- meta[, .(rail_id, rl_group)]

junc_long_list <- list()
for (i in seq_len(N_TOP_CDF)) {
  g <- cdf_genes[i]; chr <- cdf_chrs[i]; lf <- cdf_lfacts[i]
  if (is.na(chr)) next
  res <- tryCatch(readRDS(gene_res_path(chr, g)), error = function(e) NULL)
  if (is.null(res) || is.null(res$F)) next
  fi <- suppressWarnings(as.integer(sub("^k", "", lf)))
  if (is.na(fi) || fi < 1L || fi > ncol(res$F)) next

  loads     <- res$F[, fi]
  top_load  <- max(loads)
  if (!is.finite(top_load) || top_load <= 0) next
  cand_juncs <- rownames(res$F)[loads >= top_load / LOAD_RATIO]
  if (!length(cand_juncs)) next

  psi_path <- file.path(PSI_DIR, chr, g, "psi_table.tsv")
  if (!file.exists(psi_path)) next
  psi <- fread(psi_path, sep = "\t", data.table = FALSE)
  rownames(psi) <- psi$V1; psi$V1 <- NULL
  cols <- names(psi)
  rids <- as.integer(rownames(psi))
  rl_lookup_v <- meta_lookup[match(rids, meta_lookup$rail_id), rl_group]

  for (j in cand_juncs) {
    hits <- cols[endsWith(cols, paste0(":", j))]
    if (!length(hits)) next
    v  <- suppressWarnings(as.numeric(psi[[hits[1L]]]))
    ok <- is.finite(v) & !is.na(rl_lookup_v)
    if (!any(ok)) next
    v_ok  <- v[ok]; rl_ok <- rl_lookup_v[ok]; rids_ok <- rids[ok]

    mean_48 <- mean(v_ok[rl_ok == "48 bp"], na.rm = TRUE)
    if (!is.finite(mean_48) || mean_48 > MAX_PSI_48) next

    junc_long_list[[length(junc_long_list) + 1L]] <- data.frame(
      junction      = paste0(g, ":", j),
      gene          = g,
      chr           = chr,
      local_loading = loads[match(j, rownames(res$F))],
      rail_id       = rids_ok,
      psi           = v_ok,
      rl_group      = rl_ok,
      stringsAsFactors = FALSE
    )
  }
}

cdf_df <- do.call(rbind, junc_long_list)
cdf_df$rl_group <- factor(cdf_df$rl_group, levels = c("48 bp", "76 bp"))

n_juncs_kept <- length(unique(cdf_df$junction))
cat(sprintf("Candidate junctions kept (within factor %g of top + mean_PSI_48bp <= %.3g): %d\n",
            LOAD_RATIO, MAX_PSI_48, n_juncs_kept))
cat(sprintf("Total PSI observations: 48bp = %d, 76bp = %d\n",
            sum(cdf_df$rl_group == "48 bp"),
            sum(cdf_df$rl_group == "76 bp")))

## ---- chunk: per_cohort_cdf
# meta_full: rail_id, project_id, avg_len, rl_group (and a TCGA filter).
meta_full <- fread(TCGA_META)
meta_full <- meta_full[, .(
  rail_id    = as.integer(rail_id),
  project_id = gdc_cases.project.project_id,
  avg_len    = as.numeric(avg_len)
)]
meta_full <- meta_full[!is.na(rail_id) & !is.na(avg_len) & !is.na(project_id)]
meta_full[, rl_group := dplyr::case_when(
  round(avg_len) <= 50L ~ "48 bp",
  round(avg_len) >= 75L ~ "76 bp",
  TRUE                  ~ NA_character_
)]
meta_full <- meta_full[!is.na(rl_group)]
meta_full$rl_group <- factor(meta_full$rl_group, levels = c("48 bp", "76 bp"))

cohort_counts <- meta_full[, .(n_48 = sum(rl_group == "48 bp"),
                                n_76 = sum(rl_group == "76 bp")),
                            by = project_id]
qualifying_cohorts <- cohort_counts[n_48 >= MIN_SAMPLES
                                     & n_76 >= MIN_SAMPLES,
                                     project_id]
cat("Cohorts with >= ", MIN_SAMPLES,
    " samples in BOTH read-length buckets:\n", sep = "")
print(cohort_counts[project_id %in% qualifying_cohorts][order(-n_48)])

# Reuse the candidate junction set from the build_junction_set chunk
# (already has rail_id, psi, junction, rl_group).
all_psi <- merge(cdf_df[, c("rail_id", "junction", "psi")],
                  meta_full[, .(rail_id, project_id, rl_group)],
                  by = "rail_id")
all_psi <- all_psi[all_psi$project_id %in% qualifying_cohorts, ]

# Mean PSI per (cohort, junction, bucket). Drop (cohort, junction)
# pairs that don't have data in BOTH buckets so each CDF compares the
# SAME set of junctions per panel.
agg <- as.data.table(all_psi)[
  , .(mean_psi = mean(psi), n = .N),
  by = .(project_id, junction, rl_group)]

both_buckets <- agg[, .(both = uniqueN(rl_group) == 2L),
                     by = .(project_id, junction)]
agg <- merge(agg, both_buckets[both == TRUE,
                                  .(project_id, junction)],
              by = c("project_id", "junction"))

n_juncs_per_cohort <- agg[, .(n_juncs = uniqueN(junction)),
                           by = project_id]
cat("\nJunctions per cohort (after `both buckets present` filter):\n")
print(n_juncs_per_cohort)

p_cohort <- ggplot(agg, aes(x = mean_psi, colour = rl_group)) +
  stat_ecdf(geom = "step", linewidth = 0.8) +
  facet_wrap(~ project_id, ncol = 2) +
  scale_colour_d3() +
  coord_cartesian(xlim = c(0.0, 0.5), ylim = c(0.75, 1.0)) +
  labs(title    = sprintf("GSP #%d candidate junctions: mean PSI per junction, by read length",
                          GSP),
       subtitle = sprintf("within factor %g of top + mean_PSI_48bp <= %.3g; >= %d samples per bucket per cohort",
                          LOAD_RATIO, MAX_PSI_48, MIN_SAMPLES),
       x = "Mean PSI (over cohort samples in bucket)",
       y = "Cumulative fraction (over junctions)",
       colour = "Read length") +
  theme_minimal(base_size = 11) +
  theme(plot.title      = element_text(face = "bold"),
        legend.position = "bottom",
        strip.text      = element_text(face = "bold"),
        aspect.ratio    = 1)

n_cohorts <- length(unique(agg$project_id))
fig_w     <- 8
fig_h     <- ceiling(n_cohorts / 2) * 4 + 1.5
if (FALSE) ggplot2::ggsave("read_length_gsp_per_cohort_cdf.pdf", p_cohort,
                 width = fig_w, height = fig_h, device = cairo_pdf)
p_cohort

## ---- chunk: build_junction_set_76
MAX_PSI_76 <- MAX_PSI_48   # same numeric threshold, applied to the
                            # 76 bp pool instead of the 48 bp pool

junc_long_list_76 <- list()
for (i in seq_len(N_TOP_CDF)) {
  g <- cdf_genes[i]; chr <- cdf_chrs[i]; lf <- cdf_lfacts[i]
  if (is.na(chr)) next
  res <- tryCatch(readRDS(gene_res_path(chr, g)), error = function(e) NULL)
  if (is.null(res) || is.null(res$F)) next
  fi <- suppressWarnings(as.integer(sub("^k", "", lf)))
  if (is.na(fi) || fi < 1L || fi > ncol(res$F)) next

  loads     <- res$F[, fi]
  top_load  <- max(loads)
  if (!is.finite(top_load) || top_load <= 0) next
  cand_juncs <- rownames(res$F)[loads >= top_load / LOAD_RATIO]
  if (!length(cand_juncs)) next

  psi_path <- file.path(PSI_DIR, chr, g, "psi_table.tsv")
  if (!file.exists(psi_path)) next
  psi <- fread(psi_path, sep = "\t", data.table = FALSE)
  rownames(psi) <- psi$V1; psi$V1 <- NULL
  cols <- names(psi)
  rids <- as.integer(rownames(psi))
  rl_lookup_v <- meta_lookup[match(rids, meta_lookup$rail_id), rl_group]

  for (j in cand_juncs) {
    hits <- cols[endsWith(cols, paste0(":", j))]
    if (!length(hits)) next
    v  <- suppressWarnings(as.numeric(psi[[hits[1L]]]))
    ok <- is.finite(v) & !is.na(rl_lookup_v)
    if (!any(ok)) next
    v_ok  <- v[ok]; rl_ok <- rl_lookup_v[ok]; rids_ok <- rids[ok]

    mean_76 <- mean(v_ok[rl_ok == "76 bp"], na.rm = TRUE)
    if (!is.finite(mean_76) || mean_76 > MAX_PSI_76) next

    junc_long_list_76[[length(junc_long_list_76) + 1L]] <- data.frame(
      junction      = paste0(g, ":", j),
      gene          = g,
      chr           = chr,
      local_loading = loads[match(j, rownames(res$F))],
      rail_id       = rids_ok,
      psi           = v_ok,
      rl_group      = rl_ok,
      stringsAsFactors = FALSE
    )
  }
}

cdf_df_76 <- do.call(rbind, junc_long_list_76)
cdf_df_76$rl_group <- factor(cdf_df_76$rl_group, levels = c("48 bp", "76 bp"))

n_juncs_kept_76 <- length(unique(cdf_df_76$junction))
cat(sprintf("Candidate junctions kept (within factor %g of top + mean_PSI_76bp <= %.3g): %d\n",
            LOAD_RATIO, MAX_PSI_76, n_juncs_kept_76))

## ---- chunk: per_cohort_cdf_76
all_psi_76 <- merge(cdf_df_76[, c("rail_id", "junction", "psi")],
                     meta_full[, .(rail_id, project_id, rl_group)],
                     by = "rail_id")
all_psi_76 <- all_psi_76[all_psi_76$project_id %in% qualifying_cohorts, ]

agg_76 <- as.data.table(all_psi_76)[
  is.finite(psi),
  .(mean_psi = mean(psi), n = .N),
  by = .(project_id, junction, rl_group)]

both_buckets_76 <- agg_76[, .(both = uniqueN(rl_group) == 2L),
                           by = .(project_id, junction)]
agg_76 <- merge(agg_76, both_buckets_76[both == TRUE,
                                          .(project_id, junction)],
                 by = c("project_id", "junction"))

n_juncs_per_cohort_76 <- agg_76[, .(n_juncs = uniqueN(junction)),
                                  by = project_id]
cat("\nJunctions per cohort (after `both buckets present` filter):\n")
print(n_juncs_per_cohort_76)

p_cohort_76 <- ggplot(agg_76, aes(x = mean_psi, colour = rl_group)) +
  stat_ecdf(geom = "step", linewidth = 0.8) +
  facet_wrap(~ project_id, ncol = 2) +
  scale_colour_d3() +
  coord_cartesian(xlim = c(0.0, 0.5), ylim = c(0.75, 1.0)) +
  labs(title    = sprintf("GSP #%d candidate junctions: mean PSI per junction, by read length",
                          GSP),
       subtitle = sprintf("within factor %g of top + mean_PSI_76bp <= %.3g; >= %d samples per bucket per cohort",
                          LOAD_RATIO, MAX_PSI_76, MIN_SAMPLES),
       x = "Mean PSI (over cohort samples in bucket)",
       y = "Cumulative fraction (over junctions)",
       colour = "Read length") +
  theme_minimal(base_size = 11) +
  theme(plot.title      = element_text(face = "bold"),
        legend.position = "bottom",
        strip.text      = element_text(face = "bold"),
        aspect.ratio    = 1)

n_cohorts_76 <- length(unique(agg_76$project_id))
fig_w_76     <- 8
fig_h_76     <- ceiling(n_cohorts_76 / 2) * 4 + 1.5
if (FALSE) ggplot2::ggsave("read_length_gsp_per_cohort_cdf_low76.pdf", p_cohort_76,
                 width = fig_w_76, height = fig_h_76, device = cairo_pdf)
p_cohort_76

OUT <- file.path(QROOT, "figures/supplementary/cache")
fwrite(agg,    file.path(OUT, "s4_low48_per_cohort.tsv"), sep = "\t")
fwrite(agg_76, file.path(OUT, "s4_low76_per_cohort.tsv"), sep = "\t")
