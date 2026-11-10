# Figure 4f inputs:
#   cache/panel_f_high_psi_hits.tsv
#       YLPM1 20v60 leafcutter, alt 3'ss + alt 5'ss, dPSI > 0.05, FC > 2.5
#   cache/panel_f_tcga_mat.tsv
#       logit(PSI), scaled, over those junctions, for the top 13 GSP #30
#       samples + 200 project-matched controls; NA -> 0, first occurrence
#       per junction, clamp +-5.
# The control draw is seeded (set.seed(30)).
# Logic shared with scripts/phenotype_axis.R (load / classify). Run from
# figure4/:  Rscript --vanilla scripts/build_cache_panel_f_tcga.R [results_root]
source("../scripts/phenotype_axis.R")
source("../scripts/vt_paths.R")
args <- commandArgs(TRUE)
results_root <- if (length(args)) args[1] else file.path(QROOT, "leafcutter_analyses/diff_splicing/results")
TOP_N <- 13; N_CTRL <- 200

A <- pa_annotation()
fls <- pa_load_comparison(file.path(results_root, "ylpm1_20v60"), A)
# Alt 3'ss + alt 5'ss only, consistent with the rest of the paper.
high_psi_hits <- pa_assign_splice_types(fls, A, dpsi_thresh = 0.05, fc_thresh = 2.5) %>%
  filter(type %in% c("alt 3'ss", "alt 5'ss")) %>%
  mutate(end = as.numeric(end), start = as.numeric(start))
dir.create("cache", showWarnings = FALSE)
write_tsv(high_psi_hits, "cache/panel_f_high_psi_hits.tsv")
cat(sprintf("[4f tcga] heatmap junctions (alt 3'ss + 5'ss, dPSI > 0.05, FC > 2.5): %d\n",
            nrow(high_psi_hits)))

info <- fread(vt_path("YLPM1_full_plus_variants"), data.table = FALSE) %>% arrange(desc(k30))
condition <- head(info$rail_id, TOP_N)
project_props <- info %>% slice_head(n = TOP_N) %>% count(project_id, name = "count") %>%
  mutate(n_control = round(count / sum(count) * N_CTRL))
set.seed(30)
control <- info %>% filter(!rail_id %in% condition) %>%
  inner_join(project_props, by = "project_id") %>%
  group_by(project_id) %>%
  group_modify(~ dplyr::slice_sample(.x, n = min(dplyr::first(.x$n_control), nrow(.x)))) %>%
  ungroup() %>% pull(rail_id)
to_use <- c(condition, control)

running_df <- tibble(rail_id = to_use)
suppressWarnings(for (g in unique(high_psi_hits$gene_name)) {
  if (grepl("ENSG", g)) next
  chr <- unique(high_psi_hits$chrom[high_psi_hits$gene_name == g]); chr <- chr[grepl("chr", chr)]
  p <- file.path(file.path(QROOT, "leafcutter_analyses/psi_files"), "tcga", chr, g, "psi_table.tsv")
  if (!file.exists(p) || file.size(p) == 0) next
  tmp <- fread(p, showProgress = FALSE) %>% filter(V1 %in% to_use) %>%
    rename(rail_id = V1) %>% rename_with(~ paste0(g, ":", .x), -rail_id)
  running_df <- running_df %>% left_join(tmp, by = "rail_id")
})
keep_cols <- colnames(running_df)[str_split_i(colnames(running_df), ":", 3) %in%
                                    paste0(high_psi_hits$start, "-", high_psi_hits$end)]
mat <- as.matrix(running_df[, keep_cols])
mat[is.na(mat)] <- 0
m <- scale(qlogis(pmin(pmax(mat, 1e-6), 1 - 1e-6)))
m <- m[, colSums(is.na(m)) == 0]
collapsed <- sub("^([^:]+):[^:]+:([0-9]+-[0-9]+)$", "\\1:\\2", colnames(m))
keep <- !duplicated(collapsed)
m <- m[, keep]; colnames(m) <- collapsed[keep]
m[m > 5] <- 5; m[m < -5] <- -5
fwrite(as.data.frame(m, check.names = FALSE), "cache/panel_f_tcga_mat.tsv", sep = "\t")
fwrite(tibble(rail_id = running_df$rail_id), "cache/panel_f_tcga_rail_ids.tsv", sep = "\t")
cat(sprintf("[4f tcga] matrix %d samples (%d top + %d controls) x %d junctions\n",
            nrow(m), length(condition), length(control), ncol(m)))
