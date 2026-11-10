# Supplementary table: every TCGA sample used in the nine GSP differential-
# splicing comparisons (SF3B1 / NKAP / YLPM1 x 10v30, 20v60, 30v90), test and
# control, with the GDC file used. Built from the LeafCutter sample sheets that
# were actually run (results/<gene>_<comp>/leaf_pheno_data.txt), so it records
# what was analysed, not what a selection script would pick today.
#
# Columns: comparison, gsp, role, rail_id, tcga_barcode, project, GSP loading
# (normalised flashier L, as in the figures), matching covariates (mean read
# length, FFPE), GDC file_id / filename / md5, and download_source:
#   "planned"  -- in download/sample_plan.tsv
#   "earlier"  -- downloaded in earlier exploratory analyses and matched from the
#                 on-disk pool (control slots only)
# Any cwd.
if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
suppressPackageStartupMessages({ library(yaml); library(data.table); library(tidyverse) })
DS  <- dirname(file.path(QROOT, "leafcutter_analyses/diff_splicing/unified_results"))                      # .../diff_splicing
FIG <- file.path(QROOT, "figures")
source(file.path(FIG, "scripts/vt_paths.R"))
OUT <- file.path(FIG, "supplementary/out/tables"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

GSP <- c(sf3b1 = 10, nkap = 13, ylpm1 = 30)
load_tab <- list(sf3b1 = "SF3B1_full_plus_variants", nkap = "NKAP_full_plus_variants",
                 ylpm1 = "YLPM1_full_plus_variants")

sheets <- list.dirs(file.path(DS, "results"), recursive = FALSE)
sheets <- sheets[grepl("^(sf3b1|nkap|ylpm1)_[0-9]+v[0-9]+$", basename(sheets))]
used <- rbindlist(lapply(sheets, function(d) {
  p <- fread(file.path(d, "leaf_pheno_data.txt"), header = FALSE)
  data.table(comparison = basename(d), gene = sub("_.*", "", basename(d)),
             role = ifelse(p$V2 == "control", "control", "test"),
             rail_id = as.integer(sub("\\..*", "", p$V1)))
}))
used[, gsp := paste0("GSP #", GSP[gene])]

loads <- rbindlist(lapply(names(load_tab), function(g) {
  x <- fread(vt_path(load_tab[[g]]), data.table = FALSE)
  data.table(gene = g, rail_id = as.integer(x$rail_id), gsp_loading = x[[paste0("k", GSP[g])]])
}))
meta <- fread(file.path(QROOT, "inputs/tcga_v2_samples.tsv"), data.table = FALSE) %>%
  transmute(rail_id = as.integer(rail_id), project = gdc_cases.project.project_id,
            mean_read_length = avg_len,
            ffpe = ifelse(is.na(gdc_cases.samples.is_ffpe), FALSE, gdc_cases.samples.is_ffpe))
# GDC file per sample. used_bams.rail_map.tsv indexes two relabelled BAMs
# under their original rail_id; translate through bam_relabel and use the
# barcode verified on disk (aliquot_on_disk). Re-downloaded BAMs are in
# the rerun rail_map.
relab <- fread(file.path(DS, "download/bam_relabel_2026-09-25.tsv"), data.table = FALSE) %>%
  filter(action == "relabel") %>% transmute(old = as.integer(old_rail_id), new = as.integer(true_rail_id))
rmap <- fread(file.path(DS, "download/used_bams.rail_map.tsv"), data.table = FALSE) %>%
  transmute(rail_id = as.integer(rail_id), tcga_barcode = aliquot_on_disk, gdc_file_id = file_id,
            gdc_filename = filename, md5) %>%
  mutate(rail_id = ifelse(rail_id %in% relab$old, relab$new[match(rail_id, relab$old)], rail_id))
rr <- fread(file.path(DS, "download/rerun_2026-10-05.rail_map.tsv"), data.table = FALSE) %>%
  transmute(rail_id = as.integer(rail_id), tcga_barcode = barcode, gdc_file_id = file_id,
            gdc_filename = filename, md5 = NA_character_)
rmap <- bind_rows(rr, rmap %>% filter(!rail_id %in% rr$rail_id)) %>% distinct(rail_id, .keep_all = TRUE)
plan <- fread(file.path(DS, "download/sample_plan.tsv"), data.table = FALSE)

tab <- used %>% as_tibble() %>%
  left_join(loads, by = c("gene", "rail_id")) %>%
  left_join(meta, by = "rail_id") %>%
  left_join(rmap, by = "rail_id") %>%
  mutate(download_source = ifelse(rail_id %in% plan$rail_id, "planned", "earlier"),
         # ylpm1_6v18: the YLPM1 SVM's junction run (cases = the 6 P samples).
         comparison = factor(comparison, levels = c(paste0(rep(c("sf3b1", "nkap", "ylpm1"), each = 3), "_",
                                                           c("10v30", "20v60", "30v90")), "ylpm1_6v18"))) %>%
  arrange(comparison, desc(role == "test"), desc(gsp_loading)) %>%
  select(comparison, gsp, role, rail_id, tcga_barcode, project, gsp_loading,
         mean_read_length, ffpe, download_source, gdc_file_id, gdc_filename, md5)

stopifnot(!anyNA(tab$comparison), !anyNA(tab$tcga_barcode), !anyNA(tab$gdc_file_id))
fwrite(tab, file.path(OUT, "STable_diff_splicing_samples.tsv"), sep = "\t")
s <- tab %>% count(comparison, role, download_source) %>% pivot_wider(names_from = c(role, download_source),
                                                                      values_from = n, values_fill = 0)
print(as.data.frame(s))
cat(sprintf("[table] %d rows, %d unique samples -> %s\n", nrow(tab), n_distinct(tab$rail_id),
            file.path(OUT, "STable_diff_splicing_samples.tsv")))
