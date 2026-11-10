# SF3B1 top-500 alt-splice-site PSI matrix + SVM phenotype-axis scores
# (30v90, alt 3'ss, k10). Logic: ../../scripts/phenotype_axis.R (ported from
# diff_splicing/sf3b1/top500_alt3_fc_vs_loading.Rmd). Loadings from the
# quetzal_code variant table SF3B1_full_info (vt_path). Writes cache/sf3b1_top500_alt3_fc5_psi.tsv
# and cache/sf3b1_phenotype_axis_scores.tsv. Run from figure2/:
#   Rscript --vanilla scripts/build_cache_sf3b1_phenotype.R [results_root]
source("../scripts/phenotype_axis.R")
source("../scripts/vt_paths.R")
args <- commandArgs(TRUE)
info <- fread(vt_path("SF3B1_full_info"), data.table = FALSE)
pa_build("sf3b1", "30v90", c("alt 3'ss"), info, "k10", "sf3b1_top500_alt3_fc5_psi", "cache",
         results_root = if (length(args)) args[1] else
           file.path(QROOT, "leafcutter_analyses/diff_splicing/results"))
