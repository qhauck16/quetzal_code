# YLPM1 top-500 alt-splice-site PSI matrix + SVM phenotype-axis scores
# (6v18, alt 3'ss + alt 5'ss, k30). 6v18, not 10v30 (2026-11-09): its cases are only the 6 P samples
# (scaled GSP #30 >= 0.1) with junction files; 10v30's cases included 4 intermediates the SVM then scores (leakage).
# 461535, the 7th P sample, has no .junc. See diff_splicing/scripts/run_ylpm1_6v18.md. Logic: ../../scripts/phenotype_axis.R (ported from
# diff_splicing/ylpm1/top500_alt3_fc_vs_loading.Rmd). Loadings from the
# quetzal_code variant table YLPM1_full_info (vt_path). Writes cache/ylpm1_top500_alt3_5_fc5_psi.tsv
# and cache/ylpm1_phenotype_axis_scores.tsv. Run from figure4/:
#   Rscript --vanilla scripts/build_cache_ylpm1_phenotype.R [results_root]
source("../scripts/phenotype_axis.R")
source("../scripts/vt_paths.R")
args <- commandArgs(TRUE)
info <- fread(vt_path("YLPM1_full_info"), data.table = FALSE)
pa_build("ylpm1", "6v18", c("alt 3'ss", "alt 5'ss"), info, "k30", "ylpm1_top500_alt3_5_fc5_psi", "cache",
         results_root = if (length(args)) args[1] else
           file.path(QROOT, "leafcutter_analyses/diff_splicing/results"))
