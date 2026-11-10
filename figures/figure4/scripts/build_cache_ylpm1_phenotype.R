# YLPM1 top-500 alt-splice-site PSI matrix + SVM phenotype-axis scores
# (6v18, alt 3'ss + alt 5'ss, k30). Cases are only the 6 P samples (scaled GSP #30 >= 0.1)
# with junction files, so no sample the SVM scores is also a training case.
# 461535, the 7th P sample, has no .junc. Logic: ../../scripts/phenotype_axis.R.
# Loadings from the variant table YLPM1_full_info (vt_path). Writes cache/ylpm1_top500_alt3_5_fc5_psi.tsv
# and cache/ylpm1_phenotype_axis_scores.tsv. Run from figure4/:
#   Rscript --vanilla scripts/build_cache_ylpm1_phenotype.R [results_root]
source("../scripts/phenotype_axis.R")
source("../scripts/vt_paths.R")
args <- commandArgs(TRUE)
info <- fread(vt_path("YLPM1_full_info"), data.table = FALSE)
pa_build("ylpm1", "6v18", c("alt 3'ss", "alt 5'ss"), info, "k30", "ylpm1_top500_alt3_5_fc5_psi", "cache",
         results_root = if (length(args)) args[1] else
           file.path(QROOT, "leafcutter_analyses/diff_splicing/results"))
