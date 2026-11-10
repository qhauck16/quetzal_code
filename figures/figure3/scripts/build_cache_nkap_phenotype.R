# NKAP top-500 alt-splice-site PSI matrix + SVM phenotype-axis scores
# (20v60, alt 3'ss, k13). Logic: ../../scripts/phenotype_axis.R (ported from
# diff_splicing/nkap/top500_alt3_fc_vs_loading.Rmd). Loadings from the
# quetzal_code variant table NKAP_full_plus_variants (vt_path). Writes cache/nkap_top500_alt3_fc5_psi.tsv
# and cache/nkap_phenotype_axis_scores.tsv. Run from figure3/:
#   /scratch/midway2/qhauck/conda_env/rstudio-server/bin/Rscript --vanilla scripts/build_cache_nkap_phenotype.R [results_root]
source("../scripts/phenotype_axis.R")
source("../scripts/vt_paths.R")
args <- commandArgs(TRUE)
info <- fread(vt_path("NKAP_full_plus_variants"), data.table = FALSE)
pa_build("nkap", "20v60", c("alt 3'ss"), info, "k13", "nkap_top500_alt3_fc5_psi", "cache",
         results_root = if (length(args)) args[1] else
           file.path(QROOT, "leafcutter_analyses/diff_splicing/results"))
