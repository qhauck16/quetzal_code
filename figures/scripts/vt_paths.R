# Where each per-gene variant table lives (built by
# figure{2,3,4}/scripts/build_cache_*_variants.R from
# scripts/variant_tables.R). Readers call vt_path("<name>").
if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
VT_HOME <- c(SF3B1_full_info          = "figure2",
             SF3B1_SUGP1_full_info    = "figure2",
             SF3B1_full_plus_variants = "figure2",
             k13_full_info            = "figure3",
             NKAP_full_plus_variants  = "figure3",
             YLPM1_full_info          = "figure4",
             YLPM1_full_plus_variants = "figure4",
             # top-500 PSI matrices + SVM phenotype-axis scores
             # (scripts/phenotype_axis.R)
             sf3b1_top500_alt3_fc5_psi     = "figure2",
             sf3b1_phenotype_axis_scores   = "figure2",
             nkap_top500_alt3_fc5_psi      = "figure3",
             nkap_phenotype_axis_scores    = "figure3",
             ylpm1_top500_alt3_5_fc5_psi   = "figure4",
             ylpm1_phenotype_axis_scores   = "figure4")
vt_path <- function(name) {
  stopifnot(name %in% names(VT_HOME))
  file.path(file.path(QROOT, "figures"),
            VT_HOME[[name]], "cache", paste0(name, ".tsv"))
}
