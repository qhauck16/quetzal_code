# Supplementary Table 2: targeted mutation re-calling (Methods, "Additional
# mutation calling"). One row per sample sent to re-calling for GSP #10, #13 and
# #30: test samples (loading > 0.01, no known driver of that GSP's phenotype) and
# cohort-matched controls (loading < 1e-4). Mutations are those passing every
# filter (Mutect2 PASS/multiallelic on WXS; DP > 10; VAF > 10%; VEP MODERATE/HIGH;
# AlphaMissense > 0.9 for missense; protein-truncating for YLPM1), as built into the
# figure variant tables by figures/scripts/variant_tables.R.
# Inputs: mutation_recalling/unified/<gene>/{samples,wxs_manifest,rna_manifest}.tsv,
#         figures/<fig>/cache/<GENE>_full_plus_variants.tsv (vt_path()).
# Output: figures/supplementary/out/tables/STable2_mutation_recalls.tsv
if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
suppressPackageStartupMessages({ library(data.table); library(dplyr) })
source(file.path(QROOT, "figures/scripts/vt_paths.R"))   # vt_path()

SETS <- list(list(gsp = 10, dir = "sf3b1", vt = "SF3B1_full_plus_variants"),
             list(gsp = 13, dir = "nkap",  vt = "NKAP_full_plus_variants"),
             list(gsp = 30, dir = "ylpm1", vt = "YLPM1_full_plus_variants"))
UNI <- file.path(QROOT, "mutation_recalling/unified")
na_blank <- function(x) ifelse(is.na(x) | x == "NA", "", gsub(",", ", ", x))

tab <- bind_rows(lapply(SETS, function(s) {
  smp <- fread(file.path(UNI, s$dir, "samples.tsv"), data.table = FALSE)
  wxs <- fread(file.path(UNI, s$dir, "wxs_manifest.tsv"), data.table = FALSE)
  rna <- fread(file.path(UNI, s$dir, "rna_manifest.tsv"), data.table = FALSE)
  vt  <- fread(vt_path(s$vt), data.table = FALSE) %>% select(rail_id, wxs_mutation, rnaseq_mutation)
  wxs_ok <- intersect(wxs$rail_id[wxs$sample_type == "Tumor"], wxs$rail_id[wxs$sample_type == "Normal"])
  smp %>% left_join(vt, by = "rail_id") %>%
    transmute(gsp = paste0("GSP #", s$gsp), rail_id, tcga_barcode, tcga_cohort = project_id,
              gsp_loading = loading,
              group = ifelse(role == "case", "test", "control"),
              wxs_tumor_normal_analysed = ifelse(rail_id %in% wxs_ok, "yes", "no"),
              rnaseq_analysed = ifelse(rail_id %in% rna$rail_id, "yes", "no"),
              mutations_wxs = na_blank(wxs_mutation),
              mutations_rnaseq = na_blank(rnaseq_mutation))
})) %>% arrange(gsp, desc(group == "test"), desc(gsp_loading))

dir.create(file.path(QROOT, "figures/supplementary/out/tables"), showWarnings = FALSE, recursive = TRUE)
out <- file.path(QROOT, "figures/supplementary/out/tables/STable2_mutation_recalls.tsv")
fwrite(tab, out, sep = "\t", quote = FALSE)
print(tab %>% group_by(gsp, group) %>%
        summarise(n = n(), wxs = sum(wxs_tumor_normal_analysed == "yes"),
                  rna = sum(rnaseq_analysed == "yes"),
                  with_mutation = sum(mutations_wxs != "" | mutations_rnaseq != ""), .groups = "drop"))
cat("wrote", out, "\n")
