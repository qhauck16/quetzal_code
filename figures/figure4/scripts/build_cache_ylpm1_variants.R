# Variant tables for YLPM1: YLPM1_full_info, YLPM1_full_plus_variants
# Logic: ../../scripts/variant_tables.R (ported from flash_fact_char.Rmd).
# Light: reads shared_cache/gsp_loadings.tsv, the AlphaMissense-merged
# mutation table and the unified pipeline's VEP output. Run from figure4/:
#   Rscript --vanilla scripts/build_cache_ylpm1_variants.R
source("../scripts/variant_tables.R")
I <- vt_inputs()
write_variant_tables(build_ylpm1_tables(I), "cache")
