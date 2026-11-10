# Variant tables for SF3B1 (+ SUGP1): SF3B1_full_info, SF3B1_SUGP1_full_info, SF3B1_full_plus_variants
# Logic: ../../scripts/variant_tables.R (ported from flash_fact_char.Rmd).
# Light: reads shared_cache/gsp_loadings.tsv, the AlphaMissense-merged
# mutation table and the unified pipeline's VEP output. Run from figure2/:
#   Rscript --vanilla scripts/build_cache_sf3b1_variants.R
source("../scripts/variant_tables.R")
I <- vt_inputs()
write_variant_tables(build_sf3b1_tables(I), "cache")
