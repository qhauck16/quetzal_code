# Variant tables for NKAP / PRPF8 / FAM50A / SLU7: k13_full_info, NKAP_full_plus_variants
# Logic: ../../scripts/variant_tables.R (ported from flash_fact_char.Rmd).
# Light: reads shared_cache/gsp_loadings.tsv, the AlphaMissense-merged
# mutation table and the unified pipeline's VEP output. Run from figure3/:
#   Rscript --vanilla scripts/build_cache_nkap_variants.R
source("../scripts/variant_tables.R")
I <- vt_inputs()
write_variant_tables(build_nkap_tables(I), "cache")
