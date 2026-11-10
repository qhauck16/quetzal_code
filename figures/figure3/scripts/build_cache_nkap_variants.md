# build_cache_nkap_variants.R

Thin driver: `build_nkap_tables()` from `../../scripts/variant_tables.R` → `figure3/cache/`.
See `scripts/variant_tables.md` for the logic, provenance and gotchas. Run from `figure3/`
after `scripts/build_gsp_loadings.R`, and re-run whenever the unified calling pipeline's
VEP output changes. It reads the 560 MB mutation table, so use caslake rather than a login
node for comfort (~4 min).
