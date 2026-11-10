# build_cache_ylpm1_variants.R

Thin driver: `build_ylpm1_tables()` from `../../scripts/variant_tables.R` → `figure4/cache/`.
See `scripts/variant_tables.md` for the logic, provenance and gotchas. Run from `figure4/`
after `scripts/build_gsp_loadings.R`, and re-run whenever the unified calling pipeline's
VEP output changes. It reads the 560 MB mutation table, so use caslake rather than a login
node for comfort (~4 min).
