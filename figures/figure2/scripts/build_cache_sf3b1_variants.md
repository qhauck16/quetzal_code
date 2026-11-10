# build_cache_sf3b1_variants.R

Thin driver: `build_sf3b1_tables()` from `../../scripts/variant_tables.R` → `figure2/cache/`.
See `scripts/variant_tables.md` for the logic, provenance and gotchas. Run from `figure2/`
after `scripts/build_gsp_loadings.R`, and re-run whenever the unified calling pipeline's
VEP output changes. It reads the 560 MB mutation table, so use caslake rather than a login
node for comfort (~4 min).
