# cache/build_s5_heatmap_cache.R

## Purpose & context
Cache for S5. Copy of `../figure2/scripts/build_cache_panel_c.R` (main Fig 2c
GSP x covariate Pearson heatmap; `flashier/flash_fact_char.Rmd` correlation
chunks) with `N_TOP = 150` and outputs redirected to `cache/s5_*.tsv`.

## Decision rationale
A copy (not a parameterised call) so fig 2c's signed-off cache is untouched.
Upstream chain is verbatim, quirks included -- do not "clean it up".

## Gotchas & warnings
- Freshness guard: exits immediately if all three outputs are newer than the
  flashier RDS. Delete `s5_*.tsv` to force a rebuild.
- Needs tens of GB (3.5 G RDS + 536 M mutation table) -> SLURM build partition.
- The `[panel_c] mem ...` log lines over-report (units bug in `mem()`); ignore.
- Log: 9800 samples x 300 GSPs; 300 x 42 covariates correlated; writes 150 x 46
  (33 tissue / 9 technical / 4 summary), 54 rows labelled.

## User preferences
Absolute paths; `.libPaths()` at top.

## Input/output
In: `softimpute_flash_300_qc_filtered.RDS`, `tcga_v2_samples.tsv`,
`merged_grch38_with_am.tsv`, `standardized_covariates.RDS`.
Out: `s5_matrix.tsv`, `s5_cols.tsv`, `s5_rows.tsv` (rows unused by figS5). Any cwd.

## File flow schematic
```
flashier RDS ────────────────┐
merged_grch38_with_am.tsv ───┼──► build_s5_heatmap_cache.R ──► s5_{matrix,cols,rows}.tsv ──► figS5/build.R
standardized_covariates.RDS ─┤
tcga_v2_samples.tsv ─────────┘
```
