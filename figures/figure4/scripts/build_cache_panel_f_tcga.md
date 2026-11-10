# build_cache_panel_f_tcga.R

## Purpose & context
Builds the two TCGA inputs of figure 4f's GIST heatmap that previously came
from `diff_splicing/ylpm1/leafcutter_prep_and_analysis.Rmd`:
`high_psi_hits` (which junctions are drawn) and the TCGA logit-PSI matrix
whose heatmaply column clustering sets the column order of 4f.

## Decision rationale
- Reuses `pa_load_comparison` / `pa_assign_splice_types` from
  `scripts/phenotype_axis.R` with the Rmd's 4f thresholds (YLPM1 20v60,
  alt 3'ss + 5'ss, dPSI > 0.05, FC > 2.5) instead of copying them again.
- Validated on the pre-fix results: identical 561 hits, identical 512 columns.

## Gotchas & warnings
- The Rmd drew its 200 project-matched controls with an UNSEEDED
  slice_sample; here it is `set.seed(30)`. The column order can therefore
  differ slightly from the published one even on identical inputs.
- After rerunning, delete `cache/panel_f_{col_order,tcga_cluster_ids,data}.rds`
  so `panel_f/build.R` rebuilds them (it only builds when missing).
- `panel_f/build.R` now prints a message, rather than stopping, when the
  matrix is not 17 x 512 or the first column is not PTPRO.

## Input/output
In: `results/ylpm1_20v60`, `YLPM1_full_plus_variants` (vt_path), TCGA psi tables.
Out: `cache/panel_f_high_psi_hits.tsv`, `cache/panel_f_tcga_mat.tsv`
(213 samples x junctions), `cache/panel_f_tcga_rail_ids.tsv`.

## File flow
```
results/ylpm1_20v60 ──┐
YLPM1_full_plus_variants ┼──► build_cache_panel_f_tcga.R ──► cache/panel_f_high_psi_hits.tsv ──► panel_f/build.R
psi tables ───────────┘                                └► cache/panel_f_tcga_mat.tsv ──► (heatmaply col order)
```
