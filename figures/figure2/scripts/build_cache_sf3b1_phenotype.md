# build_cache_sf3b1_phenotype.R

## Purpose & context
Thin driver for `../scripts/phenotype_axis.R` (read that .md for logic and
rationale). Builds the SF3B1 phenotype-axis caches: comparison 30v90,
splice types alt 3'ss, loading column k10 from `SF3B1_full_info`. Replaces
`diff_splicing/sf3b1/top500_alt3_fc_vs_loading.Rmd`.

## Gotchas & warnings
- Run from `figure2/`: `Rscript --vanilla scripts/build_cache_sf3b1_phenotype.R [results_root]`.
- Must be re-run whenever `diff_splicing/results/sf3b1_30v90` or the loading table changes.
- Heavy-ish (reads hundreds of psi_table.tsv); run on caslake, not the login node.

## Input/output
Out: `figure2/cache/sf3b1_top500_alt3_fc5_psi.tsv`, `figure2/cache/sf3b1_phenotype_axis_scores.tsv` -> panels 2f,S13,S18.

## File flow
```
results/sf3b1_30v90 ──┐
SF3B1_full_info ─────────────┼──► build_cache_sf3b1_phenotype.R ──► cache/sf3b1_top500_alt3_fc5_psi.tsv
psi tables ─────┘                                  └► cache/sf3b1_phenotype_axis_scores.tsv
```
