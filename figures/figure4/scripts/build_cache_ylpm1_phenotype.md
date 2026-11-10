# build_cache_ylpm1_phenotype.R

## Purpose & context
Thin driver for `../scripts/phenotype_axis.R` (read that .md for logic and
rationale). Builds the YLPM1 phenotype-axis caches: comparison 10v30,
splice types alt 3'ss + alt 5'ss, loading column k30 from `YLPM1_full_info`. Replaces
`diff_splicing/ylpm1/top500_alt3_fc_vs_loading.Rmd`.

## Gotchas & warnings
- Run from `figure4/`: `Rscript --vanilla scripts/build_cache_ylpm1_phenotype.R [results_root]`.
- Must be re-run whenever `diff_splicing/results/ylpm1_10v30` or the loading table changes.
- Heavy-ish (reads hundreds of psi_table.tsv); run on caslake, not the login node.

## Input/output
Out: `figure4/cache/ylpm1_top500_alt3_5_fc5_psi.tsv`, `figure4/cache/ylpm1_phenotype_axis_scores.tsv` -> panels S13,S17,S18.

## File flow
```
results/ylpm1_10v30 ──┐
YLPM1_full_info ─────────────┼──► build_cache_ylpm1_phenotype.R ──► cache/ylpm1_top500_alt3_5_fc5_psi.tsv
psi tables ─────┘                                  └► cache/ylpm1_phenotype_axis_scores.tsv
```
