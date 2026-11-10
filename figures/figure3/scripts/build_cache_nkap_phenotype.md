# build_cache_nkap_phenotype.R

## Purpose & context
Thin driver for `../scripts/phenotype_axis.R` (read that .md for logic and
rationale). Builds the NKAP phenotype-axis caches: comparison 20v60,
splice types alt 3'ss, loading column k13 from `NKAP_full_plus_variants`. Replaces
`diff_splicing/nkap/top500_alt3_fc_vs_loading.Rmd`.

## Gotchas & warnings
- Run from `figure3/`: `Rscript --vanilla scripts/build_cache_nkap_phenotype.R [results_root]`.
- Must be re-run whenever `diff_splicing/results/nkap_20v60` or the loading table changes.
- Heavy-ish (reads hundreds of psi_table.tsv); run on caslake, not the login node.

## Input/output
Out: `figure3/cache/nkap_top500_alt3_fc5_psi.tsv`, `figure3/cache/nkap_phenotype_axis_scores.tsv` -> panels S13,S18.

## File flow
```
results/nkap_20v60 ──┐
NKAP_full_plus_variants ─────────────┼──► build_cache_nkap_phenotype.R ──► cache/nkap_top500_alt3_fc5_psi.tsv
psi tables ─────┘                                  └► cache/nkap_phenotype_axis_scores.tsv
```
