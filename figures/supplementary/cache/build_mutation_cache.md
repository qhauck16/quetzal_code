# cache/build_mutation_cache.R

## Purpose & context
One-off cache for S9 and S12. Reproduces `flashier/flash_fact_char.Rmd`'s
AlphaMissense-filtered (am_pathogenicity > 0.9 or NA) sample x gene mutation
matrix and its cohort-stratified one-sided CMH tests ("SF3B1-axis CMH",
"SUGP1 along GSP 10", "NKAP finding other mutants along the same axis",
co-occurrence heatmap), verbatim in logic.

## Decision rationale
Loadings taken from `SF3B1_full_plus_variants.tsv` (k10) and
`NKAP_full_plus_variants.tsv` (k13) (= norm_L), so the 3.7 GB flashier RDS is never loaded.

## Gotchas & warnings
- Genes tested only if mutated in > 2 top samples.
- Log shows many-to-many join warnings on `tcga_metadata` (rail_id/sampleId);
  inherited from the Rmd, left as-is.
- Log: `mut_mat: 9970 19585`; `S9 SUGP1: 2.72e-83`; S12 sig genes
  FAM50A, PRPF8, SLU7, KAT5, SBNO2.
- Header's sbatch line uses bare `Rscript`; use the conda env's absolute path.

## User preferences
Absolute paths; `.libPaths()` at top.

## Input/output
In: `tcga/tcga_v2_samples.tsv`, `flashier/just_for_fun/mutation_merge/merged_grch38_with_am.tsv`
(536 M), `{SF3B1,NKAP}_full_plus_variants.tsv`. Out: `s9_cmh_sf3b1_excluded.tsv`,
`s9_top100_gsp10.tsv`, `s12_cmh_nkap_excluded.tsv`, `s12_cooccurrence.tsv`.

Run (any cwd): `sbatch -p build --account=pi-yangili1 --mem=24G --wrap "/scratch/midway2/qhauck/conda_env/rstudio-server/bin/Rscript --vanilla cache/build_mutation_cache.R"`

## File flow schematic
```
merged_grch38_with_am.tsv ─────────┐                          ┌─► s9_*.tsv  ──► figS10/build.R
tcga_v2_samples.tsv ───────────────┼──► build_mutation_cache.R ┤
{SF3B1,NKAP}_full_plus_variants ───┘                          └─► s12_*.tsv ──► figS12/build.R
```
