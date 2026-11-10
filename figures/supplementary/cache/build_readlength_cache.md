# cache/build_readlength_cache.R

## Purpose & context
Cache for S4 (GSP #3 read-length junctions). The paths / top5 / meta /
build_junction_set / per_cohort_cdf / build_junction_set_76 / per_cohort_cdf_76
chunks of `flashier/just_for_fun/nmd_investigation/covariate_gsp_figures.Rmd`, verbatim.

## Decision rationale
Verbatim so the junction sets match the Rmd; its `ggsave` calls are disabled
and only the per-cohort mean-PSI tables are written.

## Gotchas & warnings
- `setwd()`s itself to `just_for_fun/nmd_investigation` (FLASH_RDS is relative
  there); runnable from any cwd.
- Loads the 3.7 GB flashier RDS -> SLURM build partition.
- Driving junctions: top 200 (gene, local_factor) groups on GSP #3, within a
  factor of 2 of the group's top junction; mean PSI <= 0.001 in the 48 (or 76) bp bucket.
- Log: 282 junctions (48 bp set), 45 (76 bp set); cohorts with >= 50 samples in
  both buckets: COAD (354/192), UCEC (209/380), READ (105/72).

## User preferences
`.libPaths()` at top; absolute output dir.

## Input/output
In: `softimpute_flash_300_qc_filtered.RDS`, `tcga_v2_samples.tsv`,
`260116_filters` per-gene fits, leaf2_try `psi_table.tsv`s.
Out: `s4_low48_per_cohort.tsv`, `s4_low76_per_cohort.tsv`.

## File flow schematic
```
flashier RDS (3.7 GB) ─┐
per-gene res.RDS ──────┼──► build_readlength_cache.R ──► s4_low{48,76}_per_cohort.tsv ──► figS4/build.R
psi_table.tsv ─────────┘
```
