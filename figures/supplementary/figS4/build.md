# figS4/build.R -- Supplementary Figure S4

## Purpose & context
Junctions driving GSP #3 (the read-length GSP). Left: driving junctions with
mean PSI <= 0.1% in 48 bp samples; right: same for 76 bp. Per cohort, ECDF of
per-junction mean PSI in each read-length bucket. Original:
`flashier/just_for_fun/nmd_investigation/covariate_gsp_figures.Rmd`
(`per_cohort_cdf`, `per_cohort_cdf_76`).

## Decision rationale
All heavy work (3.7 GB flashier RDS, per-gene fits) is in
`cache/build_readlength_cache.R`; this script only plots the two TSVs.

## Gotchas & warnings
- Layout: two rows (48 bp-low, 76 bp-low) x three cohorts, collected legend,
  6.5 x 5.6 in. Replaced two side-by-side 2x2 facet blocks that each had an
  empty fourth slot; the overall "Junctions driving GSP #3" title was dropped
  (it is the figure caption).
- View is zoomed with `coord_cartesian(xlim 0-0.5, ylim 0.75-1)`, not scale
  limits, so the ECDF is computed on all data.
- Only cohorts with >= 50 samples in both buckets appear (COAD, READ, UCEC per the cache log).
- Header comment says "<= 50 / >= 75 bp"; the legend/titles say 48 / 76 bp (the cache's `rl_group` labels).

## User preferences
tab10 blue/orange for 48/76 bp; square panels; legend collected at bottom.

## Input/output
In: `cache/s4_low48_per_cohort.tsv`, `cache/s4_low76_per_cohort.tsv`
(project_id, rl_group, mean_psi, ...). Out: `out/figS4.pdf` (8.5 x 5.2).

## File flow schematic
```
cache/build_readlength_cache.R ──► cache/s4_low48_per_cohort.tsv ─┬──► figS4/build.R ──► out/figS4.pdf
                               └─► cache/s4_low76_per_cohort.tsv ─┘
```
