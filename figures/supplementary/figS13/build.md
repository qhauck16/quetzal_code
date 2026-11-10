# figS13/build.R -- Supplementary Figure S13

## Purpose & context
GSP #13 and #30 characterisation. (a, b) pies of splice types of the driving
unannotated junctions of the top 250 single-gene groups on GSP #13 / #30
(`flashier/spec_flash_fact.Rmd`, "plotting information about splice site
selection"). (c, d) per-sample mean PSI of the top 500 LeafCutter junctions as
fold change over the global mean, boxplots by log-loading bin (0.5-wide bins,
(-Inf,-8] ... (-0.5,0]) -- NKAP 20v60 alt 3'ss and YLPM1 10v30 alt 3'ss + alt 5'ss
(`diff_splicing/{nkap,ylpm1}/top500_alt3_fc_vs_loading.Rmd`).

## Decision rationale
- Pies use the shared `../scripts/pie_chart.R`; drawn pie radius measured off
  the render and iterated to a fixed point, as for the main-figure pies.
- YLPM1 uses alt 3'ss AND alt 5'ss (project rule for YLPM1).

## Gotchas & warnings
- `fact_{13,30}_top_junctions.tsv` has one row per snaptron_id.
- Empty bins kept (`scale_x_discrete(drop = FALSE)`).

## User preferences
Pie numbers on the chart; angled x labels.

## Input/output
In: `$FLASHIER/fact_{13,30}_top_junctions.tsv`,
`diff_splicing/nkap/nkap_top500_alt3_fc5_psi.tsv`,
`diff_splicing/ylpm1/ylpm1_top500_alt3_5_fc5_psi.tsv`,
`{NKAP,YLPM1}_full_plus_variants.tsv` (log k13 / k30). Out: `out/figS13.pdf` (8.5 x 7.2).

## File flow schematic
```
fact_{13,30}_top_junctions.tsv ───────────┐
{nkap,ylpm1}_top500_*_psi.tsv ────────────┼──► figS13/build.R ──► out/figS13.pdf
{NKAP,YLPM1}_full_plus_variants.tsv ──────┘
```

## Polish notes (2026-09)
- Panel letters: `plot_annotation(theme=)` does not style patchwork tags
  (each tag takes its own panel's theme), so they rendered plain and at a
  different size on the wrapped pies vs the boxplots. Fixed with
  `& theme(plot.tag = ...)` after plot_annotation_publication().
