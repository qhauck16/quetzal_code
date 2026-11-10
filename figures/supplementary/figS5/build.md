# figS5/build.R -- Supplementary Figure S5

## Purpose & context
Extended GSP x covariate heatmap: all 150 top GSPs (rows) x 33 TCGA cohorts +
technical covariates + 4 summary columns, Pearson r. Extends main Figure 2c
(`flashier/flash_fact_char.Rmd` correlation-heatmap chunks) from 60 to 150 GSPs.

## Decision rationale
Matrix from `cache/build_s5_heatmap_cache.R` (fig 2c's cache builder with
N_TOP = 150), so the statistic is identical to fig 2c.

## Gotchas & warnings
- GSP row labels at 0.8 x SIZE_BASE (6.4 pt): 0.7 was below the 6 pt floor.
- **Known deviation:** uses the SAME +/-0.4 colour scale (squished) and group-strip
  colours as main Fig 2c. The published S5 used +/-0.6 and different strip colours.
- **Known deviation:** tissue columns ordered alphabetically; rows in GSP order
  (top = GSP 1), not clustered/tissue-ordered as in 2c.
- Technical/summary columns ordered by `display_order` from `s5_cols.tsv`.
- Row ticks every 10 GSPs so readers can find a GSP.

## User preferences
`heatmap_diverging()` palette from theme.R; angled x labels; absolute-inch page.

## Input/output
In: `cache/s5_matrix.tsv` (gsp + one column per covariate), `cache/s5_cols.tsv`
(column, label, group, display_order). `s5_rows.tsv` is written but unused here.
Out: `out/figS5.pdf` (8.5 x 9).

## File flow schematic
```
cache/build_s5_heatmap_cache.R ──► cache/s5_matrix.tsv ─┬──► figS5/build.R ──► out/figS5.pdf
                               └─► cache/s5_cols.tsv ───┘
```

**2026-10-05:** colour scale now +-0.6 (same as main Fig 2c, which moved from +-0.4; also matches the published S5).
