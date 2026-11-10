# figS8/build.R -- Supplementary Figure S8

## Purpose & context
LeafCutter summary: counts of splice-event types (alt 3'ss, alt 5'ss,
canonical, skipped exon, other) among differentially spliced junctions for
SF3B1 / NKAP / YLPM1 and each gene's cohort-proportion-matched control run.
Facets: rows = dPSI > {0.01, 0.05} x FC > {1, 2.5, 5}; columns = 10v30 / 20v60 / 30v90.
Original: `flashier/diff_splicing/unified_results/leafcutter_meta_analysis.Rmd`.

## Decision rationale
Counts built by `cache/build_leafcutter_counts.R` (the Rmd's code verbatim);
this script only plots `s8_counts.tsv` (540 rows).

## Gotchas & warnings
- Colours = house `PAL$splice_type` (matches fig 2e) + canonical = Set2 yellow
  (#FFD92F). The published S8 used plain Set2 order, which made canonical the
  lavender that fig 2 uses for "other" (changed 2026-09-24 polish pass).
- Free y per row facet; total per bar group printed under the bars at 6 pt,
  inside the band opened by the 11% lower y expansion (previously printed on
  the axis line). A plot caption explains the totals.
- Threshold strips horizontal (angle 0), not vertical.

## User preferences
Totals printed on the chart, with thousands separators.

## Input/output
In: `cache/s8_counts.tsv` (type, n, gene, condition real/control, comparison,
dpsi, fc). Out: `out/figS8.pdf` (8.5 x 11).

## File flow schematic
```
diff_splicing/{results,control_results}/ ──► cache/build_leafcutter_counts.R ──► cache/s8_counts.tsv ──► figS8/build.R ──► out/figS8.pdf
```

**2026-09-24 (Quinn review):** removed the explanatory caption under the plot; per-group totals under the bars kept.
