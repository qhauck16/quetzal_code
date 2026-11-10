# figS9/build.R -- Supplementary Figure S9

## Purpose & context
Overlap of splice-factor GSP driving junctions with LeafCutter results. 3x3
grid: rows GSP #10 (SF3B1, alt 3'ss, n = 168) / #13 (NKAP, alt 3'ss, 203) /
#30 (YLPM1, alt 3'/5'ss, 188); columns LeafCutter runs top-N vs 3N matched
controls (10v30 / 20v60 / 30v90). Each cell = `overlap_heatmap()` with one row
(the GSP driving junctions), up/down columns, and a binomial test of up vs down
among the hits (`binom.test(up, up + down, p = 0.5)`). Original: `flashier/spec_flash_fact.Rmd`
("flashier-vs-leafcutter GSP comparison panels").

## Decision rationale
Counts and p come from `cache/build_overlap_cache.R`. Until 2026-11-10 each panel also
had an RPM/PSI-matched MANE-intron row and a Fisher test of GSP vs matched; Quinn dropped
it because canonical MANE introns are not a like-for-like baseline for unannotated
driving junctions. Backup of the old cache: `cache/.bak_s9_2026-11-10/`.

## Gotchas & warnings
- Row labels only in the left column: repeated per panel they squeezed tiles to slivers.
- Panel order = order of `panel` in the cache TSV (9 panels).
- Fisher p is on each panel's title line ("**SF3B1 10v30**  p = 3.7 x 10^-43",
  real superscript via markdown), NOT in overlap_heatmap()'s caption: a
  caption sits below the x title, so the bottom row's p hung lower than the
  others. Caption is nulled locally; common.R untouched.
- gridtext drops `&nbsp;`; the title gap is a literal em space (\u2003).
- One "LeafCutter direction" via `plot_layout(axis_titles = "collect")`.

## User preferences
Gene symbols italic in panel titles (`gene_md`); one collected legend.

## Input/output
In: `cache/s7_leafcutter_cells.tsv` (panel, set, direction, n_overlap, n_set,
binom_p, fact). Out: `out/figS9.pdf` (8.5 x 4.6).

## File flow schematic
```
cache/build_overlap_cache.R ──► cache/s7_leafcutter_cells.tsv ──► figS9/build.R ──► out/figS9.pdf
```
