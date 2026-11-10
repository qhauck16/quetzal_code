# figS6/build.R -- Supplementary Figure S6

## Purpose & context
GSP #10 vs Alsafadi et al. (2016) SF3B1-MUT vs WT uveal-melanoma junctions.
(a) 2x2 exact-coordinate overlap heatmap: top-250 GSP #10 alt 3'ss driving
unannotated junctions (n = 168) and RPM/PSI-matched controls vs Alsafadi
up/down; Fisher's exact on the upregulated column. (b) Venn of Alsafadi
significant upregulated alt 3'ss vs the GSP #10 set.
Original: `flashier/spec_flash_fact.Rmd`.

## Decision rationale
All set logic lives in `cache/build_overlap_cache.R`; this script only draws.
Heatmap via `overlap_heatmap()` in common.R (shared with S7).

## Gotchas & warnings
- Fisher p is the second line of panel a's x title; the common.R caption
  drifted out under the legend. No coord_fixed (square tiles were too narrow
  for "upregulated"/"downregulated" at a 3 in page height).
- Tags need `& theme(plot.tag = bold)` (see S2 build.md).
- Venn is hand-drawn (`ggforce::geom_circle`); titles are three lines anchored
  at the OUTER circle edges (hjust 0/1) -- two-line versions collided.
- Cache log values: GSP up 84/168, matched up 0/168, Fisher p = 8.0e-32;
  Venn only-Alsafadi 1219 / shared 84 / only-GSP 84.

## User preferences
Numbers on the chart (tiles and Venn regions), not in a legend.

## Input/output
In: `cache/s6_alsafadi_cells.tsv` (panel, set, direction, n_overlap, n_set,
fisher_p), `cache/s6_venn.tsv` (only_alsafadi, shared, only_gsp).
Out: `out/figS6.pdf` (8.5 x 3.4).

## File flow schematic
```
cache/build_overlap_cache.R ──► cache/s6_alsafadi_cells.tsv ─┬──► figS6/build.R ──► out/figS6.pdf
                            └─► cache/s6_venn.tsv ───────────┘
```
