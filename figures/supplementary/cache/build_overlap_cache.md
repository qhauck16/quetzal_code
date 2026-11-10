# cache/build_overlap_cache.R

## Purpose & context
Cache for S6 (Alsafadi 2016 overlap + Venn) and S7 (GSP vs LeafCutter).
Logic from `flashier/spec_flash_fact.Rmd` ("alsafadi GSP10 overlap",
"matched background junctions by RPM and PSI", "4-cell overlap heatmap",
"venn diagram", "flashier-vs-leafcutter GSP comparison panels").

## Decision rationale
Top-250 driving-junction sets are read from `fact_{10,13,30}_top_junctions.tsv`
(the Rmd's own saved `build_gsp_top()` + `assign_types()` output) instead of
rebuilding from the 3.7 GB RDS and 250 per-gene fits. Set sizes match the
published panels: GSP #10 alt 3'ss n = 168, #13 = 203, #30 alt 3'/5'ss = 188.

## Gotchas & warnings
- Controls matched on RPM and PSI from `unified_results/shared/all_junction_rpm_and_psi.tsv`.
- Log (S6): GSP up 84/168, matched up 0; Fisher p 8.0e-32; Venn 1219/84/84.
- Log (S7) overlaps: SF3B1 112/125/128 of 168; NKAP 141/158/155 of 203;
  YLPM1 62/73/71 of 188 (10v30/20v60/30v90).
- Light enough to run on a login node (no RDS).

## User preferences
Absolute paths

## Input/output
In: `hg38_granges.RDS`, `fact_{10,13,30}_top_junctions.tsv`,
`all_junction_rpm_and_psi.tsv`, `Alsafadi2016_supp1_hg38.tsv`,
`diff_splicing/results/<GENE>_<comp>/leafcutter_ds_{effect_sizes,cluster_significance}.txt`.
Out: `s6_alsafadi_cells.tsv`, `s6_venn.tsv`, `s7_leafcutter_cells.tsv`. Run from any cwd.

## File flow schematic
```
fact_{10,13,30}_top_junctions.tsv ─┐                           ┌─► s6_alsafadi_cells.tsv, s6_venn.tsv ──► figS6
Alsafadi2016_supp1_hg38.tsv ───────┼──► build_overlap_cache.R ──┤
all_junction_rpm_and_psi.tsv ──────┤                           └─► s7_leafcutter_cells.tsv ──► figS9
diff_splicing/results/*/ ──────────┘
```
