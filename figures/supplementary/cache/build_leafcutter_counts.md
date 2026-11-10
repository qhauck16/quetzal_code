# cache/build_leafcutter_counts.R

## Purpose & context
Cache for S8. The config / mane / helpers / load_leaf_sigs / counts chunks of
`flashier/diff_splicing/unified_results/leafcutter_meta_analysis.Rmd`, verbatim.
Counts splice types among DS junctions per gene x real/control x comparison x
(dpsi, fc) threshold.

## Decision rationale
Verbatim copy (not a re-implementation) so the counts equal the Rmd's.

## Gotchas & warnings
- `setwd()`s itself to `.../diff_splicing/unified_results` so the Rmd's relative
  paths (`../results`, `../control_results`, `../gencode_v46_introns.tsv`,
  `../MANE.GRCh38.v1.4.summary.txt`) resolve -- runnable from any cwd.
- Output path is absolute.
- Log: "50 or more warnings" (inherited); output tibble 540 rows.

## User preferences
Absolute output path

## Input/output
In: LeafCutter `leafcutter_ds_effect_sizes.txt` / `leafcutter_ds_cluster_significance.txt`
under `results/` and `control_results/<gene>_<comparison>`, introns, MANE.
Out: `s8_counts.tsv` (type, n, gene, condition, comparison, dpsi, fc).

## File flow schematic
```
diff_splicing/results/*/ ─────────┐
diff_splicing/control_results/*/ ─┼──► build_leafcutter_counts.R ──► s8_counts.tsv ──► figS8/build.R
gencode_v46_introns.tsv, MANE ────┘
```
