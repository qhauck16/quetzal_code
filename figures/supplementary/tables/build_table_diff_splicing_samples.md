# build_table_diff_splicing_samples.R

## Purpose & context
Supplementary table listing every TCGA sample used in the nine GSP differential-splicing
comparisons (SF3B1, NKAP and YLPM1, each at 10v30, 20v60 and 30v90), tests and controls, with
the exact GDC file analysed. Quinn chose (2026-10-08) to document the sample set rather than
rerun the comparisons.

## Decision rationale
- **Why document instead of rerun.** 302 of the 540 control slots (138 unique BAMs) were matched
  from BAMs downloaded in earlier exploratory analyses, not from the 2026-04-20 planned pool
  (`download/sample_plan.tsv`). Those controls can't be regenerated from code. Listing them
  makes the analysis auditable.
- **Built from what actually ran.** The source is the LeafCutter sample sheets
  (`results/<comp>/leaf_pheno_data.txt`), not from re-running any selection.
- **`download_source` column.** Marks each sample "planned" or "earlier". All tests are planned.

## Gotchas & warnings
- **Relabelled BAMs.** Two BAMs were relabelled in the 2026-09-25 aliquot fix, and
  `used_bams.rail_map.tsv` indexes them under their OLD rail_id:
  - 111276 is now 443609
  - 116871 is now 443194

  They are translated through `bam_relabel_2026-09-25.tsv`, and the barcode comes from
  `aliquot_on_disk`, which was verified against the BAM's SM tag.
- **Rerun BAMs.** BAMs re-downloaded for the 2026-10-05 rerun come from
  `rerun_2026-10-05.rail_map.tsv`; md5 is not recorded there, so it is NA for those rows.
- **Validated (2026-10-08):**
  - all 720 barcodes equal `tcga_v2_samples` barcodes;
  - no sample is missing a loading;
  - every control has loading < 1e-4 on its GSP.
- **Not included:** the matched-vs-matched null runs in `control_results/`.

## Input/output
- **In:** `diff_splicing/results/*/leaf_pheno_data.txt`, `download/{used_bams.rail_map,
  rerun_2026-10-05.rail_map, bam_relabel_2026-09-25, sample_plan}.tsv`, the variant tables
  (loadings, via `vt_path`), and `tcga_v2_samples`.
- **Out:** `supplementary/out/tables/STable_diff_splicing_samples.tsv`, one row per
  (comparison, sample).

## File flow
```
results/*/leaf_pheno_data.txt ─┐
download/*rail_map*, plan ─────┼─► build_table_diff_splicing_samples.R ─► out/tables/STable_diff_splicing_samples.tsv
variant tables, tcga metadata ─┘
```
