# phenotype_axis.R

## Purpose & context
Shared library that turns a LeafCutter differential-splicing comparison into
(1) a top-500 alt-splice-site PSI matrix over all TCGA samples and (2) a
per-sample linear-SVM "phenotype axis" score. Ported 2026-10-06 from
`tcga/genome_wide/flashier/diff_splicing/{sf3b1,nkap,ylpm1}/top500_alt3_fc_vs_loading.Rmd`
so that quetzal_code is self-contained for release. Called only through the
three per-figure drivers (`figure{2,3,4}/scripts/build_cache_*_phenotype.R`).

## Decision rationale
- The three Rmds were near-identical copies, so they became one function set
  parameterised by gene / comparison / splice types / loading column.
- Logic was copied verbatim (classification, logit + scale, LiblineaR type 2,
  cost search toward ~10% FDR among training samples). Validated against
  the pre-fix leafcutter results (`diff_splicing/results.pre_fix_2026-10-05`):
  PSI matrices identical, scores within ~1e-14, groups identical (all 3 genes).
- One fix: the YLPM1 Rmd trained on a stale April PSI file
  (`ylpm1_top500_alt3_fc5_psi.tsv`); here the SVM always uses the matrix
  built in the same call. Same 467 junctions at the time, so no change.

## Gotchas & warnings
- `results_root` defaults to `diff_splicing/results` (derived from
  `paths.yaml: gencode_introns`). Pass another root (e.g. a `.pre_fix` backup)
  as the driver's first arg to rebuild against it.
- P/N training sets: loading >= 0.1 is P, <= 1e-4 is N; everything else is
  "intermediate" and is scored but not trained on.
- The cost search is a deterministic grid + refinement; LiblineaR is
  deterministic for type 2, so reruns are byte-stable.
- `select/filter/rename` are rebound to dplyr at the top (GenomicRanges
  masking elsewhere in the pipeline).

## User preferences
Quinn wants all figure-generating code inside quetzal_code; outputs land in
`figureN/cache/` and readers find them via `vt_path()` (scripts/vt_paths.R).

## Input/output
In: `results/<gene>_<comp>/leafcutter_ds_{cluster_significance,effect_sizes}.txt`,
GENCODE introns + MANE summary (paths.yaml), PSI tables
`<psi_root>/tcga/<chr>/<gene>/psi_table.tsv`, loading table (`*_full_info` /
`*_full_plus_variants` via vt_path).
Out: `cache/<psi_name>.tsv` (rail_id x `gene:cluster:start-end` PSI) and
`cache/<gene>_phenotype_axis_scores.tsv` (rail_id, k, group, score).
Readers: 2f, S13 (PSI matrices); S17, S18 (scores).

## File flow
```
results/<gene>_<comp>/ ──┐
gencode/MANE introns ────┼──► pa_build() ──► cache/<psi_name>.tsv ──► 2f, S13
psi_table.tsv (TCGA) ────┤                 └► cache/<gene>_phenotype_axis_scores.tsv ──► S17, S18
loading table (vt_path) ─┘
```

> 2026-11-10: `paths.yaml` and `scripts/paths.R` were removed. Paths are now
> `file.path(QROOT, "<folder>/<file>")`, with `QROOT` the quetzal_code root
> (found via the `.quetzal_root` marker); references to paths.yaml keys above are historical.
