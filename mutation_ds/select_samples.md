# select_samples.R

## Purpose & context
Supplementary note (IDH1 / SPOP). Builds the mutant-vs-wild-type sample sheets for the
LeafCutter differential-splicing comparisons. The DS junctions are then tested for overlap
with the associated GSP: IDH1 R132 in LGG against GSP #1, and SPOP residues in PRAD against
GSP #12. Quinn asked for this on 2026-10-08, starting from scratch rather than from the old
`diff_splicing/{idh1,spop}` runs.

## Decision rationale
- **1:3 matched controls**, the same design as the GSP comparisons: MatchIt nearest
  neighbour on mean read length + FFPE, within one cohort.
- **IDH1:** 20 mutants drawn at random (seed 123) from 406 R132 carriers.
- **SPOP:** every carrier of each residue: F133 (20), W131 (14), F102 (8), Y87 (7).
  PRAD has 439 SPOP-wild-type samples, so 1:3 is feasible for all four. Quinn decided
  against resampling.
- **Controls** have mutation calls but no call anywhere in the gene.
- **Same patient allowed.** Several tumour samples of one patient (e.g. 01 + 06) may all
  be used, for consistency with the first paper (Quinn, 2026-10-08).
- **One RNA aliquot per tumour sample.** 01A is preferred over 01B, then a rail_id in
  flashier; this avoids counting one tumour twice.

## Gotchas & warnings
- Read lengths are all 48 and FFPE is 0 for every selected sample, so matching is
  effectively random among controls; the seed makes it reproducible.
- `have_junc` marks samples that already have a `.junc` in `diff_splicing/junc_files`
  (same regtools settings, corrected by the 2026-09-25 aliquot fix). Only the rest are
  downloaded; all TCGA BAMs were deleted on 2026-09-25.

## Input/output
- **In:** `paths.yaml` (`tcga_metadata`, `mutation_am_merged`, `unified_results_dir`) and
  `figures/figure2/cache/SF3B1_full_info.tsv` (the set of flashier rail_ids).
- **Out:** `samples/<comparison>.tsv` (role test/control) and `samples/download_list.tsv`.

## File flow
```
tcga_v2_samples ───┐
AM-merged calls ───┼─► select_samples.R ─► samples/<comp>.tsv ─► Snakefile
diff_splicing/junc ┘                    └─► samples/download_list.tsv ─► build_manifest.py
```

> 2026-11-10: `paths.yaml` and `scripts/paths.R` were removed. Paths are now
> `file.path(QROOT, "<folder>/<file>")`, with `QROOT` the quetzal_code root
> (found via the `.quetzal_root` marker); references to paths.yaml keys above are historical.
