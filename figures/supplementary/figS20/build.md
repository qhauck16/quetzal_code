# figS19/build.R -- Supplementary Figure S19

## Purpose & context
SVM method for estimating phenotype prevalence: SVM phenotype-axis score vs
log GSP loading for GSP #10 (SF3B1), #13 (NKAP), #30 (YLPM1). Labelled
positives/negatives trained the boundary (dotted blue line, score = 0);
intermediates were scored by it. Original: `diff_splicing/svm_plots.Rmd`.

## Decision rationale
Reads the saved `*_phenotype_axis_scores.tsv` directly; no SVM refit, no cache.

## Gotchas & warnings
- Non-finite log loadings (loading 0) are dropped.
- Positives are drawn last (on top) via `arrange(group)` on a factor N < intermediate < P.

## User preferences
Positives #D62728, intermediates #F0A400, negatives grey65.

## Input/output
In: `$FLASHIER/diff_splicing/{sf3b1,nkap,ylpm1}/<gene>_phenotype_axis_scores.tsv`
(rail_id, k<fact>, group N/intermediate/P, score). Out: `out/figS19.pdf` (8.5 x 3.3).

## File flow schematic
```
diff_splicing/sf3b1/sf3b1_phenotype_axis_scores.tsv ─┐
diff_splicing/nkap/nkap_phenotype_axis_scores.tsv ───┼──► figS19/build.R ──► out/figS19.pdf
diff_splicing/ylpm1/ylpm1_phenotype_axis_scores.tsv ─┘
```

## Polish notes (2026-09)
- y title only on the left panel (scales stay per-panel); legend pulled up
  to the axes (legend.box.spacing 4 pt) -- the collected bottom legend
  otherwise left a blank band under the x titles.

## Panel b: dropped 2026-10-26
Quinn: S19 is cited only from the Methods, so the near-boundary C* panel was removed. The test (8 of 50 vs
216 of 9,696, one-sided Fisher P = 1.6e-5) is in `figure3/paper_stats.Rmd`.
