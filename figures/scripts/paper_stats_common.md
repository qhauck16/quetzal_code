# paper_stats_common.R  (+ figure{1,2,3,4}/paper_stats.Rmd)

## Purpose & context
Every number stated in the manuscript text is recomputed by a `paper_stats.Rmd` in the folder
of the figure it relates to most closely. Claims with no figure of their own are grouped with the
nearest figure: gene-level/cancer-type numbers go to figure1, GSP #10 to figure2, GSP #13 to
figure3, and GSP #30 plus the cross-GSP stringent counts to figure4.

This helper holds the paths, the input switch and `claim()`. `claim()` records each number next
to the value printed in the text, and each Rmd ends by writing `paper_stats_<inputs>.tsv` with a
CHANGED flag. Built 2026-10-07 after the diff-splicing, mutation-recall and chr10/11/18 PSI reruns.

## Decision rationale
- **Ported, not reinvented.** Definitions come from the analysis Rmds that produced the text:
  `diff_splicing/paper_stats.Rmd`, `{sf3b1,nkap,ylpm1}/top500_alt3_fc_vs_loading.Rmd`,
  `unified_results/ylpm1_stats.Rmd` and `plotting_genes/gene_cancer_type_hits.Rmd`.
- **Two input modes.** `params$inputs = "pre_rerun"` points at the backed-up inputs:
  - `unified_results/.pre_fix_2026-10-05`;
  - the Rmd-era phenotype scores in `diff_splicing/<gene>/`;
  - for figure1, today's per-gene tree with the 2,114 rerun genes swapped for
    `gene_level_poisson/pre_chr18fix_2026-10-06/`.

  In that mode each Rmd must reproduce the text. That proves the port is faithful, so any
  difference in "current" mode is a real data change.

## Gotchas & warnings
- **Re-calls inside `gene_mutation`.** `gene_mutation` in the variant tables already contains the
  WXS / RNA-seq re-calls. "Annotated" (PanCanAtlas/GDC) means `gene_mutation` minus the entries
  that also appear in `wxs_mutation` / `rnaseq_mutation`. Counting all three columns double-counts.
- **Phenotype-positive** = P or intermediate with SVM score > 0. P samples scored <= 0 are excluded.
- **Mutation variant tables** have no pre-rerun copy (the old ones were overwritten 2026-10-05).
  Mutation-based claims run on current calls in both modes.
- **Unmatched definition: "851 disrupted introns".** No definition reproduces it. Every test intron
  has CoLa-seq branchpoint data (1,863 pre-rerun, 1,822 now).
- **figure1's re-walk is slow.** It walks ~13.7k gene dirs and caches to
  `figure1/cache/paper_stats_glp_<inputs>.rds`; delete the cache to refresh.

## User preferences
Quinn wants every paper number generated from code inside quetzal_code, grouped by figure.

## Input/output
- **In:** `paths.yaml` (`unified_results_dir`, `cosmic_tier1`, `mutation_am_merged`,
  `gene_level_poisson_root`, `tcga_metadata`), figure caches, and variant tables via `vt_path()`.
- **Out:** `figureN/paper_stats_{current,pre_rerun}.{html,tsv}`.

## File flow
```
paths.yaml + figure caches + vt_path tables ─┐
unified_results (+ .pre_fix backup) ─────────┼──► figureN/paper_stats.Rmd ──► paper_stats_<inputs>.tsv / .html
gene_level_poisson per-gene outputs ─────────┘        (sources scripts/paper_stats_common.R)
```

> 2026-11-10: `paths.yaml` and `scripts/paths.R` were removed. Paths are now
> `file.path(QROOT, "<folder>/<file>")`, with `QROOT` the quetzal_code root
> (found via the `.quetzal_root` marker); references to paths.yaml keys above are historical.
