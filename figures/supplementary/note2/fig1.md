# note2/fig1.R  (+ note2/build_cache.R) -- Supplementary Note 2, Figure 1 (was Supplementary Fig. S8)

## Purpose & context
Supplementary figure linking GSP #13 to NKAP-related neurodevelopmental syndrome. The two TCGA
tumours carrying the syndrome substitutions show the C* cryptic 3'ss phenotype, each ranking
first of all samples in its cohort:
- p.R330C: rail 886325, TCGA-READ
- p.R330H: rail 857257, TCGA-LAML

Ported 2026-10-08 from `tcga/genome_wide/flashier/diff_splicing/nkap_rare_disease.Rmd`, using the
same event set (below). Quinn asked for this as its own supplementary figure, backing a short
Results paragraph.

## Decision rationale
- **Every cohort sample, not just the median.** The original plotted each carrier against its
  cohort's per-event median only. The port draws every cohort sample (grey) so the carrier's
  extremity is visible, and adds two tests:
  - empirical rank of the carrier's mean PSI within its cohort, P = rank / (n + 1);
  - one-sided paired Wilcoxon of carrier vs cohort-median PSI across events (the original's
    comparison, made formal).
- **Cohort membership from `tcga_v2_samples`, not the flashier sample set.** The R330H tumour was
  removed by the pre-flashier QC (one of the 7.9%), so it has no GSP #13 score; its evidence here
  is direct PSI.
- **`coord_cartesian(xlim)` instead of `xlim()`.** The original's `xlim()` dropped data before
  `stat_ecdf`, which slightly distorted its curves.

## Gotchas & warnings
- **Events:** NKAP 20v60 df_final at dPSI > 0.05, FC > 2.5; alt 3'ss with
  `distance_to_cryptic` in [-8, 20]. That is 469 events, 445 of them with PSI.
- **PSI column matching:** tries start and start-1, as the original did.
- **The empirical P has a floor.** With rank 1 of ~170-180 it is ~0.006, the minimum attainable;
  the Wilcoxon P carries the effect size.
- **Cohort size:** READ has 177 samples in metadata but 167 with PSI for these events.
- **Cohort references are conservative:** other C*-mutant samples are not removed from them.

## User preferences
Raw data + own ggplot; Quinn's theme tokens; carrier in red; gene names italic.

## Input/output
- **In:** `paths.yaml` (`unified_results_dir`, `psi_root`, `tcga_metadata`).
- **Cache out:** `note2/cache/psi_long.tsv` (rail_id, event, psi, cohort) and `note2/cache/carriers.tsv`.
- **Figure out:** `note2/cache/stats.tsv` (read by `figure3/paper_stats.Rmd`) and `out/note2/SN2_Fig1.pdf`.

## File flow
```
nkap 20v60 df_final ──┐
psi_table.tsv (TCGA) ─┼─► note2/build_cache.R ─► psi_long.tsv ─► note2/fig1.R ─► out/note2/SN2_Fig1.pdf
tcga_v2_samples ──────┘                                                  └─► stats.tsv ─► figure3/paper_stats.Rmd
```

> 2026-11-10: `paths.yaml` and `scripts/paths.R` were removed. Paths are now
> `file.path(QROOT, "<folder>/<file>")`, with `QROOT` the quetzal_code root
> (found via the `.quetzal_root` marker); references to paths.yaml keys above are historical.
