# figS2/build.R -- Supplementary Figure S2

## Purpose & context
Quetzal's recovery of LeafCutter differentially-spliced junctions (as main
Figure 1b) over |dPSI| cutoffs, repeated for fold-change requirements
FC >= 1, 1.05, 1.1, 1.25, 1.5, 2 (rows) x KICH vs LUAD / KICH vs TGCT /
TGCT vs LUAD (columns); LeafCutter DS vs matched-control junctions.
Original: `tcga/model_eval/quetzal_benchmarking.Rmd` via `figure1/compose.R build_panel_b()`.

## Decision rationale
- Reuses figure 1 panel b's `recovery_cache` (`*_fls_test.RDS`, `*_fls_matched.RDS`,
  `*_matched_{test,background}.tsv`); only the FC loop (`run_recovery_fc`) is new.
- FC >= 1 is the plain direction test (success = correct sign of change).

## Gotchas & warnings
- FC row strips are horizontal (`strip.text.y angle = 0`).
- **Known deviation:** at matched FC >= 1 the first two points of KICH vs LUAD
  diverge slightly from the published panel.
- Gene -> chromosome via gencode GRanges from `../../paths.yaml`; genes without
  a match are dropped.

## User preferences
PAL$leafcutter / PAL$matched colours shared with fig 1b; percent y axis.

## Input/output
In: `figure1/panel_b/recovery_cache/*`, gencode granges. Out: `out/figS2.pdf` (6.5 x 8).

## File flow schematic
```
figure1/panel_b/recovery_cache/*_fls_{test,matched}.RDS ─┐
                     *_matched_{test,background}.tsv ────┼──► figS2/build.R ──► out/figS2.pdf
paths.yaml ─► gencode_granges ───────────────────────────┘
```

**2026-09-24 (Quinn review):** major gridlines back on (grey90) — fine to omit in main figures, odd without them in this supplement.

> 2026-11-10: `paths.yaml` and `scripts/paths.R` were removed. Paths are now
> `file.path(QROOT, "<folder>/<file>")`, with `QROOT` the quetzal_code root
> (found via the `.quetzal_root` marker); references to paths.yaml keys above are historical.
