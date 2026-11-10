# common.R

## Purpose & context
Shared setup sourced by every `figS<N>/build.R` (`source("common.R")`). Loads
the tidyverse/patchwork stack, sources the main-figure design tokens
(`../theme.R`, `../scripts/coverage_plot.R`, `../scripts/pie_chart.R`), defines
path constants (`TCGA`, `FLASHIER`, `REF_SUPP` = published supplementary PDFs,
`OUT_DIR` = `supplementary/out/`) and two helpers:
- `save_supp(p, name, w, h)` -- cairo_pdf to `out/<name>.pdf`, inches.
- `overlap_heatmap(cells, gsp_label, ...)` -- the 2x2 exact-overlap heatmap
  (GSP set / RPM-PSI matched rows x up/down columns) used by S6 and S7.

## Decision rationale
- One source for fonts/palette/`gene_md()` shared with figures 1-4, so the
  supplement cannot drift from the main figures (same rule as `../figure2/compose.md`).
- Always `cairo_pdf`: every text measurement in the project is calibrated to it.

## Gotchas & warnings
- Relative `source("common.R")` + `cache/...` reads in every build.R => figure
  scripts MUST be run with cwd = `supplementary/` (build_all.sh cd's there).
- Re-binds `select`/`rename`/`filter` to dplyr (GenomicRanges/S4Vectors mask them).
- Editing `../theme.R` or `../scripts/*.R` also changes the main figures.

## User preferences
Gene symbols italic via `gene_md()`/`gene_pm()`; x tick labels 45 deg
(`angled_x_labels()`); page sizes in absolute inches

## Input/output
In: `../theme.R`, `../scripts/{coverage_plot,pie_chart}.R`. Out: nothing itself;
`save_supp()` writes `out/figS<N>.pdf`.

## File flow schematic
```
../theme.R ─────────────┐
../scripts/*.R ─────────┼──► common.R ──► figS<N>/build.R ──► save_supp() ──► out/figS<N>.pdf
                        │
FLASHIER / TCGA paths ──┘
```
