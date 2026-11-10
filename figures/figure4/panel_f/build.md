# figure4/panel_f/build.R

## Purpose & context

Figure 4 panel f: the GIST validation heatmap. 17 paired tumour/normal GIST
RNA-seq samples x the same 512 junctions the TCGA YLPM1 heatmap (panel-less,
upstream) clusters, as scaled logit(PSI), with a row dendrogram on the right,
three annotation bars, and a colourbar.

The biological point is the row clustering: it splits 11 YLPM1-wildtype from
6 YLPM1-LoF samples cleanly, with the LoF block uniformly shifted positive
(row means -0.50..+0.02 for the 11 NO rows, +0.29..+0.83 for the 6 YES rows,
and the shift holds across all eight 64-column blocks). That is what
"recapitulates the TCGA differential splicing" means here.

Source: `<FLASHIER>/diff_splicing/ylpm1/leafcutter_prep_and_analysis.Rmd`,
chunk `additional GIST heatmaps` (lines ~504-634).

## Decision rationale

**Rebuilt natively rather than insetting the original.** The Rmd draws this
with `heatmaply()`, i.e. a plotly HTML widget; the published
`Figure4-YLPM1_improved.pdf` is a Safari print of the RStudio Viewer. A widget
is not a grob and cannot enter a patchwork, and the printed raster has browser
type in it. So: `geom_raster` for the matrix, `geom_tile` strips for the bars,
`dendextend::as.ggdend()` segments for the dendrogram.

**Everything lives in ONE ggplot and ONE coordinate system.** Heatmap, strips,
dendrogram, leader lines and labels are all placed in the same x/y data units.
This is the whole reason the dendrogram lines up row-for-row: leaf *i* is drawn
at `y = N_ROW + 1 - i`, which is by construction the centre of heatmap row *i*.
Composing three aligned sub-plots with patchwork would have made that alignment
a layout coincidence instead of an identity. Measured at 400 dpi: every leaf
stem is within 0.6 px of its strip-row centre.

**heatmaply's clustering is reproduced, not re-invented.** `heatmapr()` does
`dist()` (euclidean) -> `hclust(method = "complete")` -> `as.dendrogram()` ->
`dendextend::seriate_dendrogram(method = "OLO")`, rows only. Those four calls
are repeated verbatim, so no plotly is needed at render time. The published
17-row order is reproduced exactly (`stopifnot` on the 11-NO/6-YES split, plus
a `warning()` naming the order if it ever drifts).

**plotly is called exactly once, offline.** Only to recover the TCGA heatmap's
own column order, which is cached to `cache/panel_f_col_order.rds`. After that
the figure has no heatmaply/plotly dependency.

**The annotation palette is hardcoded.** heatmaply melts all three bars into a
single factor and spreads one palette over the union of their 18 levels, so the
bars *share* a palette (individual "87" and "NO" come from the same ramp).
Three independent per-bar `scale_fill_*` calls cannot reproduce that, and the
result would not match the published colours. The 18-entry map is literal.

## Gotchas & warnings

- **TCGA clusters come from the current PSI tables** (`new_psi_files`, paths.yaml `psi_root`).
  Until 2026-10-28 this panel read the older Dec-2025 `leaf2_try/psi_files` set by mistake.

- **The PSI denominator is the TCGA LeafCutter cluster, not the GIST one.**
  The Rmd drops GIST's own `cluster_id` and picks the cluster up by joining to
  the TCGA `psi_table.tsv` column headers. This looks incidental and is not:
  grouping by GIST's clusters gives a matrix with the same shape and almost the
  same value range (-3.761..3.881 vs -3.762..3.881) but **permutes rows 2-9 of
  the dendrogram** (87N/73T swap, and the 17T/14T/19T/7T/35T/65T block
  reverses). The 11/6 split survives either way, so this will not announce
  itself — only the exact published row order catches it. That is why
  `cache/panel_f_tcga_cluster_ids.rds` (the psi_table headers, 89,760 rows)
  exists at all.
- **`as.vector()` of the samples x junctions matrix walks DOWN columns**, so
  the sample index varies fastest: `x = rep(1:N_COL, each = N_ROW)` and
  `y = rep(..., times = N_COL)`. Getting this backwards (or slipping a `t()`
  in) still renders a plausible heatmap-looking texture — it is the transpose
  smeared across the panel, and it cost one build iteration to spot by eye.
  Three `stopifnot` spot-checks now pin corner and interior cells.
- **`heatmapr(...)$matrix$data` returns rows REVERSED** (plotly draws
  bottom-up). The displayed top-to-bottom order is plain
  `order.dendrogram(dend)`, *not* its reverse.
- **`guide_colourbar` multiplies `legend.key.height` by 5** for a vertical bar
  (ggplot2 >= 3.5). `unit(12.0 / 5, "mm")` is the drawn height over 5, not a
  typo. At the un-divided value the 2 / 0 / -2 tick labels overlapped.
- **The leader geometry is derived from font metrics**, via `grid::grobWidth`
  on each label at `SIZE_BASE` inside a throwaway cairo device, then converted
  to data units with `F_PANEL_W_IN` / `F_PANEL_H_IN` (measured off a 400 dpi
  render). Re-measure those two if `ROW_H[4]` or `COL_W$fg` change in
  `compose.R`.
- **Labels are LEFT-aligned, so the longest one has the shortest leader.**
  Widening the pad does not help — it shifts the whole block. `F_MIN_RUN`
  sets the block back by a minimum horizontal run that every leader inherits;
  without it "YLPM1 mutant status" gets a bare riser that reads as
  disconnected from its bar.
- **`panel_cached()` serialises the panel**, so no `aes()` may close over a
  local scalar. Leader endpoints are carried as columns of `f_lead_df` /
  `f_lab_df`. Verified: a fresh build and a cached build produce a
  byte-identical 400 dpi crop.
- The main colour ramp spans the **full** data range (-3.76..3.88), matching
  `heatmaply(scale = "none")`. The 2 / 0 / -2 in the legend are tick *breaks*,
  not limits.
- No row labels and no junction x-ticks: the published figure deletes both.
- `compose.R` runs panels concurrently with other agents in this repo; a
  `wrap_dims()` / `is.grob()` error mid-compose usually means another agent
  rewrote `panels/<x>.rds` during the read. Just re-run.

## User preferences

- Quinn's tokens only: `FONT_FAMILY`, `SIZE_BASE`, `SIZE_TITLE`,
  `theme_publication()`. No hand-set point sizes.
- Stretch the plot, don't shrink the type: the panel fills its whole 3.0 in
  cell. The published reference leaves roughly a third of the cell empty
  because it was a browser print at browser type sizes; reproducing that
  whitespace was not worth it.
- Labels sit on the object with leader lines, not in a legend.

## Input/output

Inputs (all read-only, under `<FLASHIER>` =
`/project2/yangili1/qhauck/spliceNMF/tcga/genome_wide/flashier`):

| Path | Shape |
|---|---|
| `y_gene_mut/paired_tumor_normal/ylpm1_perind_numers.counts.gz` | 322,837 x 18, `V1` = `chr:start:end:cluster` |
| `y_gene_mut/YLPM1_mutant_status.txt` | 73 GIST cases |
| `y_gene_mut/accession_to_experiment.txt` | HRA005970 run table |
| `y_gene_mut/paired_tumor_normal/HRS*/` | 17 accessions present |
| `diff_splicing/ylpm1/high_psi_hits.tsv` | 561 x 20 |
| `diff_splicing/ylpm1/tcga_heatmap_mat_to_use.tsv` | 213 x 512 |
| `leafcutter/leaf2_try/new_psi_files/tcga/{chr}/{gene}/psi_table.tsv` (paths.yaml `psi_root`) | headers only |

Caches written under `figure4/cache/` (~450 KB total, all rebuilt on demand):

- `panel_f_col_order.rds` — 512 junction ids, the TCGA column order
- `panel_f_tcga_cluster_ids.rds` — 89,760 TCGA `cluster_id:start-end` headers
- `panel_f_data.rds` — `list(mat = 17 x 512 numeric, pheno = 17-row data.frame)`

## File flow schematic

```
tcga_heatmap_mat_to_use.tsv ──► heatmaply+plotly (ONCE) ──► cache/panel_f_col_order.rds ──┐
                                                                                           │
new_psi_files/tcga/*/*/psi_table.tsv (headers) ──► cache/panel_f_tcga_cluster_ids.rds ──┐      │
                                                                                    ├──────┤
ylpm1_perind_numers.counts.gz ──► PSI per TCGA cluster ──────────────────────────────┤      │
high_psi_hits.tsv ───────────────────────────────────────────────────────────────────┘      │
YLPM1_mutant_status.txt ──┐                                                                 │
accession_to_experiment.txt ├──► phenotype bars ────────────────────────────────────────────┤
paired_tumor_normal/HRS*/ ─┘                                                                │
                                                                                            ▼
                                                                          cache/panel_f_data.rds
                                                                                            │
                                    dist -> hclust(complete) -> OLO seriation ◄──────────────┤
                                                                                            ▼
                                                        panel_f (wrap_elements) ──► compose.R
                                                                                 ──► out/figure4.pdf
```

## Verification

```
cd quetzal_code/figures/figure4
Rscript --vanilla compose.R f
pdftoppm -r 400 -png -singlefile -x 0 -y 3300 -W 1800 -H 1200 out/figure4.pdf /tmp/crop
```

Checked at 400 dpi: dendrogram leaf stems within 0.6 px of their strip-row
centres; f/g titles and tags share baselines; no text overlaps; colourbar
labels clear of each other; fresh and cached builds byte-identical.

> 2026-11-10: `paths.yaml` and `scripts/paths.R` were removed. Paths are now
> `file.path(QROOT, "<folder>/<file>")`, with `QROOT` the quetzal_code root
> (found via the `.quetzal_root` marker); references to paths.yaml keys above are historical.
