# compose.R — Figure 1 assembly

Single-entry pipeline that regenerates every R-native panel of Figure 1
with a shared theme and lays it out as one PDF. Reads its knobs from
`theme.R`, cached inputs from each `panel_*/output/*/gene_data.RDS` (so
no SLURM rerun is needed), and turns external PDFs (illustrator
schematic + pygenome sashimi + NF1 transcript model) into aspect-ratio-
preserving raster grobs via `pdftocairo` and png::readPNG with auto-crop.

## Regenerate

```
Rscript compose.R
```

Cold: ~30 s (pdftocairo renders + tsne on the structure_plot). Warm:
<10 s (raster cache hits, gene_data.RDS reused).

Writes `out/figure1.pdf`.

## What tweaks propagate

Every knob in `theme.R` re-flows through the ggplots on the next compose
run:
- `FONT_FAMILY` — swap font (any name cairo can find)
- `SIZE_*`      — text sizes (base / axis / strip / title / legend / panel label)
- `PAL$*`       — palette (leafcutter/matched, TCGA/GTEx accents, factor arcs)
- `GAP_GEOM`    — `"chevron"` (default here — cleaner) or `"arrow"` (the
  original noisy density). Threaded via
  `options(figure1.gap_geom = GAP_GEOM)` into the ggbio
  `geom_alignment()` calls in every gene-track panel.

Sub-panel labels a–e come from `plot_annotation_publication()` in
theme.R. Each panel is wrapped in `wrap_elements()` before joining so
patchwork tags only the five top-level panels (not every sub-plot).

## Adding a new panel or tweaking layout

Edit only compose.R. The four panel bodies (1b, 1c, 1d, 1e) each build
their sub-plots, apply `theme_publication()` / `theme_track()`, and
join with patchwork. Widths / heights of the outer `plot_layout()` at
the bottom of the file control page proportions.

## External PDF handling

`gg_from_pdf(path)` rasterizes at 300 dpi to `raster_cache/`, autocrops
whitespace, and returns a ggplot with `coord_equal()` so patchwork
sizes the outer cell without distorting the image. Cache is
path-derived (basename collisions on `sashimi_simple_coverage.pdf`
across panels would otherwise overwrite each other -- guarded).

To force re-rasterize a specific PDF: `rm raster_cache/<full_path>.png`.

## Cache invalidation

- Change theme.R → recompose. (No cache to invalidate.)
- Change an external PDF (sashimi rebuild, schematic update) → delete
  its `raster_cache/*.png` and recompose.
- Change a gene_data.RDS (rare) → recompose.

## Gotchas

- **Relative paths in cached gene_data.RDS**: the plotting_genes
  pipeline stored per-file paths relative to its own dir. `gd_repath()`
  rewrites them to the absolute paths in `../../paths.yaml` at readRDS
  time so `plot_loading_violins` / `plot_psi_cdfs_full_static` can
  still find TCGA/GTEx metadata + PSI tables.
- **Palette injection**: `LC$palet` / `LDE$palet` are assigned into the
  sourced env of each lib.R BEFORE calling their plot functions --
  this overrides their lib-internal defaults with `PAL$factor`.
- **arrow.rate warning**: benign ggbio warning when we pass
  `gap.geom = "chevron"`; the chevron geom ignores those args.

> 2026-11-10: `paths.yaml` and `scripts/paths.R` were removed. Paths are now
> `file.path(QROOT, "<folder>/<file>")`, with `QROOT` the quetzal_code root
> (found via the `.quetzal_root` marker); references to paths.yaml keys above are historical.
