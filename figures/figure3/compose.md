# Figure 3 (NKAP / GSP #13) — compose.R

## Purpose & context

Regenerates `Figure3-NKAP_improved.pdf` from raw data instead of
hand-assembly in Illustrator, using the same shared toolkit as figures
1 and 2: `../theme.R` (typography + palette), `../scripts/raster_utils.R`,
`../scripts/coverage_plot.R` (native sashimi/coverage/gene-model geoms),
`../scripts/figure_qa.R` (layout measurement) and
`../scripts/title_align.R` (title-row alignment measurement).

Reference: `tcga/genome_wide/flashier/paper_figures/Figure3-NKAP_improved.pdf`
(letter, 8.5 x 11 in).

## Panels

| Panel | Content | Source |
|---|---|---|
| a | Somatic mutations in NKAP/PRPF8/FAM50A/SLU7 enriched on GSP #13 — CDF + labelled zoom | `flashier/flash_fact_char.Rmd` |
| b | Mutations cluster in specific regions of NKAP and PRPF8 — two domain tracks with stacked lollipops coloured by TCGA cohort | bespoke Illustrator; mutation data in `NKAP_full_plus_variants` |
| c | Four mis-splicing loci (SENP2, PRKDC, EZH2, BRD9): per locus a zoom callout, two sashimi tracks (High GSP #13 / Matched Controls), a gene model with a purple trapezoid | `flashier/pygenome/NKAP` |
| d | Cryptic 3'ss displacement histogram + event-type pie | `diff_splicing/unified_results` |
| e | Expression ECDF by consequence class (frameshift / in-frame vs RPM-matched controls) | `diff_splicing/combined_expression_analysis.Rmd` |
| f | Structure screenshot — raster, embedded AS IS | extracted from the reference PDF |

Panel a is the figure-2 panel-d pattern (CDF + labelled zoom). Panel d
is the figure-2 panel-e pattern (histogram + pie inset). Both should
reuse those panels' solved problems rather than rediscovering them —
read `../figure2/compose.md` first.

## Layout

`ROW_H` is in **inches** and `PAGE_H = sum(ROW_H)`; see
`../figure2/compose.md` for why relative heights are a trap. Four rows:
`(a|b) / c / (d|e) / f`.

Panel c is the tall one and is itself a 2x2 grid of loci, each locus a
stack of: zoom callout, two coverage tracks, gene model.

## Gotchas

- **Panel c needs a new capability**: a zoom callout box above the
  coverage, connected to a boxed sub-region of the track by leader
  lines, annotated with the shift ("+9 nt", "-5 nt"). `coverage_plot.R`
  does not have this yet. The user made this zoom by hand in the
  reference; it has to be derived here.
- **Panel f is a screenshot** — do not attempt to redraw it. Extract
  losslessly from the reference PDF and re-embed at its true aspect
  ratio via `gg_from_pdf(preserve_aspect = TRUE)`, the same way
  figure 2's panel a handles its vector art. Getting `preserve_aspect`
  wrong stretches it.
- **NKAP's factor index must be verified, never assumed.** Memory
  records k13 for an older RDS and explicitly flags it as stale for the
  current one.

## File flow schematic

```
flash_fact_char.Rmd ─────────► cache/panel_a_*.tsv ──┐
NKAP_full_plus_variants ─────► cache/panel_b_*.tsv ──┤
pygenome/NKAP/*.bw + .links ─────────────────────────┼──► compose.R ──► out/figure3.pdf
unified_results/ ────────────► cache/panel_d_*.tsv ──┤
combined_expression_analysis ─► cache/panel_e_*.tsv ─┤
Figure3-NKAP_improved.pdf ───► panel_f raster ───────┘
```

## 2026-10-05 relayout (Quinn): new full-width panel f
- f = `panel_f/Improved_structure.pdf` (overview + orange/blue zooms), on its OWN
  bottom row spanning the page; title "Mutations spatially cluster" re-set in R.
  Rasterised with ghostscript at 600 dpi and ink-cropped (`raster_cache/
  panel_f_improved.png`; rebuilt when the PDF is newer). Its legend's Arial is
  embedded, so it stays Arial. The old reference-crop build is in
  `panel_f/.build.R.pre_2026-10-05`.
- Rows: a|b 2.95, c 3.35, d|e 2.75 (was d|e|f 3.05), f derived from the art's
  aspect at full width + title band (2.12) -> page 11.17 in. Rows 1-2 unchanged.
- `COL_W$def` -> `COL_W$de = c(0.47, 0.53)`.
- Re-measured for the wider cells: d `D_PIE_BOX` left 0.04 (pie was on the y
  axis) and `D_PIE_R_IN` 0.499 -> 0.645 (fixed point, saturation scan);
  e `E_IN_PER_X/Y` 1.0967/1.8400 -> 1.6688/1.6025 (y from tick spacing, x from
  the drawn axis length).

## Zoom connectors (2026-10-05)
The two lines from the CDF zoom box to the zoom panel end at constants in the
CDF's data units (`*_CALL_DX/_TOP/_BOT`). These are MEASURED, not tuned: run
`../scripts/measure_zoom_connectors.py out/<fig>.pdf x0 y0 x1 y1` (crop = page
fractions around the CDF + zoom) after any layout change to that row, and paste
the DX / TOP / BOT it prints. Hand-tuned values had drifted off the corners.
