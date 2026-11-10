# figS15/build.R

## Purpose & context
Supplementary Fig. S15: sashimi-style coverage plots for the top 15 NKAP-GSP
cryptic alt-3'ss events (rank order DUS1L ... ELAC2), 5 x 3 grid on an
8.5 x 11 in page, no panel letters. Native-ggplot redraw of
`tcga/genome_wide/flashier/NKAP_detailed/custom_plot.R` + `assemble_metaplot.R`
(published: `paper_figures/all_supplementary_figures/nkap_extended_sashimi_plots_supp.pdf`).

## Decision rationale
- Data come from `cache/s14/*.rds` (built once by `cache/build_s14_cache.R`,
  which needs bigwigs/GRanges); this script only draws, so it runs in seconds.
- Each gene = patchwork stack (cond track / ctrl track / zoom trapezoid / gene
  model); 15 stacks go through `wrap_plots(ncol = 3)`. All vector.
- Arcs: sin-shaped rise from acceptor coverage toward the off-panel donor,
  linewidth `max(0.2, sqrt(PSI)/11)`. Unlike the original, the two arcs exit
  the panel edge at staggered heights (ARC_END, ARC_END - ARC_GAP) instead of
  converging on one point, so each PSI label can sit clearly on its own arc.
- PSI labels placed by `place_labels()`: grid search over arc fractions
  (0.15-0.85, nearest 0.5 first) rejecting boxes that hit each other, the
  coverage bars, the panel top, or (strict pass) the other arc.
- Window highlight on the gene model widened to >= 1.5% of the gene span so it
  stays visible on long genes (NUMB, CTNNAL1).
- Published version's second subtitle line ("alt-3ss @ chr:pos") dropped;
  title is "GENE | -8 nt | CDS" only.

## Gotchas
- `TRACK_W_IN`/`TRACK_H_IN` (label size -> data units) were measured off the
  render for 8.5 x 11; change the page/HEIGHTS and re-measure.
- A "relaxed label placement: GENE grp" message = strict pass failed for that
  track (label may sit over the other arc). ARC_GAP was raised 0.30 -> 0.40
  (2026-09 polish) so the strict pass now succeeds on all 30 tracks; at 0.30
  seven cond tracks fell back to relaxed. TRACK_W_IN/TRACK_H_IN re-measured
  then (2.28 x 0.50 in) and still match.
- PSI "0.0" (NADK/AAAS controls) is real (~0.02%), not missing.
- y max tick = ymax rounded to 2 sig figs (>= data max), shared by both tracks.

## User preferences
Italic gene names (gene_md, fontface italic), plain-weight short y titles,
two-tick capped y axis, white-backed inline arc labels "%.1f", angled x ticks,
FONT_FAMILY/SIZE_* tokens, no rasterizing, grow page rather than shrink type.

## Input / output
In: `cache/s14/<GENE>.rds` (list: gene, chr, strand, shift, region, acceptors,
donor, win, cov[pos,cov,group], psi[junction,acceptor,cond,ctrl], exons, codons).
Out: `out/figS15.pdf` (cairo_pdf).

## File flow
```
NKAP_detailed/{events_meta.tsv, samples} ─┐
recount3 bigwigs, leafcutter PSI, MANE ───┼─► cache/build_s14_cache.R ─► cache/s14/<GENE>.rds
                                          │                                   │
common.R + ../theme.R ────────────────────┴──────────────► figS15/build.R ◄──┘
                                                               └─► out/figS15.pdf
```

**2026-09-24 (Quinn review):** start/stop-codon star key under the grid; group y titles shown on the left column only (other columns transparent, so track widths — and TRACK_W_IN — are unchanged).
