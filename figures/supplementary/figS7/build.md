# figS7/build.R -- Supplementary Figure S7

## Purpose & context
SF3B1 p.Y765C mutants: the GSP #10 per-class cumulative-percentage CDF of main
Figure 2d (left half), with the two p.Y765C carriers labelled.
Original: `tcga/genome_wide/flashier/cdf_supplements.Rmd` (chunk "SF3B1 Y765C").

## Decision rationale
Small input table read directly; no cache needed.

## Gotchas & warnings
- The two p.Y765C labels are placed by hand (geom_segment + white-backed
  geom_label, offsets in `lx`/`ly`), not ggrepel: repel put the BRCA label on
  the green/blue curves. BRCA label sits between the steep rise and the legend.
- `stopifnot(length(y765c) == 2)` -- the caption says "the two samples"; the
  build fails loudly if the variant table changes.
- Label placement is hand-tuned above the red curve's 2% bottom run, clear of
  the hotspot rise at the right; recheck if data change.

## User preferences
Labels on the curve, not in a legend; x = log(GSP #10).

## Input/output
In: `$FLASHIER/SF3B1_full_plus_variants.tsv`. Out: `out/figS7.pdf` (5.5 x 4.2).

## File flow schematic
```
flashier/SF3B1_full_plus_variants.tsv ──► figS7/build.R ──► out/figS7.pdf
```
