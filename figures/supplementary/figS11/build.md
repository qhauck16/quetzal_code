# figS11/build.R -- Supplementary Figure S11

## Purpose & context
Other NKAP mutants: ECDF of GSP #13 (as main Figure 3a), any NKAP variant vs
every other sample, labelling every NKAP mutant affected at amino acid <= 329.
Original: `tcga/genome_wide/flashier/cdf_supplements.Rmd` (chunk "NKAP aa <= 329").

## Decision rationale
Mutation parsing copied verbatim from that chunk: it reads all three mutation
columns (TCGA call, WXS re-call, RNA-seq re-call).

## Gotchas & warnings
- Mutants with log GSP #13 < -9 are labelled in a hand-built column at
  x = -9.3 (y evenly spaced 3-44%, same order as their points so leaders do
  not cross); ggrepel stacked them on both curves. Only the 3 high-loading
  mutants still use ggrepel (separate layer, disjoint region).
- Depends on the current `NKAP_full_plus_variants.tsv`; the WXS re-call (see
  project memory on the SLU7-exclusion re-call) can change which mutants appear.

## User preferences
Variant colour #8E44AD vs grey70; labels on the chart.

## Input/output
In: `$FLASHIER/NKAP_full_plus_variants.tsv`. Out: `out/figS11.pdf` (5.5 x 4.2).

## File flow schematic
```
flashier/NKAP_full_plus_variants.tsv ──► figS11/build.R ──► out/figS11.pdf
```

**2026-09-24 (Quinn review):** samples that also carry an NKAP R402 change are not labelled (R402Q+L313M sat high because of R402). L313V (separate sample, no R402) stays labelled.
