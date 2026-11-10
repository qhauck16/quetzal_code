# figS10/build.R -- Supplementary Figure S10

## Purpose & context
Enrichment of SUGP1 mutants on GSP #10: samples ranked 40-100 on GSP #10,
coloured by SF3B1 / SUGP1 status, SUGP1 carriers labelled with protein
change; annotated with the cohort-stratified CMH p_adj for SUGP1 among the
top-100 non-SF3B1 samples (Bonferroni over all genes).
Original: `flashier/flash_fact_char.Rmd` ("SUGP1 along GSP 10").

## Decision rationale
Mutation matrix + CMH computed once in `cache/build_mutation_cache.R`
(AlphaMissense > 0.9 or NA filter, as in the Rmd); here only plotting.

## Gotchas & warnings
- p_adj rendered as plotmath so the exponent is a real superscript.
- Cache log: SUGP1 p_adj = 2.7e-83. Annotation reads "*SUGP1* p_adj = ...".
- Repel direction "both": with "y", rank 94's leader ran through rank 95's label.

## User preferences
Carrier labels on the points (ggrepel), not in a legend.

## Input/output
In: `cache/s9_top100_gsp10.tsv`, `cache/s9_cmh_sf3b1_excluded.tsv`.
Out: `out/figS10.pdf` (6.5 x 4).

## File flow schematic
```
cache/build_mutation_cache.R ──► cache/s9_top100_gsp10.tsv ────────┬──► figS10/build.R ──► out/figS10.pdf
                             └─► cache/s9_cmh_sf3b1_excluded.tsv ──┘
```
