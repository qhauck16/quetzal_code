# figS12/build.R -- Supplementary Figure S12

## Purpose & context
Mutations enriched among non-NKAP-mutant samples on GSP #13: pairwise
co-occurrence heatmap of AlphaMissense-filtered mutations among the top 50
non-NKAP samples, for genes passing cohort-stratified CMH (Bonferroni
p_adj < 0.05). Diagonal = samples mutated in that gene (outlined); off-diagonal
= mutated in both. Original: `flashier/flash_fact_char.Rmd`
("NKAP sig-gene pairwise co-occurrence heatmap").

## Decision rationale
Built from `cache/build_mutation_cache.R` outputs; genes already in p_adj order there.

## Gotchas & warnings
- p_adj shown as "*p*<sub>adj</sub> = 2.2 x 10^-15" (markdown superscript, 2 s.f.).
- Legend is inside the empty lower-left triangle (horizontal bar); page
  4.6 x 3.75 in to drop the old side column / blank bands.
Cache log significant genes: FAM50A, PRPF8, SLU7, KAT5, SBNO2.

## User preferences
Counts printed on the tiles.

## Input/output
In: `cache/s12_cooccurrence.tsv` (row, col, count), `cache/s12_cmh_nkap_excluded.tsv`.
Out: `out/figS12.pdf` (4.6 x 3.75).

## File flow schematic
```
cache/build_mutation_cache.R ──► cache/s12_cooccurrence.tsv ──────┬──► figS12/build.R ──► out/figS12.pdf
                             └─► cache/s12_cmh_nkap_excluded.tsv ─┘
```
