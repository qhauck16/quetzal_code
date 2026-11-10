# figS17/build.R

**Purpose & context** — Supplementary Figure S17: intron coverage over the
top 100 YLPM1-axis (GSP #30) disrupted introns, top-15 GSP #30 samples vs
project-matched controls. Port of
`tcga/genome_wide/flashier/diff_splicing/ylpm1/plot_intron_retention_heatmap.py --n_genes 100`
(published PDF: `paper_figures/all_supplementary_figures/ylpm1_intron_meta_supp.pdf`).

**Decision rationale** — Data half is the original Python run unchanged
(cache/build_s16_cache.py imports its functions) so sample picks, seed, filters
and normalisation are identical; only the drawing moved to ggplot for the
shared fonts. 100 rows = row count measured on the published PDF. Row order
(centre of mass of the top-15 trace, 3'-most on top) reproduces the published
order exactly (SLC25A15, TTLL7, FAM131C, ...).

**Gotchas & warnings**
- Minus-strand traces are reversed in the cache, so every row reads 5' → 3'.
- The matched-control count comes from the cache log (`n_ctrl 201`): per-project
  rounding gives 201, the published title said 200. The title uses the real n.
- Every layer's `column` must be a factor with the same levels, or facet_wrap
  sorts panels alphabetically and controls appear first.

**User preferences** — italic gene names, FONT_FAMILY/SIZE_* tokens, vector output.

**Input/output** — `cache/s16_traces.tsv` (dpsi_rank, gene, strand, type,
deltapsi, com, column top|ctrl, bin 1..200, ratio) + `cache/build_s16_cache.log`
→ `out/figS17.pdf` (8.5 × 11 in; page stretched from 10 in so the 100 gene
labels can be set at the 6 pt floor, LABEL_PT -- they were 4.6 pt at 10 in).

```
YLPM1_full_plus_variants.tsv ─┐
df_final (ylpm1 20v60) ───────┼─► cache/build_s16_cache.py ─► s16_traces.tsv ─► figS17/build.R ─► out/figS17.pdf
TCGA bigwigs ─────────────────┘   (python3.12 env, 10 procs)
```

**2026-09-24 (Quinn review):** no row baseline — it made enrichment look weak. Zeroed-bin runs are NA so the ribbon breaks; bins touching coverage kept so islands taper.
