# figS14/build.R -- Supplementary Figure S14

## Purpose & context
Gene-expression effect of mis-splicing in genes driving GSP #13 across cutoffs:
for every (dPSI > 0.01/0.02/0.05/0.10) x (FC > 1/2.5/5) cell, CDFs of mean
scaled RPM in the top 20 GSP #13 samples for genes with a frameshifting vs
in-frame alt 3'ss (cryptic site -20..0 nt from the annotated 3'ss), one-sided
Wilcoxon (frameshifted < in-frame). Same as main Figure 3e, extended.
Original: `diff_splicing/combined_expression_analysis.Rmd`,
`plot_faceted_cdf("NKAP", nkap_df, "frameshift", dist_window = c(-20, 0), require_alt3 = TRUE)`.

## Decision rationale
Expression comes from figure 3 panel e's cache (`figure3/cache/panel_e_scaled.tsv`,
the 13.4 GB matrix pre-sliced to the same 20 samples), so the full matrix is never read.

## Gotchas & warnings
- The LeafCutter table is the `20v60/dpsi0.01_fc2.5` instance, so every junction
  already has FC > 2.5: the FC > 1 row necessarily equals the FC > 2.5 row --
  exactly as in the original figure. Not a bug.
- Genes in both classes are counted as frameshifted only (in-frame = setdiff).
- `fread(..., header = TRUE)` is required: autodetect treats the header as data.

## User preferences
Frameshifted #B22222 vs in-frame #4682B4 (figure 3 panel e E_PAL, so S15 reads as "same as 3e").

## Input/output
In: `figure3/cache/panel_e_scaled.tsv` (gene_name x top-20 samples),
`diff_splicing/unified_results/nkap/20v60/dpsi0.01_fc2.5/df_final.tsv`.
Out: `out/figS14.pdf` (8.5 x 5.8).

## File flow schematic
```
figure3/cache/panel_e_scaled.tsv ─────────────────┬──► figS14/build.R ──► out/figS14.pdf
unified_results/nkap/20v60/dpsi0.01_fc2.5/df_final.tsv ─┘
```
