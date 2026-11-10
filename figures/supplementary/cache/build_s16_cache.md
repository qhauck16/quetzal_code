# cache/build_s16_cache.py

## Purpose & context
Cache for S16 (YLPM1 intron-coverage meta-plot / heatmap). Runs the data half
of `tcga/genome_wide/flashier/diff_splicing/ylpm1/plot_intron_retention_heatmap.py`
unchanged (imported as module `H`: same sample picks, seed, Coding/UTR filter,
per-gene dedup, bigwig normalisation) for the top N_MAX = 100 genes by
delta-PSI, and writes binned intron/exon ratio traces as a long TSV.

## Decision rationale
Importing the original script rather than copying it guarantees identical
sample selection and normalisation. `figS17/build.R` picks the top N and
orders rows by centre of mass of the top trace, as the original does.

## Gotchas & warnings
- Run with the python3.12 env (pyBigWig), from any directory: it imports
  vendored copies of `plot_intron_retention_heatmap.py` / `plot_top25_intron_panels.py`
  from `s16_lib/` (copied 2026-10-06 from diff_splicing/ylpm1; paths made absolute,
  INFO_FILE repointed to `figure4/cache/YLPM1_full_plus_variants.tsv`).
- Controls: per-cohort counts use largest-remainder allocation (2026-10-26) so they sum to exactly
  200. Independent rounding of 15ths of 200 gave 201.
- 10-process Pool; bigwig reads dominate (~20 s/intron serial).
- Ratios below `H.RATIO_THRESH` set to 0; minus-strand traces reversed (5' left).
- At time of writing (2026-09-24 12:52) the log showed only sample picking
  (`col1 15 ctrl 201 bw 216`, control pool k30 < 1e-4: 9448) and `s16_traces.tsv`
  did not yet exist -- script was still being developed/run.

## User preferences
N_MAX = 100 to match the published S16 (`--n_genes 100`).

## Input/output
In: `H.INFO_FILE`, `H.TCGA_META`, `H.DF_FINAL`, `H.LEAF2_CLASS`, TCGA bigwigs (`H.BW_DIR`).
Out: `cache/s16_traces.tsv` (dpsi_rank, gene, strand, type, deltapsi, com, column top/ctrl, bin, ratio).

## File flow schematic
```
plot_intron_retention_heatmap.py (H) ─┐
df_final / leaf2 class ───────────────┼──► build_s16_cache.py ──► s16_traces.tsv ──► figS17/build.R
TCGA bigwigs, metadata ───────────────┘
```
