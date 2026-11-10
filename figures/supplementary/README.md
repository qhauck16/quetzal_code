# Supplementary figures S1-S20

Each figure is `figS<N>/build.R`, sourcing `common.R`, writing `out/figS<N>.pdf`.
Run from this directory:

    Rscript --vanilla figS5/build.R
    ./build_all.sh          # all, numeric order, pass/fail summary, logs in out/logs/

Caches in `cache/` must exist first; see the header of `build_all.sh` for each
builder's working-directory and SLURM requirements. Cache file names (`s7_*`, `s9_*`,
`s12_*`, `s14`, `s16`, `cache/s2/`) do not always match the figure number they feed.

| Fig | Content | Cache / direct inputs |
|---|---|---|
| S1 | Single-gene Quetzal examples (C2, MVK, NHEJ1; MAP2K2, KLHDC4) | showcase-gene config + per-gene RDS/PSI |
| S2 | Recovery vs FC threshold | `figure1/panel_b/recovery_cache` |
| S3 | TCF12 / SRC / SETD2 cancer-gene mis-splicing | `cache/s2/` (`s2/run.sh`) |
| S4 | GSP #3 read-length junctions | `s4_low{48,76}_per_cohort.tsv` (`build_readlength_cache.R`) |
| S5 | 150-GSP x covariate heatmap | `s5_{matrix,cols}.tsv` (`build_s5_heatmap_cache.R`) |
| S6 | GSP #10 vs Alsafadi 2016 | `s6_alsafadi_cells.tsv`, `s6_venn.tsv` (`build_overlap_cache.R`) |
| S7 | SF3B1 p.Y765C CDF | `SF3B1_full_plus_variants` variant table |
| S8 | LeafCutter splice-type summary | `s8_counts.tsv` (`build_leafcutter_counts.R`) |
| S9 | GSP vs LeafCutter overlaps (3x3) | `s7_leafcutter_cells.tsv` (`build_overlap_cache.R`) |
| S10 | SUGP1 on GSP #10 | `s9_*.tsv` (`build_mutation_cache.R`) |
| S11 | NKAP aa<=329 mutants CDF | `NKAP_full_plus_variants` variant table |
| S12 | Non-NKAP mutant co-occurrence | `s12_*.tsv` (`build_mutation_cache.R`) |
| S13 | GSP #13/#30 splice types + PSI FC by loading bin | `fact_{13,30}_top_junctions.tsv`, top-500 PSI matrices |
| S14 | Expression: frameshift vs in-frame across cutoffs | `figure3/cache/panel_e_scaled.tsv`; NKAP 20v60 `df_final.tsv` |
| S15 | NKAP extended sashimi (top 15 events) | `cache/s14/<GENE>.rds` (`build_s14_cache.R`) |
| S16 | YLPM1 splice-acceptor variant (TCGA-CJ-4638, rank 2 on GSP #30): RNA coverage, reads, exons 9-12 with sequence callout | `figS16/bams/` from `figS16/fetch.sh` (GDC BAM slicing); five-panel version with tumour WGS = `figS16/build_full.R` |
| S17 | YLPM1 intron-coverage heatmap | `cache/s16_traces.tsv` (`build_s16_cache.py`) |
| S18 | GSP #30 vs YLPM1 expression | `cache/ylpm1_expression_row.tsv` (`build_ylpm1_expression_row.sh`); YLPM1 variant table; phenotype-axis scores |
| S19 | YLPM1 LOH profile | PanCanAtlas ABSOLUTE segments; phenotype-axis scores |
| S20 | SVM phenotype-axis scores | `{sf3b1,nkap,ylpm1}` phenotype-axis scores |

Supplementary Notes 1 and 2 live in `note1/` and `note2/`.
