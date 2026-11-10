# Supplementary figures S1-S20

Each figure is `figS<N>/build.R` (companion `build.md`), sourcing `common.R`, writing
`out/figS<N>.pdf`. Run from this directory:

    /scratch/midway2/qhauck/conda_env/rstudio-server/bin/Rscript --vanilla figS5/build.R
    ./build_all.sh          # all, numeric order, pass/fail summary, logs in out/logs/

Caches in `cache/` must exist first; see the header of `build_all.sh` for each
builder's working-directory and SLURM requirements, and `cache/*.md`.
Paths below: `FL` = `tcga/genome_wide/flashier`, `DS` = `FL/diff_splicing`.

| Fig | Content | Cache / direct inputs | Original source (tcga/ tree) | Status (2026-09-24) |
|---|---|---|---|---|
| S1 | Single-gene Quetzal examples (C2, MVK, NHEJ1; MAP2K2, KLHDC4) | none; `showcase_genes/supplementary_plots/config.yaml` + per-gene RDS/PSI | `gene_level_poisson/showcase_genes/supplementary_plots` (Snakefile + scripts) | built; extra PRAD/KIRP 50-bp groups shown (present in data) (was S1) |
| S2 | Recovery vs FC threshold | `figure1/panel_b/recovery_cache` | `model_eval/quetzal_benchmarking.Rmd` | built; KICH vs LUAD matched FC>=1 first two points differ slightly (was S3) |
| S3 | TCF12 / SRC / SETD2 cancer-gene mis-splicing | `cache/s2/` (`s2/run.sh`) | `figure1/panel_d_e/scripts/make_gene_outputs.R` (fig 1d code) | built (was S2) |
| S4 | GSP #3 read-length junctions | `s4_low{48,76}_per_cohort.tsv` (`build_readlength_cache.R`) | `FL/just_for_fun/nmd_investigation/covariate_gsp_figures.Rmd` | built (was S4) |
| S5 | 150-GSP x covariate heatmap | `s5_{matrix,cols}.tsv` (`build_s5_heatmap_cache.R`) | `FL/flash_fact_char.Rmd` (via fig 2c cache script) | built; +/-0.4 scale + fig 2c colours (published +/-0.6), alphabetical tissues, GSP-order rows (was S5) |
| S6 | GSP #10 vs Alsafadi 2016 | `s6_alsafadi_cells.tsv`, `s6_venn.tsv` (`build_overlap_cache.R`) | `FL/spec_flash_fact.Rmd` | built (was S6) |
| S7 | SF3B1 p.Y765C CDF | direct: `FL/SF3B1_full_plus_variants.tsv` | `FL/cdf_supplements.Rmd` | built (was S10) |
| S8 | LeafCutter splice-type summary | `s8_counts.tsv` (`build_leafcutter_counts.R`) | `DS/unified_results/leafcutter_meta_analysis.Rmd` | built (was S8) |
| S9 | GSP vs LeafCutter overlaps (3x3) | `s7_leafcutter_cells.tsv` (`build_overlap_cache.R`) | `FL/spec_flash_fact.Rmd` | built; counts/p match published (was S7) |
| S10 | SUGP1 on GSP #10 | `s9_*.tsv` (`build_mutation_cache.R`) | `FL/flash_fact_char.Rmd` | built (was S9) |
| S11 | NKAP aa<=329 mutants CDF | direct: `FL/NKAP_full_plus_variants.tsv` | `FL/cdf_supplements.Rmd` | built (was S11) |
| S12 | Non-NKAP mutant co-occurrence | `s12_*.tsv` (`build_mutation_cache.R`) | `FL/flash_fact_char.Rmd` | built (was S12) |
| S13 | GSP #13/#30 splice types + PSI FC by loading bin | direct: `FL/fact_{13,30}_top_junctions.tsv`, `DS/{nkap,ylpm1}/*top500*_psi.tsv` | `FL/spec_flash_fact.Rmd`; `DS/{nkap,ylpm1}/top500_alt3_fc_vs_loading.Rmd` | built (was S13) |
| S14 | Expression: frameshift vs in-frame across cutoffs | `figure3/cache/panel_e_scaled.tsv`; `DS/unified_results/nkap/20v60/dpsi0.01_fc2.5/df_final.tsv` | `DS/combined_expression_analysis.Rmd` | built; FC>1 row = FC>2.5 row by construction (was S15) |
| S15 | NKAP extended sashimi (top 15 events) | `cache/s14/<GENE>.rds` (`build_s14_cache.R`) | `FL/NKAP_detailed/custom_plot.R` + `assemble_metaplot.R` | done (was S14) |
| S16 | YLPM1 splice-acceptor variant (TCGA-CJ-4638, rank 2 on GSP #30): RNA coverage + equal-height arcs, squished reads, exons 9-12 with sequence callout (MTIF2 slide style) | `figS16/bams/` from `figS16/fetch.sh` (GDC BAM slicing) | new 2026-10-26; five-panel version with tumour WGS = `figS16/build_full.R` | built (new) |
| S17 | YLPM1 intron-coverage heatmap | `cache/s16_traces.tsv` (`build_s16_cache.py`) | `DS/ylpm1/plot_intron_retention_heatmap.py` | done (was S16) |
| S18 | GSP #30 vs YLPM1 expression | `cache/ylpm1_expression_row.tsv` (`build_ylpm1_expression_row.sh`); `FL/YLPM1_full_plus_variants.tsv`; `DS/ylpm1/ylpm1_phenotype_axis_scores.tsv` | `FL/ylpm1_lof_investigation.Rmd`; `DS/ylpm1/top500_alt3_fc_vs_loading.Rmd` | built; Phen- n=105 vs published 106; X1783_splice counted LoF (was S17) |
| S20 | SVM phenotype-axis scores | direct: `DS/{sf3b1,nkap,ylpm1}/*_phenotype_axis_scores.tsv` | `DS/svm_plots.Rmd` | built (was S18) |

"built" = `out/figS<N>.pdf` exists as of 2026-09-24 12:52; published references are
in `FL/paper_figures/all_supplementary_figures` (`REF_SUPP` in common.R; OLD numbering).

**Renumbered 2026-10-26** to match the order of citation in the manuscript draft. "(was S<n>)" gives the old number;
cache file names (`s9_*`, `s12_*`, `s14`, `s16`, `cache/s2/`) keep their OLD numbers. Order follows the draft (2026-10-26): S8 = LeafCutter splice-type summary (cited at the 30v90 GSP #10
benchmark), S9 = GSP vs LeafCutter overlaps, S19 = YLPM1 LOH profile (2026-11-09), S20 = SVM phenotype axis (cited from the Methods, so last).
