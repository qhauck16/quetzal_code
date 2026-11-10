# run_ylpm1_6v18.sbatch

**Purpose & context** — One-off LeafCutter run (`-m 5` clustering + leafcutter_ds) for `results/ylpm1_6v18`, used ONLY to pick junctions for the YLPM1 / GSP #30 SVM phenotype axis (2026-11-09). Every other YLPM1 analysis still uses 10v30/20v60/30v90.

**Decision rationale** — 10v30 cases include 3 intermediates (scaled GSP #30 < 0.1: 865440, 866724, 434773) that the SVM later scores, so the junctions it picks are leaked from samples it classifies. Only 7 samples are SVM "P" (scaled GSP #30 ≥ 0.1); 461535 (rank 4) has no BAM/.junc (absent from the April GDC manifest; GDC has it now, file 9769a8a6-…, but the token had expired). Running 7v21 made the selector backfill 865440, which brings the leakage back, so Quinn chose 6v18 = the 6 P samples with junctions (882990, 869984, 441042, 870567, 939375, 440985) + 18 matched controls.

**Gotchas** — Generated from `rerun_m5.sbatch` by sed (array removed, `R=` fixed); it reuses that script's conda python + leafcutter_ds env. Selection was made by running `select_and_match_bam.R --n_test 6 --control_ratio 3` directly, not via Snakemake, so the Snakefile's `config.yaml` comparisons don't list 6v18. The run isn't in `unified_results`.

**User preferences** — No leakage between junction selection and the samples the SVM scores.

**Input/output** — In: `results/ylpm1_6v18/{juncfiles.txt,leaf_pheno_data.txt}`. Out: `leafcutter_perind_numers.counts.gz`, `leafcutter_ds_{cluster_significance,effect_sizes}.txt`, read by `quetzal_code/figures/scripts/phenotype_axis.R::pa_build("ylpm1", "6v18", ...)`.

```
select_and_match_bam.R ──► results/ylpm1_6v18/{juncfiles,leaf_pheno_data}
junc_files/*.junc ──┐
                    ├──► run_ylpm1_6v18.sbatch ──► leafcutter_ds_*.txt ──► phenotype_axis.R (SVM)
```
