# Figure 1b — NMF-factor recovery of LeafCutter-called junctions (FC ≥ 1.5)

Reproduces the three FC ≥ 1.5 recovery-curve panels shown in Figure 1b:
TCGA-KICH vs TCGA-LUAD, TCGA-TGCT vs TCGA-LUAD, TCGA-KICH vs TCGA-TGCT.

**Source:** `tcga/model_eval/quetzal_benchmarking.Rmd` (only the FC = 1.5
chunks and their setup; other FC thresholds are shipped separately with
the supplementary figures).

**Method:**
1. Load LeafCutter DS output for each pairwise TCGA comparison; keep
   junctions with p.adjust ≤ 0.05 and |ΔPSI| > 0.1.
2. Match those "test" junctions 1:1 to a set of "control" junctions
   (Mahalanobis-nearest on median PSI + median RPM, via `MatchIt`).
3. For each test/control junction, check whether any Poisson-NMF
   factor in the per-gene fit has a top loading at that junction AND
   discriminates the two cohorts in the correct direction with a
   between-cohort loading-ratio ≥ 1.5.
4. Sweep across |ΔPSI| cutoffs from 0.5 down to 0.1 and record the
   fraction of test / control junctions recovered.

**Cache:** `recovery_cache/` holds the MatchIt outputs + the
factor-loading summaries for each comparison. Regenerating them takes
~1 hour per comparison (Poisson-NMF file reads dominate). The cache is
pre-populated with the paper's exact numbers -- delete a
`{comp}_*.RDS`/`*.tsv` pair to force a rebuild.

**Inputs (via `paths.yaml`):**
- `gencode_granges` — intron reference (junction → gene mapping).
- `tcga_metadata` — TCGA rail_id ↔ cohort map.
- `benchmark_leafcutter_dir` — three `{comp}_ds_cluster_significance.txt`
  + `{comp}_ds_effect_sizes.txt` pairs.
- `gene_level_poisson_root` — per-gene Poisson NMF fits laid out as
  `{chr}/{gene}/res.RDS`.
- `aggregated_psi_medians`, `aggregated_junction_rpms`,
  `all_junction_rpm_and_psi` — background pool for MatchIt.

**Output:** `out/{KICH_vs_LUAD,KICH_vs_TGCT,TGCT_vs_LUAD}_recovery_fc1.5.pdf`.

**Run:** `Rscript panel_b_recovery_fc1.5.R` (uses the cache).

> 2026-11-10: `paths.yaml` and `scripts/paths.R` were removed. Paths are now
> `file.path(QROOT, "<folder>/<file>")`, with `QROOT` the quetzal_code root
> (found via the `.quetzal_root` marker); references to paths.yaml keys above are historical.
