# Figure 1c — BDH1 cohort composition + coverage sashimi

Reproduces the two Figure 1c sub-panels:
- **structure_plot.pdf** — per-sample stacked bar of `k` (Poisson NMF)
  factor usage in each of TCGA-LIHC, SKCM, LGG, COAD, UVM.
- **regions/default/simple_coverage.pdf** + **sashimi_simple_coverage.pdf**
  — BDH1 coverage tracks across the same five cohorts, with `k4`-unique
  arcs drawn on top.

**Source:** copied verbatim from `tcga/gene_level_poisson/showcase_genes/`
(Snakefile, scripts/, envs/, config.yaml). One showcase entry is
active: `BDH1_k4`.

**Pipeline:**
```
rule gene_outputs      -> junctions_classified.tsv, junction_arcs_mane.pdf,
                          structure_plot.pdf, gene_data.RDS
rule prep_region       -> regions/default/{junctions.bed, sample_bigwigs.tsv,
                          groups_file.txt, region.txt, mane_gene_model.pdf}
rule pygenome_tracks   -> regions/default/{simple_coverage.pdf,
                          sashimi_simple_coverage.pdf}
```

**Inputs (via panel_c/config.yaml):**
- Poisson-NMF fit for BDH1: `<tcga_chr_dir>/chr3/BDH1/{res.RDS,de_res.RDS,results.tsv,dat.tsv}`.
- MANE Select summary + hg38 GRanges + TCGA metadata + Leaf2 junction
  classifications.
- Recount3 bigWigs for the five TCGA cohorts under `tcga_bw_dir`.
- `~/software/GenometracksByGenotype/` install for the coverage
  aggregator + bed12 track.

**Output already present:** `output/BDH1_k4/` is copied over so users
can view the exact paper figure files without rebuilding. Delete the
directory and re-run snakemake to regenerate.

**Run:**
```
snakemake --executor slurm --default-resources slurm_partition=broadwl slurm_account=pi-yangili1 --jobs 20 --use-conda
```

**Gotcha:** the `paths:` block in `config.yaml` must match the equivalent
entries in `quetzal_code/paths.yaml`. If you edit paths.yaml for your
own environment, edit this file's `paths:` block to match.

**Gotcha:** `fastTopics::structure_plot` orders samples with an unseeded
t-SNE embedding (`loadings_order = "embed"`, the default) and subsamples
with `n = 5000`. Without a seed the panel reshuffles on every build —
a rerun moved 6034 px at 110 dpi while every other panel in figure 1 was
byte-identical, which makes a before/after figure diff useless for
catching real regressions. `plot_structure_static()` in `scripts/lib.R`
now calls `set.seed(20260922)` immediately before it; two consecutive
`figure1/compose.R` runs are byte-identical. Keep the seed.

> 2026-11-10: `paths.yaml` and `scripts/paths.R` were removed. Paths are now
> `file.path(QROOT, "<folder>/<file>")`, with `QROOT` the quetzal_code root
> (found via the `.quetzal_root` marker); references to paths.yaml keys above are historical.
