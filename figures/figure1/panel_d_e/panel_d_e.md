# Figures 1d + 1e — EZH2 (LAML) and NF1 (LUSC) showcase panels

Both panels come out of the same pipeline (copied verbatim from
`tcga/gene_level_poisson/plotting_genes/`). It is narrowed here to just
the two samples the figure needs — `EZH2_LAML` and `NF1_LUSC` — via a
two-row `gene_level_poisson_sig_results.tsv`.

## Figure 1d — EZH2 in TCGA-LAML

Files in `output/EZH2_LAML/EZH2/`:
- `loading_violin.pdf` — EZH2 k6 loadings, TCGA-LAML vs the pooled Other-TCGA.
- `psi_cdfs_tcga_vs_gtex.pdf` — five per-junction PSI empirical CDFs
  (TCGA-LAML vs GTEx-BLOOD) — the five k6-unique junctions.
- `junction_arcs_mane.pdf` — gene track with the 5 arcs.
- Coverage + sashimi in `regions/default/`.

## Figure 1e — NF1 in TCGA-LUSC

Files in `output/NF1_LUSC/NF1/`:
- `loading_violin.pdf` — NF1 k1 loading violin, TCGA-LUSC vs Other TCGA.
- `regions/narrow/simple_coverage.pdf` + `sashimi_simple_coverage.pdf`
  — the tight NF1 sashimi (chr17:31,248,000-31,259,000) matched to
  GTEx-LUNG.
- `regions/default/...` — wider-view versions.

`NF1_transcript_models/plot_nf1_transcripts.R` renders the isoform track
overlaid on the sashimi (`NF1_MANE.pdf`, `NF1-216.pdf`,
`NF1_MANE_and_216.pdf`).

## Pipeline

Same three-rule pattern as panel_c:
1. `gene_outputs`  — build the violin, CDFs, arcs, `gene_data.RDS`.
2. `prep_region`   — resolve the coverage window + write bed/samples/
   groups/mane files.
3. `pygenome_tracks` — call the GenometracksByGenotype aggregator to
   render coverage + sashimi PDFs.

## Inputs (via panel_d_e/config.yaml)

- Poisson-NMF fits: `<tcga_chr_dir>/{chr}/{gene}/{res.RDS, de_res.RDS,
  results.tsv, dat.tsv}`.
- PSI tables under `psi_dir` (LeafCutter output).
- TCGA + GTEx bigWigs under `tcga_bw_dir` / `gtex_bw_dir`.
- `gene_level_poisson_sig_results.tsv` in this directory has been
  pre-filtered to the EZH2_LAML + NF1_LUSC rows.

## Output already present

`output/EZH2_LAML/` and `output/NF1_LUSC/` are copied over from the
original pipeline. Delete them and re-run snakemake to regenerate.

`NF1_transcript_models/{NF1_MANE.pdf, NF1-216.pdf, NF1_MANE_and_216.pdf}`
are the transcript overlay renders; re-run `Rscript
NF1_transcript_models/plot_nf1_transcripts.R` to regenerate.

## Run

```
snakemake --executor slurm --default-resources slurm_partition=broadwl slurm_account=pi-yangili1 --jobs 10 --use-conda
```

## Gotcha

The `paths:` block in `config.yaml` must match the equivalent entries in
`quetzal_code/paths.yaml`. If you edit paths.yaml for your own
environment, edit this file's `paths:` block to match.

> 2026-11-10: `paths.yaml` and `scripts/paths.R` were removed. Paths are now
> `file.path(QROOT, "<folder>/<file>")`, with `QROOT` the quetzal_code root
> (found via the `.quetzal_root` marker); references to paths.yaml keys above are historical.
