# cache/s2/run.sh

## Purpose & context
Builds the S2 inputs: runs figure 1 panel_d_e's `scripts/make_gene_outputs.R`
for TCF12 chr15 k10 (TCGA-LGG and TCGA-GBM), SRC chr20 k8 (TCGA-COAD) and
SETD2 chr3 k2 (TCGA-LAML) -- the same generator as main Figure 1d (EZH2).

## Decision rationale
Reusing fig 1d's generator + its `config.yaml` keeps S2 byte-compatible with
the fig 1d drawing code (`figure1/panel_d_e/scripts/lib.R`) used by figS3/build.R.

## Gotchas & warnings
- `cd`s to `figure1/panel_d_e` itself (config and script paths are relative there); `set -e`.
- Output dir `<GENE>_<COHORT without TCGA->/<GENE>/`.
- Log has benign fread "column names vs columns" warnings on psi_table.tsv (unlabeled
  rail_id column) and removed non-finite rows from stat_ecdf.

## User preferences
Conda-env Rscript by absolute path.

## Input/output
In: per-gene fits / PSI via `figure1/panel_d_e/config.yaml`. Out: `s2/{TCF12_LGG,TCF12_GBM,SRC_COAD,SETD2_LAML}/<GENE>/`
(`gene_data.RDS`, `junctions_classified.tsv`, plus the generator's own plots).

## File flow schematic
```
figure1/panel_d_e/config.yaml ─┐
make_gene_outputs.R ───────────┴──► run.sh ──► s2/<GENE>_<COHORT>/<GENE>/gene_data.RDS ──► figS3/build.R
                                              └─► junctions_classified.tsv ──────────────┘
```
