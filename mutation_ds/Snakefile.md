# Snakefile (mutation_ds)

## Purpose & context
Runs LeafCutter DS for each sample sheet in `samples/`. The steps and settings are copied
from `tcga/genome_wide/flashier/diff_splicing/Snakefile`, so results are comparable to the
GSP comparisons.

## Decision rationale
- **Runs on midway2.** The `leafcutter_ds` conda env and `~/software/leafcutter` exist
  only there, and rebuilding the env is fragile.
- **Existing envs by absolute path.** The regtools and python2.7 binaries come from
  `diff_splicing/.snakemake/conda/`. There is no `conda:` directive, so nothing gets
  rebuilt and `--use-conda` is not needed.
- **Existing `.junc` files are reused.** Samples that have one in `diff_splicing/junc_files`
  get a symlink to it; no BAM is needed for them.

## Gotchas & warnings
- A dry run fails with MissingInput until every BAM in `download/rail_map.tsv` is
  downloaded.
- `leaf_pheno_data.txt` lists controls first: leafcutter_ds.R encodes groups by first appearance, so this makes dPSI mutant − control. The first run (2026-10-08) was mutant-first; overlap.R corrects the sign.

## Run (midway2, from this dir)
`snakemake --executor slurm --default-resources slurm_partition=broadwl slurm_account=pi-yangili1 --jobs 90`

## File flow
```
samples/<comp>.tsv ─► setup ─► results/<comp>/{leaf_pheno_data,juncfiles}.txt ─┐
download/bams ─► junc (regtools / symlink) ─► junc_files/<rail>.junc ───────────┼─► leafcutter_cluster ─► leafcutter_ds ─► results/<comp>/leafcutter_ds_*.txt
```
