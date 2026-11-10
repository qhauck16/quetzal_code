# build_manifest.py

## Purpose & context
Builds the GDC download manifest for the mutation_ds samples that have no `.junc` yet.

## Decision rationale
- **Exact aliquot match.** Each sample is matched to GDC's aliquot ID using its full
  barcode, restricted to `STAR 2-Pass Genome` BAMs. The script exits if any barcode has
  zero or several matches.
- **Why this matters.** The old diff-splicing download joined on the 12-character patient
  ID, which filed some patients' BAMs under the wrong aliquot.

## Gotchas & warnings
- Uses GDC's public metadata API, so no token is needed here.
- The download itself is `download/download_array.sbatch`, a copy of
  `tcga/.../gdc_parallel`: 12 size-balanced caslake shards, with the token passed as
  `TOKEN=<path>` through `--export`. Check completeness with
  `cd download && bash verify.sh [--md5]`.

## Input/output
- **In:** `samples/download_list.tsv`.
- **Out:** `download/manifest.tsv` (gdc-client columns + `dir`), `download/rail_map.tsv`
  (rail_id → file), and `download/shards/`.

## File flow
```
samples/download_list.tsv ─► build_manifest.py ─► download/manifest.tsv ─► shards ─► download_array.sbatch ─► download/bams/<file_id>/
                                               └─► download/rail_map.tsv ─► Snakefile
```
