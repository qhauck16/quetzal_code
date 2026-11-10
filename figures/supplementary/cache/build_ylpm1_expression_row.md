# build_ylpm1_expression_row.sh

**Purpose & context** — S17 needs YLPM1 expression per sample. The source
matrix (`tcga/genome_wide/flashier/tcga_gene_expression/combined_scaled_rpm.tsv`,
~7.7 GB, first column `gene_name`, then one column per rail_id) is too big
to read inside a figure build, so this pulls the header + the YLPM1 row.

**Decision rationale** — awk streams the file (no memory cost); the output
keeps the original layout so figS18/build.R reads it with fread(header = TRUE).

**Gotchas** — Always read the output (and source) with `header = TRUE`:
fread's header autodetect misfires on these all-numeric-column matrices.
Streaming 7.7 GB takes a few minutes; fine on a login node, or wrap in sbatch.

**Input/output** — combined_scaled_rpm.tsv → `cache/ylpm1_expression_row.tsv` (2 lines).

```
combined_scaled_rpm.tsv ──► build_ylpm1_expression_row.sh ──► ylpm1_expression_row.tsv ──► figS18/build.R
```
