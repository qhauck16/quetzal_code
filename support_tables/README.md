# support_tables

| Output | Code | Reads |
|---|---|---|
| `standardized_covariates.RDS` (+ `.tsv`) | `covariate_table.R` | `../inputs/tcga_v2_samples.tsv`, `../inputs/analyte.tsv` |
| `aggregated_psi_medians.tsv` | `support_table_generation.Rmd`, chunk "aggregated psi medians" | `../leafcutter_analyses/psi_files/tcga/` |
| `aggregated_junction_rpms.tsv` | `support_table_generation.Rmd`, chunk "junction rpm medians" | `../leafcutter_analyses/junction_tables/`, `../inputs/tcga_v2_samples.tsv` |

Run from this directory, after `../leafcutter_analyses/` steps 1-3. The two
aggregated tables feed `../leafcutter_analyses/diff_splicing/unified_results/`.
