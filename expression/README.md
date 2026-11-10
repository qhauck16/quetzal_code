# expression

TCGA gene expression used for Fig. 4g and the YLPM1 neighbouring-gene check.

1. `bash download_gene_sums.sh`: recount3 TCGA gene counts, 33 cohorts (`tcga_counts/`).
2. `Rscript gen_count_file.R`: per cohort, gene counts -> RPM (GENCODE v46 gene names, first
   Ensembl ID per name), then z-scored per gene within cohort; genes present in every cohort.
   Writes `combined_rpm.tsv` and `combined_scaled_rpm.tsv` (genes x samples, rail_id columns).

Reads `../inputs/gencode_v46_granges.RDS` and `../inputs/tcga_v2_samples.tsv`.
