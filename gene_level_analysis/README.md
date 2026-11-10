# gene_level_analysis

Aggregation of the per-gene, per-cohort tests and the Supplementary Fig. 1 panels.

The tests themselves (`cancer_specific_factors.R`) run as the second rule of
`../quetzal_runs/gene_level/updated_Snakefile`, which writes `results.tsv`,
`junction_results.tsv` and `n_bb_tests.tsv` next to each gene's fits.

| Code | Produces |
|---|---|
| `gene_cancer_type_hits.Rmd` (delete `gene_cancer_type_hits_cache.RDS` to force a rebuild) | `gene_level_poisson_sig_results.tsv`, `junction_level_sig_upregulated.tsv`, `tsg_aberrant_splicing_candidates.tsv` |
| `showcase_plots/` | Supp. Fig. 1 per-gene panels |

Needs `../inputs/cosmic_tier_1.tsv`. On midway, `bash make_links.sh` links the
manuscript's aggregated tables.
