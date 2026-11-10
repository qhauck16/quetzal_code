# quetzal_runs

`gene_level/updated_Snakefile` is the Snakefile the manuscript's gene-level run used,
unchanged except for input paths: rule `run_gene_plots_and_objs` fits each gene
(`gene_plots_and_objs.R`, Quetzal v0.1.1's gene-level script) and rule
`run_cancer_specific_factors` runs the per-cohort tests (`cancer_specific_factors.R`).
Run it from `gene_level/`:
`snakemake -s updated_Snakefile --use-conda --jobs N` (add your cluster executor).
The genome-wide outputs were made with Quetzal v0.1.1 as below.

Outputs of Quetzal v0.1.1 (github.com/qhauck16/Quetzal, tag `v0.1.1`), laid out
exactly as Quetzal's own `output/` directory, so a re-run can be dropped in place:

```bash
git clone https://github.com/qhauck16/Quetzal && cd Quetzal && git checkout v0.1.1
# inputs: ../inputs/snaptron/tcga/<chr>/snaptron_output/ -> Quetzal/data/all_genes/<chr>/snaptron_output/,
#         ../inputs/tcga_v2_samples.tsv -> Quetzal/data/tcga_v2_samples.tsv
cd scripts/gene_level  && snakemake --use-conda --jobs N        # -> output/gene_level/
cd ../genome_wide      && bash setting_up_snakemake.sh          # -> output/genome_wide/<chr>/FastTopics_output/
Rscript fasttopics_to_flashier.R                                 # -> output/genome_wide/softimpute_features.tsv,
                                                                 #    softimpute_flash_300_qc_filtered.RDS
# then: cp -r ../../output/{gene_level,genome_wide} <this directory>/
```

| Path | What |
|---|---|
| `gene_level/<chr>/<gene>/` | per-gene fits (`res.RDS`, `de_res.RDS`, `dat.tsv`, `whole_factor.html`) and per-cohort tests (`results.tsv`, `junction_results.tsv`, `n_bb_tests.tsv`) |
| `genome_wide/<chr>/FastTopics_output/<gene>/res.RDS` | per-gene fits for the genome-wide factorization |
| `genome_wide/softimpute_features.tsv` | merged sample x unannotated-factor matrix (31,242 features) |
| `genome_wide/softimpute_flash_300_qc_filtered.RDS` | the flashier object; GSP #k = factor k |

The NMF fits are multi-threaded, so a re-run reproduces the inputs and settings but not
bit-identical factors; the figures use the manuscript's objects. On midway,
`bash make_links.sh` links them.
