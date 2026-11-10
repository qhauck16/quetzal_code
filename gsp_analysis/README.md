# gsp_analysis

Characterisation of the genome-wide splicing programs (GSPs) from the Quetzal v0.1.1
flashier object (`../quetzal_runs/genome_wide/softimpute_flash_300_qc_filtered.RDS`).
Each script is the chunk(s) of the original analysis notebook that wrote the file, with
paths updated; each was re-run and its output is byte-identical to the manuscript's.

| Script (run from this directory) | Writes | From notebook |
|---|---|---|
| `top_junctions.R` | `fact_{10,13,30}_top_junctions.tsv` (top unannotated junction per top-loaded feature, with splice type) | `spec_flash_fact.Rmd` |
| `alsafadi_to_hg38.R` | `Alsafadi2016_supp1_hg38.tsv` (Alsafadi et al. 2016 SF3B1 junctions, lifted to hg38) | `spec_flash_fact.Rmd` |
| `cmh_enrichment.R` | `cmh_mutation_enrichment_am0_9.tsv` (GSP x mutated gene, CMH by cohort; Fig. 2d) | `flash_fact_char.Rmd` |

Each takes an optional output directory (default: here).
