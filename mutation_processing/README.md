# mutation_processing

TCGA somatic mutations from cBioPortal, merged on GRCh38 and annotated with
AlphaMissense. Run from this directory, in order:

| Step | Code | Produces |
|---|---|---|
| 1. Download + compile | `get_tcga_mutations.Rmd` (cBioPortal via `cbioportalR`: every `*_tcga_pan_can_atlas` and `*_tcga_gdc` study) | `mutations/{pan_can_atlas,gdc}/<cohort>.tsv`, `mutations/{pan_can_atlas,gdc}/compiled.tsv` |
| 2. Lift to GRCh38 + merge | `liftover_and_merge.R` (UCSC `liftOver`; set `LIFTOVER` if it is not on PATH) | `merged_grch38.tsv` |
| 3. AlphaMissense subsets | `make_am_tables.sh` | `missense_keys.tsv`, `am_hits.tsv`, `am_refs/<GENE>_alphamissense.tsv` |
| 4. Join AlphaMissense | `join_alphamissense.R` | `merged_grch38_with_am.tsv` |

Inputs: `../inputs/alphamissense/` (AlphaMissense hg38, hg19->hg38 chain).

The manuscript's download is from 2025-09-17. cBioPortal's GDC studies are updated
over time, so a fresh step 1 can differ. Steps 1 (compilation), 3 and their outputs
were re-run and are byte-identical to the manuscript's files. On midway,
`bash make_links.sh` links the manuscript's copies.
