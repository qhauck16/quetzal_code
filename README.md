# quetzal_code

Everything needed to regenerate the analyses and figures of the Quetzal manuscript.
Each folder holds its code and writes its outputs inside itself, so the folder
layout is the map of every file the figures read.

| Folder | What | Run order |
|---|---|---|
| `inputs/` | public inputs: `fetch_inputs.sh` (single files), `fetch_snaptron.sh`, `fetch_bigwigs.sh` | 1 |
| `quetzal_runs/` | Quetzal v0.1.1 gene-level and genome-wide outputs (run Quetzal, drop `output/` here) | 2 |
| `expression/` | recount3 TCGA gene expression | 2 |
| `mutation_processing/` | cBioPortal mutations, GRCh38 merge, AlphaMissense | 2 |
| `leafcutter_analyses/` | junction tables, PSI, LeafCutter2 classes, differential splicing | 3 |
| `support_tables/` | covariates, aggregated PSI / RPM | 3 (after leafcutter steps 1-3) |
| `gene_level_analysis/` | per-cohort tests on gene-level fits; Supp. Fig. 1 panels | 4 |
| `gsp_analysis/` | GSP characterisation (top junctions, CMH enrichment) | 4 |
| `variant_recalling/` | WXS / RNA-seq re-calling, YLPM1 biallelic analysis | 4 |
| `coverage_tracks/` | genome-browser track inputs | 4 |
| `benchmarking/` | Quetzal vs LeafCutter (Fig. 1b) | 4 |
| `gist/` | GIST cohort (HRA005970), Fig. 4f | 4 |
| `figures/` | figure code; each `figure<N>/compose.R`, `supplementary/build_all.sh` | 5 |

Figure scripts locate this directory by walking up to the marker file
`.quetzal_root` and build paths as `file.path(QROOT, "<folder>/<file>")`, so run
them from inside it. The analysis pipelines use paths relative to their own folder. Several steps need controlled-access TCGA data (GDC/dbGaP); their
outputs are deposited instead. On midway, each folder's `make_links.sh` links the
manuscript's existing outputs.
