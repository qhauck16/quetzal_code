# gist

Public GIST cohort (Xie et al.; GSA-Human HRA005970): paired tumour/normal RNA-seq,
used for Fig. 4f.

| Step | Code | Produces |
|---|---|---|
| 1. Sample choice | `processing_file_ids.Rmd` (from `YLPM1_mutant_status.txt`, `sample_to_accession.txt`, `accession_to_experiment.txt`, `HRA005970.xlsx`) | `GIST_files_to_download.tsv` (18 runs) |
| 2. Download | `paired_tumor_normal/download_fqgz.sh` (CNCB GSA-Human, `HRA005970_files.tsv`) | `paired_tumor_normal/HRS*/` FASTQs |
| 3. Trim, align, junctions, cluster | `paired_tumor_normal/align_Snakefile`: concatenate lanes, fastp, STAR 2.7.5c (GDC GRCh38.d1.vd1 + GENCODE v36 index), regtools (`-a 8 -m 50 -M 500000`), LeafCutter clustering (`-m 0`) | `paired_tumor_normal/<HRS>_f1_combined.fqAligned.sortedByCoord.out.bam`, `junc_files/`, `ylpm1_perind_numers.counts.gz` |

Run step 3 from `paired_tumor_normal/`: `snakemake -s align_Snakefile --use-conda --jobs N`
(set `LEAFCUTTER_DIR`; STAR index from `../inputs/fetch_inputs.sh --only star_index`).
No differential-splicing test is run within GIST. Fig. 4f scores each GIST sample's PSI at
the junctions found in the TCGA YLPM1 run (`../figures/figure4/scripts/build_cache_panel_f_tcga.R`:
`leafcutter_analyses/diff_splicing/unified_results/ylpm1/20v60`, alt 3'ss + alt 5'ss) and
clusters them (`../figures/figure4/panel_f/build.R`), reading `ylpm1_perind_numers.counts.gz`,
`YLPM1_mutant_status.txt` and `accession_to_experiment.txt`. On midway, `bash make_links.sh` links every processed file.
