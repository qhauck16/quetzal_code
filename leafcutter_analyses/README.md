# leafcutter_analyses

Code that produces every LeafCutter-derived file the figures read. Each step
writes its outputs inside this directory, at the paths the figure code reads, so
it finds them whether you re-ran the step or downloaded its outputs.

Inputs come from `../inputs/`, the GSP loadings from the Quetzal v0.1.1 run
(`../quetzal_runs/`), and `../support_tables/`.

## Steps, in run order

| Step | Code | Produces | Needs |
|---|---|---|---|
| 1. Per-gene junction tables | `junction_tables/Snakefile` (`make_dat.R`) | `junction_tables/chr*/<gene>/dat.tsv` | Snaptron TCGA tables, TCGA metadata, GENCODE v46 |
| 1b. GTEx junction tables | `junction_tables_gtex/setting_up_snakemake.sh` (`lf_Snakefile`, `gtex_LF_saving.R`) | `junction_tables_gtex/chr*/FastTopics_output/<gene>/dat.tsv` | Snaptron GTEx tables, GTEx metadata, GENCODE v46 (no factorization: counts only, all samples, `min_samples_per_junc = 0`) |
| 2. LeafCutter2 clustering | `leaf2_clustering/regtools_leafcutter_cluster.sh` | `leaf2_clustering/reg.*.unclassified.gz` | 330 regtools `.junc` files (`leaf2_clustering/junc_files/`), GENCODE v46 GTF, `LEAFCUTTER2_DIR` |
| 3. PSI tables | `psi/Snakefile` (`reg_dat_to_PSI.R`) | `psi_files/{tcga,gtex}/chr*/<gene>/psi_table.tsv` | steps 1, 1b, 2 |
| 4. Junction classification | `leaf2_classification/run_classifier_from_junc_tables.sh` | `leaf2_classification/clustering/Leaf2_junction_classifications.txt` | step 1, GENCODE v46 GTF, GRCh38 FASTA |
| 5. Differential splicing | `diff_splicing/Snakefile` (+ `config.yaml`, `scripts/`) | `diff_splicing/results/`, `control_results/`, `junc_files/` | GSP loadings (`scripts/extract_loadings_from_rds.R`), TCGA RNA-seq BAMs (controlled access), `LEAFCUTTER_DIR` |
| 6. Junction characterisation | `diff_splicing/unified_results/Snakefile` | `diff_splicing/unified_results/{nkap,sf3b1,ylpm1,shared}/` | step 5, `diff_splicing/gencode_v46_introns.tsv` (`build_gencode_v46_introns.R`), MANE, CoLa-seq branchpoints, `../support_tables/` |

Step 5 notes:
- Clustering uses `leafcutter_cluster_regtools.py -m 5`.
- `config.yaml` `extra_comparisons` adds `results/ylpm1_6v18`, used only to choose the YLPM1 SVM junctions. Its cases are the 6 samples with scaled GSP #30 ≥ 0.1 that have junction files.
- `download/bam_inventory.tsv` is the final BAM inventory (sample-mislabel fix already applied); re-running sample selection on it reproduces the current runs exactly (checked for nkap_20v60, sf3b1_30v90, ylpm1_20v60). The other `download/` tables record which GDC files were used and are read by the supplementary sample table.

Software: leafcutter (github.com/davidaknowles/leafcutter) and leafcutter2
(github.com/cfbuenabadn/leafcutter2, commit 62cccc6). Conda envs are in each step's `envs/`.
