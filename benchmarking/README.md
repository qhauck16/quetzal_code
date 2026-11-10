# benchmarking

Quetzal vs LeafCutter on three TCGA cohorts (Fig. 1b, Supplementary Fig. 2).

| Step | Code | Produces |
|---|---|---|
| 1. Samples | `choosing_ids.R`: 20 primary tumours each from LUAD, TGCT, KICH (`set.seed(123)`) | sample list |
| 2. Download | GDC manifest/sample sheet (`gdc_manifest.2026-01-05.173855.txt`, `gdc_sample_sheet.2026-01-05.tsv`); `barcodes_to_manifest.R` -> `manifest_tumor_rnaseq.tsv`; `renaming_downloaded_bams.R` -> `bams/<cohort>/<rail_id>.bam` | BAMs (controlled access) |
| 3. Junctions | `leafcutter/regtools_extract.sh` (conda env with regtools) | `leafcutter/junc_files/` |
| 4. Clustering | `leafcutter/clustering.sh` (`leafcutter_cluster_regtools.py -m 50`, Python 2.7 env) | `leafcutter/all_samples_perind_numers.counts.gz` |
| 5. Differential splicing | `leafcutter/diff_splicing.sh` (`leafcutter_ds`, env `leafcutter_ds`) | `leafcutter/<A>_vs_<B>_ds_{cluster_significance,effect_sizes}.txt` |

Set `LEAFCUTTER_DIR` to a clone of github.com/davidaknowles/leafcutter. The Quetzal
side of the comparison uses the main gene-level run (`../quetzal_runs/gene_level/`).
