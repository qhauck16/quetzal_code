# variant_recalling

Somatic re-calling of the splicing-factor genes in TCGA WXS and RNA-seq, plus the
YLPM1 biallelic-loss analysis. All BAM-level steps need controlled-access TCGA data
(GDC/dbGaP).

| Step | Code | Produces |
|---|---|---|
| 1. Target exons | `make_target_beds.sh` | `targets/{nkap,ylpm1,sf3b1}.target_exons.merged.bed` (GENCODE v46 exons of NKAP/SLU7/FAM50A/PRPF8, YLPM1, SF3B1/SUGP1) |
| 2. Per-gene re-calling | `unified/Snakefile` (+ `config.yaml`, `scripts/`, `envs/`): sample selection, WXS Mutect2 (low-threshold) + FilterMutectCalls, RNA-seq `bcftools mpileup | call`, VEP (RefSeq, release 115) | `unified/<gene>/vep/`, merged calls (read by `../figures/scripts/variant_tables.R`) |
| 3. YLPM1 copy number | `ylpm1_biallelic/ylpm1_absolute_cn.R` | `ylpm1_biallelic/ylpm1_absolute_cn.tsv` (PanCanAtlas ABSOLUTE allelic state at YLPM1) |
| 4. YLPM1 locus re-calls | `ylpm1_biallelic/{top13_recall,lof_low_gsp30}/`: `build_manifest.py` -> `slice.sh` (GDC BAM slicing; set `GDC_TOKEN_FILE`) -> `call.sbatch` (same callers as step 2) -> `summarize.py` | `recall_calls.tsv`; top 13 on GSP #30 and the 91 YLPM1 LoF carriers with GSP #30 < 1e-4 |

Inputs: `../inputs/` (GDC GRCh38.d1.vd1 reference, VEP cache, GENCODE v46 GTF, TCGA
metadata, PanCanAtlas ABSOLUTE), `../mutation_processing/` (compiled mutations),
`../quetzal_runs/` (GSP loadings), `../expression/` (neighbouring-gene check,
`top13_recall/neighbor_expr.R`). Tools: GATK 4.6.2.0, bcftools/samtools, Ensembl VEP 115
(conda envs in `unified/envs/`). Step 1 was re-run and is byte-identical to the
manuscript's BEDs. On midway, `bash make_links.sh` links the manuscript's outputs.
