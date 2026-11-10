#!/bin/bash
#SBATCH --job-name=lc_m5
#SBATCH --account=pi-yangili1
#SBATCH --partition=caslake
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=06:00:00
#SBATCH --array=1-18
#SBATCH --output=logs/rerun_m5/%A_%a.log
# Re-cluster (leafcutter_cluster_regtools -m 5, was 10) and re-run leafcutter_ds for one
# existing run dir, reusing its juncfiles.txt / leaf_pheno_data.txt (no re-selection).
# Mirrors Snakefile rules leafcutter_cluster + leafcutter_ds; run on midway3 because the
# envs live on /project2. Old outputs: {results,control_results}.pre_m5_2026-10-28/.
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"   # submit from leafcutter_analyses/diff_splicing/
R=$(sed -n "${SLURM_ARRAY_TASK_ID}p" logs/rerun_m5/runs.txt)
PY=${PY:-python}   # python from envs/python2.7.yml (leafcutter_cluster_regtools.py is python 2)
LC=${LEAFCUTTER_DIR:?set LEAFCUTTER_DIR to a leafcutter clone}
E=${LEAFCUTTER_DS_ENV:?set LEAFCUTTER_DS_ENV to the leafcutter_ds conda env prefix}
P=$(basename $R)
echo "[$(date)] $R"
cd $R
rm -f ${P}_pooled ${P}_refined ${P}_sortedlibs ${P}_perind*.counts.gz *.junc.${P}.sorted.gz leafcutter_perind_numers.counts.gz leafcutter_ds_*.txt
$PY $LC/clustering/leafcutter_cluster_regtools.py -j $PWD/juncfiles.txt -m 5 -o $P -l 500000 -p 0.0001
zcat ${P}_perind_numers.counts.gz | awk 'NR==1 {for(i=1;i<=NF;i++) $i=$i".junc"}1' | gzip > leafcutter_perind_numers.counts.gz
echo "[$(date)] clustered: $(zcat leafcutter_perind_numers.counts.gz | wc -l) introns"
PATH=$E/bin:$PATH $E/bin/Rscript $LC/scripts/leafcutter_ds.R --num_threads 4 leafcutter_perind_numers.counts.gz leaf_pheno_data.txt -o leafcutter_ds
echo "[$(date)] done: $(wc -l < leafcutter_ds_cluster_significance.txt) clusters"
