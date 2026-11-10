#!/bin/bash

######################################################################
# @author      : qhauck (qhauck@midway2-login2.rcc.local)
# @file        : diff_splicing
# @created     : Wednesday Jan 07, 2026 11:49:12 CST
#
# @description : 
######################################################################

set -euo pipefail

# ================================
# Configuration
# ================================
CONDA_ENV="leafcutter_ds"
COUNTS="all_samples_perind_numers.counts.gz"
BAM_DIR="../bams"
THREADS=4

# Optional exon annotation (uncomment and set if desired)
# EXON_FILE="gencode.v31.exons.txt"
# EXON_ARG="--exon_file ${EXON_FILE}"

SCRIPT_DIR="${LEAFCUTTER_DIR:?set LEAFCUTTER_DIR to a leafcutter clone}/scripts"
OUTDIR="leafcutter_ds_results"
mkdir -p "${OUTDIR}"

# ================================
# Activate conda environment
# ================================
eval "$(conda shell.bash hook)"
conda activate "${CONDA_ENV}"

# ================================
# Helper: make groups file
# ================================
make_groups_file () {
    local g1=$1
    local g2=$2
    local outfile=$3

    rm -f "${outfile}"

    for bam in ${BAM_DIR}/${g1}/*.bam; do
        echo "$(basename "${bam}") ${g1}" >> "${outfile}"
    done

    for bam in ${BAM_DIR}/${g2}/*.bam; do
        echo "$(basename "${bam}") ${g2}" >> "${outfile}"
    done
}

# ================================
# Pairwise comparisons
# ================================

# 1) TCGA-KICH vs TCGA-TGCT
GROUPS_FILE="${OUTDIR}/groups_KICH_vs_TGCT.txt"
make_groups_file "TCGA-KICH" "TCGA-TGCT" "${GROUPS_FILE}"

Rscript ${SCRIPT_DIR}/leafcutter_ds.R \
    --num_threads ${THREADS} \
    ${COUNTS} \
    ${GROUPS_FILE} \
    > ${OUTDIR}/KICH_vs_TGCT.log 2>&1
mv leafcutter_ds_cluster_significance.txt KICH_vs_TGCT_ds_cluster_significance.txt
mv leafcutter_ds_effect_sizes.txt KICH_vs_TGCT_ds_effect_sizes.txt

# 2) TCGA-KICH vs TCGA-LUAD
GROUPS_FILE="${OUTDIR}/groups_KICH_vs_LUAD.txt"
make_groups_file "TCGA-KICH" "TCGA-LUAD" "${GROUPS_FILE}"

Rscript ${SCRIPT_DIR}/leafcutter_ds.R \
    --num_threads ${THREADS} \
    ${COUNTS} \
    ${GROUPS_FILE} \
    > ${OUTDIR}/KICH_vs_LUAD.log 2>&1
mv leafcutter_ds_cluster_significance.txt KICH_vs_LUAD_ds_cluster_significance.txt
mv leafcutter_ds_effect_sizes.txt KICH_vs_LUAD_ds_effect_sizes.txt

# 3) TCGA-TGCT vs TCGA-LUAD
GROUPS_FILE="${OUTDIR}/groups_TGCT_vs_LUAD.txt"
make_groups_file "TCGA-TGCT" "TCGA-LUAD" "${GROUPS_FILE}"

Rscript ${SCRIPT_DIR}/leafcutter_ds.R \
    --num_threads ${THREADS} \
    ${COUNTS} \
    ${GROUPS_FILE} \
    > ${OUTDIR}/TGCT_vs_LUAD.log 2>&1
mv leafcutter_ds_cluster_significance.txt TGCT_vs_LUAD_ds_cluster_significance.txt
mv leafcutter_ds_effect_sizes.txt TGCT_vs_LUAD_ds_effect_sizes.txt


echo "All LeafCutter differential splicing analyses completed."

