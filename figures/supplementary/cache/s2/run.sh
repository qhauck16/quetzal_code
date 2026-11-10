#!/bin/bash
# Supplementary S2 inputs: figure 1 panel_d_e's make_gene_outputs.R, run for
# the S2 genes (TCF12 k10 in LGG and GBM, SRC k8 in COAD, SETD2 k2 in LAML).
set -e
cd "$(dirname "$0")/../../../figure1/panel_d_e"
for spec in "TCF12 chr15 k10 TCGA-LGG" "TCF12 chr15 k10 TCGA-GBM" "SRC chr20 k8 TCGA-COAD" "SETD2 chr3 k2 TCGA-LAML"; do
  set -- $spec
  /scratch/midway2/qhauck/conda_env/rstudio-server/bin/Rscript --vanilla scripts/make_gene_outputs.R --gene $1 --chr $2 --factor $3 --cancer_type $4 --config config.yaml --out_dir ../../supplementary/cache/s2/${1}_${4#TCGA-}/$1
done
