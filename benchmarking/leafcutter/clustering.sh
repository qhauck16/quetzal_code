#!/usr/bin/env sh

######################################################################
# @author      : qhauck (qhauck@midway2-login2.rcc.local)
# @file        : clustering
# @created     : Wednesday Jan 07, 2026 10:41:23 CST
#
# @description : 
######################################################################

conda deactivate
conda activate python2.7

python ${LEAFCUTTER_DIR:?set LEAFCUTTER_DIR to a leafcutter clone}/clustering/leafcutter_cluster_regtools.py -j test_juncfiles.txt -m 50 -o all_samples -l 500000

