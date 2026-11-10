#!/usr/bin/env sh

######################################################################
# @file        : setting_up_snakemake
# @created     : Thursday Jan 23, 2025 14:47:54 CST
#
# @description : 
######################################################################

# Run from leafcutter_analyses/junction_tables_gtex/. One Snakemake job per chromosome
# (each builds chr<N>/FastTopics_output/<gene>/dat.tsv). Set your SLURM partition/account.
SLURM_PARTITION=${SLURM_PARTITION:?set SLURM_PARTITION}; SLURM_ACCOUNT=${SLURM_ACCOUNT:?set SLURM_ACCOUNT}
for d in ../../inputs/snaptron/gtex/chr*/
do
    chr=$(basename $d)
    mkdir -p $chr

    cp lf_Snakefile $chr/
    cd $chr
    snakemake -s lf_Snakefile --executor slurm --default-resources slurm_partition=$SLURM_PARTITION slurm_account=$SLURM_ACCOUNT --jobs 90 --use-conda --rerun-incomplete --keep-going --retries 2
    cd ..
done
