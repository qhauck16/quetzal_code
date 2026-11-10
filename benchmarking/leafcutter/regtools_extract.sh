#!/bin/sh
#SBATCH --job-name=regtools_junc_extract
#SBATCH --output=regtools_junc_extract.%j.out
#SBATCH --error=regtools_junc_extract.%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50G
#SBATCH --time=24:00:00

######################################################################
# @file        : regtools_extract
# @created     : Tuesday Jan 06, 2026 21:30:30 CST
#
# @description : 
######################################################################
conda deactivate
conda activate rna-seq
for bamfile in ../bams/*/*.bam; do
    echo Converting $bamfile to $bamfile.junc
    base=$(basename $bamfile)
    echo $base
    regtools junctions extract -a 8 -m 50 -M 500000 $bamfile -s XS -o junc_files/$base.junc
    echo junc_files/$base.junc >> test_juncfiles.txt
done


