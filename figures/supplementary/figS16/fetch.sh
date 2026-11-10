#!/bin/bash
# Slice the reads for the YLPM1 splice-site figure from GDC (BAM slicing API):
#   TCGA-CJ-4638 (KIRC), rank 2 on GSP #30, YLPM1 X1783_splice,
#   chr14:74812626 A>G = intron position -2 of MANE exon 11's acceptor.
#   RNA-seq STAR 2-Pass Genome (01A-02R-1325-07): whole YLPM1 locus.
#   Tumour WGS (01A-01D-A97E-36; 295 GB file, so slicing only): +-300 bp.
# The token is passed by path inside the curl header and never printed.
set -euo pipefail
cd "$(dirname "$0")"
T=${GDC_TOKEN_FILE:?set GDC_TOKEN_FILE to your GDC token file (controlled access)}
S=${SAMTOOLS:-samtools}
mkdir -p bams
slice() {  # file_id region out
  curl -s --fail -o "bams/$3.raw.bam" -H "X-Auth-Token: $(cat $T)" \
    "https://api.gdc.cancer.gov/slicing/view/$1?region=$2"
  $S sort -o "bams/$3.bam" "bams/$3.raw.bam" && rm "bams/$3.raw.bam" && $S index "bams/$3.bam"
  echo "$3: $($S view -c bams/$3.bam) reads, SM=$($S view -H bams/$3.bam | grep -m1 '^@RG' | tr '\t' '\n' | grep '^SM:' | cut -c4-)"
}
slice 5d10fb20-d7c0-4b83-9d2d-d6c722d15154 chr14:74760000-74840000 rna_YLPM1
slice e7d6c16f-1094-4ee8-9c23-eab8fe90d16b chr14:74812326-74812926 wgs_tumour_site
