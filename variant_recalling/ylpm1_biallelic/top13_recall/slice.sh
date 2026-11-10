#!/bin/bash
# Slice the YLPM1 locus out of every BAM in manifest.tsv via the GDC BAM slicing API -> bams/<rail>_<role>.bam
# The token is read from GDC_TOKEN_FILE inside the curl header and never printed.
cd "$(dirname "$0")"
T=${GDC_TOKEN_FILE:?set GDC_TOKEN_FILE to your GDC token file (controlled access)}
S=${SAMTOOLS:-samtools}
REG=chr14:74760000-74862000
mkdir -p bams
tail -n +2 manifest.tsv | while IFS=$'\t' read rank rail role sample fid fname size; do
  o=bams/${rail}_${role}
  [ -s $o.bam.bai ] && continue
  curl -s --fail --retry 3 -o $o.raw.bam -H "X-Auth-Token: $(cat $T)" "https://api.gdc.cancer.gov/slicing/view/$fid?region=$REG" || { echo "FAIL $o"; continue; }
  $S sort -o $o.bam $o.raw.bam && rm $o.raw.bam && $S index $o.bam
  echo "$o $($S view -c $o.bam) reads SM=$($S view -H $o.bam | grep -m1 '^@RG' | tr '\t' '\n' | grep '^SM:' | cut -c4-)"
done
echo SLICED
