#!/bin/bash
#SBATCH --job-name=ylpm1_lowlof
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --array=1-91%30
#SBATCH --output=logs/call_%a.log
# Re-call YLPM1 in the 91 LoF mutants with GSP #30 < 1e-4 on GDC slices (slice.sh), mirroring
# wxs_mut_calling/unified/Snakefile exactly: paired Mutect2 (low-threshold flags) + FilterMutectCalls on
# WXS (and WGS where it exists), bcftools mpileup|call on RNA, VEP (--refseq, same cache) on both.
# Target = ../../targets/ylpm1.target_exons.merged.bed, reference GRCh38.d1.vd1.fa (as unified).
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"   # submit from this directory: sbatch -A <account> -p <partition> call.sbatch
# Tools: activate the envs in ../../unified/envs/ (gatk.yml: java; bcftools.yml; vep.yml) or put them on PATH.
GATK=${GATK:-gatk}   # GATK 4.6.2.0
BT=$(dirname "$(command -v bcftools)")
VEPE=$(dirname "$(command -v vep)")
REF=../../../inputs/GRCh38.d1.vd1.fa; BED=../../targets/ylpm1.target_exons.merged.bed; CACHE=../../../inputs/vep_cache/
RAIL=$(awk -F'\t' -v n=$SLURM_ARRAY_TASK_ID 'NR == n + 1 {print $2}' samples.tsv)
mkdir -p mutect2 mpileup vep
for ES in wxs wgs; do
  t=bams/${RAIL}_${ES}_tumor.bam; n=bams/${RAIL}_${ES}_normal.bam
  [ -s $t ] && [ -s $n ] || continue
  NSM=$($BT/samtools view -H $n | grep -m1 '^@RG' | tr '\t' '\n' | grep '^SM:' | cut -c4-)
  $GATK Mutect2 -R $REF -I $t -I $n -normal "$NSM" -L $BED \
    --max-reads-per-alignment-start 0 --disable-adaptive-pruning true --pruning-lod-threshold -100 \
    --initial-tumor-lod 0 --tumor-lod-to-emit 0 --active-probability-threshold 0.001 \
    -O mutect2/${RAIL}_${ES}.raw.vcf
  $GATK FilterMutectCalls -R $REF -V mutect2/${RAIL}_${ES}.raw.vcf -O mutect2/${RAIL}_${ES}.filtered.vcf
  $VEPE/vep -i mutect2/${RAIL}_${ES}.filtered.vcf -o vep/${RAIL}_${ES}.vep.vcf --cache --offline --dir_cache $CACHE \
    --refseq --use_given_ref --format vcf --force_overwrite --vcf
done
r=bams/${RAIL}_rna_tumor.bam
if [ -s $r ]; then
  $BT/bcftools mpileup -R $BED -f $REF -Q 20 -q 20 $r | $BT/bcftools call -mv -Ov -o mpileup/${RAIL}_rna.vcf
  $VEPE/vep -i mpileup/${RAIL}_rna.vcf -o vep/${RAIL}_rna.vep.vcf --cache --offline --dir_cache $CACHE \
    --refseq --use_given_ref --format vcf --force_overwrite --vcf
fi
echo "DONE $RAIL"
