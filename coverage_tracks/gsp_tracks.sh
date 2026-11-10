#!/usr/bin/env bash
# Second half of the Fig. 4b track preparation (after gsp_tracks.R): the commands
# prep_pygenometrack() prints, for KDM4C and TRIM38. Needs bgzip/tabix and
# GenometracksByGenotype (github.com/bfairkun/GenometracksByGenotype; set GTBG to the clone).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
: "${GTBG:?set GTBG to a GenometracksByGenotype clone}"
OUT=${1:-fig4b_tracks}
for spec in "TRIM38 chr6:25,973,077-25,984,559" "KDM4C chr9:7,011,426-7,014,212"; do
  set -- $spec; g=$1; region=$2; d=$OUT/$g
  (head -n 1 $d/junctions.bed && tail -n +2 $d/junctions.bed | sort -k1,1 -k2,2n) | bgzip -c > $d/junctions.bed.gz
  tabix -f -p bed $d/junctions.bed.gz
  python $GTBG/AggregateBigwigsForPlotting.py --BigwigList $d/${g}_sample_bigwigs.tsv --BigwigListType KeyFile \
    --OutputPrefix $d/ --Region $region --GroupSettingsFile $d/${g}_groups_file.txt \
    --Bed12GenesToIni $GTBG/PremadeTracks/gencode.v26.FromGTEx.genes.bed12.gz \
    --Normalization WholeGenome --BedfileForSashimiLinks $d/junctions.bed.gz
done
