#!/usr/bin/env bash
# Midway only: link the manuscript's re-calling outputs (per-gene calls + VEP, YLPM1 re-calls, ABSOLUTE CN).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
W=${TCGA:-/project2/yangili1/qhauck/spliceNMF/tcga}/genome_wide/flashier/wxs_mut_calling
link() { [ -e "$1" ] || { echo "MISSING $1" >&2; return 1; }
         if [ -e "$2" ] && [ ! -L "$2" ]; then echo "SKIP (real file) $2" >&2; return 0; fi
         mkdir -p "$(dirname "$2")"; ln -sfn "$1" "$2"; }
for g in nkap ylpm1 sf3b1; do link $W/unified/$g unified/$g; done
link $W/ylpm1_absolute_cn.tsv ylpm1_biallelic/ylpm1_absolute_cn.tsv
for d in top13_recall lof_low_gsp30; do
  for x in manifest.tsv bams mutect2 mpileup vep recall_calls.tsv; do [ -e $W/$d/$x ] && link $W/$d/$x ylpm1_biallelic/$d/$x; done
done
echo "variant_recalling links OK"
