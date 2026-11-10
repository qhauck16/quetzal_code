#!/usr/bin/env bash
# Midway only: link the manuscript's GSP-analysis outputs the figures read.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
F=${TCGA:-/project2/yangili1/qhauck/spliceNMF/tcga}/genome_wide/flashier
for f in "$F"/fact_*_top_junctions*.tsv "$F"/cmh_mutation_enrichment_am0_*.tsv "$F"/Alsafadi2016_supp1_hg38.tsv; do
  b=$(basename "$f"); if [ -e "$b" ] && [ ! -L "$b" ]; then echo "SKIP (real file) $b" >&2; continue; fi
  ln -sfn "$f" "$b"
done
echo "gsp_analysis links OK"
