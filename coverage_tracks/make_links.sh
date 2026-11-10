#!/usr/bin/env bash
# Midway only: link the manuscript's prepared coverage-track folders.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
F=${TCGA:-/project2/yangili1/qhauck/spliceNMF/tcga}/genome_wide/flashier
link() { if [ -e "$2" ] && [ ! -L "$2" ]; then echo "SKIP (real file) $2" >&2; return 0; fi; ln -sfn "$1" "$2"; }
for g in NKAP YLPM1 SF3B1 neuronal_gsp; do [ -e "$F/pygenome/$g" ] && link "$F/pygenome/$g" "pygenome/$g"; done
for g in KDM4C TRIM38; do
  mkdir -p fig4b_tracks/$g
  for f in ${g}_sample_bigwigs.tsv ${g}_groups_file.txt junctions.bed junctions.bed.gz junctions.bed.gz.tbi \
           High_k30-3-..bw Matched_Controls-3-..bw MeanPSI.High_k30-3-..links MeanPSI.Matched_Controls-3-..links; do
    link "$F/GSP_pygenome/YLPM1/$g/$f" "fig4b_tracks/$g/$f"
  done
done
echo "coverage_tracks links OK"
