#!/usr/bin/env bash
# Midway only: link the manuscript's expression matrices and recount3 downloads.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
X=${TCGA:-/project2/yangili1/qhauck/spliceNMF/tcga}/genome_wide/flashier/tcga_gene_expression
for f in combined_rpm.tsv combined_scaled_rpm.tsv tcga_counts; do
  if [ -e "$f" ] && [ ! -L "$f" ]; then echo "SKIP (real file) $f" >&2; continue; fi
  ln -sfn "$X/$f" "$f"
done
echo "expression links OK"
