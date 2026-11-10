#!/usr/bin/env bash
# Midway only: link the manuscript's aggregated gene-level tables (gene_cancer_type_hits.Rmd outputs).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
P=${TCGA:-/project2/yangili1/qhauck/spliceNMF/tcga}/gene_level_poisson/plotting_genes
for f in gene_level_poisson_sig_results.tsv junction_level_sig_upregulated.tsv tsg_aberrant_splicing_candidates.tsv; do
  if [ -e "$f" ] && [ ! -L "$f" ]; then echo "SKIP (real file) $f" >&2; continue; fi
  ln -sfn "$P/$f" "$f"
done
echo "gene_level_analysis links OK"
