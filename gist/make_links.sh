#!/usr/bin/env bash
# Midway only: link every downloaded / processed GIST file at the path the pipeline writes it to
# (paired_tumor_normal/: HRS*/ FASTQ + trimmed reads, STAR BAMs, junc_files/, LeafCutter counts).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/paired_tumor_normal"
P=${TCGA:-/project2/yangili1/qhauck/spliceNMF/tcga}/genome_wide/flashier/y_gene_mut/paired_tumor_normal
for f in "$P"/*; do
  b=$(basename "$f")
  case $b in download_fqgz.sh|align.sh|align_Snakefile|leaf_Snakefile|envs) continue ;; esac
  if [ -e "$b" ] && [ ! -L "$b" ]; then echo "SKIP (real file) $b" >&2; continue; fi
  ln -sfn "$f" "$b"
done
echo "gist links OK"
