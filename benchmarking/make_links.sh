#!/usr/bin/env bash
# Midway only: link the manuscript's benchmark LeafCutter outputs.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
L=${TCGA:-/project2/yangili1/qhauck/spliceNMF/tcga}/model_eval/leafcutter
for f in "$L"/*_ds_cluster_significance.txt "$L"/*_ds_effect_sizes.txt "$L"/all_samples_perind_numers.counts.gz "$L"/junc_files "$L"/leafcutter_ds_results; do
  [ -e "$f" ] || continue; h=leafcutter/$(basename "$f")
  if [ -e "$h" ] && [ ! -L "$h" ]; then echo "SKIP (real file) $h" >&2; continue; fi
  ln -sfn "$f" "$h"
done
echo "benchmarking links OK"
