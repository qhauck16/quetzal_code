#!/usr/bin/env bash
# Midway only: point this folder's outputs at the manuscript's copies under spliceNMF/tcga/.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
T=${TCGA:-/project2/yangili1/qhauck/spliceNMF/tcga}; M=$T/genome_wide/flashier/just_for_fun/mutation_merge
link() { [ -e "$1" ] || { echo "MISSING $1" >&2; return 1; }
         if [ -e "$2" ] && [ ! -L "$2" ]; then echo "SKIP (real file) $2" >&2; return 0; fi
         mkdir -p "$(dirname "$2")"; ln -sfn "$1" "$2"; }
link $T/tcga_mutation_data/pan_can_atlas mutations/pan_can_atlas
link $T/tcga_mutation_data/gdc           mutations/gdc
for f in merged_grch38.tsv missense_keys.tsv am_hits.tsv merged_grch38_with_am.tsv; do link $M/$f $f; done
for f in $M/refs/*_alphamissense.tsv; do link "$f" am_refs/$(basename "$f"); done
echo "mutation_processing links OK"
