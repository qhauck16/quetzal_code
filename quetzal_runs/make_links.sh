#!/usr/bin/env bash
# Midway only: point quetzal_runs/ at the manuscript's Quetzal v0.1.1 outputs under
# spliceNMF/tcga/ (layout = Quetzal's own output/ directory).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
T=${TCGA:-/project2/yangili1/qhauck/spliceNMF/tcga}
link() { [ -e "$1" ] || { echo "MISSING $1" >&2; return 1; }
         if [ -e "$2" ] && [ ! -L "$2" ]; then echo "SKIP (real file) $2" >&2; return 0; fi
         mkdir -p "$(dirname "$2")"; ln -sfn "$1" "$2"; }
# gene_level/<chr>/<gene>/{res.RDS,de_res.RDS,dat.tsv,whole_factor.html} (+ results.tsv etc.
# written there by ../gene_level_analysis/)
for d in "$T"/gene_level_poisson/chr*/; do link "${d%/}" "gene_level/$(basename "$d")"; done
# genome_wide/<chr>/FastTopics_output/<gene>/res.RDS, merged features, flashier object
for d in "$T"/genome_wide/260116_filters/chr*/; do link "${d%/}/FastTopics_output" "genome_wide/$(basename "$d")/FastTopics_output"; done
link "$T/genome_wide/flashier/260116_multinom_tcga_just_unn_just_factors_of_interest.tsv" genome_wide/softimpute_features.tsv
link "$T/genome_wide/flashier/just_for_fun/nmd_investigation/softimpute_flash_300_qc_filtered.RDS" genome_wide/softimpute_flash_300_qc_filtered.RDS
echo "quetzal_runs links OK"
