#!/usr/bin/env bash
# Midway only: point the outputs of this folder at the manuscript's copies under spliceNMF/tcga/.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
T=${TCGA:-/project2/yangili1/qhauck/spliceNMF/tcga}
link() { [ -e "$1" ] || { echo "MISSING $1" >&2; return 1; }
         if [ -e "$2" ] && [ ! -L "$2" ]; then echo "SKIP (real file) $2" >&2; return 0; fi
         ln -sfn "$1" "$2"; }
link $T/standardized_covariates.RDS                          standardized_covariates.RDS
link $T/standardized_covariates.tsv                          standardized_covariates.tsv
link $T/genome_wide/flashier/aggregated_psi_medians.tsv       aggregated_psi_medians.tsv
link $T/genome_wide/flashier/aggregated_junction_rpms.tsv     aggregated_junction_rpms.tsv
echo "support_tables links OK"
