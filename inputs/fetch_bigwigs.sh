#!/usr/bin/env bash
# recount3 per-sample coverage bigwigs used for the coverage tracks. URL lists are the
# recount3 study-explorer exports: lists/tcga_bigwigs.csv (33 TCGA cohorts) and
# lists/gtex_bigwigs.csv (GTEx tissues). ~25,000 files, several TB.
#
#   bash fetch_bigwigs.sh            # both
#   bash fetch_bigwigs.sh tcga       # one source
#
# Saved as bigwigs/<src>/<study>/<rail_id>.bw (skips files already present).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
for src in ${1:-tcga gtex}; do
  tail -n +2 lists/${src}_bigwigs.csv | tr -d '"' | while IFS=',' read -r rail ext study url; do
    out=bigwigs/$src/$study/$rail.bw
    [ -s "$out" ] && continue
    mkdir -p "$(dirname "$out")"
    curl -fsSL --retry 3 -o "$out.part" "$url" && mv "$out.part" "$out"
  done
done
echo "done"
