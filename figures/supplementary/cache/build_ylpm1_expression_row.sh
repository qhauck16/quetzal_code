#!/usr/bin/env bash
# Cache for S17: header + YLPM1 row of the 7.7 GB TCGA expression matrix
# (gene x rail_id scaled RPM), so figS18/build.R never has to fread it.
set -euo pipefail
SRC="$(dirname "$0")/../../../expression/combined_scaled_rpm.tsv"   # quetzal_code/expression/
OUT="$(dirname "$0")/ylpm1_expression_row.tsv"
awk -F'\t' 'NR == 1 || $1 == "YLPM1"' "$SRC" > "$OUT.tmp" && mv "$OUT.tmp" "$OUT"
