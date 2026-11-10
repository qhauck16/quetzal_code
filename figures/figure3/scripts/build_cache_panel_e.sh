#!/bin/bash
# Figure 3 panel e: slice the two huge expression matrices down to the
# 20 columns the panel actually uses.
#
# combined_scaled_rpm.tsv is 7.7 GB and combined_rpm.tsv is 5.7 GB, both
# with 11,174 columns. The analysis only ever reads the top-20 NKAP
# samples (by k13) -- every row, twenty columns. fread()ing 13.4 GB to
# use 0.2% of it is what makes this panel look expensive; it is not.
#
# Writes cache/panel_e_scaled.tsv and cache/panel_e_raw.tsv, both
# gene_name + the 20 rail_id columns, in the ORDER given by the ids
# file (so downstream code can rely on the column order).
#
#   bash scripts/build_cache_panel_e.sh
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$HERE"
RSCRIPT=Rscript
EXPR="$HERE/../../expression"   # quetzal_code/expression/
mkdir -p cache

# 1. Top-20 rail_ids by k13, derived from the variant table rather than
#    hardcoded.
$RSCRIPT --vanilla -e '
suppressPackageStartupMessages(library(data.table))
source("../scripts/vt_paths.R")
f <- vt_path("NKAP_full_plus_variants")
d <- fread(f, data.table = FALSE)
stopifnot("k13" %in% names(d))
d <- d[order(-d$k13), ]
writeLines(as.character(d$rail_id[1:20]), "cache/panel_e_top_ids.txt")
cat("top-20 k13:", paste(round(d$k13[1:20], 4), collapse = " "), "\n")'

echo "[panel e] ids:" $(tr '\n' ' ' < cache/panel_e_top_ids.txt)

# 2. Slice each matrix. awk reads the header once, records the column
#    index of gene_name and of each wanted id, then emits only those.
slice () {
  local src="$1" dst="$2"
  echo "[panel e] slicing $(basename "$src") -> $dst"
  awk -F'\t' -v OFS='\t' -v idfile="cache/panel_e_top_ids.txt" '
    BEGIN { n = 0; while ((getline id < idfile) > 0) if (id != "") want[++n] = id }
    NR == 1 {
      for (c = 1; c <= NF; c++) col[$c] = c
      idx[1] = col["gene_name"]
      if (!idx[1]) { print "gene_name column not found" > "/dev/stderr"; exit 1 }
      out = $idx[1]
      for (i = 1; i <= n; i++) {
        if (!(want[i] in col)) { print "missing id: " want[i] > "/dev/stderr"; exit 1 }
        idx[i+1] = col[want[i]]
        out = out OFS want[i]
      }
      print out
      next
    }
    {
      line = $idx[1]
      for (i = 2; i <= n + 1; i++) line = line OFS $idx[i]
      print line
    }
  ' "$src" > "$dst"
  echo "[panel e]   $(wc -l < "$dst") rows, $(du -h "$dst" | cut -f1)"
}

slice "$EXPR/combined_scaled_rpm.tsv" cache/panel_e_scaled.tsv
slice "$EXPR/combined_rpm.tsv"        cache/panel_e_raw.tsv
echo "[panel e] done"
