#!/usr/bin/env bash
# recount3 TCGA gene-level counts (gene_sums, GENCODE v26 "G026"), one file per cohort,
# into tcga_counts/ (URLs: tcga_recount3_urls.txt). Then: Rscript gen_count_file.R
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
mkdir -p tcga_counts/rpm_files
while read -r url; do
  out=tcga_counts/$(basename "$url")
  if [ ! -s "$out" ]; then curl -fsSL --retry 3 -o "$out.part" "$url" && mv "$out.part" "$out"; fi
done < tcga_recount3_urls.txt
echo "done"
