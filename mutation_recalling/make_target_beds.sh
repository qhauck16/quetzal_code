#!/usr/bin/env bash
# Target exon BEDs for the WXS / RNA-seq re-calling (unified/config.yaml: target_beds).
# Same commands as the manuscript's {nkap,ylpm1,sf3b1}/making_bed.sh: GENCODE v46 exons
# of each gene set, sorted and merged with bedtools. The awk patterns are kept verbatim.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
GTF=../inputs/gencode.v46.annotation.gtf
mkdir -p targets
for spec in "nkap:NKAP|SLU7|FAM50A|PRPF8" "ylpm1:YLPM1" "sf3b1:SF3B1|SUGP1"; do
  name=${spec%%:*}; GENES=${spec#*:}
  awk -v genes="$GENES" '
      $3=="exon" && $0 ~ "gene_name \""genes"\"" {
          OFS="\t";
          print $1, $4-1, $5
      }
  ' $GTF | bedtools sort -i - | bedtools merge > targets/$name.target_exons.merged.bed
done
echo "done"
