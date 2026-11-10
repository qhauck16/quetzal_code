#!/usr/bin/env bash
# AlphaMissense subsets used downstream (run from mutation_processing/, after
# liftover_and_merge.R and before join_alphamissense.R):
#   missense_keys.tsv              unique (chr, pos, ref, alt) of missense SNVs in merged_grch38.tsv,
#                                  in first-seen order
#   am_hits.tsv                    AlphaMissense hg38 rows at those loci (AlphaMissense file order)
#   am_refs/<GENE>_alphamissense.tsv  every AlphaMissense row for one protein (UniProt ID), with header;
#                                  read by ../figures/scripts/variant_tables.R
# These were first made with one-off shell commands; this script reproduces them
# byte-for-byte from the same inputs.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
AM=../inputs/alphamissense/AlphaMissense_hg38.tsv.gz

# merged_grch38.tsv columns: 12 startPosition, 14 referenceAllele, 15 variantAllele,
# 17 mutationType, 21 chr
awk -F'\t' 'NR>1 && $17=="Missense_Mutation" && length($14)==1 && length($15)==1 {
              k=$21"\t"$12"\t"$14"\t"$15; if (!(k in s)) { s[k]=1; print k } }' \
    merged_grch38.tsv > missense_keys.tsv

zcat "$AM" | awk -F'\t' 'NR==FNR { k[$1"\t"$2"\t"$3"\t"$4]=1; next } /^#/ { next }
                         { c=$1; sub(/^chr/, "", c); if ((c"\t"$2"\t"$3"\t"$4) in k) print }' \
    missense_keys.tsv - > am_hits.tsv

# gene -> UniProt accession for the per-protein tables
mkdir -p am_refs
zcat "$AM" | awk -F'\t' -v OFS='\t' '
  BEGIN { split("FAM50A:Q14320 NKAP:Q8N5F7 PRPF8:Q6P2Q9 SF3B1:O75533 SLU7:O95391 SUGP1:Q8IWZ8 YLPM1:P49750", g, " ")
          for (i in g) { split(g[i], x, ":"); gene[x[2]] = x[1] } }
  /^#CHROM/ { sub(/^#/, ""); hdr = $0; next }
  /^#/ { next }
  ($6 in gene) { f = "am_refs/" gene[$6] "_alphamissense.tsv"
                 if (!(f in started)) { print hdr > f; started[f] = 1 }
                 print > f }'
echo "done"
