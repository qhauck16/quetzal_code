#!/usr/bin/env bash
# Per-gene Snaptron junction tables for TCGA (tcgav2) and GTEx (gtexv2), queried
# exactly as for the manuscript: one `qs --region <GENE> --datasrc <src>` call per
# protein-coding gene in hugo_all_pc_genes.tsv (chr1-22, X).
#
#   REQUIRES the Snaptron command-line client `qs` from
#   https://github.com/ChristopherWilks/snaptron-experiments (manuscript: commit d03371d),
#   run under Python 2.7 (conda env: envs/snaptron.yml). Point QS at it:
#
#     git clone https://github.com/ChristopherWilks/snaptron-experiments
#     conda env create -f envs/snaptron.yml && conda activate python2.7
#     QS=$PWD/snaptron-experiments/qs bash fetch_snaptron.sh
#
# Writes snaptron/{tcga,gtex}/chr<N>/snaptron_output/<GENE>_snaptron.tsv (skips files
# already present; ~19.5k genes per source). hugo_all_pc_genes.tsv comes from
# fetch_inputs.sh (shipped with Quetzal v0.1.1).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
: "${QS:?set QS to the snaptron-experiments qs client (see header)}"
[ -s hugo_all_pc_genes.tsv ] || { echo "run fetch_inputs.sh --only quetzal_data first" >&2; exit 1; }

# gene -> chromosome from the HGNC cytoband, as tcga/all_genes/hugo_genes_by_chr.R
awk -F'\t' 'NR>1 { c=$2; sub(/[pq].*/, "", c); if (c == "X" || (c ~ /^[0-9]+$/ && c+0 >= 1 && c+0 <= 22)) print "chr" c "\t" $1 }' \
    hugo_all_pc_genes.tsv > snaptron_genes.tsv

for src in tcga gtex; do
  while IFS=$'\t' read -r chr gene; do
    out=snaptron/$src/$chr/snaptron_output/${gene}_snaptron.tsv
    [ -e "$out" ] && continue
    mkdir -p "$(dirname "$out")"
    "$QS" --region "$gene" --datasrc ${src}v2 > "$out" || true
  done < snaptron_genes.tsv
done
echo "done"
