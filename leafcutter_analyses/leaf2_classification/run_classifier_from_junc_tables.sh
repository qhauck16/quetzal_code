#!/usr/bin/env bash
# Run inside the leafcutter2 conda env (../leaf2_clustering/envs/leafcutter2_env.yml).
set -euo pipefail

######################################################################
# Pulls junctions from junction_tables_for_psi/chr*/gene/dat.tsv,
# formats them into leafcutter2 perind style, then runs
# ForwardSpliceJunctionClassifier.py
#
# Skips MAN1B1 (hangs the classifier).
######################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
JUNC_DIR="${SCRIPT_DIR}/../junction_tables"
GTF="${SCRIPT_DIR}/../../inputs/gencode.v46.annotation.gtf"
GENOME="${SCRIPT_DIR}/../../inputs/GRCh38.primary_assembly.genome.fa"
OUTFILE="${SCRIPT_DIR}/junc_table_all_junctions.tsv.gz"
RUNDIR="${SCRIPT_DIR}"
SKIP_GENES="MAN1B1"

###############################################
# 1. Build gene -> strand map from GTF
###############################################
echo "Building gene-strand map from GTF..."
declare -A GENE_STRAND
while IFS=$'\t' read -r chr _ type start end _ strand _ info; do
    if [[ "$type" != "gene" ]]; then continue; fi
    if [[ "$chr" == "chrY" ]]; then continue; fi
    gname=$(echo "$info" | grep -oP 'gene_name "\K[^"]+')
    GENE_STRAND["$gname"]="$strand"
done < <(grep -P '^\w+\t\w+\tgene\t' "$GTF")
echo "  Found ${#GENE_STRAND[@]} genes in GTF."

###############################################
# 2. Extract junctions from dat.tsv headers
###############################################
echo "Extracting junctions from junction_tables_for_psi..."
clu_counter=0
tmpfile=$(mktemp)

for chr_dir in "$JUNC_DIR"/chr*; do
    chrom=$(basename "$chr_dir")
    [[ -d "$chr_dir" ]] || continue

    for gene_dir in "$chr_dir"/*/; do
        gene=$(basename "$gene_dir")
        dat="$gene_dir/dat.tsv"
        [[ -f "$dat" ]] || continue

        # Skip problematic genes
        if [[ "$gene" == "$SKIP_GENES" ]]; then
            echo "  Skipping $gene"
            continue
        fi

        # Look up strand from GTF
        strand="${GENE_STRAND[$gene]:-+}"

        clu_counter=$((clu_counter + 1))

        if (( clu_counter % 20 == 0 )); then
            echo "  [$clu_counter] Processing $gene ($chrom)"
        fi

        # Parse header: tokens are "start-end" (quoted, 1-based snaptron intron coords)
        # start-1 converts intron_start to last exon base (0-based donor)
        # end stays as-is: classifier adds 1 internally to get first base of next exon (1-based acceptor)
        head -n 1 "$dat" \
        | tr ' ' '\n' \
        | sed 's/"//g' \
        | awk -F '-' -v chrom="$chrom" -v clu="$clu_counter" -v strand="$strand" '
            NF==2 && $1 ~ /^[0-9]+$/ && $2 ~ /^[0-9]+$/ {
                print chrom ":" ($1-1) ":" $2 ":clu_" clu "_" strand
            }
        ' >> "$tmpfile"
    done
done

# Deduplicate and sort
sort -u "$tmpfile" > "${tmpfile}.sorted"
n_juncs=$(wc -l < "${tmpfile}.sorted")
echo "  Extracted $n_juncs unique junctions."

###############################################
# 3. Compress into perind-style file
###############################################
echo "Writing $OUTFILE ..."
# Add a header line so the classifier skips it, then junction lines
(echo "chrom dummy"; cat "${tmpfile}.sorted") | gzip > "$OUTFILE"
rm -f "$tmpfile" "${tmpfile}.sorted"

###############################################
# 4. Run ForwardSpliceJunctionClassifier
###############################################
echo "Running ForwardSpliceJunctionClassifier.py ..."
mkdir -p "${RUNDIR}/clustering"

python3 "${SCRIPT_DIR}/ForwardSpliceJunctionClassifier.py" \
    -c "$OUTFILE" \
    -A "$GTF" \
    -G "$GENOME" \
    -r "$RUNDIR" \
    -o Leaf2 \
    --max_juncs 5000

echo "Done! Output in ${RUNDIR}/clustering/"
