#!/bin/bash
# Download transcriptomic fq.gz files from CNCB for GIST samples
# Skips files that already exist. Sequential downloads to avoid throttling.

BASE_URL="https://download.cncb.ac.cn/gsa-human/HRA005970"
DEST_DIR="$(cd "$(dirname "$0")" && pwd)"
GIST_FILE="$(dirname "$DEST_DIR")/GIST_files_to_download.tsv"
FILES_FILE="$(dirname "$DEST_DIR")/HRA005970_files.tsv"

# Build download list: only transcriptomic fq.gz, skip existing
awk -F'\t' '
  NR==FNR && FNR>1 { hrs[$1]=1; next }
  FNR>1 && ($7 in hrs) && $2 ~ /fq\.gz$/ && $5 ~ /transcriptomic/ { print $7, $1, $2 }
' "$GIST_FILE" "$FILES_FILE" | while read -r hrs hrr filename; do
  outdir="$DEST_DIR/$hrs"
  outfile="$outdir/$filename"
  if [ -f "$outfile" ] && [ -s "$outfile" ]; then
    continue
  fi
  # Skip if decompressed .fq already exists
  fqfile="$outdir/${filename%.gz}"
  if [ -f "$fqfile" ] && [ -s "$fqfile" ]; then
    continue
  fi
  mkdir -p "$outdir"
  echo "$BASE_URL/$hrr/$filename $outfile"
done > /tmp/cncb_download_list.txt

TOTAL=$(wc -l < /tmp/cncb_download_list.txt)
echo "$(date): Files to download: $TOTAL"

if [ "$TOTAL" -eq 0 ]; then
  echo "Nothing to download — all files already exist."
  exit 0
fi

# Download sequentially with retries
COUNT=0
while IFS=' ' read -r url outfile; do
  COUNT=$((COUNT + 1))
  echo "$(date): [$COUNT/$TOTAL] Downloading $(basename "$outfile")..."
  for attempt in 1 2 3; do
    wget -c --tries=3 --read-timeout=600 --connect-timeout=30 -O "$outfile" "$url" 2>&1
    if [ $? -eq 0 ] && [ -s "$outfile" ]; then
      echo "$(date): [$COUNT/$TOTAL] Done: $(basename "$outfile")"
      break
    else
      echo "$(date): Attempt $attempt failed for $(basename "$outfile"), retrying in 10s..."
      sleep 10
    fi
  done
done < /tmp/cncb_download_list.txt

echo "$(date): All downloads complete."
