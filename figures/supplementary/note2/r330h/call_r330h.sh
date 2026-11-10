#!/bin/bash
# Re-call NKAP p.R330H in TCGA-AB-2932 (LAML; rail 857257) from the GDC BAMs:
# tumour WXS (03A), matched normal WXS (11A) and tumour RNA-seq (03A).
# NKAP is on the minus strand; p.R330 codon = chrX:119930101-119930103 (hg38;
# R330C in the call table is chrX:119930101 G>A), so R330H (CGC>CAC, codon
# base 2) is chrX:119930100 C>T on the + strand. Counts bases with
# samtools mpileup (MAPQ >= 20, BASEQ >= 20, no BAQ), plus R330C's site.
# Run after the download (manifest.tsv) finishes; verifies md5 first.
set -euo pipefail
cd "$(dirname "$0")"
S=${SAMTOOLS:-samtools}
while IFS=$'\t' read -r id fn md5 size state; do
  [ "$id" = id ] && continue
  f=bams/$id/$fn
  [ "$(stat -c %s "$f")" = "$size" ] || { echo "SIZE MISMATCH $f"; exit 1; }
  [ "$(md5sum "$f" | cut -c1-32)" = "$md5" ] || { echo "MD5 MISMATCH $f"; exit 1; }
  [ -s "$f.bai" ] || $S index -@ 4 "$f"
  sm=$($S view -H "$f" | grep -m1 '^@RG' | tr '\t' '\n' | grep '^SM:' | cut -c4-)
  for pos in 119930100 119930101; do
    $S mpileup -r chrX:$pos-$pos -q 20 -Q 20 -B -d 100000 "$f" 2>/dev/null | \
      awk -v sm="$sm" -v pos="$pos" '{b=toupper($5); gsub(/\^./,"",b); gsub(/\$/,"",b);
        n=split("A C G T",B," "); out=""; for(i=1;i<=n;i++){c=gsub(B[i],"",b); out=out" "B[i]"="c}
        r=$5; ref=gsub(/[.,]/,"",r); print sm, "chrX:"pos, "ref("$3")="ref, out, "depth="$4}'
  done
done < manifest.tsv | tee r330h_pileup.txt
