#!/bin/bash
# Check every manifest file is present with the exact GDC size (and md5 with --md5).
# Run from gdc_parallel/:  bash verify.sh [--md5]
ok=0; bad=0
while IFS=$'\t' read -r id fn md5 size state dir; do
  [ "$id" = id ] && continue
  f="$dir/$id/$fn"
  if [ ! -s "$f" ] || [ "$(stat -c %s "$f")" != "$size" ]; then echo "MISSING/SHORT $f"; bad=$((bad+1)); continue; fi
  if [ "${1:-}" = --md5 ] && [ "$(md5sum "$f" | cut -c1-32)" != "$md5" ]; then echo "MD5 BAD $f"; bad=$((bad+1)); continue; fi
  ok=$((ok+1))
done < manifest.tsv
echo "complete: $ok   problems: $bad"
