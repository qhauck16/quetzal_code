#!/bin/bash
#SBATCH --job-name=gdc_dl
#SBATCH --partition=caslake
#SBATCH --account=pi-yangili1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=logs/gdc_dl_%A_%a.log
# One array task per shard (sbatch --array=0-23 --export=TOKEN=... ).
# Re-submittable: files already present at full size are skipped, and
# gdc-client resumes partial downloads.
set -uo pipefail
GDC=/project2/yangili1/qhauck/software/gdc-client/venv/bin/gdc-client
S=shards/shard_${SLURM_ARRAY_TASK_ID}.tsv
[ -s "$S" ] || { echo "no shard $S"; exit 0; }
[ -r "$TOKEN" ] || { echo "token not readable: $TOKEN"; exit 1; }
tail -n +2 "$S" | cut -f6 | sort -u | while read -r dir; do
  mkdir -p "$dir"
  sub=$(mktemp -p "${TMPDIR:-/tmp}" sub.XXXX.tsv)
  head -1 "$S" | cut -f1-5 > "$sub"
  tail -n +2 "$S" | awk -F'\t' -v d="$dir" '$6==d' | while IFS=$'\t' read -r id fn md5 size state d; do
    f="$dir/$id/$fn"
    if [ -s "$f" ] && [ "$(stat -c %s "$f")" = "$size" ]; then echo "skip $id (complete)"; else
      printf "%s\t%s\t%s\t%s\t%s\n" "$id" "$fn" "$md5" "$size" "$state" >> "$sub"; fi
  done
  n=$(( $(wc -l < "$sub") - 1 ))
  echo "[$(date +%T)] $dir: $n file(s) to fetch"
  [ $n -gt 0 ] && "$GDC" download -m "$sub" -t "$TOKEN" -d "$dir" -n "$SLURM_CPUS_PER_TASK" \
      --retry-amount 20 --wait-time 10 --no-annotations --no-related-files < /dev/null
  rm -f "$sub"
done
echo "[$(date +%T)] shard ${SLURM_ARRAY_TASK_ID} done"
