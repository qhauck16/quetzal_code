#!/usr/bin/env bash
# Run pyGenomeTracks for one region directory.
# Usage: run_pygenome.sh <region_dir> <aggregator_path> <bed12_path>
#
# region_dir must already contain (from prep_pygenome.R):
#   junctions.bed, sample_bigwigs.tsv, groups_file.txt, region.txt
set -euo pipefail

DIR=$1
AGGREGATOR=$2
BED12=$3

# Snakemake activates the conda env but the slurm worker can land in a
# context where `python` still resolves to a system binary that doesn't have
# pyBigWig. Pin python (and the rest of the env binaries) to the activated
# env explicitly when CONDA_PREFIX is set.
if [[ -n "${CONDA_PREFIX:-}" && -x "$CONDA_PREFIX/bin/python" ]]; then
  export PATH="$CONDA_PREFIX/bin:$PATH"
  PY="$CONDA_PREFIX/bin/python"
else
  PY="$(command -v python)"
fi
echo "[run_pygenome] python: $PY"

REGION=$(cat "$DIR/region.txt")

# 1. bgzip + tabix the junctions bed
(head -n 1 "$DIR/junctions.bed" && \
 tail -n +2 "$DIR/junctions.bed" | sort -k1,1 -k2,2n) \
  | bgzip -c > "$DIR/junctions.bed.gz"
tabix -f -p bed "$DIR/junctions.bed.gz"

# 2. Build tracks.ini via the GenometracksByGenotype aggregator
"$PY" "$AGGREGATOR" \
  --BigwigList         "$DIR/sample_bigwigs.tsv" \
  --BigwigListType     KeyFile \
  --OutputPrefix       "$DIR/" \
  --Region             "$REGION" \
  --GroupSettingsFile  "$DIR/groups_file.txt" \
  --Bed12GenesToIni    "$BED12" \
  --Normalization      WholeGenome \
  --BedfileForSashimiLinks "$DIR/junctions.bed.gz"

# 3. Full-coverage PDF
pyGenomeTracks --region "$REGION" \
  --tracks "$DIR/tracks.ini" \
  -o "$DIR/simple_coverage.pdf"

# 4. Build sashimi-only ini and render it
SCRIPT_DIR=$(dirname "$0")
"$PY" "$SCRIPT_DIR/make_sashimi_ini.py" \
  "$DIR/tracks.ini" "$DIR/sashimi_tracks.ini"

pyGenomeTracks --region "$REGION" \
  --tracks "$DIR/sashimi_tracks.ini" \
  -o "$DIR/sashimi_simple_coverage.pdf"
