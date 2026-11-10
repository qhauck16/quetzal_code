#!/bin/bash
set -euo pipefail
# LEAFCUTTER2_DIR: clone of github.com/cfbuenabadn/leafcutter2 at commit 62cccc6; run inside its
# conda env (envs/leafcutter2_env.yml).
# junc_files.txt lists the 330 per-sample regtools .junc files (junc_files/).
cd "$(dirname "${BASH_SOURCE[0]}")"
: "${LEAFCUTTER2_DIR:?set LEAFCUTTER2_DIR to a leafcutter2 clone}"
PY=${PY:-python}
$PY --version
$PY "$LEAFCUTTER2_DIR"/scripts/leafcutter2.py -j junc_files.txt -m 330 -C -o reg --minreadstd 0 -A ../../inputs/gencode.v46.annotation.gtf -L -T
echo CLUSTERDONE
