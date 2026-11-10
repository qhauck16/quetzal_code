#!/bin/bash
# Build every supplementary figure (figS<N>/build.R -> out/figS<N>.pdf) in
# numeric order. Continues past failures; prints a pass/fail summary.
#
#   cd quetzal_code/figures/supplementary && ./build_all.sh
#
# Every build.R sources "common.R" and reads "cache/..." by RELATIVE path,
# so this must run from supplementary/ (the script cd's there itself).
#
# PREREQUISITE: the caches in cache/ must exist first. The figure scripts
# only read them. Builders, with their working-directory requirement
# (heavy ones as SLURM jobs):
#   sbatch --mem=<N>G --wrap "<command>"
#
#   cache/build_mutation_cache.R     S10, S12   any cwd (absolute paths); --mem=24G
#   cache/build_overlap_cache.R      S6, S9    any cwd (absolute paths)
#   cache/build_leafcutter_counts.R  S8        any cwd (setwd()s itself)
#   cache/build_readlength_cache.R   S4        any cwd (setwd()s itself);
#                                              loads the 3.7 GB flashier RDS -> SLURM
#   cache/build_s5_heatmap_cache.R   S5        any cwd (absolute paths); flashier RDS +
#                                              536 M mutation table -> SLURM, tens of GB
#   cache/build_s14_cache.R          S15       any cwd (absolute paths); bigwig reads -> SLURM
#   cache/build_s16_cache.py         S17       python3.12 env (pyBigWig); imports cache/s16_lib/
#   cache/s2/run.sh                  S3        any cwd (cd's itself to figure1/panel_d_e)
#   note2/: build_cache.R + fig1.R (Supplementary Note 2), run from supplementary/
#
#   cache/build_ylpm1_expression_row.sh  (any dir; streams the 7.7 GB matrix) -> ylpm1_expression_row.tsv (S18)
# Also read directly (not cached here): figure1/panel_b/recovery_cache (S2),
# figure3/cache/panel_e_scaled.tsv (S14), and variant/phenotype tables (S7, S11, S13, S18, S20).
#   figS16 (YLPM1 splice variant): run figS16/fetch.sh first (GDC slicing; needs the token path).
RSCRIPT=Rscript
cd "$(dirname "$(readlink -f "$0")")" || exit 1
mkdir -p out/logs

pass=(); fail=()
for n in $(seq 1 20); do
  f="figS$n/build.R"
  if [ ! -f "$f" ]; then echo "[skip] $f missing"; fail+=("S$n(missing)"); continue; fi
  echo "[build] $f"
  if "$RSCRIPT" --vanilla "$f" > "out/logs/figS$n.log" 2>&1; then
    pass+=("S$n")
  else
    echo "  FAILED -- see out/logs/figS$n.log"; tail -5 "out/logs/figS$n.log" | sed 's/^/    /'
    fail+=("S$n")
  fi
done

echo
echo "PASS (${#pass[@]}): ${pass[*]}"
echo "FAIL (${#fail[@]}): ${fail[*]}"
[ ${#fail[@]} -eq 0 ]
