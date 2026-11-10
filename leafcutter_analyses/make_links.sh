#!/usr/bin/env bash
# Midway only. Points every output path of leafcutter_analyses/ at the result
# that already exists under the original spliceNMF/tcga/ tree, so the figure code
# reads the paper's actual files without re-running anything. Anyone else either
# runs the pipelines here (outputs land at these same paths) or unpacks the
# deposited outputs into them.
#
#   bash make_links.sh            # create / refresh links (idempotent)
#
# Never overwrites a real file or directory: a path that exists and is not a
# symlink is reported and skipped.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
TCGA=${TCGA:-/project2/yangili1/qhauck/spliceNMF/tcga}
DS=$TCGA/genome_wide/flashier/diff_splicing
L2=$TCGA/leafcutter/leaf2_try

link() {  # link <target> <path-here>
  local tgt=$1 here=$2
  [ -e "$tgt" ] || { echo "MISSING target: $tgt" >&2; return 1; }
  if [ -e "$here" ] && [ ! -L "$here" ]; then echo "SKIP (real file): $here" >&2; return 0; fi
  mkdir -p "$(dirname "$here")"; ln -sfn "$tgt" "$here"
}

# junction_tables/: per-gene dat.tsv from Snaptron (junction_tables/Snakefile)
for d in "$TCGA"/junction_tables_for_psi/chr*/; do link "${d%/}" "junction_tables/$(basename "$d")"; done

# junction_tables_gtex/: GTEx per-gene dat.tsv (junction_tables_gtex/setting_up_snakemake.sh)
GT=${GTEX:-/project2/yangili1/qhauck/spliceNMF/gtex}/genome_wide/260120_filters
for d in "$GT"/chr*/; do link "${d%/}/FastTopics_output" "junction_tables_gtex/$(basename "$d")/FastTopics_output"; done

# leaf2_clustering/: 330 regtools .junc inputs + leafcutter2 cluster outputs
link "$TCGA/leafcutter/junc_files"                   leaf2_clustering/junc_files
link "$L2/reg.junction_counts.unclassified.gz"        leaf2_clustering/reg.junction_counts.unclassified.gz
link "$L2/reg.cluster_ratios.unclassified.gz"         leaf2_clustering/reg.cluster_ratios.unclassified.gz

# psi_files/: per-gene PSI tables (psi/Snakefile)
link "$L2/new_psi_files"                              psi_files

# leaf2_classification/: junction classifications (run_classifier_from_junc_tables.sh)
link "$TCGA/leafcutter/leafcutter2/clustering"        leaf2_classification/clustering
link "$TCGA/leafcutter/leafcutter2/junc_table_all_junctions.tsv.gz" leaf2_classification/junc_table_all_junctions.tsv.gz

# diff_splicing/: regtools .junc cache, LeafCutter runs, loadings, intron table
for x in junc_files results control_results download flashier_loadings.tsv gencode_v46_introns.tsv; do
  link "$DS/$x" "diff_splicing/$x"
done

# diff_splicing/unified_results/: per-gene threshold sweeps + shared background tables
for x in nkap sf3b1 ylpm1 shared; do link "$DS/unified_results/$x" "diff_splicing/unified_results/$x"; done

echo "leafcutter_analyses links OK"
