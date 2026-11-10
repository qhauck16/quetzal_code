#!/usr/bin/env bash
# Fetch the single-file public inputs (sample tables, GENCODE, MANE, AlphaMissense,
# PanCanAtlas ABSOLUTE, Quetzal v0.1.1 data files) from their original sources into
# this directory, at the paths ../paths.yaml and the analysis folders expect.
#
#   bash fetch_inputs.sh                 # all single-file inputs
#   bash fetch_inputs.sh --only <step>   # one of: metadata quetzal_data gencode mane alphamissense
#                                        #   gdc_reference vep_cache star_index alsafadi pancanatlas
#   bash fetch_inputs.sh --midway        # midway only: symlink the copies already on disk
#                                        #   (also links the Snaptron tables, bigwigs and mutation downloads)
#
# Separate scripts here: fetch_snaptron.sh (per-gene Snaptron junction tables) and
# fetch_bigwigs.sh (recount3 coverage bigwigs). The cBioPortal mutation download is the
# first step of ../mutation_processing/.
#
# Two inputs cannot be downloaded by script and are placed by hand:
#   cosmic_tier_1.tsv          COSMIC Cancer Gene Census, Tier 1 (export of 2025-11-03; needs a
#                              COSMIC account: https://cancer.sanger.ac.uk/census)
#   cola_seq_branchpoints.csv  CoLa-seq human branchpoints (Zeng et al., bioRxiv
#                              10.1101/2021.10.18.464728)
#
# Needs: curl, gzip.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
IN=$PWD

MIDWAY=0; ONLY=""
while [ $# -gt 0 ]; do
  case $1 in
    --midway)  MIDWAY=1 ;;
    --only)    ONLY=$2; shift ;;
    *) echo "unknown option: $1" >&2; exit 1 ;;
  esac; shift
done
want() { [ -z "$ONLY" ] || [ "$ONLY" = "$1" ]; }
get()  { [ -s "$2" ] && { echo "  have $2"; return; }; mkdir -p "$(dirname "$2")"; echo "  get  $2"
         curl -fsSL --retry 3 -o "$2.part" "$1" && mv "$2.part" "$2"; }
getgz(){ [ -s "$2" ] && { echo "  have $2"; return; }; get "$1" "$2.gz"; gunzip "$2.gz"; }

QUETZAL_RAW=https://raw.githubusercontent.com/qhauck16/Quetzal/v0.1.1/data
GENCODE=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46

# ---------------------------------------------------------------------------
if [ $MIDWAY = 1 ]; then
  S=${SPLICENMF:-/project2/yangili1/qhauck/spliceNMF}; T=$S/tcga
  lnk() { [ -e "$1" ] || { echo "MISSING $1" >&2; return 1; }
          if [ -e "$2" ] && [ ! -L "$2" ]; then echo "SKIP (real file) $2" >&2; return 0; fi
          mkdir -p "$(dirname "$2")"; ln -sfn "$1" "$2"; }
  lnk $T/tcga_v2_samples.tsv                                  tcga_v2_samples.tsv
  lnk $S/gtexv2_metadata.tsv                                   gtexv2_metadata.tsv
  for d in $T/all_genes/chr*/;   do lnk "${d%/}" snaptron/tcga/$(basename "$d"); done
  for d in $S/gtex/all_genes/chr*/; do lnk "${d%/}" snaptron/gtex/$(basename "$d"); done
  lnk $T/leafcutter/leafcutter2/gencode.v46.annotation.gtf     gencode.v46.annotation.gtf
  lnk $T/leafcutter/leafcutter2/genome.fa                      GRCh38.primary_assembly.genome.fa
  lnk $T/leafcutter/leafcutter2/genome.fa.fai                  GRCh38.primary_assembly.genome.fa.fai
  W=$T/genome_wide/flashier/wxs_mut_calling
  lnk $W/GRCh38.d1.vd1.fa                                      GRCh38.d1.vd1.fa
  lnk $W/GRCh38.d1.vd1.fa.fai                                  GRCh38.d1.vd1.fa.fai
  lnk $W/GRCh38.d1.vd1.dict                                    GRCh38.d1.vd1.dict
  lnk /project2/yangili1/qhauck/conda/.vep                     vep_cache
  lnk $T/reference_genome/STAR_index/star-2.7.5c_GRCh38.d1.vd1_gencode.v36 star_index/star-2.7.5c_GRCh38.d1.vd1_gencode.v36
  lnk "$T/genome_wide/flashier/41467_2016_BFncomms10615_MOESM961_ESM (1).xlsx" Alsafadi2016_supp1.xlsx
  lnk $T/hg38_granges.RDS                                      gencode_v46_granges.RDS
  lnk $T/all_genes/hugo_all_pc_genes.tsv                       hugo_all_pc_genes.tsv
  lnk $T/gene_level_poisson/analyte.tsv                        analyte.tsv
  for f in $T/genome_wide/flashier/just_for_fun/nmd_investigation/chr*_productive_unproductive.tsv; do
    lnk "$f" productive_unproductive/$(basename "$f"); done
  lnk $T/genome_wide/flashier/diff_splicing/MANE.GRCh38.v1.4.summary.txt MANE.GRCh38.v1.4.summary.txt
  lnk $T/genome_wide/flashier/cola_seq_branchpoints.csv        cola_seq_branchpoints.csv
  lnk $T/gene_level_poisson/cosmic_tier_1_11_3_2025.tsv        cosmic_tier_1.tsv
  M=$T/genome_wide/flashier/just_for_fun/mutation_merge/refs
  lnk $M/AlphaMissense_hg38.tsv.gz                             alphamissense/AlphaMissense_hg38.tsv.gz
  lnk $M/hg19ToHg38.over.chain.gz                              alphamissense/hg19ToHg38.over.chain.gz
  for f in TCGA_mastercalls.abs_segtabs.fixed.txt TCGA_mastercalls.abs_tables_JSedit.fixed.txt; do
    lnk $T/pancanatlas/$f pancanatlas/$f; done
  lnk $T/recount3_bigwigs/tcga_bws                             bigwigs/tcga
  lnk $T/recount3_bigwigs/gtex_bws                             bigwigs/gtex
  echo "inputs linked (midway)"; exit 0
fi

# ---------------------------------------------------------------------------
if want metadata; then echo "== Snaptron sample metadata"
  get https://snaptron.cs.jhu.edu/data/tcgav2/samples.tsv tcga_v2_samples.tsv
  get https://snaptron.cs.jhu.edu/data/gtexv2/samples.tsv gtexv2_metadata.tsv
fi

if want quetzal_data; then echo "== Files shipped with Quetzal v0.1.1"
  get $QUETZAL_RAW/hg38_granges.RDS       gencode_v46_granges.RDS
  get $QUETZAL_RAW/hugo_all_pc_genes.tsv  hugo_all_pc_genes.tsv
  get $QUETZAL_RAW/analyte.tsv            analyte.tsv
  for c in $(seq 1 22) X; do
    get $QUETZAL_RAW/productive_unproductive/chr${c}_productive_unproductive.tsv \
        productive_unproductive/chr${c}_productive_unproductive.tsv
  done
fi

if want gencode; then echo "== GENCODE release 46 (GTF, GRCh38 primary assembly)"
  getgz $GENCODE/gencode.v46.annotation.gtf.gz             gencode.v46.annotation.gtf
  getgz $GENCODE/GRCh38.primary_assembly.genome.fa.gz      GRCh38.primary_assembly.genome.fa
fi

if want mane; then echo "== MANE v1.4"
  getgz https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.4/MANE.GRCh38.v1.4.summary.txt.gz \
        MANE.GRCh38.v1.4.summary.txt
fi

if want alphamissense; then echo "== AlphaMissense + hg19->hg38 chain"
  get https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz alphamissense/AlphaMissense_hg38.tsv.gz
  get https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz alphamissense/hg19ToHg38.over.chain.gz
fi

if want gdc_reference; then echo "== GDC GRCh38.d1.vd1 reference (WXS / RNA-seq re-calling)"
  [ -s GRCh38.d1.vd1.fa ] || { get https://api.gdc.cancer.gov/data/254f697d-310d-4d7d-a27b-27fbf767a834 GRCh38.d1.vd1.fa.tar.gz
                               tar -xzf GRCh38.d1.vd1.fa.tar.gz && rm GRCh38.d1.vd1.fa.tar.gz; }
fi

if want vep_cache; then echo "== Ensembl VEP RefSeq cache, release 115, GRCh38 (~25 GB)"
  [ -d vep_cache/homo_sapiens_refseq ] || { get https://ftp.ensembl.org/pub/release-115/variation/indexed_vep_cache/homo_sapiens_refseq_vep_115_GRCh38.tar.gz vep_cache/homo_sapiens_refseq_vep_115_GRCh38.tar.gz
                                            tar -xzf vep_cache/homo_sapiens_refseq_vep_115_GRCh38.tar.gz -C vep_cache; }
fi

if want star_index; then echo "== GDC STAR index, star-2.7.5c GRCh38.d1.vd1 + GENCODE v36 (~25 GB; GIST alignment)"
  [ -d star_index/star-2.7.5c_GRCh38.d1.vd1_gencode.v36 ] || { get https://api.gdc.cancer.gov/data/c0008693-0583-4eac-bd5c-583070763893 star_index/star-2.7.5c_GRCh38.d1.vd1_gencode.v36.tgz
                                                             tar -xzf star_index/star-2.7.5c_GRCh38.d1.vd1_gencode.v36.tgz -C star_index; }
fi

if want alsafadi; then echo "== Alsafadi et al. 2016 (Nat. Commun. 10.1038/ncomms10615) Supplementary Data 1"
  get 'https://static-content.springer.com/esm/art%3A10.1038%2Fncomms10615/MediaObjects/41467_2016_BFncomms10615_MOESM961_ESM.xlsx' Alsafadi2016_supp1.xlsx
fi

if want pancanatlas; then echo "== PanCanAtlas ABSOLUTE (GDC open access; hg19)"
  get https://api.gdc.cancer.gov/data/0f4f5701-7b61-41ae-bda9-2805d1ca9781 pancanatlas/TCGA_mastercalls.abs_segtabs.fixed.txt
  get https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5 pancanatlas/TCGA_mastercalls.abs_tables_JSedit.fixed.txt
fi

for f in cosmic_tier_1.tsv cola_seq_branchpoints.csv; do
  [ -e "$f" ] || echo "NOTE: place $f by hand (see the header of this script)"
done
echo "done"
