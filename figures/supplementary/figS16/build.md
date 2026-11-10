# figS16/ (fetch.sh + build.R + build_full.R) -- Supplementary Fig. 16: YLPM1 splice-acceptor variant

## Purpose & context
An IGV-style figure of the YLPM1 splice-site variant in TCGA-CJ-4638 (TCGA-KIRC; rail 869984;
rank 2 on GSP #30). It replaces Quinn's IGV screenshots.

- **The variant:** chr14:74812626 A>G (`X1783_splice` in the call table; 50 alt / 35 ref in
  WXS). It is the A of the AG acceptor of MANE exon 11 (ENST00000325680.12; exon 11 starts at
  74812628), so the acceptor becomes GG.
- **The consequence in RNA:**
  - the exon 10→11 junction falls to 2 reads (flanking junctions 9–14);
  - intron 10 is partly retained (mean depth 2.8 vs 0.2 in intron 11; flanking exons ~10.5);
  - all 9 RNA reads covering the site carry G.
- **Tumour WGS:** 33 G / 40 A at the site.

## Two versions
- **`build.R`** (main, `figS16.pdf`) is one figure in the style of
  `/project/yangili1/qhauck/kinnex/fig_zoom/MTIF2/fig_MTIF2_slide.pdf`, as Quinn asked on
  2026-10-26:
  - no titles or panel letters; track labels on the left;
  - grey tracks and exons, with no colour on exons or regions;
  - red only for the variant;
  - junction arcs all the same height (thickness ∝ reads, count on top);
  - squished reads;
  - a callout from exon 11's acceptor to the reference / variant sequence
    (intron lowercase | exon uppercase).
- **`build_full.R`** (`figS16_full.pdf`) is the earlier five-panel version,
  including the tumour WGS reads (33 G / 40 A) and the whole-gene track.

## Decision rationale
- **Reads come from GDC BAM slices, not full BAMs.** The tumour WGS file is 295 GB.
  `fetch.sh` uses `/slicing/view` with the token file's path inside the curl header; the
  contents are never printed. Quinn approved this use on 2026-10-26.
- **Panels of the five-panel version (build_full.R), top to bottom:**
  - (a) RNA coverage with junction arcs (thickness ∝ reads; exon 10→11 in red);
  - (b) squished RNA reads, variant base coloured;
  - (c) MANE exons 9–12 with the variant;
  - (d) tumour WGS ±40 bp with the reference row and the AG acceptor boxed; only non-reference
    bases are coloured, as in IGV;
  - (e) the whole YLPM1 MANE gene, with the window boxed.

  This is the order Quinn asked for: RNA on top, then the zoom and WGS, then the gene track.
- **Base colours:** A green, C blue, G orange, T red.

## Gotchas & warnings
- **The sequence callout is drawn in genomic data units,** with the letters spread across
  the window inside the gene panel, so the callout lines are plain segments and no
  cross-panel drawing is needed.
- **Use ASCII `>` for the intron chevrons.** The `›` glyph renders as a slash in Nimbus.
- **The reads panel must keep clipping on.** With `clip = "off"`, reads ran past the window.
- **Window views use `coord_cartesian()`.** `scale_x(limits=)` silently dropped reads and arcs
  crossing the window edges.
- **Panel a needs explicit y headroom,** otherwise the arc counts are clipped.
- **Per-read base:** `stackStringsFromGAlignments()` returns only reads overlapping the
  site, so reads are first filtered to those whose aligned blocks (not splice gaps) cover it.
- **Bioconductor masks dplyr verbs** (`count`, `filter`, `select`, …), so they are re-bound
  after `library()`.
- **Reference sequence:** `/project2/yangili1/yangili/chRNA/GRCh38.primary_assembly.genome.fa`
  (has `.fai`).
- **No matched normal is shown.** Only the tumour WGS was sliced; the normal WGS (11A,
  b57022b7…) could be added the same way.

## Input/output
- **In:** `bams/rna_YLPM1.bam` and `bams/wgs_tumour_site.bam` from `fetch.sh`,
  `tcga/hg38_granges.RDS`, and the GRCh38 fasta.
- **Out:** `supplementary/out/figS16/figS16.pdf` (build.R) and
  `figS16_full.pdf` (build_full.R). Run from `supplementary/`.

```
GDC slicing API ─► fetch.sh ─► bams/rna_YLPM1.bam ────┐
                            └► bams/wgs_tumour_site.bam ┼─► build.R ─► out/figS16/figS16.pdf
hg38_granges.RDS, GRCh38 fasta ─────────────────────────┘
```
