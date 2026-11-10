# figS3/build.R -- Supplementary Figure S3

## Purpose & context
Cancer-gene mis-splicing from single-gene fits: (a) TCF12 k10 in TCGA-LGG and
TCGA-GBM, (b) SRC k8 in TCGA-COAD, (c) SETD2 k2 in TCGA-LAML. Each panel:
loading violin (cohort vs all other cohorts, Wilcoxon), one PSI CDF per DE
junction (cohort / matched GTEx / other TCGA / other GTEx) and the gene model
with DE junctions as numbered arcs (red `*` = non-coding junction).
Drawing code is main Figure 1d's (EZH2; `../figure1/compose.R`,
`../figure1/panel_d_e/scripts/lib.R`).

## Decision rationale
- Reuses fig 1d's `make_gene_outputs.R` output (`gene_data.RDS`) via
  `cache/s2/run.sh` and fig 1d's `lib.R` (`plot_loading_violins`) so S2 and
  fig 1d look identical.
- `load_gd()` overrides the RDS's stored metadata/PSI paths with `../../paths.yaml`
  (`QP`), since paths baked into the RDS may be stale.

## Gotchas & warnings
- The violin's ggpubr bracket (GeomSignif layer) is removed and redrawn with
  the exact Wilcoxon p (`fmt_p`, size `P_SIZE`, same as CDF p-values); ggpubr
  printed the "p < 2.22e-16" floor at ~10 pt. Fig 1d still uses ggpubr.
- Rows use fixed panel widths in inches (violin 1.35, CDF 1.45) plus a null
  spacer at the right, so CDFs sit next to the violin; CDFs no longer use
  aspect.ratio = 1 (that centred them in an over-wide cell).
- Panel tags need `& theme(plot.tag = bold)`: `plot_annotation_publication()`'s
  theme does not reach the tags of nested wrap_elements (rendered regular).
- Needs `cache/s2/<GENE>_<COHORT>/<GENE>/{gene_data.RDS,junctions_classified.tsv}`
  -- run `cache/s2/run.sh` first.
- p-value shown only if the junction coefficient > 0, else "ns".
- Arc-number x positions are de-overlapped by a simple 20-iteration push loop.
- Non-coding lookup shifts junction start by +1 to match junction naming.

## User preferences
Thin arcs scaled by factor loading; numbers on the arcs; PAL colours from theme.R.

## Input/output
In: cache/s2 gene_data.RDS + junctions_classified.tsv; `../../paths.yaml`
(TCGA/GTEx metadata, psi_root); other-cohort `psi_table.tsv`. Out: `out/figS3.pdf` (8.5 x 11).

## File flow schematic
```
cache/s2/run.sh ──► cache/s2/*/gene_data.RDS ─┐
                    junctions_classified.tsv ─┼──► figS3/build.R ──► out/figS3.pdf
figure1/panel_d_e/scripts/lib.R ──────────────┤
paths.yaml ─► psi_table.tsv, metadata ────────┘
```

**2026-09-24 (Quinn review):** violin+CDF rows centred between flexible spacers over the full-width gene tracks; panel c violin ~square (VIOLIN_C_W 2.1 in, vertically centred).

> 2026-11-10: `paths.yaml` and `scripts/paths.R` were removed. Paths are now
> `file.path(QROOT, "<folder>/<file>")`, with `QROOT` the quetzal_code root
> (found via the `.quetzal_root` marker); references to paths.yaml keys above are historical.
