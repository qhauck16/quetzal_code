# figS1/build.R -- Supplementary Figure S1

## Purpose & context
Representative single-gene Quetzal fits. Rows: C2 (k10, prostate), MVK (k6,
kidney), NHEJ1 (k2, colon/rectum) -- structure plot + CDFs of the factor's DE
junctions (TCGA cohort of interest vs matched GTEx tissue); then MAP2K2 (k5)
and KLHDC4 (k6) -- structure plot + read-length boxplot (<=50 vs >=75 bp) of
the top DE junction. Ports `tcga/gene_level_poisson/showcase_genes/supplementary_plots`
(Snakefile + `scripts/make_{structure_plot,cdf_plot,psi_box}.R`, `config.yaml`)
as functions returning ggplots, restyled to the paper theme.

## Decision rationale
- Reads that pipeline's `config.yaml` (paths, thresholds, cohorts, factor per
  gene) instead of duplicating it, so gene choices stay in one place.
- DE junctions = `svalue < sval_thresh & postmean > lfc_thresh & F > abs_thresh`
  (same thresholds as the Snakefile scripts).

## Gotchas & warnings
- Right-hand plots are `wrap_elements(full = ...)`: patchwork otherwise aligns
  them to the structure plot's tall angled tick-label row, leaving a big gap
  between the CDF/box axis and its x title/legend. `free()` crashed
  (patchwork 1.3.0 `simplify_gt.free_table` error) on these rows.
- Read-length structure plots use `gap = 220` (cohort rows 100): tiny groups
  (PRAD-50bp) otherwise put 45-deg labels touching.
- CDF strips are plain-weight and the right column is 1.35/3.65 of the width so
  9-digit coordinate strips fit; facet spacing 1 line so "1.0"/"0.0" don't touch.
- **Known deviation:** read-length structure plots show extra PRAD/KIRP 50-bp
  groups. They are present in the data (the cohort list includes PRAD/KIRP and
  groups are `<cohort>-<len>bp`); the published S1 did not show them.
- `fastTopics::structure_plot` orders samples stochastically -- `set.seed(1)`
  before each call keeps renders diffable.
- No cache: reads per-gene `res.RDS`, `de_res.RDS`, `psi_table.tsv` live (slow-ish, fine on login).

## User preferences
d3 category10/20 topic colours; gene names italic in titles (`gene_md`); angled x labels.

## Input/output
In: `config.yaml` paths -> TCGA/GTEx metadata, `<chr>/<gene>/{res,de_res}.RDS`,
TCGA and GTEx `psi_table.tsv`. Out: `out/figS1.pdf` (8.5 x 11 in).

## File flow schematic
```
supplementary_plots/config.yaml ─┐
<chr>/<gene>/res.RDS, de_res.RDS ─┼──► figS1/build.R ──► out/figS1.pdf
TCGA + GTEx psi_table.tsv ────────┤
TCGA / GTEx metadata ─────────────┘
```
