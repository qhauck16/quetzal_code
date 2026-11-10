# build_gsp_top.R

## Purpose & context
Builds the top-250 driving unannotated junctions, with splice types, for any GSP. It is a
verbatim port of the `spec_flash_fact.Rmd` chunk that wrote the published
`fact_{10,13,30}_top_junctions.tsv`. Needed for GSP #1 (IDH1), GSP #12 (SPOP) and
GSP #40 (exploratory) in the IDH1/SPOP supplementary note.

## Decision rationale
- **Ported, not reused.** `tcga/.../flashier/fact_1_top_junctions.tsv` is stale: 745 rows,
  from April 2026, built with different settings. Every GSP is rebuilt here from the
  current RDS.
- **Validated.** Running GSP 10 reproduces the published `fact_10_top_junctions.tsv`
  exactly.

## Gotchas & warnings
- Loads the 3.7 GB flashier RDS, so run it on `amd` with 64 GB (~10 min for 4 GSPs).
- The many-to-many join warnings also appear in the original Rmd, and duplicated rows are
  removed by its `n() < 2` filter.
- GSP #40 was added because SPOP mutants separate on k40 more strongly (Wilcoxon) than on
  k12. Its median loading is 0 in both arms.

## Input/output
- **In:** the flashier RDS, `260116_filters/<chr>/FastTopics_output/<gene>/res.RDS`, the
  snaptron tables, and the gencode/MANE introns.
- **Out:** `gsp_top/fact_<N>_top_junctions.tsv`.

```
flashier RDS + per-gene res.RDS + snaptron ─► build_gsp_top.R <N...> ─► gsp_top/fact_<N>_top_junctions.tsv ─► overlap.R
```
