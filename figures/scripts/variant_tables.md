# variant_tables.R

**Purpose & context** — Builds the per-gene tables every NKAP / SF3B1 / YLPM1 figure
reads: the GSP loading per TCGA sample, catalogued mutations, and self-called
WXS / RNA-seq variants. Ported 2026-10-05 from
`tcga/genome_wide/flashier/flash_fact_char.Rmd` so quetzal_code is self-contained.
The functions mirror these Rmd chunks:

| Function | Rmd chunk |
|---|---|
| `vt_inputs()` | "loading re-used files", "correlating flashier factors with mutation status", "per-gene AlphaMissense lookup" |
| `process_vep_vcf()` | "Function to process a single VEP VCF" (verbatim) |
| `build_sf3b1_tables()` | "SF3B1 cdf", "SUGP1 CDFs", "adding variants to SF3B1 df" |
| `build_nkap_tables()` | "processing variants" |
| `build_ylpm1_tables()` | "YLPM1 parsing vcf" |

**Decision rationale** — The Rmd's three `extract_*_variants()` differed only in their
transcripts, region filters and PTC-vs-missense mode, so they became one
`extract_called_variants()` with per-gene arguments. Everything else is the Rmd's code
unchanged, quirks included (e.g. SF3B1 tests `gene_mutation != 'NA'`, a string), so the
outputs reproduce the Rmd's exactly. Verified 2026-10-05: all 7 tables were cell-for-cell
identical to the Rmd's outputs.

**Gotchas & warnings**
- Loadings are `flashier::ldf(out, type = "i")$L`, cached in
  `shared_cache/gsp_loadings.tsv` by `build_gsp_loadings.R`. They are NOT `out$L_pm`.
- Read outputs with `read_variant_table()` or fread. A plain `readr::read_tsv` guesses
  the ~all-NA mutation columns as logical and silently NAs every call.
- The YLPM1 caller keeps PTC consequences only, so missense calls never enter.
- The calling cohort excludes samples with catalogued variants (`exclude_mut_genes`), so
  `wxs_mutation` / `rnaseq_mutation` are filled for only a handful of samples.
- If the Rmd logic changes, port the change here and re-check against its outputs.

**User preferences** — Quinn: quetzal_code must contain everything needed to regenerate
the paper's analyses; readers point at these outputs, never at `flashier/`.

**Input/output** — `paths.yaml` keys `tcga_metadata`, `mutation_am_merged`,
`alphamissense_refs_dir`, `unified_calls_dir`, plus `shared_cache/gsp_loadings.tsv`.
Outputs go to `figure{2,3,4}/cache/*.tsv`, located via `vt_path()` (`vt_paths.R`).

```
flashier RDS ─► build_gsp_loadings.R ─► shared_cache/gsp_loadings.tsv ─┐
merged_grch38_with_am.tsv ─────────────────────────────────────────────┤
refs/<GENE>_alphamissense.tsv ─────────────────────────────────────────┼─► variant_tables.R
wxs_mut_calling/unified/<gene>/vep/{wxs,rna}/*.vep.vcf ────────────────┘      │
   figure2/scripts/build_cache_sf3b1_variants.R ─► figure2/cache/SF3B1_{full_info,SUGP1_full_info,full_plus_variants}.tsv
   figure3/scripts/build_cache_nkap_variants.R  ─► figure3/cache/{k13_full_info,NKAP_full_plus_variants}.tsv
   figure4/scripts/build_cache_ylpm1_variants.R ─► figure4/cache/YLPM1_{full_info,full_plus_variants}.tsv
        └─► readers via vt_path(): fig 2d, 3a/b/e, 4a/g, S9-S13, S17
```

> 2026-11-10: `paths.yaml` and `scripts/paths.R` were removed. Paths are now
> `file.path(QROOT, "<folder>/<file>")`, with `QROOT` the quetzal_code root
> (found via the `.quetzal_root` marker); references to paths.yaml keys above are historical.
