# build_gsp_loadings.R

**Purpose** — Caches `flashier::ldf(out, type = "i")$L` from the 3.7 GB flashier RDS
(9,800 samples × 300 GSPs) as `shared_cache/gsp_loadings.tsv`, so the variant-table
builders don't each load the RDS. This is `norm_L` in flash_fact_char.Rmd.

**Gotchas** — Use SLURM with an absolute `-o` path (a path relative to `--chdir`
silently fails the job). Not the same as `diff_splicing/flashier_loadings.tsv`,
which is `out$L_pm`. Re-run whenever `paths.yaml: flashier_rds` changes.

```
flashier_rds ─► build_gsp_loadings.R (caslake, ~1 min, 48 GB) ─► shared_cache/gsp_loadings.tsv ─► variant_tables.R
```

> 2026-11-10: `paths.yaml` and `scripts/paths.R` were removed. Paths are now
> `file.path(QROOT, "<folder>/<file>")`, with `QROOT` the quetzal_code root
> (found via the `.quetzal_root` marker); references to paths.yaml keys above are historical.
