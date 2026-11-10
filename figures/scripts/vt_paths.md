# vt_paths.R

Defines `vt_path("<table>")`, the single source of truth for where each per-gene
variant table lives (`figure{2,3,4}/cache/<table>.tsv`). Sourced by `theme.R`
(so every figure has it) and directly by the data-only cache builders
(`figure4/scripts/build_cache_panel_g.R`, `figure3/scripts/build_cache_panel_e.sh`,
`supplementary/cache/build_mutation_cache.R`). To add or move a table, edit `VT_HOME`.

> 2026-11-10: `paths.yaml` and `scripts/paths.R` were removed. Paths are now
> `file.path(QROOT, "<folder>/<file>")`, with `QROOT` the quetzal_code root
> (found via the `.quetzal_root` marker); references to paths.yaml keys above are historical.
