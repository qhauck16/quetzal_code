# figure4/panel_g — YLPM1 expression ECDF + bootstrap ribbon

## Purpose & context

Panel g of Figure 4. It is the "so what" panel for YLPM1: genes that pick
up a Leafcutter2-**unproductive** junction in the top-20 GSP #30 samples
are shifted *down* in expression relative to what the same genes do in
matched random sample groups. That is the NMD-mediated
downregulation claim in the title.

Reproduces
`flashier/diff_splicing/expression_CDFs/publication/YLPM1_leaf2_psi0.01_fc2.5_bootstrap.pdf`
(published as panel g of `flashier/paper_figures/Figure4-YLPM1_improved.pdf`),
re-drawn with the shared `../theme.R` tokens instead of `theme_classic(14)`.

The analysis is transcribed from
`flashier/diff_splicing/combined_expression_analysis.Rmd`
(`build_set_pair` / `stratified_random` / `bootstrap_random_cdfs` /
`plot_single_panel_cdf`), `pub_specs` entry 5 (line 887):

```
gene_label "YLPM1", criterion "leaf2", dist_window c(-Inf, Inf),
require_alt3 FALSE, include_complement FALSE, include_bootstrap TRUE,
x range c(-1.3, 1.3)
```

## The analysis

Constants: `N_TOP = 20`, `N_BOOTSTRAP = 100`, `N_GRID = 200`,
`seed_base = 1000 * N_TOP + 1 = 20001`.

1. **Target gene set.** `df_final.tsv` from
   `unified_results/ylpm1/20v60/dpsi0.01_fc2.5/`, joined to
   `Leaf2_junction_classifications.txt` on
   `Intron_coord = chrom:(start-1)-end`. Hits are `deltapsi > 0.01 &
   test/control > 2.5`; the target set is the unique `gene_name` of hits
   with `!Coding & !UTR`, intersected with the expression matrix's rows.
   `require_alt3 = FALSE`, so — unlike NKAP/SF3B1 — there is **no**
   alt-3'ss filter: YLPM1 hits both alt 3'ss and alt 5'ss.
   **604 genes**, all 604 present in the matrix.
2. **Red curve.** `rowMeans` of `combined_scaled_rpm.tsv` over the 20
   rail_ids with the largest `k30` in `YLPM1_full_plus_variants.tsv`,
   ECDF'd.
3. **Bootstrap.** 100 iterations; each draws a stratified random sample
   from the *whole* rail_id pool matching the top-20's per-`project_id`
   counts, then takes the same 604 genes' per-gene means. The ribbon is
   `rowMeans` and the 2.5 % / 97.5 % quantiles of the 100 ECDFs on a
   shared 200-point x grid spanning the pooled range.
4. **p-value.** one-sided `wilcox.test(red_means, pooled_boot_means,
   alternative = "less")` → **p = 1.4e-05**. Matches the reference PDF.
5. **Gene labels.** 10 highlights (ATM, WRN, FH, STAG2, SPOP, LATS1,
   MSH6, TP53BP1, XPC, PBRM1); **all 10** are in the target set. Each
   sits at its own percentile within the RED curve,
   `ecdf_y = rank(x, ties.method = "first") / n()`.

Verified against the reference render: n = 604 and p = 1.4e-05 both agree
(the reference's legend block prints "Unproductive: 604" twice — once for
the curve, once for the bootstrap ribbon, same gene set).

## Caching strategy

`combined_scaled_rpm.tsv` is **7.7 GB** (34,665 genes x 11,173 rail_ids).
Figure 3 panel e could slice it to 20 columns; here that is impossible,
because the bootstrap samples from the whole column pool. So the slice
goes the other way: **all columns, only the 604 target rows.**

`scripts/build_cache_panel_g.R` does everything expensive:

1. computes the candidate gene list from `df_final` + Leaf2 (no big file
   touched), writes `cache/panel_g_cand_genes.txt`;
2. shells out to a single `awk` pass that streams the 7.7 GB matrix and
   emits the header plus the matching rows → `cache/panel_g_expr.tsv`
   (**134 MB**). The awk match *is* the Rmd's
   `intersect(rownames(expr_mat))`, so the sliced file's row names are
   the final target set;
3. runs the top-20 selection, the red curve, the 100-iteration bootstrap
   and the Wilcoxon in R;
4. writes the small artifacts `compose.R` actually reads.

`panel_g/build.R` reads only those and never sees the big matrix. Whole
build is a couple of seconds.

Re-run the cache builder only when the upstream data changes:

```
cd quetzal_code/figures/figure4
Rscript --vanilla scripts/build_cache_panel_g.R
```

It short-circuits the awk pass if `cache/panel_g_expr.tsv` already
exists; delete that file to force a re-slice.

## Gotchas & warnings

- **`start - 1` before building `Intron_coord`.** `df_final` is 1-based,
  Leaf2's `Intron_coord` is not. Without the offset the join returns all
  `NA` and the `!Coding & !UTR` filter silently yields an empty set.
- **Never put `limits =` on the x scale.** A scale limit DROPS rows
  before `stat_ecdf` runs, so the drawn curve would be the ECDF of a
  truncated sample (the red means run from -1.49 to +1.39, outside the
  ±1.3 view). The view is set with `coord_cartesian(xlim =, ylim =)`.
- **Series names travel as DATA COLUMNS (`red$grp`, `ribbon$grp`).**
  `panel_cached()` round-trips the panel through `saveRDS`; an `aes()`
  that referenced a scalar from the enclosing `local()` would resolve to
  nothing after `readRDS` and the whole figure would render blank. Same
  for the p-value string, which is carried in `ann$label`.
- **ggrepel is blind to the curves and the ribbon.** The step curve, the
  dashed mean and the grey band are all fed into both repel layers as
  zero-text obstacle rows, along with a tiled rectangle standing in for
  the legend. Without them the label column walks straight onto the
  ribbon.
- **Two repel layers, disjoint y ranges.** The corridor layer
  (`ecdf_y >= 0.10`) is capped just under the legend block; PBRM1
  (`ecdf_y = 0.02`) gets its own band in the bottom-left corner, because
  at that height the red curve's tail has already run out past x = -1.4
  and there is no shared corridor. Disjoint ranges also guarantee the two
  layers cannot collide with each other.
- **`max.time = 1e6`** — it is a wall-clock budget, not an iteration
  count. A small value makes `seed` useless and the panel stops being
  reproducible run to run.
- **The x axis title is wrapped onto two lines.** On one line it is
  ~4.0 in against a 4.07 in cell, so it overran the page edge and ran
  into panel f.
- **The long legend entry is wrapped in the series NAME**
  (`G_BOOT`), not by a `labels =` argument, so nothing downstream has to
  re-map it. Unwrapped it is ~2.4 in wide and there is no 2.4 x 0.3 in
  empty pocket in this panel.
- **`G_IN_PER_X` / `G_IN_PER_Y` are MEASURED** off the composed render
  (1.125 / 1.610 in per data unit at row 4's right cell). The legend
  obstacle rectangle is derived from them, so re-measure if `ROW_H` or
  `COL_W$fg` change.
- The reference labels the genes in **bold**; here they are **italic**
  per the figure-wide gene-name convention (same as figure 3 panel e).

## Input/output

```
Leaf2_junction_classifications.txt ─┐
ylpm1/20v60/dpsi0.01_fc2.5/         ├─► scripts/build_cache_panel_g.R
  df_final.tsv ─────────────────────┤     │  (awk row-slice + bootstrap)
YLPM1_full_plus_variants.tsv ───────┤     │
combined_scaled_rpm.tsv (7.7 GB) ───┘     ▼
                                    cache/panel_g_expr.tsv (134 MB)
                                    cache/panel_g_red.tsv        (604 rows)
                                    cache/panel_g_ribbon.tsv     (200 rows)
                                    cache/panel_g_highlights.tsv (10 rows)
                                    cache/panel_g_top_ids.txt    (20 ids)
                                    cache/panel_g_meta.rds       (n, p-value)
                                              │
                                              ▼
                                    panel_g/build.R ──► compose.R ──► out/figure4.pdf
```

`panel_g_meta.rds` carries `n_target`, `n_top`, `p_value`, `p_label`,
`n_highlight`, `missing_highlight`, `n_boot`; `build.R` prints the first
three on every rebuild so a silent data change is visible in the log.

## 2026-11-09: gene set restricted to alt 3'ss + alt 5'ss
- **Change (Quinn):** `build_cache_panel_g.R` now keeps only the alt 3'ss / alt 5'ss hits (`TYPES`), so the set
  matches the Fig. 4f heatmap's YLPM1 event classes. All other filters are unchanged: ΔPSI > 0.01, FC > 2.5 and
  LeafCutter2-unproductive (not Coding, not UTR), 20v60, after the -m 5 re-run.
- **Result:** 537 genes, p = 6.0e-6 (one-sided Wilcoxon). The all-types set was 583 genes, p = 9.2e-6.
- **Highlights:** 9 of 10 genes are labelled; ATM drops out because its unproductive junction is not an alt-ss event.
- **Previous cache:** `cache/pre_alt35_2026-11-09/`.
