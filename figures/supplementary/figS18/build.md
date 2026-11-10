# figS18/build.R -- Supplementary Figure S18

## Purpose & context
GSP #30 vs YLPM1 expression. (a) top 100 samples on GSP #30: log loading vs
YLPM1 expression (TCGA cohort-scaled RPM), YLPM1 mutants labelled
(`ylpm1_lof_investigation.Rmd`). (b) YLPM1 LoF samples split by SVM phenotype
status (Phen+ score > 0 vs Phen-). (c) YLPM1 expression in Phen+, the 25
phenotype-negative intermediates nearest the boundary, and the rest
(`diff_splicing/ylpm1/top500_alt3_fc_vs_loading.Rmd`, last two chunks).
Wilcoxon rank-sum p-values.

## Decision rationale
LoF = protein change matching `*`, `fs` or `splice` -- the rule panel b's
source and main Figure 4a use, so all three panels and fig 4 agree.

## Gotchas & warnings
- **Known deviation (a):** X1783_splice is classified LoF here; the original
  panel a used the table's `plot_color` column, which left it "other".
- **Known deviation (b):** Phen- n = 105 vs published 106 --
  `YLPM1_full_plus_variants.tsv` was modified after the figure was made.
- Prints `[S17] b p = ... (n pos/neg); c p = ...` -- check it in the build log.
- `cache/ylpm1_expression_row.tsv` (one YLPM1 row, columns = rail_ids) comes from
  `cache/build_ylpm1_expression_row.sh` (header + YLPM1 row of combined_scaled_rpm.tsv).

## User preferences
Mutant labels on the points (ggrepel); gene symbols italic. The repel layer
gets ALL 100 points (unlabelled ones label "") so labels avoid the grey
cloud too; M1180Cfs*27 gets a hand nudge because the solver otherwise left
it on its own point. b/c titles wrapped to two lines (b ran into c tag);
b/c x ticks angled 45 deg. Bold tags set with `& theme(plot.tag=...)`
because plot_annotation(theme=) does not reach patchwork tags.

## Input/output
In: `$FLASHIER/YLPM1_full_plus_variants.tsv` (k30, proteinChange),
`cache/ylpm1_expression_row.tsv`, `diff_splicing/ylpm1/ylpm1_phenotype_axis_scores.tsv`.
Out: `out/figS18.pdf` (8.5 x 9; widened from 7 x 8.5 to give panel a labels room).

## File flow schematic
```
YLPM1_full_plus_variants.tsv ───────┐
cache/ylpm1_expression_row.tsv ─────┼──► figS18/build.R ──► out/figS18.pdf
ylpm1_phenotype_axis_scores.tsv ────┘
```

**2026-09-24 (Quinn review):** panel c x labels 'Next 25 intermediates' / 'Other samples'.
