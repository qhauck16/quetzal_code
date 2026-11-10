# overlap.R

## Purpose & context
The statistic for the IDH1/SPOP note: are a GSP's top driving alt-ss junctions enriched
among the junctions LeafCutter calls as up in mutants (vs matched wild-type)? It is the S7
test, applied to:
- IDH1 R132 vs GSP #1;
- each SPOP residue group vs GSP #12;
- the SPOP groups vs GSP #40 (exploratory).

## Decision rationale
- `build_sets` / `cells_for` are copied unchanged from `figures/supplementary/cache/
  build_overlap_cache.R` (S7), so the statistic is identical:
  - **background:** junctions matched 1:1 on the MANE intron's median PSI + RPM;
  - **test:** Fisher's exact test on overlap with LeafCutter-up junctions
    (p.adjust < 0.05, dPSI > 0).
- Reported for alt 3'ss, alt 5'ss and the two combined. Skipped-exon and "other" junctions
  have no MANE partner intron, so they can't be matched (the same limit as in S7).

## Gotchas & warnings
- **LeafCutter's dPSI sign follows the order of the phenotype file.** `leafcutter_ds.R`
  encodes the first group it sees as 0 and reports dPSI = group 1 − group 0.
  - The first IDH1/SPOP run listed mutants first, so its dPSI is control − mutant.
  - `load_leaf()` reads the first row of `leaf_pheno_data.txt` and flips the sign, so the
    output is always mutant − control.
  - The Snakefile now writes controls first.
  - The published GSP runs (`diff_splicing/results/`) are control-first and correct;
    `control_results/` null runs are test-first, but they have no meaningful direction.
- Uses `set.seed(1)`, as in S7.
- **Fold-change filter (`min_fc`), added 2026-10-09.** It asks whether the GSP junctions
  are newly used or already-used junctions nudged slightly.
  - FC = mutant PSI / control PSI (control / mutant for down), using LeafCutter's PSI
    columns, which are named by group. This matches the diff-splicing `df_final`
    definition.
  - Four filters run: none, |dPSI| > 0.05, FC ≥ 2.5, and both.
  - Result: under FC ≥ 2.5 the matched background drops to 0 everywhere. IDH1/SPOP keep
    significant GSP overlap (e.g. IDH1 11/161, P = 8e-4; F133 16/160, P = 2e-5). So their
    GSP junctions are low-PSI junctions switched on ≥ 2.5× but by < 0.05 PSI.
- **Effect-size sensitivity (`min_dpsi`).** S7 has no dPSI floor: a junction counts if its
  cluster has p.adjust < 0.05, with any dPSI sign. Quinn asked (2026-10-09) for a
  |dPSI| > 0.05 version, so every comparison is run at `min_dpsi` = 0 and 0.05.
  - The nine S7 GSP comparisons are included.
  - Their `min_dpsi` = 0 rows reproduce `s7_leafcutter_cells.tsv` exactly (all 36 cells,
    identical p-values), which validates the port.
  - Those S7 sets come from `tcga/.../flashier/fact_{13,30}_top_junctions.tsv`, the
    fallback when a set isn't in `gsp_top/`.

## Input/output
- **In:** `gsp_top/`, `results/<comp>/leafcutter_ds_*.txt`, `leaf_pheno_data.txt`, and
  `all_junction_rpm_and_psi.tsv`.
- **Out:** `out/overlap_cells.tsv`. Run on `amd` with 48 GB; Mahalanobis matching over
  the full junction pool is heavy.

```
gsp_top/fact_N ─┐
results/<comp> ─┼─► overlap.R ─► out/overlap_cells.tsv
rpm_psi table ──┘
```
