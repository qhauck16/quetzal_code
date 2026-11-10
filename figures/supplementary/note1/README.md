# Supplementary Note 1: IDH1 / SPOP

## Purpose & context
This note asks whether non-spliceosome mutations
associated with a GSP (IDH1 R132 in LGG with GSP #1; SPOP hotspots in PRAD with GSP #12)
reproduce that GSP's splicing in direct within-cohort mutant-vs-wild-type comparisons.
This folder groups every Note 1 panel. Outputs go to `out/note1/`.

## Scripts
Run all of them from `supplementary/`, like the figS builds.

| Script | Output | Content |
|---|---|---|
| `cdf_loadings.R` | `note1_gsp_cdf.pdf`, `note1_cdf_stats.tsv` | GSP loading CDFs, all solid lines. (a) IDH1: LGG mutant and wild-type, plus R132 carriers and IDH1 wild-type pooled over non-brain cohorts. (b) SPOP: PRAD mutant and wild-type, plus SPOP wild-type of the other cohorts, with the 3 hotspot carriers outside PRAD as labelled points on that curve. (c) GSP #13 within TCGA-UCEC: NKAP C-term (300–413) / PRPF8 RNase H (1570–1900) variant carriers vs other UCEC vs other cohorts. (d) GSP #10 within TCGA-UVM: SF3B1 hotspot vs other UVM vs other cohorts. c/d show that for spliceosome genes the GSP follows the mutation, not the tissue. Per-cohort stats for every cohort go to the tsv |
| `cdf_scale_compare.R` | `note1_cdf_scale_compare.pdf` | SF3B1, NKAP and YLPM1 (pan-cancer, the exact classes of Fig 2d/3a/4a), IDH1 (LGG) and SPOP (PRAD) CDFs on ONE shared x range (log GSP loading). Shows how small the IDH1/SPOP shift is next to the spliceosome genes |
| `build_psi_cache.R` (SLURM, ~10 min) | `note1/cache/psi_group_means.tsv` | per-junction mean PSI of the GSP #1 / #12 top driving junctions (`mutation_ds/gsp_top`, all types) in brain / non-brain / IDH1 R132 brain / other brain and PRAD / non-PRAD / SPOP-mutant PRAD / other PRAD |
| `psi_cdf.R` | `note1_psi_cdf.pdf`, `note1_psi_contrast_stats.tsv` | Tissue (covariate) vs mutation (splicing): CDFs of mean PSI per junction, and per-junction differences split additively (tissue = wild-type tissue − other tissues; mutation = mutant − wild-type tissue) |
| `signature_plots.R` | `note1_signature.pdf` | (a/b) LeafCutter signature score (mean PSI up − mean PSI down; adj. P < 0.05, \|dPSI\| > 0.05) in HELD-OUT LGG / PRAD tumours, mutant vs wild-type, one-sided Wilcoxon. (c/d) PCA of the signature junctions' scaled logit PSI (clipped to 0.001–0.999; junctions with ≤ 20% NA, NA → median) over each whole cohort. Inputs: `mutation_ds/out/signature_scores.tsv` and `signature_psi_<comp>.rds` from `mutation_ds/signature_check.R` |
| `dpsi_overlap.R` | `note1_dpsi_overlap.pdf` | LeafCutter dPSI of GSP vs matched junctions, and the overlap test with and without \|dPSI\| > 0.05 |

## Decision rationale
- **The CDFs use the whole cohort,** not the 20 v 60 LeafCutter selection. Samples are
  flashier tumours with mutation calls, and wild-type means no call anywhere in the gene.
- **Loadings** come from `shared_cache/gsp_loadings.tsv`, the same normalised L as the
  main figures. They are NOT `L_pm` loadings, which are on
  another scale.
- **Plotted like Fig 3a / 4a**: x = natural-log GSP loading, labelled
  `log(GSP #N)`, y = cumulative percentage, step curves.
- **All mutants are pooled into one curve** per gene. Mutant = any call
  in the gene.
- **No per-cohort facet panels**; per-cohort numbers are in `note1_cdf_stats.tsv`.
- **SPOP outside PRAD is 3 labelled points, not a curve:** DLBC F102I, LGG W131C,
  UCEC W131R. That's too few for a CDF; the evidence is weak and is presented as such.
- **Test:** one-sided Wilcoxon, mutant > wild-type.

## Gotchas & warnings
- **Signature PCA:** LGG separates by IDH1 status along PC1 (11.7% of variance). PRAD
  does not cluster by SPOP (PC1 6.8%; mutants only slightly shifted). The held-out
  score difference in PRAD is real but small relative to other variation.
- **Tissue vs mutation (psi_cdf.R).** Panels a/d use all brain / all PRAD. The difference panels use wild-type tissue, because 412 of the 660 brain
  tumours are R132.
  - GSP #1: tissue median dPSI +0.006, mutation +0.005. The two are comparable, with
    93% and 86% of junctions up.
  - GSP #12: tissue +0.009, SPOP +0.0015 (65% up). GSP #12 is mostly prostate
    tissue, with a small SPOP shift on top.
- **The PRAD SPOP group (n = 56) is bigger than the LeafCutter groups combined.** It pools
  every residue and keeps duplicate aliquots. The CDFs keep every flashier rail_id, including duplicate aliquots, while
  `mutation_ds/select_samples.R` dedupes aliquots and draws from the metadata, not the
  fit. One selected F102 sample has no loading.
- **No "Other IDH1" group in LGG.** Every IDH1-mutant LGG with a call is R132.
- **Outside the home cohorts, only GBM IDH1 was significant** (7 mutants, all R132,
  P = 0.004; excluded as brain). CHOL's 5 R132 tumours are not, and SPOP outside PRAD (almost no hotspots)
  shows nothing.
- **IDH1 "other cohorts" means non-brain**: GBM is dropped from panel
  a's other-cohort curves. GSP #1 is brain-associated; with GBM in,
  the pooled R132 curve gives P = 0.005, and without it P = 0.15 (40 R132 vs 8,498).
  SPOP has only 3 hotspot carriers outside PRAD.
- **IDH1 wild-type LGG also loads highly on GSP #1** (median 0.43 vs 0.63). GSP #1 is
  broadly glioma-associated, and R132 shifts it further.
- **`dpsi_overlap.R` needs `mutation_ds/out/junction_dpsi.tsv`,** which `overlap.R`
  writes (run it on amd, ~15 min).

## Input/output
```
shared_cache/gsp_loadings.tsv ─┐
AM-merged mutation calls ──────┼─► cdf_loadings.R ─► out/note1/note1_gsp_cdf.pdf (+ stats tsv)
tcga_v2_samples ───────────────┘
mutation_ds/out/{junction_dpsi,overlap_cells}.tsv ─► dpsi_overlap.R ─► out/note1/note1_dpsi_overlap.pdf
```
