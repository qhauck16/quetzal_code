# cache/build_s14_cache.R

## Purpose & context
Cache for S14 (NKAP extended sashimi / coverage plots). Sections 1-5 (data half)
of `tcga/genome_wide/flashier/NKAP_detailed/custom_plot.R`, looped over the
top-15 events. The plotting half is redrawn in `figS15/build.R`. One RDS per
gene in `cache/s14/`, plus `s14/events.tsv` and `s14/samples.tsv`.

Since 2026-10-06 the event list (`NKAP_detailed/pick_events.R`) and the sample
list (`pygenome/scripts/select_samples.R`) are computed here instead of read
from static files, so S14 follows the diff-splicing results automatically.
Both ports were validated: on the pre-fix df_final the 15 events are identical
to the old events_meta.tsv, and the samples are the same 15 + 150 rail_ids
as the old NKAP/samples.tsv (k13 from the quetzal_code NKAP table, seed 123).

## Decision rationale
Bigwig coverage extraction is slow; caching it lets figS15/build.R iterate on
layout without touching bigwigs.

## Gotchas & warnings
- Any cwd (absolute paths). Bigwig reads are slow -> SLURM.
- figS15/build.R takes its gene order from `s14/events.tsv`; on the
  2026-10-06 results the same 15 genes as published, slightly re-ranked.
- Stale `s14/<GENE>.rds` from genes that drop out of the top 15 are not deleted
  automatically, but are ignored (build.R only reads genes in events.tsv).
- Window = 100 nt intron + 50 nt exon around the annotated acceptor (strand-aware).
- Missing/unreadable bigwigs become NA rows, not errors.
- Log prints per-gene PSI (condition vs control) for all 15 genes; 15 warnings.

## User preferences
Raw bigwigs + own ggplot drawing (no pygenometracks / ggbio).

## Input/output
In: TCGA recount3 bigwigs, leaf2_try `psi_table.tsv`, `hg38_granges.RDS`, MANE,
`unified_results/nkap/20v60/dpsi0.05_fc2.5/df_final.tsv`,
`figure3/cache/NKAP_full_plus_variants.tsv` (k13), `tcga_v2_samples.tsv`.
Out: `s14/<GENE>.rds` = list(gene, chr, strand, rank, shift, region, annot_acc,
cryp_acc, donor, win, cov, psi, exons, codons, ...).

## File flow schematic
```
df_final (nkap 20v60) ────┐
NKAP_full_plus_variants ──┤                        ┌► s14/events.tsv, s14/samples.tsv
tcga bigwigs ─────────────┼──► build_s14_cache.R ──┤
psi_table.tsv ────────────┤                        └► s14/<GENE>.rds (x15) ──► figS15/build.R
gencode, MANE ────────────┘
```
