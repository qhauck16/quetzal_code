# coverage_tracks

Coverage / sashimi track inputs for the genome-browser panels (Figs. 2-4, Supp.).

| Code | Produces |
|---|---|
| `pygenome/Snakefile` (+ `config.yaml`, `regions.tsv`, `scripts/`, `envs/`) | `pygenome/{NKAP,YLPM1,SF3B1,neuronal_gsp}/<gene>/`: per-group averaged bigwigs, MeanPSI links, junction BEDs |
| `gsp_tracks.R` then `gsp_tracks.sh` | `fig4b_tracks/{KDM4C,TRIM38}/`: the same for the Fig. 4b loci (top-13 GSP #30 vs matched controls) |

Both average recount3 bigwigs (`../inputs/bigwigs/`) per sample group with
GenometracksByGenotype's `AggregateBigwigsForPlotting.py` (github.com/bfairkun/GenometracksByGenotype;
path in `pygenome/config.yaml`, or `GTBG` for `gsp_tracks.sh`). The figures draw the tracks
themselves from these files; pyGenomeTracks is not used. On midway,
`bash make_links.sh` links the manuscript's track folders.
