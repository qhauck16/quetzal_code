# inputs

Public data every analysis and figure reads. Three scripts download it from the original sources into this directory, at the
paths the analysis folders and figure code expect:

- `fetch_inputs.sh`: single-file inputs (sample tables, GENCODE, MANE, AlphaMissense, PanCanAtlas, Quetzal v0.1.1 data).
- `fetch_snaptron.sh`: per-gene Snaptron junction tables. Requires the Snaptron client `qs` (snaptron-experiments) under Python 2.7; see its header.
- `fetch_bigwigs.sh`: recount3 bigwigs (several TB).

| File / directory | Source |
|---|---|
| `tcga_v2_samples.tsv`, `gtexv2_metadata.tsv` | Snaptron sample tables (tcgav2, gtexv2) |
| `snaptron/{tcga,gtex}/chr*/snaptron_output/<gene>_snaptron.tsv` | Snaptron per-gene junction queries (`fetch_snaptron.sh`: `qs --region <gene>`, snaptron-experiments d03371d), protein-coding genes from `hugo_all_pc_genes.tsv` |
| `gencode_v46_granges.RDS`, `hugo_all_pc_genes.tsv`, `analyte.tsv` (RIN), `productive_unproductive/` | shipped with Quetzal v0.1.1 (`data/`) |
| `gencode.v46.annotation.gtf`, `GRCh38.primary_assembly.genome.fa` | GENCODE release 46 |
| `MANE.GRCh38.v1.4.summary.txt` | NCBI MANE v1.4 |
| `GRCh38.d1.vd1.fa` (+ `.fai`, `.dict`) | GDC GRCh38.d1.vd1 reference (WXS / RNA-seq re-calling) |
| `vep_cache/` | Ensembl VEP RefSeq cache, release 115, GRCh38 |
| `star_index/` | GDC STAR 2.7.5c index, GRCh38.d1.vd1 + GENCODE v36 (GIST alignment) |
| `alphamissense/` | AlphaMissense hg38 predictions; UCSC hg19→hg38 liftOver chain |
| `Alsafadi2016_supp1.xlsx` | Alsafadi et al. 2016 Supplementary Data 1 (SF3B1-mutant uveal melanoma junctions, hg19) |
| `pancanatlas/` | PanCanAtlas ABSOLUTE purity/ploidy and allelic segments (GDC, hg19) |
| `bigwigs/{tcga,gtex}/` | recount3 coverage bigwigs, URLs in `lists/` (`fetch_bigwigs.sh`) |
| `cosmic_tier_1.tsv` | COSMIC Cancer Gene Census Tier 1 (account required; place by hand) |
| `cola_seq_branchpoints.csv` | CoLa-seq human branchpoints (place by hand) |

Checked byte-identical to
the manuscript's copies: both sample tables, the Quetzal v0.1.1 data files, MANE
and both PanCanAtlas files; the derived gene-to-chromosome split matches the
original exactly.

On midway, `bash fetch_inputs.sh --midway` symlinks the copies already on disk instead.
