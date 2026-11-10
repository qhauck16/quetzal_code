# figure4/panel_e/build.R

## Purpose & context

Panel e of Figure 4 (YLPM1 / GSP #30): *"Tumor suppressor genes and
frequently mutated oncogenes are mis-spliced in samples with high
loadings upon GSP #30."* A full-page-width 2x2 grid of loci —
**TP53BP1** (top-left), **PBRM1** (top-right), **MSH6** (bottom-left),
**ASXL1** (bottom-right) — each drawn as

```
High GSP #30 coverage + junction arcs
Matched Controls coverage + junction arcs
lavender taper (trapezoid_gg)
full MANE gene model (transcript_zoom_gg)
genomic coordinate row (bespoke, same x scale)
```

No zoom callouts (figure 3 panel c's top row has them; this panel's
reference does not).

Reference being reproduced:
`tcga/genome_wide/flashier/paper_figures/Figure4-YLPM1_improved.pdf`,
row 3. `compose.R` sources this file and expects a `panel_e` object; it
is wrapped in `panel_cached("e", ...)`.

Everything is drawn from the raw bigwigs and `.links`. **No
pyGenomeTracks, no ggbio, no re-used `*_MANE_gene_model.pdf`** — Quinn
has rejected raster/pre-rendered gene models repeatedly, and each of
these four directories contains a tempting `<GENE>_MANE_gene_model.pdf`
that must be left alone.

## Decision rationale

Same three structural decisions as `../panel_b/build.md`, for the same
reasons — read that file first; only the differences are repeated here.

**Arc ends sit on the BASELINE at the junction coordinate.**
`sashimi_track_gg` is called with `min_psi = Inf` so its own arc engine
(which pins ends to the local coverage top) draws nothing; it is used
for the coverage and the shared theming, and the arcs are a local
`geom_path` layer. X is always the raw `.links` coordinate.

**Apex tier by SPAN, longest highest**, band `A_LO = 0.10`,
`A_HI = 0.78` of `y_max`. The band's WIDTH is what lets a PSI number
sit in the *gap between two arcs* instead of being lifted above the
whole bundle. ASXL1 stacks four junctions off one acceptor
(32,366,383), so the tier gap is `(A_HI - A_LO)/3 * y_max` and it has
to exceed `clr + clr_a` ≈ 0.205 × `y_max`. At the 0.16–0.62 band first
tried the gap was 0.153 and "4.8" was lifted clear above the 74.0 arc,
where it read as a second label for 74.0. `A_HI` is then capped by the
label ceiling: `A_HI <= 1 - 2 * clr / y_max`.

**PSI labels solved in two dimensions** by `place_psi()` — identical to
panel b's. ASXL1 is the locus that needs it.

**No orange here.** The "unannotated junction driving GSP #30"
annotation belongs to panel b only; these arcs are all black.

**Lavender taper `#ccccff`**, matching figure 3 panel c rather than the
reference's pink. Change it in both sashimi panels or neither.

**Captions on the right column only**, as in the reference, which is
why the two columns get different `MARG` values — otherwise the
right-hand tracks would be as wide as the left-hand ones and the
captions would run off the page.

**The italic gene name's SIDE is derived, not listed per locus.**
`transcript_zoom_gg()` anchors the name on the transcript's own
`tx_lo` / `tx_hi` edge and draws it INSIDE the gene row — the same row
the lavender highlight rect spans from `y = 0` to `y = 2`. So the name
lands on the taper's foot whenever the zoom window reaches that end of
the gene. `make_gene_rows()` picks an end the window does not reach,
testing with the label's real drawn width converted to bp
(`lab_bp`, Helvetica caps ≈ 0.62 em — the same estimator this file
already uses for the coordinate labels), plus a 0.07 in visible gap.
**LEFT is preferred when both ends are free**, so the four loci read
consistently and the one exception is forced by the data rather than
chosen. It resolves to left / left / left for TP53BP1, PBRM1 and MSH6
and **right for ASXL1**, whose window starts AT `tx_lo`. The chosen
side is printed at build time.

## Gotchas & warnings

- **The published "High GSP #13 (15)" label is HALF stale, and the
  stale half is the GSP number.** Every `sashimi_tracks.ini` in these
  four directories reads `title = High_k30 (15)` and
  `Matched_Controls (150)`, and the bigwigs are `High_k30-*`. So the
  sample counts 15/150 are right and `#13` is wrong — it is figure 3's
  NKAP factor index left in the YLPM1 figure. Corrected to
  **High GSP #30 (15)**. (Panel b's labels were stale on *both* counts;
  see `../panel_b/build.md`.)
- **Ignore every `High_k65-*` file.** k65 is a stale factor key for
  YLPM1; the current key is k30. TP53BP1, MSH6 and ASXL1 all carry k65
  leftovers; PBRM1 does not.
- **ASXL1's region starts 201 bp BEFORE its MANE transcript does**
  (32,358,130 vs 32,358,331). `make_gene_rows()` clamps the window with
  `win <- c(max(wide[1], tx_lo), min(wide[2], tx_hi))` so the taper's
  feet and the highlight rect stay on the gene model instead of
  spilling into the left gutter under `clip = "off"`. The same fact is
  why ASXL1 is the locus whose gene name has to go RIGHT: its window
  begins at `tx_lo`, so a left-anchored "ASXL1" sits inside the
  lavender highlight. It did, in the published render.
- **`label_side` is deliberately NOT a field of `LOCI`.** It is not a
  property of the locus, and hand-setting it per gene is exactly what
  put ASXL1's name on its own taper while PBRM1's sat clear. If a
  locus or a region ever changes, the side follows automatically.
- **The two y-extreme labels are pushed `AX_IN` lines INWARD at both
  ends.** Consecutive tracks butt flush (the gutters carry no
  top/bottom margin), so the only thing separating this track's "0"
  from the next track's "27.3" is `2 * AX_IN` lines. At the 0.15 first
  used that gap was 0.027 in and the pair read as one stacked number;
  0.40 gives 0.074 in and still costs only ~7% of the track's drawn
  height. The labels sit in the left gutter (`clip = "off"`), so
  moving them inward cannot collide with the coverage.
- **`panel_cached()` round-trips through `saveRDS`.** Any `aes()` that
  references a scalar from the enclosing environment renders blank
  after `readRDS`. Everything here is a data column; keep it that way.
  See `../figure2/compose.md`.
- **`E_TITLE_H` changes nothing about the layout** — patchwork heights
  are relative. Its only job is to make `TRACK_H_IN` the track's true
  drawn height, which the PSI solver converts point sizes through. It
  was MEASURED: at 400 dpi panel e's cell runs y = 1805..3186 px and
  the first track's top edge is at 1900, giving a 3.215 in content band
  and `TRACK_H_IN = 0.536 in`. Re-measure after any change to the title
  or the row heights.
- **Two clearances, not one** (`clr` off the label's own arc, the
  smaller `clr_a` from every other arc), and the ceiling is the label's
  TOP edge (`y_max - clr`), not its centre. Collapsing either back
  causes overprinting or digits drawn into the track above.
- **`place_psi()` warns when a label has no clear slot.** That warning
  is the regression alarm — a silent fallback puts the number at its
  arc's apex whether or not that collides.
- `PBRM1`'s MANE transcript is `ENST00000707071.1`; it resolves, but it
  is a recent ID, so if the GENCODE GRanges is ever refreshed re-run
  the `stopifnot(length(me) > 0)` check rather than assuming.
- The arc/label machinery is deliberately DUPLICATED from
  `panel_b/build.R`. `../scripts/coverage_plot.R` is shared with
  figures 1–3 and must not be edited for one panel's benefit. Keep the
  two copies in step by hand.

## User preferences

- Stretch the plot rather than shrink the type; `SZ_PSI` is
  `SIZE_BASE * 0.80`, and the arc band was widened rather than the
  numbers made smaller.
- Thin curved arcs scaled by `sqrt(PSI)`
  (`0.25 + 0.60 * sqrt(psi/100)` pt), never filled shapes.
- Labels on the object: every PSI number sits on its own arc.
- Only `FONT_FAMILY` and the `SIZE_*` tokens from `../theme.R`.
- No pyGenomeTracks, no ggbio, no pre-rendered gene-model PDFs.

## Input/output

Reads, per locus, from
`<FLASHIER>/pygenome/YLPM1/<GENE>/<chrom>_<start>-<end>/` — the
directory name IS the region:

| locus | region | ini `max_value` |
|---|---|---|
| TP53BP1 | chr15:43,469,677–43,474,976 | 6.9 |
| PBRM1   | chr3:52,633,846–52,642,490  | 5.0 |
| MSH6    | chr2:47,800,733–47,803,847  | 27.3 |
| ASXL1   | chr20:32,358,130–32,366,548 | 10.5 |

(printed to one decimal, as the reference does; both tracks of a locus
share it so the cohorts are comparable).

Files used per directory: `High_k30-3-..bw`,
`Matched_Controls-3-..bw`, `MeanPSI.High_k30-3-..links`,
`MeanPSI.Matched_Controls-3-..links`, `sashimi_tracks.ini`
(`max_value` + titles, read by hand).

Arcs surviving the `s1 >= start & e2 <= end` filter, verified against
the reference:

| locus | high PSI | control PSI |
|---|---|---|
| TP53BP1 | 69.3, 12.9, 7.81, 1.54 | 91.2, 1.5 |
| PBRM1   | 83.8, 5.12, 1.19       | 99.0 |
| MSH6    | 85.4, 11.5             | 99.4 |
| ASXL1   | 74.0, 16.9, 4.8, 1.28  | 85.1, 5.2, 4.93 |

Gene models come from `tcga/hg38_granges.RDS` filtered to the MANE
transcript listed in
`<FLASHIER>/diff_splicing/MANE.GRCh38.v1.4.summary.txt`, cached once
into `figure4/cache/panel_e_exons.rds`. Delete that file to rebuild it.

Writes nothing but the `panel_e` object; `compose.R` serialises it to
`figure4/panels/e.rds` and renders `figure4/out/figure4.pdf`.

## File flow schematic

```
High_k30-3-..bw ─────────┐
Matched_Controls-3-..bw ─┼─► bigwig_binned() ──► coverage area
                         │
MeanPSI.*.links ─────────┼─► read_links() ──► arc_frame() ──► geom_path arcs
                         │                              └──► place_psi() ──► PSI labels
sashimi_tracks.ini ──────┘   (max_value, sample counts: read by hand)

hg38_granges.RDS ──┐
MANE summary.txt ──┴─► cache/panel_e_exons.rds ──► transcript_zoom_gg() + trapezoid_gg()

                       all of the above ──► panel_e ──► compose.R ──► out/figure4.pdf
```

## Running it

```
cd /project2/yangili1/qhauck/spliceNMF/quetzal_code/figures/figure4
/scratch/midway2/qhauck/conda_env/rstudio-server/bin/Rscript --vanilla compose.R e
pdftoppm -r 400 -png -singlefile -x 0 -y 1800 -W 3400 -H 1400 out/figure4.pdf /tmp/crop
```

then look at `/tmp/crop.png`. Panel e is row 3, full width, y 1805–3186
at 400 dpi. Crop the gene-model rows at 600 dpi to check the gene-name
gutters (`-x 130 -y 3540 -W 1300 -H 330` for TP53BP1, `-y 4400` for
MSH6, `-x 1550` for the right column).

A row-wise blank-band scan of the render finds no blank run over
0.12 in inside the panel; the two that exist are the inter-row gutter
between the coordinate labels and the next stack, and the panel's own
bottom padding.
