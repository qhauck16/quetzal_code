# figure3/panel_c/build.R — four mis-splicing loci

## Purpose & context

Panel c of figure 3 (NKAP / GSP #13): a 2x2 grid of loci — SENP2 and
PRKDC on top (each with a zoom callout), EZH2 and BRD9 below. Every
quadrant is a vertical patchwork stack drawn entirely from raw bigwigs
plus MeanPSI `.links` with `../scripts/coverage_plot.R`. No
pyGenomeTracks, no ggbio, no raster crops — the reference
(`paper_figures/Figure3-NKAP_improved.pdf`) was assembled in
Illustrator and this rebuilds it natively so the type stays vector and
matches the rest of the figure.

`compose.R` sources this file and expects it to leave `panel_c`
defined. Layout above the panel (row heights, column widths) is
decided only in `compose.R`; everything inside the quadrant stack is
decided here.

## Decision rationale

**`sashimi_track_gg`'s arc engine is deliberately half-bypassed.**
`show_psi = FALSE`, and `apex_range` is *back-solved* rather than
passed as a plain tier spec. Two reasons:

1. **`anchor_bp = 0`.** The helper's 200 bp default pinned SENP2's arc
   end to the top of the exon *next door* (0.70 of `y_max`) instead of
   to the coverage at the splice site, leaving no room for an arc. X is
   never affected by this — only the endpoint height.
2. **PSI labels are placed here.** The helper's own solver only moves a
   label vertically, at its own arc's midpoint. These junction pairs
   share an endpoint, so both numbers want the same x, and in a track
   five text-heights tall there is no vertical slot that clears both
   arcs *and* the other label. They overprinted.

**Arc heights are relative to their own endpoints, never absolute.**
This is the bug that made Quinn reject the figure. The old code chose
each arc's height at its midpoint as a fixed fraction of `y_max`
(`ARC_MID <- c(0.26, 0.08)`), as if the arcs sprang from the baseline.
They do not — they spring from the coverage top at the splice site,
and because the *cryptic* site of each pair lies inside the exon
plateau, those tops reach 0.81 of `y_max` in PRKDC's controls. Every
chosen midpoint was therefore *below* its own endpoints,
`sin(t·pi) · (apex − max(y1,y2))` went negative, and each "arc" was
drawn drooping to the baseline and back up. The PSI labels, stacked
from `max(curve_mid)` upward, then floated in empty air well above the
drooping lines.

Now: each arc gets a target PEAK of `max(y1,y2) + ARC_RISE`, tiered
upward by at least `ARC_SEP`, capped by `ARC_CEIL`, and the sine
amplitude that puts the drawn peak exactly there is found with
`uniroot`. The peak is solved for, not assumed to be at the midpoint —
with a steep chord (PRKDC controls run 0.81 → 0.19) the true maximum
sits at t ≈ 0.35 and is a long way above the midpoint height.

**Arcs stack by endpoint height, not by junction length.** "Shortest
junction highest" is `sashimi_track_gg`'s convention and is right when
arcs start at the baseline. Here the *longer* junction of each pair is
the one anchored deeper inside the exon, so its chord dominates its
partner's at both ends. Forcing the short arc on top would mean
lifting it clear over the long arc's endpoint — out through the track
ceiling — or letting the two cross, which is exactly what the old
constants did at PRKDC. Verified: in all eight tracks the
higher-endpoint arc dominates, so the stack never crosses
(`arc_min_gap()` prints the signed minimum each build; all positive).

**PSI labels ride on their own arc and are separated HORIZONTALLY.**
`place_psi()` offsets each number one clearance from its own curve —
above it where that fits, otherwise below into the empty bowl — and
slides it along the arc until the box clears every other arc, the
coverage, the obstacles and the labels already placed. Horizontal room
is the cheap axis: an arc spans hundreds of bp and a number is about
50 of them, so two labels at different t do not interact at all.
Vertical room is what there isn't. The resulting rule reads
consistently across all eight tracks: *upper arc's number sits above
it, lower arc's number sits below it.*

**The zoom callouts carry an ARC EDGE, not a curved leader.** Quinn:
*"the zoom in 3c should have the edge of an arc and a label for the nt,
not a curvy label for the nt"*. `zoom_callout_gg`'s own annotation is a
decorative `geom_curve` hairline sweeping from the "+9 nt" text down to
the coverage; that is switched off (`nt_label = NULL`, `nt_site =
NULL`) and `zoom_arc_layers()` in this file draws the real thing
instead:

* **both junctions' arc FEET**, each anchored at its own splice-site
  coordinate on the zoomed coverage and rising away in the direction
  the junction actually runs (left for SENP2's acceptors, right for
  PRKDC's donors), in the same black stroke as the main track's arcs
  and at the same PSI-derived width `0.35 + LW_SCALE·sqrt(psi/100)`.
  The box is then literally the main track's arc pair magnified, and
  the reader can *see* the 9 nt / 6 nt offset as the horizontal
  distance between the two feet landing on the two coverage steps. In
  the control boxes the weights also tell the story on their own —
  SENP2's canonical foot is 96 % PSI (thick) and its cryptic foot
  1.38 % (hairline).
* **a plain, straight, horizontal nt label** in the headroom beside the
  feet. No leader at all: the feet already terminate on the sites, and
  the text abuts them.

Why one foot was rejected: the whole point of the zoom is that there
are two 3'ss 9 nt apart, and with a single foot drawn the "+9 nt"
string has no visible referent. Why the foot is curved at all: over the
56 bp the zoom shows, the *true* arc is 0.03 of `y_max` from straight,
so a faithful mapping would be a plain diagonal line — and the box's y
axis is the coverage rescaled to its own maximum anyway, a different
scale from the track's, so there is no faithful mapping to make. The
shape is therefore a visual convention: the last quarter of a sine arc,
`y = anchor + A·sin(s·U·pi/2)/sin(U·pi/2)`, which is exactly
`sashimi_track_gg`'s `sin(t·pi)` restricted to `t in [0.5, 1]`.
`ZOOM_ARC_U = 0.8` stops it short of the apex so the stroke is still
climbing where it meets the frame — it reads as an arc that continues
outside the box, which is true.

**The inline "-11 nt" / "-5 nt" is drawn once per locus, on the High
track.** The reference puts it just left of the splice site, over the
exon's descending shoulder. That still works for the High tracks here,
but not for the controls: `y_max` is the pipeline's ini `max_value`
and is driven by the *controls*, whose exon plateau fills 0.91 of it,
while half a line of `SZ_NT` is 0.11 of it. There is no air above a
control exon at all and the text printed **on** the grey coverage.
Note the reference has the same problem and solves it the same way —
EZH2's control track carries no "-11 nt" either.

Rejected alternatives for that label: the intron pocket right of the
site (it is a wedge that is shallow at the site, so the label had to
drift a tenth of the window inward and ended up hugging the underside
of an arc, reading like a second PSI number); a white-boxed label over
the coverage (punches a hole in a coverage track, which reads as a
dip); white text on the coverage (2.8:1 contrast at 6.8 pt).

## Gotchas & warnings

- **patchwork's `heights` govern PANELS, not cells.** `TRK_GAP_PT` (the
  5 pt breathing room under each coverage track, which stops this
  track's y-axis "0" overprinting the next track's y_max number) is a
  plot MARGIN, so it is spent *before* the rows are scaled: `C_SCALE`
  subtracts all four gaps from the budget, and `TRACK_H_IN` is then the
  plain `H_TRK · C_SCALE`. The same rule bites the zoom leaders — a
  `src_drop` of `H_TRK / H_CALL` lands on the High track's BASELINE,
  not on the Controls track's top edge; the gap has to be added
  separately as `(TRK_GAP_PT / 72) / CALL_H_IN` (see `ctrl_drop`).
- **Each zoom inset leads to its OWN track's source box**, which means
  the Controls inset's two dashed leaders cross the whole High track
  diagonally. Two consequences: (1) the inset nearer the source box in
  x must be the HIGH one, or the two pairs of leaders cross — that is
  why `box_hi`/`box_ct` are swapped between SENP2 (source at the right)
  and PRKDC (source at the left); (2) a PSI number solved without
  knowing about those leaders gets struck through (PRKDC's "22.1" was),
  so `leader_obst()` hands `place_psi()` a stack of thin boxes
  approximating each leader's corridor. The leaders still cross the
  High track's arcs — unavoidable for a line that has to reach a box
  two rows down — but they no longer touch any text.
- **`ZOOM_INNER0` is a COPY of a hardcoded number in
  `coverage_plot.R`.** `zoom_callout_gg` computes the box interior as
  `inner_lo = by_lo + 0.035·h`, `inner_hi = by_lo + cov_frac·h`, and
  0.035 is not a parameter. `zoom_feet()` has to reproduce the box
  geometry exactly in order to land the arc feet on the coverage, so
  `ZOOM_BOX_Y`, `ZOOM_COVF` and `ZOOM_NBINS` are *passed* to the helper
  (they cannot drift) but `ZOOM_INNER0` is mirrored by hand. If the
  feet ever float above or sink into the coverage, check that number
  first. `coverage_plot.R` is shared with the other figures, so this
  file must not "fix" it by editing the helper.
- **The feet must not cross.** They end on the same box edge, so the
  one whose site is further from that edge is the flatter of the pair;
  their END heights are additionally fanned by `ZOOM_ARC_SEP` in anchor
  order, because in the SENP2 controls the coverage step they spring
  from is only 0.08 of the box and without the fan the two strokes read
  as one doubled line. `zoom_feet()` measures the signed minimum gap
  between the drawn curves and `stop()`s if it goes negative; the gap
  is printed every build (currently 0.130 of a box height in all four
  boxes).
- **PRKDC's callouts are cramped and that is structural.** Its zoom
  window is `chr8:47858859-47858999` and the two donors sit at 86.4 %
  and 90.7 % across it, so the feet only have 9–14 % of the box width
  to climb in. Shifting `zoom_start` right would buy room but would
  also shrink the source box drawn on the tracks (it is *derived* from
  the zoom coordinates), and that box is approved. Left as is.
- **The nt label's y is shared across a locus's two boxes.**
  `zoom_arc_layers()` takes the max of the per-box floors. Solved
  independently, PRKDC's High label sat at 0.80 of its box and the
  Controls label at 0.73, and the two "-6 nt" read as misaligned.
- **`ANCHOR_BP` and the x of an arc end are different things.** The
  arcs' x must stay exactly on the `.links` coordinates — that has been
  the most painful recurring bug in this project. Nothing here touches
  x; `anchor_bp` only chooses which coverage value sets the endpoint
  *height*.
- **`apex_range` is only two numbers and the helper INTERPOLATES them
  linearly** over its own short→long ordering. With `n <= 2` that
  reproduces the solved peaks exactly; with three or more arcs it will
  not, and `arc_geom()` warns. If a locus ever gains a third junction,
  the fix is per-arc apexes in `coverage_plot.R`, not a fudge here.
- **`min_arc_rise` in `coverage_plot.R` must stay 0 for this panel.**
  It exists for exactly this sag, but it is a *floor* applied on top of
  whatever `apex_range` says — it would silently raise a back-solved
  apex and put the arcs through the labels that were placed for the
  un-raised ones.
- **Every text clearance is derived from `TRACK_H_IN`/`TRACK_W_IN`,**
  which are derived from `ROW_H[2]` in `compose.R`. Point sizes become
  data units through the track's *drawn* size. Change any `H_*` here or
  `ROW_H[2]` there and the solver re-solves correctly — but hardcoding a
  data-unit offset anywhere would go stale invisibly.
- **A ceiling test must not double-count the pad.** `ceil_y` is
  `y_max − (hh + 0.5·pad)`, i.e. the highest label *centre* whose
  glyphs stay in the track. Using the collision half-height there
  instead cost ~0.06 of `y_max` and the top label in four of the eight
  tracks found no slot at all and fell through to the fallback.
- **`for (t in ...) for (side in ...)` with one `break` is a trap.** In
  R that leaves only the inner loop, so the search returned the LAST
  acceptable slot rather than the first and every number drifted out to
  t = 0.78. The candidates are one flat list (`PSI_CAND`) for that
  reason.
- **Another agent's `compose.R` can overwrite `out/figure3.pdf` while
  you are looking at it.** Panels are cached to `panels/*.rds`, so a
  concurrent run that started before your rebuild will re-emit the PDF
  with the *previous* panel c. Copy the PDF aside immediately after
  rendering and check a known label position against the build log
  before trusting a measurement.
- **`panel_cached()` serialises with `saveRDS`.** An `aes()` that
  closes over a local scalar renders BLANK after `readRDS`. Every layer
  here carries its values as data-frame columns.
- **The PSI numbers are bare 3-significant-digit values** as the
  reference prints them, formatted one at a time — `format()` on the
  whole vector pads to a common width ("92.40").
- **Corrections to the reference that are deliberate:** SENP2 is chr3
  (the reference's "chr2" is a typo) and PRKDC's `y_max` is 9.433 (the
  reference's "5.4" is a hand-retype). Track captions are
  "High GSP #13 (15)" / "Matched Controls (150)".

## User preferences

- Stretch the plot rather than shrink the type. `H_TRK` is what it is
  because the PSI numbers set the floor on it; everything else in the
  stack was squeezed to pay for that. Do not reach for a smaller
  `SZ_PSI` when something collides.
- No label may overlap another label, a line or the coverage, and no
  line may cross text. The build log prints every label's position and
  every arc gap so this is checkable without opening the PDF — but
  check the PDF anyway, at 400+ dpi.
- Percentages/values belong on the object, not in a legend.

## Input/output

Reads, per locus, from
`flashier/pygenome/NKAP/<GENE>/<chrom>_<start>-<end>/`:

| file | what |
|---|---|
| `High_k13-3-..bw` | High GSP #13 coverage (15 samples) |
| `Matched_Controls-3-..bw` | matched control coverage (150 samples) |
| `MeanPSI.High_k13-3-..links` | `c1 s1 e1 c2 s2 e2 psi`, High |
| `MeanPSI.Matched_Controls-3-..links` | same, controls |

`High_k20-*` in those directories is a STALE factor key — do not use
it. Gene models come from `tcga/hg38_granges.RDS` filtered to the MANE
transcript and cached as `cache/panel_c_exons.rds` (the source RDS is
26 MB, so the cache is what the build actually reads).

Output: the `panel_c` object, cached by `compose.R` to
`panels/c.rds` and composed into `out/figure3.pdf`.

## File flow schematic

```
High_k13-3-..bw ─────────────┐
Matched_Controls-3-..bw ─────┤
MeanPSI.High_k13.links ──────┼──► arc_geom() ──► peaks + curves ──┐
MeanPSI.Matched_Controls ────┘        │                           │
                                      ├──► place_psi() ───────────┤
                                      └──► nt_place() ────────────┤
                                                                  │
zoom-region bigwigs ──► zoom_callout_gg(nt_label = NULL) ──────────┤
                    └─► zoom_feet() ──► arc edges + plain nt ──────┤
                                                                  ▼
hg38_granges.RDS ──► cache/panel_c_exons.rds ──► make_gene_rows() ─┼──► panel_c
MANE.GRCh38.v1.4.summary.txt ──┘                                   │
                                                                   ▼
                                                        panels/c.rds ──► out/figure3.pdf
```

## Running it

```
/scratch/midway2/qhauck/conda_env/rstudio-server/bin/Rscript --vanilla compose.R c
cp out/figure3.pdf /tmp/mine.pdf        # out/ is shared; crop a COPY
pdftoppm -r 400 -png -singlefile -x 0 -y 1180 -W 3400 -H 1360 /tmp/mine.pdf /tmp/c
```

The four zoom boxes, at 800 dpi, are

```
pdftoppm -r 800 -png -singlefile -x  720 -y 2600 -W 1000 -H 300 /tmp/mine.pdf /tmp/b1
pdftoppm -r 800 -png -singlefile -x 1810 -y 2600 -W 1010 -H 300 /tmp/mine.pdf /tmp/b2
pdftoppm -r 800 -png -singlefile -x 3590 -y 2600 -W 1010 -H 300 /tmp/mine.pdf /tmp/b3
pdftoppm -r 800 -png -singlefile -x 4690 -y 2600 -W 1010 -H 300 /tmp/mine.pdf /tmp/b4
```

and a cached-vs-fresh render must be **0 differing pixels** (`compose.R c`
then `compose.R`, crop both, compare) — that is the `panel_cached()`
`saveRDS` trap.
