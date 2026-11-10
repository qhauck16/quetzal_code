# figure4/panel_b/build.R

## Purpose & context

Panel b of Figure 4 (YLPM1 / GSP #30): *"Top junctions of GSP #30
capture usage of deep intronic splice sites."* Two loci stacked
vertically — **KDM4C** above, **TRIM38** below — each drawn as

```
High GSP #30 coverage + junction arcs
Matched Controls coverage + junction arcs
lavender taper + neck (taper_gg, LOCAL to this file)
full MANE gene model (transcript_zoom_gg)
genomic coordinate row (bespoke, same x scale)
```

plus a legend line for the orange "unannotated junction driving
GSP #30" arcs.

Reference being reproduced:
`tcga/genome_wide/flashier/paper_figures/Figure4-YLPM1_improved.pdf`,
row 1 right. `compose.R` sources this file and expects a `panel_b`
object; it is wrapped in `panel_cached("b", ...)`.

Everything is drawn from the raw bigwigs and `.links`. **No
pyGenomeTracks, no ggbio, no re-used `*_MANE_gene_model.pdf`** — Quinn
has rejected raster/pre-rendered gene models repeatedly.

Start from `figure3/panel_c/build.R`; this panel is the same idea
minus the zoom callouts, plus the orange annotation.

**`../scripts/coverage_plot.R` is NOT modified by this panel** and must
not be — figures 1–3 depend on it. The two things panel b needs that
it does not provide (a padded x window, a tapered-with-neck
trapezoid) are done locally: the x scale is *replaced* on the returned
ggplot, and `taper_gg()` lives in this file.

## Decision rationale

### The orange set comes out of the factorisation, not the geometry

Quinn, verbatim: *"TRIM38 orange junctions should be based upon the
TRIM38 junctions that have high intensities in the TRIM38 poisson NMF
junction group in the top 100 junction groups that make up GSP #30 —
same for the orange junctions of the other gene in that figure."*

The previous version derived the orange set geometrically (smallest
positive intronic gap with both facing ends unannotated) and got
TRIM38 **wrong**: it coloured 6.91 **and 22.7**. `build_gsp_tsv()`
replaces that with the real chain, cached in
`figure4/cache/panel_b_gsp30_juncs.tsv`:

1. **GSP #30's top junction groups** are
   `flashier/fact_30_top_junctions.tsv`. `spec_flash_fact.Rmd`'s
   chunk *"top 250 GSPs 10 13 30 unannotated junctions + splice
   types"* writes one row per feature of
   `rownames(norm_F)[order(norm_F[, 30], decreasing = TRUE)][1:250]`,
   **in rank order**, with
   `norm_F = t(diag(ldf$D) %*% t(ldf$F))` of the 3.7 GB
   `softimpute_flash_300_qc_filtered.RDS`. That loop can `next` past a
   feature, which would shift every later row's rank — the file has all
   250 rows, so nothing was skipped and row *i* is rank *i*. The code
   **asserts** that (`GSP_ROWS_EXPECTED`) rather than assuming it.
   → **KDM4C rank 40, TRIM38 rank 55.** TRIM38 has a *second* junction
   group at rank **114**, outside the top 100.
2. **Which `kN`** each group is: the tsv drops the `.kN` suffix but
   keeps the group's argmax UNANNOTATED junction, which identifies the
   factor uniquely. Recompute the per-factor argmax over the gene's
   unannotated junctions (`snaptron_output/{GENE}_snaptron.tsv`,
   `annotated == 0`) and take the column that reproduces it.
   → **KDM4C.k31, TRIM38.k21.**
3. **Which junctions are "high intensity"**: exactly
   `poisson_NMF_junction_plots.Rmd`'s own rule — every unannotated
   junction of the group within 2× of the group's top loading,
   `loadings[loadings >= top_val / 2]` on `poisson2multinom(res)$F`.

| locus | group | rank | orange junctions | loading | group top |
|---|---|---|---|---|---|
| KDM4C  | k31 | 40 | 11.2 (chr9:7,013,007-7,013,787)     | 0.1121 | 0.1217 |
| KDM4C  | k31 | 40 | 10.7 (chr9:7,011,879-7,012,910)     | 0.1013 | " |
| TRIM38 | k21 | 55 | 1.27 (chr6:25,973,285-25,978,534)   | 0.4254 | 0.4254 |
| TRIM38 | k21 | 55 | 6.91 (chr6:25,973,285-25,981,800)   | 0.2841 | " |

KDM4C reproduces the required 10.7 / 11.2 exactly (its group's own top
junction, 0.1217 at chr9:7,172,713-7,174,552, is 160 kb outside the
drawn window, so both in-window junctions clear top/2 = 0.0609
comfortably). TRIM38 reproduces the **6.91** Quinn confirmed, and adds
**1.27**, which is k21's single highest-loading junction — it cannot be
excluded by any reading of the rule, and it is the same kind of object
(an unannotated deep-intronic acceptor at 25,978,534). The two
junctions the old heuristic got wrong are now correctly black:
**22.7 is k24** and **8.38 is k22, which is rank 114** — outside the
top 100.

`res$F` rownames are **1-based**; `.links` is BED, so
`s1 = start - 1`, `e2 = end`. Delete the cache tsv to re-derive; it
needs `fastTopics` and reads two per-gene `res.RDS` (3 and 5.6 MB). It
does **not** need the 3.7 GB RDS.

### Arcs — rewritten 2026-09-23

Quinn rejected the previous version on two counts, verbatim:

1. *"ensure junctions in the sashimi plots map to coverage height, not
   to the bottom of the figure"*
2. *"you should align junctions in the sashimi plots to reduce overlap,
   notice how in my figures they can go under other junctions, or flip
   completely — fix that, they look like crap in 4b"*

#### Feet on the coverage

The old curve was `y = apex * sin(pi t)^p`, which is **zero at t = 0
and t = 1**. Both feet therefore sat on the track baseline no matter
where the coverage was — that is defect #1, exactly. The curve is now
the same parametric form `sashimi_track_gg()` and `figure3/panel_c`
use:

```
y(t) = (1 - t) * y1 + t * y2 + sin(t pi)^ARC_P * amp
```

with `y1`/`y2` the coverage **at** `s1` and **at** `e2`.
`ANCHOR_BP = 0`, as in figure 3 panel c: the foot takes the value of
the bin the junction coordinate falls in, which is a bin that is
*drawn*, so the foot sits ON the coverage silhouette and can never
float above it. Widening the window walks up into the exon plateau
(KDM4C's 66.2 right foot is 0.13 of `y_max` at the site and 0.66 fifty
bp in) and the foot would then hang a tenth of an inch over the fill.

The asymmetry this produces is real data, not a bug: `.links` is BED,
so `s1` is the last **exonic** base of the donor side and `e2` is the
last **intronic** base before the acceptor exon. KDM4C's 66.2 runs
0.53 → 0.13 for that reason.

**X is untouched.** Every foot is still drawn at the raw `.links`
coordinate. Re-verified on the render at **1600 dpi**, calibrating x
off the red coverage outline (which spans exactly `wide`) and y off
its baseline row: all four orange feet and every black foot not
overdrawn by an orange one land **within 0.7 px = 0.0004 in** of the
computed (junction coordinate, coverage height).

#### `ARC_P = 0.50`

Not cosmetic — it is what makes requirement 2 solvable. With a plain
sine the arc leaves each foot *linearly*, so a long arc is nearly flat
near its own feet and a short arc sharing that foot climbs straight out
through it. Raising the long arc does not help, because the sine's tail
stays flat: dominating TRIM38's 22.7 over its last 1.3 kb needed
amplitude 0.96, i.e. a peak of **1.09 of `y_max`**, off the track. At
`ARC_P = 0.5` the bump's derivative at the feet diverges, the tails are
steep, and the same constraint is met at amplitude 0.51.

0.5 is also right on its own merits: `sin(t pi)^0.5` is within 0.03 of
the half-ellipse `sqrt(4 t (1-t))` that pyGenomeTracks draws,
everywhere. The arcs read as sashimi arcs, land crisply ON their splice
sites instead of melting into the coverage, and carry a flat top that
the contour-style numbers can sit on.

#### Heights are SOLVED, not tiered

A plain tier ladder — tier *k* gets a peak `(k-1) * ARC_SEP` above tier
1 — is what `figure3/panel_c` does, and it works there because no
quadrant has more than two arcs. It **fails here, and not marginally**:
with the feet at coverage height an arc is LOW near its own feet, so a
short arc nested against the end of a long one pops out through it even
though its peak is far lower. TRIM38's controls are the clean example —
1.38 and 94.2 share the right foot at 0.20 of `y_max`, and over 1.38's
1266 bp the 9878 bp arc only climbs to 0.27 while 1.38's own peak was
0.30. It crossed by 0.074 of `y_max`.

So `arc_frame()` solves arc by arc, bottom-up in the **containment
order** (a linear extension of "a nested arc sits below the arc that
contains it", Kahn's algorithm with `(endpoint height, span)` as the
tie-break). Arc *i*'s amplitude is the smallest that satisfies all of

* `peak_i >= max(y1, y2) + ARC_RISE` — it reads as an arc, and it can
  never be *below* its own feet, so it can never "flip completely";
* `peak_i >= peak_j + ARC_SEP` for every already-placed *j* that
  contains *i* or is contained by it;
* `min` over their shared x of `y_i - y_j >= ARC_GAP`, same *j* —
  **dominance**, which is requirement 2.

Each is monotone in the amplitude, so the binding one is found by
bisection, and the drawn peak for a chosen target is back-solved with
`uniroot` (with a steep chord the maximum is a long way from t = 0.5 —
KDM4C's 10.7 runs 0.53 → 0.01 and peaks at t ≈ 0.34).

Settled values: `ARC_RISE 0.08`, `ARC_SEP 0.085`, `ARC_GAP 0.03`,
`ARC_TRIM 0.10`, and `ARC_CEIL` **derived** (see below). Found by
sweeping the four of them plus `ARC_P` over all four tracks and keeping
only the combinations with zero avoidable crossings.

**`ARC_CEIL` is derived, not chosen.** It is exactly `place_psi()`'s
own ceiling for a label *centre* (`lab_ceil_frac()`, shared by both),
so the top arc's apex is always a legal place for that arc's number.
This matters more than it sounds: at a round `0.93` the top arc spent
most of its length above the label ceiling, its number could only sit
near a foot, and at KDM4C that is inside the bundle where four arcs
leave one splice site together — the break then had two lines entering
and two leaving and named neither.

**Disjoint arcs are given no constraints at all.** That is what keeps
the budget affordable: KDM4C's 11.2 and 10.7 do not overlap in x, so
11.2 costs the four arcs above it nothing.

**`ARC_TRIM` (0.10 of the shared span, dropped at each end).** These
junction sets share endpoints — four of KDM4C's five start at
7,011,879 and three of TRIM38's seven end at 25,983,163 — and at a
shared foot the gap is 0 by construction, so including it would make
every dominance test unsatisfiable.

#### Which crossings remain, and why

`arc_min_gap()` prints the **signed minimum gap** between every
solvable pair on every build and `make_track()` **stops the build** if
it is not positive. Current state:

| track | fan | min signed gap | crossings |
|---|---|---|---|
| KDM4C high   | up   | **+0.0205** of `y_max` (1.92 over 4.87) | none |
| KDM4C ctrl   | up   | **+0.0179** (2.0 over 93.5)             | none |
| TRIM38 high  | up   | **+0.0300** (1.62 over 22.7)            | **2**, both forced |
| TRIM38 high  | down | **+0.0300** (6.91 over 1.27)            | none |
| TRIM38 ctrl  | up   | **+0.0300** (94.2 over 1.38)            | none |

TRIM38's high track was **4** forced crossings before the two orange
arcs were flipped below the baseline (2026-09-23); 1.62's partial
overlaps with 1.27 and 6.91 left the upper fan with them. The two
that remain are 1.62 against 2.13 and 8.38.

+0.0205 of `y_max` is 0.0063 in = 0.45 pt of white at the tightest
point, which is near a shared foot; everywhere else the tiers are
0.085–0.12 of `y_max` apart.

The two that remain are both **TRIM38's 1.62** (25,974,930–25,983,163)
against 2.13 and 8.38. Two arcs whose intervals *partially*
overlap — neither nested nor disjoint — must cross at some x whatever
their heights; that is in the junctions, not in the drawing. Nothing
else crosses anything.

A **nested** pair can also be forced, which is new with
coverage-anchored feet: if the inner arc's far foot sits below the
outer arc's curve there and they share the near foot, the inner arc can
only stay under the outer one if the outer one is steep enough at the
shared foot, and that lower bound can exceed the track. The solver
tries, and when it cannot it leaves the constraint unsatisfied and
reports the pair rather than spending the whole ceiling on an amplitude
that still does not clear. (Doing the latter — the first draft — put
five arcs on the ceiling and made *everything* above them cross.) With
`ARC_P = 0.5` no nested pair is currently forced.

#### TRIM38's orange arcs bow DOWN — 2026-09-23

Quinn, verbatim: *"in 4b the yellow junctions in TRIM38 should curve
underneath the gene body, its too crowded on top - you can lower the
height of the 94.2 curve on the controls to allow space for the yellow
junctions"*.

So TRIM38's high track now carries **two fans**: five black arcs above
the coverage and 1.27 / 6.91 below it. Nothing about the junctions
changed — same coordinates, same PSI, same orange set out of
`cache/panel_b_gsp30_juncs.tsv`, same `ANCHOR_BP = 0` feet on the
coverage. Only the sign of the bump changed.

**It is one solver, run mirrored.** `arc_frame(dir = -1)` negates
`y1`/`y2` and solves exactly as before, so `ARC_RISE`, `ARC_SEP`,
`ARC_GAP`, the containment order and the dominance bisection all keep
their meaning (deeper = "above"). `dir` rides on the frame and
`arc_curve()` / `arc_y_at()` map the answer back to the page; the
build log prints the down fan's feet as negatives because that is the
space it was solved in. A second implementation of the solver would
have been a second thing to keep in step with `panel_e`.

**The band is bought as row height, not by squashing the coverage.**
The high row is `H_TRK * (1 + DROP_F)` and its `coord` ylim is
`(-DROP_F * y_max, y_max)`, so all four tracks still draw their
coverage at exactly `TRACK_H_IN`. That is not cosmetic: the two
cohorts share a `y_max` and the panel's whole claim is high-vs-control,
so a 30% shorter TRIM38 track would be a real misreading. `B_SCALE`
carries the correction, so the whole stack pays ~6% and no derived
constant goes stale.

**`DROP_F = 0.52` is a floor, not a taste.** 6.91 has to dominate 1.27
from their shared left foot, which costs it a trough at 0.363 of
`y_max`; `DOWN_CEIL` must clear that plus a label half-height, i.e.
0.484. At 0.45 (tried) the down fan reported `6.91 under 1.27
[-0.006, ceiling]`.

**`DOWN_MIN` is the mirror-world floor.** `ARC_RISE` is measured from
an arc's own feet, and both of these have a right foot nearly on the
baseline, so 1.27 bottomed out at 0.145 — a hair too shallow for its
own number to sit fully under the baseline, and it was pushed off the
curve (`off = -0.172`). `DOWN_MIN` puts every down arc a label
half-height plus 2.4 pt below zero.

**The 94.2 arc CANNOT be lowered** — the one thing Quinn suggested
that does not work, so do not try it again. 94.2 contains 1.38 *and
shares its right foot*, so all that keeps 1.38 under it there is how
steeply 94.2 leaves that foot, which is its amplitude. 0.702 of
`y_max` is already the smallest apex that clears 1.38 by `ARC_GAP`;
capping it at 0.50 (tried, `CTRL_CEIL`) turned a solved pair into
`94.2 under 1.38 [-0.040, ceiling]`, i.e. the exact defect the arc
solver exists to prevent. The room comes from the taller high row
instead, and the rendered white band between the 6.91 trough and the
94.2 apex is ~0.13 in.

### The y axis is a rule with two ticks

Quinn, for every sashimi in the paper: *"the numerical axis should
have a line indicating it is an axis with two ticks for the numbers,
and they are currently misaligned, the 0 sits entirely above the
bottom of the track, it should sit beneath/even with it"*.

A 0.3 pt grey15 rule at the panel's left edge spanning **0 .. y_max**
(never the drop band — that is not on the scale), a 2.5 pt tick at
each end, and the numbers 1.5 pt further left. Chrome deliberately
identical to the `y_axis` block `sashimi_track_gg()` grew for figures
1–3; `SZ_AXIS == SZ_PSI` here, which is that function's own default
size.

- **The shared one is switched OFF** (`y_axis = FALSE`, passed only if
  the formal exists, so this file does not break if that work is
  reverted). It draws its rule at `start`, which in this panel is
  161 bp *inside* the window — `PAD_FRAC` exists precisely to keep
  TRIM38's first exon off the frame, and an axis there would stand in
  the coverage. Leaving both on printed two rules and two sets of
  numbers.
- **The rule, the ticks and the numbers are outside the drawn
  window**, so the x SCALE has to keep `AX_GUT` of gutter and `coord`
  sets the view. A scale limit drops data before it is drawn
  (`../figure2/compose.md`), which is why the old numbers were placed
  at `x = s` with `hjust = 1.25` instead.
- **The "0" is centred on its tick** — that is the fix Quinn asked
  for. **The top number hangs from its tick** (`vjust = 1`) rather
  than being centred: centred, half its glyphs sit above the panel,
  in the cell of the track above.
- **`AX_BOT_PT` (3.2 pt) is a gutter under KDM4C's high track only.**
  A centred "0" hangs half a text height below its own baseline, and
  at that one boundary the page row belongs to the control track's
  "6.1" — the two printed on top of each other. It is taken off
  `CONTENT_H` and handed back as row height so the coverage height
  does not change. The other three boundaries need nothing: TRIM38's
  high track has the whole drop band under its baseline, and a control
  track's "0" hangs over the lavender taper, which is a tint, not ink.

### PSI labels are CONTOUR-STYLE

Every number sits in a **break in the arc it names**, the way a
contour map labels an isoline. Labels are solved first; then each
arc's curve is split wherever it enters a label box and only the
surviving runs are drawn.

Why not a clear gap between two arcs: **it does not fit and cannot be
made to.** TRIM38 carries seven tiers across 0.48 × `y_max` of a
0.33 in track — 0.026 in apart — while the digits are 0.062 in tall.
Even KDM4C, with five, is 0.085 of `y_max` between tiers against a
label half-height of 0.111. A white glyph halo does not rescue it
either — an outline only thickens the glyphs, so a near-horizontal arc
still shows through the counters and between the digits and renders as
a **strike through the number**. That was tried and looked worse;
cutting the arc is the fix.

Nothing white is painted, so the PDF stays vector-clean and the
invariant holds exactly: no line ever touches a glyph.

**TWO break sizes — this is the whole trick for a crowded track.** The
arc the number *names* is cut generously (`bw`/`bh`); every other arc
the box touches is cut as tightly as "no line may touch a glyph"
allows (`bw2`/`bh2` = the glyph box + 0.5 pt), a hairline notch the eye
reads straight through. That is the contour-map convention: only the
labelled isoline breaks. Cutting every arc the same width is what made
KDM4C's "1.92" and "4.87" unreadable — **two** lines stopped on the
left of each number and **two** started on the right, and nothing said
which one the number belonged to. (They are only 0.085 of `y_max`
apart and the label's half-height is 0.111, so a neighbour is *always*
inside the box; it cannot be placed away from one.)

**Every arc's landing is protected.** The last `FOOT_KEEP = 5 pt` of
drawn arc at *both* ends of *every* arc is off limits to every label
box. Without it TRIM38's "22.7" — 0.7 kb of label on a 1.3 kb arc that
shares its right splice site with 52.4 and 1.62 — ate its own arc
whole, and moved off the arc it then sat over the last few points of
52.4 and 1.62 and both simply **stopped in mid-air short of the exon**.
That is requirement 1 breaking via the label solver rather than via the
geometry, and it is the subtlest failure in this panel.

The search is **scored, not first-fit**. Hard constraints: inside the
panel, under the ceiling, clear of the coverage, clear of every
protected landing, and no overlap with an already-placed label (boxes
plus 6 pt of horizontal and 1 pt of vertical air). Among the rest it
minimises

```
CUT_W * (other arcs the box would cut)
  + 3.0 * (label pushed off its own arc)
  + 3.0 * min(on-page slope of its own arc under the box, 2)
  + 0.9 * |t - 0.5|
```

The **slope** term is load-bearing: a break in a near-horizontal line
reads as a labelled line, a break in a steep one reads as a line that
stopped. `CUT_W` stays **low (1.5)** now that a neighbour's break is a
hairline notch rather than a hole — at 3.0 (tried) TRIM38's 2.13 and
8.38 left their own arcs to chase positions that cut nothing, and with
a cheap notch there is nothing worth chasing.

`arc_layer()` also drops any surviving run shorter than `min_pt = 3`
on the page: the sliver of arc left between two adjacent numbers reads
as a dash ("66.2 – 4.87"), not as a line.

**Placement order is MOST-CONSTRAINED FIRST**: increasing
`(arc span / label width)`. The number that barely fits on its own arc
chooses before the ones with a whole track to slide along. Placed last,
TRIM38's "22.7" found no legal slot at all and fell through to the
fallback; placed by tier order or by peak order it landed on 52.4's and
1.62's landings. First, it takes the one slot it has and the roomy arcs
work around it.

**A label may leave its arc by HALF a push.** `sgn` runs
`0, ±0.5, ±1` and costs `3.0 * |sgn| + 1.2 * (sgn != 0)`, so a seat on
the curve always wins and a number only steps off when it has none.
It was `0, ±1` until the drop band shortened every track by ~6%, and
TRIM38's **1.38** — a 25 pt arc carrying a 15 pt number between two
5 pt protected landings — then had no legal seat and jumped a whole
label height clear of its curve, reading as an unattached number. At
half a push the box still cuts the arc it names, so attribution
survives. Current state: fourteen of the sixteen numbers sit exactly
on their own arc (`off = 0.000`); **1.38 `+0.089`** and TRIM38's
**22.7 `-0.089`** take the half step, both still inside their own
break. The fixed 1.2 term is what keeps 22.7 from buying its way off
a good seat to dodge one hairline notch.

### The taper has a NECK

`trapezoid_gg()` in `../scripts/coverage_plot.R` runs the taper to the
bottom of its row, so lavender covers the whole row at every x the
transcript spans. The italic gene name that `transcript_zoom_gg()`
draws is ~0.076 in tall in a ~0.12 in row, so it **always** spills up
out of the gene row into the taper row — and landed on the lavender at
both loci. It cannot be made to fit: the free band above the exons
inside the gene row is only ~0.035 in.

`taper_gg()` (local) reaches the window width at `NECK_F = 0.55` of the
row height and then drops straight down as a neck. Below `NECK_F` the
row is empty everywhere outside the window, which is exactly where the
gene name sits (`label_side = "left"` at both loci; the window starts
well right of `tx_lo` at both). `make_gene_rows()` **checks that
budget and warns** if a size change breaks it.

The neck is also what fixed Quinn's defect #1. TRIM38's window is 42%
of its transcript, so a full-height trapezoid is a ~2° sliver whose two
sides read as unexplained risers rather than one funnel; KDM4C's
window is 1.2% of its transcript, which is why KDM4C's looked right and
TRIM38's did not. With the neck both read as the same object.

`H_TRAP` was also raised from 0.075 to 0.155 (0.25 of `H_TRK`, the same
ratio figure 3 panel c uses).

### `PAD_FRAC` — the clipped left peak

The bigwigs hold data for **exactly** the plotted region and nothing
outside it (`rtracklayer::import` range: chr6:25,973,078-25,984,559,
chr9:7,011,427-7,014,212). TRIM38's first exon block therefore starts
one base inside the region and, drawn hard against the frame, read as
a coverage peak sliced off by it (Quinn's defect #3). **There is no
data to un-clip** — widening the window would only paint zeros — so the
fix is to let the peak sit inside the frame instead of on it:
`PAD_FRAC = 0.014` of the region is added as empty gutter at both ends
of the DRAWN x range.

`view` (drawn) and `wide` (junction selection) are therefore different.
The junction filter uses `wide`, deliberately: the pad exists to
un-crowd the frame, not to pull TRIM38's 63.3 — a 96 bp intron sitting
one base outside — into the figure as an eighth arc. The taper's bottom
corners and the gene model's highlight use `view`, so everything that
has to line up still does.

### Other

**Lavender taper, not the reference's pink.** `#ccccff`, the same
colour figure 3 panel c uses. Change both or neither.

**The region comes from the Rmd, not the ini.** `GSP_pygenome`'s
`sashimi_tracks.ini` records `max_value` and the track titles but NOT
the plotted region. The `--Region` strings are in
`flashier/poisson_NMF_junction_plots.Rmd`
(`prep_pygenometrack('TRIM38', 'chr6', '+', 'chr6:25,973,077-25,984,559')`,
`prep_pygenometrack('KDM4C', 'chr9', '+', 'chr9:7,011,426-7,014,212')`),
and they agree with the bigwigs' own extents to the base.

## Gotchas & warnings

- **The published cohort labels are stale on BOTH counts.** The
  reference reads "High GSP #13 (15)" / "Matched Controls (150)". This
  is the YLPM1 figure, the bigwigs are `High_k30-*`, and the
  `GSP_pygenome` inis say `High_k30 (13)` / `Matched_Controls (130)`.
  Panel e's pygenome run really is 15/150 — panel b's numbers were
  copied from it. Corrected here to **High GSP #30 (13)** /
  **Matched Controls (130)**.
- **Ignore every `High_k65-*` file.** k65 is a stale factor key for
  YLPM1; the current key is k30.
- **`panel_cached()` round-trips through `saveRDS`.** Any `aes()` that
  references a scalar from the enclosing environment resolves to
  nothing after `readRDS` and the whole panel renders blank. Arc
  colour, the path group id and every label box travel as DATA COLUMNS.
  **Verified**: a cached-only compose and a fresh build of panel b
  render pixel-identically (0 differing pixels at 200 dpi).
- **`/` between two column patchworks flattens them.** `build_locus()`
  returns a flat LIST of rows and the whole panel is one
  `Reduce(`/`, rows)` with one `heights` vector.
- **`B_TITLE_H` changes nothing about the layout.** patchwork heights
  are relative. Its only job is to make `TRACK_H_IN` the track's true
  drawn height, which is what the label solver converts point sizes
  through. Re-measure after any change to the title or row heights.
- **`TRACK_W_IN` is ~1.8% optimistic** — measured 1215 px of drawn
  panel at 400 dpi against a nominal 1237. It only scales the
  point→data-x conversion, so it is harmless, but do not treat it as
  exact if you start placing things at the panel edge. (It is ~3.09 in,
  not the ~1.45 in you get if you forget `COL_W$ab` — panel b is the
  *narrow* column of a 1.15 : 1 split of 8.5 in.)
- **`H_LEG` was raised 0.115 → 0.16 to stop the legend line being
  clipped.** At 0.115 the row drew 0.097 in while the 6.8 pt caption
  needs 0.088 in ascender-to-descender, and with the text centred at
  0.33 of the row the descenders of "junction driving" fell out of
  panel b's cell and were cut by the row below. The caption now sits at
  0.50 with the glyph on its baseline, and `legend_row` **warns** if
  the budget breaks again. Cost: 1.8% of the track height.
- **Only ever LOOK for a clipped descender — a build log will not tell
  you.** ggplot's `clip = "off"` does not stop patchwork's cell from
  clipping.
- `bigwig_binned()` uses 2000 integer-width bins: 1 bp at KDM4C
  (2786 bp region), 5 bp at TRIM38 (11 482 bp). Trust the returned bin
  centres; never re-derive an x position from them.
- **`place_psi()` warns** when it cannot find any legal slot and falls
  back to the arc's midpoint. **The fallback ignores every hard
  constraint**, including the protected landings — if it fires, expect
  an arc to stop in mid-air. It is currently silent for all four
  tracks, and the build log's `off=` column (0.000 for all sixteen
  numbers) is the second half of the same alarm.
- **`make_track()` STOPS the build** if any solvable pair of arcs
  crosses. Forced crossings (partial overlap, or a nested pair the
  ceiling cannot separate) are listed, not fatal. The check runs
  **once per fan** (`up` / `down`) and the log line names which, with
  the crossing COUNT. Read it every build — that is the regression
  check for Quinn's requirement 2.
- **A down fan's `y1`/`y2`/`amp`/`peak` are MIRRORED**, in the build
  log as well as on the frame. `dir` is the only thing that knows;
  anything that needs a page y must go through `arc_curve()` /
  `arc_y_at()`, never read `y1` directly.
- **`bh` and `bh2` are SCALARS, `bw`/`bw2` are VECTORS.** `bh[i]` for
  `i > 1` is silently `NA` and the first `if()` that sees it dies with
  "missing value where TRUE/FALSE needed".
- `make_gene_rows()` warns if the gene name no longer fits in the clear
  strip left of the neck, or spills past it vertically.
- The arc/label machinery is deliberately DUPLICATED in
  `panel_e/build.R`. `../scripts/coverage_plot.R` is shared with
  figures 1–3 and must not be edited for one panel's benefit. Keep them
  in step by hand — panel e has NOT been given the contour labels.

## User preferences

- Stretch the plot rather than shrink the type.
- Thin curved arcs scaled by `sqrt(PSI)`
  (`0.25 + 0.60 * sqrt(psi/100)` pt), never filled shapes.
- Labels on the object: every PSI number sits in a break in its own
  arc, and orange numbers are drawn in the arc's colour.
- No label touches an arc, another label, or the track above/below.
- Arcs must land on the COVERAGE at their own splice site, never on the
  track baseline, and must never bow downward or pass under another
  arc. Both are checked numerically on every build.
- Only `FONT_FAMILY` and the `SIZE_*` tokens from `../theme.R`.
- No pyGenomeTracks, no ggbio, no pre-rendered gene-model PDFs.

## Input/output

Reads, per locus, from
`<FLASHIER>/GSP_pygenome/YLPM1/{KDM4C,TRIM38}/`:

| file | used for |
|---|---|
| `High_k30-3-..bw` | High-cohort coverage |
| `Matched_Controls-3-..bw` | control coverage |
| `MeanPSI.High_k30-3-..links` | high arcs (7 cols: c1 s1 e1 c2 s2 e2 psi) |
| `MeanPSI.Matched_Controls-3-..links` | control arcs |
| `sashimi_tracks.ini` | `max_value` (the shared y scale) and the track titles/sample counts |

plus, for the orange set,
`<FLASHIER>/fact_30_top_junctions.tsv`,
`tcga/genome_wide/260116_filters/chr{N}/FastTopics_output/{GENE}/res.RDS`
and `tcga/all_genes/chr{N}/snaptron_output/{GENE}_snaptron.tsv`,
distilled once into `figure4/cache/panel_b_gsp30_juncs.tsv`;
and, for the gene models, `tcga/hg38_granges.RDS` +
`<FLASHIER>/diff_splicing/MANE.GRCh38.v1.4.summary.txt`, distilled once
into `figure4/cache/panel_b_exons.rds`. Delete either cache to rebuild
it.

Writes nothing but the `panel_b` object; `compose.R` serialises it to
`figure4/panels/b.rds` and renders `figure4/out/figure4.pdf`.

y `max_value` verified from the inis: **KDM4C 6.1**, **TRIM38 45.9**.

## File flow schematic

```
High_k30-3-..bw ─────────┐
Matched_Controls-3-..bw ─┼─► bigwig_binned() ──► coverage area
                         │
MeanPSI.*.links ─────────┴─► read_links() ──► arc_frame() ─┬─► place_psi() ──► PSI labels
                                                           └─► arc_layer()  ──► arcs, cut at the labels

fact_30_top_junctions.tsv ─┐   (GSP #30 top-100 groups, in rank order)
{GENE}/res.RDS ────────────┼─► build_gsp_tsv() ──► cache/panel_b_gsp30_juncs.tsv ──► orange_keys()
{GENE}_snaptron.tsv ───────┘   (poisson2multinom F, unannotated, >= top/2)

hg38_granges.RDS ──┐
MANE summary.txt ──┴─► cache/panel_b_exons.rds ──► transcript_zoom_gg() (gene model)

                      LOCI$view ──► taper_gg() (trapezoid + neck) ──► gene-model highlight

                       all of the above ──► panel_b ──► compose.R ──► out/figure4.pdf
```

## Running it

```
cd /project2/yangili1/qhauck/spliceNMF/quetzal_code/figures/figure4
/scratch/midway2/qhauck/conda_env/rstudio-server/bin/Rscript --vanilla compose.R b
pdftoppm -r 400 -png -singlefile -x 1750 -y 10 -W 1680 -H 1070 out/figure4.pdf /tmp/b
```

then look at `/tmp/b.png`. For per-track detail, 900 dpi:

```
pdftoppm -r 900 -png -singlefile -x 4230 -y  330 -W 3050 -H 350 out/figure4.pdf /tmp/kdm_hi
pdftoppm -r 900 -png -singlefile -x 4230 -y  630 -W 3050 -H 330 out/figure4.pdf /tmp/kdm_ct
pdftoppm -r 900 -png -singlefile -x 4230 -y 1290 -W 3050 -H 350 out/figure4.pdf /tmp/trim_hi
pdftoppm -r 900 -png -singlefile -x 4230 -y 1610 -W 3050 -H 340 out/figure4.pdf /tmp/trim_ct
```

To re-verify the arc feet, render the whole width of one track at
1600 dpi, calibrate x off the red coverage outline (it spans exactly
`wide`) and y off its baseline row, and look for ink at each computed
(junction coordinate, coverage height):

```
pdftoppm -r 1600 -png -singlefile -x 6700 -y  600 -W 6900 -H 600 out/figure4.pdf /tmp/foot
```

If other agents are composing figure 4 at the same time they will
truncate `out/figure4.pdf` under you mid-`pdftoppm` ("Document stream
is empty"). Render to your own path instead — `compose.R` only needs
its `out_pdf` swapped.
