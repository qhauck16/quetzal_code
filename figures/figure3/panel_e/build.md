# figure3/panel_e/build.R

## Purpose & context

Figure 3's panel e: "Mis-splicing leads to widespread effects on gene
expression". Four empirical CDFs of mean scaled RPM over the top-20
NKAP (GSP #13) samples — frameshifting cryptic-3'ss genes (n = 971),
in-frame cryptic-3'ss genes (n = 866), and an RPM-matched control set
for each — with sixteen DNA-repair / chromatin genes called out on the
frameshift curve and a one-sided Wilcoxon p.

Transcribed from `diff_splicing/combined_expression_analysis.Rmd`
(`build_set_pair` / `match_by_rpm` / `plot_single_panel_cdf`). Sourced
by `../compose.R`, which supplies `theme.R`'s tokens, the tidyverse,
ggrepel, `CACHE` and `panel_cached()`.

**The analysis is settled and verified.** The curves, the four series
and their wordings, n = 971 / n = 866, p = 4.2e-07 and the sixteen
drawn gene names are all checked against the source Rmd. Everything
below the data section is layout; do not "fix" the arithmetic.

## Decision rationale

### Type sizes — and the cell geometry they depend on

Everything here is set at the shared tokens: callouts at `SIZE_BASE`
(8 pt), legend at `SIZE_LEGEND` (8 pt). **Nothing below ever changes
those. They are the LAST constants to touch, not the first.**

**The dependency order, and it is not optional:**

1. re-measure `E_IN_PER_X` / `E_IN_PER_Y` and the drawn panel size off
   the composed render;
2. re-derive `E_PITCH` (it is a formula, not a number);
3. re-measure the legend block and choose its wrap;
4. re-cut `E_SPLIT` and re-derive the three bands and the three tier
   columns;
5. only then, if it is still impossible, consider `E_LAB_PT`.

History of the geometry (all at 400 dpi):

| `ROW_H[3]` | cell w | in/x-unit | in/y-unit | drawn panel |
|---|---|---|---|---|
| 2.45 | 3.40 | 1.4718 | 1.3937 | 3.40 x 1.72 in |
| 2.95 | 3.40 | 1.4718 | 1.7600 | 3.40 x 2.17 in |
| 3.05 | 3.23 | **1.0967** | **1.8400** | **2.535 x 2.267 in** |

The 2.45 -> 2.95 step was a *height* problem and Quinn's rule ("stretch
the plot rather than shrink the type") fixed it by growing the row; the
7 pt build was abandoned.

The 3.40 -> 3.23 step, when figure 3 went from four rows to three, was
a *width* problem, and it is much nastier: the panel lost 0.865 in of
drawn width, so `E_IN_PER_X` fell 25 % and **everything measured in
x-units got 34 % wider in data units at a stroke.** The legend block
went from 1.027 to 1.381 x-units and printed straight through CTC1,
CHFR, BRD2, BRCA1, ARID1A and RBBP8. Note that the row also got
*taller*, so the vertical budget actually improved — the breakage was
entirely horizontal, which is why "shrink the type" would not have
helped (see the legend section: at 7 pt the old wrap is still 1.24
x-units wide and still has no room beside it).

The budget, in the units that decide it (400 dpi, panel = 907 px tall
at 736 px per y-unit, 1014 px wide at 438.7 px per x-unit):

| item | at 8 pt |
|---|---|
| cap-height ink of one name | 33 px |
| `E_PITCH` (name + 0.030 in gutter) | 43.9 px = 0.1097 in = 0.0596 y-units |
| 15 names stacked | 658 px = 1.65 in |
| drawn panel height | 907 px = 2.267 in = 1.232 y-units |
| legend block | 1.032 x 0.983 in = 0.941 x 0.534 units |
| a 6-character italic name | 168 px = 0.420 in = 0.383 x-units |

`E_PITCH` is now **derived** — `(cap height + E_GUTTER_IN) /
E_IN_PER_Y` — rather than stated as a y-unit constant. Stating it is
exactly what went stale: 0.0625 y-units was 0.110 in at
`E_IN_PER_Y = 1.76` and would silently have become 0.115 in here.
`E_CHAR_EM = 0.630` (the italic label's em-per-character) is measured:
MAD1L1 renders 168 px wide at 8 pt. It is *not* the 0.56 used for
upright text elsewhere in this figure.

### Why three ladders rather than one column

The reference for this panel is `../../figure4/panel_g/build.R`, which
puts its ten gene names in one tidy right-aligned column on the free
flank of the curve with short leaders. The same idiom is used here,
but panel e has **sixteen** names, and one column of fifteen (DDB1 is
a separate corner case) is 1.65 in of stack against a curve whose free
flank is only ~1.1 in tall at any one x. So the names split into
**three tiers that step right as they climb** (`E_COL_A/B/C` =
-0.47 / -0.13 / 0.00), following the curve's left flank. Each column
is the largest value that keeps every one of its labels clear of the
frameshift curve at the label's box BOTTOM, and (for tier C) clear of
the legend on its left.

The tiers do **not** buy vertical room — they are y-disjoint, so the
panel still has to hold all fifteen boxes — they buy *horizontal
proximity*, which is what keeps the leaders short and nearly
horizontal. Bands are sized by label count (`E_BAND_A/B/C`), not by
where the points sit, and the slack is spent as the two inter-tier
gutters.

**The split is set by the legend, not by aesthetics.** Tiers A and B
sit in the legend's x-shadow (their left edges are at -0.79 and -0.45,
the legend ends at -0.452), so they have to fit *below* it: y
0.070..0.548, which is 0.478 y-units and holds **7** names. Tier C's
column is to the right of the legend, so it takes the other **8**.
The split is therefore `E_SPLIT = c(0.35, 0.55)` giving **5 / 2 / 8**,
where the previous, wider cell allowed 5 / 4 / 6 at `c(0.35, 0.70)`.

Both cuts still land on natural gaps in the label percentiles, which
the build prints every run:

```
DDB1 0.001 | MAD2L2 0.132 | PALB2 0.237 | TOPBP1 0.266 | FANCD2 0.292 |
BRD9 0.313 || XRCC5 0.413 | EZH2 0.504 || MAD1L1 0.579 | NHEJ1 0.629 |
RBBP8 0.766 | ARID1A 0.838 | BRCA1 0.846 | BRD2 0.919 | CHFR 0.953 |
CTC1 0.972
```

BRD9 0.313 | XRCC5 0.413 and EZH2 0.504 | MAD1L1 0.579 are the two
widest gaps below 0.70.

### Why the ladders are not drawn by ggrepel

ggrepel's equilibrium for a column packed against its band is "boxes
just touching", and `box.padding` is the only lever on the gutter.
Measured at 400 dpi on the old 1.72 in panel:

| box.padding | gutter between names |
|---|---|
| 0.06 | 2 px (0.005 in) |
| 0.11 | 4 px |
| 0.20 | 8 px, tiers at capacity |
| 0.22 | **boxes stop fitting the band; repel overlaps them silently** — all four of tier B printed as one block |

So the three ladders are laid out by hand (`ladder()`): even pitch,
point order preserved, stack centred on the tier's own points and
clamped into the band; `geom_segment` + `geom_text` instead of
`geom_text_repel`. That gives an exact, reproducible gutter and costs
nothing in robustness, because the clearance is audited (below). DDB1
is still a repel layer — it is one label in an empty corner and its
point is to the *left* of its text, which the ladder helper does not
model.

`ladder()` needs `diff(band) >= n * E_PITCH`. If a band is exactly at
capacity the clamp has no freedom and the stack is pinned to the band
rather than to its points, which lengthens that tier's leaders and can
push its lowest box down into the tier below's leaders. Keep ~0.01
y-units of slack in each band; the current ones have it.

### Why no leader crosses text (by construction, not by luck)

* Within a tier: every label's right edge is at `x_col` and every
  leader leaves from `x_col + E_SEG_GAP` heading **right**, so a
  leader cannot enter any box in its own tier. Point order is
  preserved, so leaders within a tier cannot cross each other either.
* Between tiers: a tier's leaders live between its labels' rows and
  its points' rows. At the current placement those spans are
  A 0.117–0.355, B 0.413–0.526, C 0.579–1.047, while the box extents
  are A 0.095–0.377, B 0.423–0.526, C 0.608–1.069. Each tier's highest
  leader clears the next tier's lowest box bottom by 0.046–0.082
  y-units (0.085–0.151 in). Tier C's lowest leader runs *down* to
  MAD1L1's point at y 0.579, below tier C's own lowest box, but tier
  B's boxes end at x = -0.13 and tier C's leaders start at x = +0.016,
  so they cannot meet.
* Against the curves: the frameshift CDF is the leftmost of the four,
  so a label is clear of all four iff that curve is right of the
  column at the label's row. The build prints the margin for all
  fifteen and `stopifnot`s it above 0.04 units. **Check the box's
  BOTTOM edge, not its centre** — the curve climbs left-to-right, so
  the bottom row is the tight one. Auditing the centre is exactly what
  once left MAD2L2's last glyph 0.007 in off the red curve.
* Against the legend: a label whose box overlaps the legend's rows
  must start right of the legend's right edge. That check is new — it
  is the one that was missing when the cell narrowed — and it is
  `stopifnot`'d above 0.03 x-units. It is printed in the same table.

### Legend — the wrap is load-bearing

Figure 4 panel g's legend geometry — inside, top left, no box,
`legend.key.height` 2.8 mm, `legend.key.width` 4.6 mm,
`legend.key.spacing.y` 0.4 mm, `lineheight` 0.85, text at
`SIZE_LEGEND` — with four entries instead of two.

**The wording is fixed; only the wrapping is a layout choice, and in a
2.535 in panel the wrapping decides whether the panel works at all.**
A legend is as wide as its widest LINE across ALL entries, so
re-wrapping only the two long "RPM-matched control genes …" entries
does literally nothing: `"Frameshifted (n = 971)"` is the widest line.
Measured as ink bounding boxes off 400 dpi renders of this legend at
this cell size:

| wrap | block |
|---|---|
| RPM entries 2-line, n entries 1-line (the old one) | 1.515 x 0.642 in |
| n entries 2-line, RPM entries 2-line | 1.407 x 0.812 in |
| n entries 1-line, RPM entries 3-line | 1.515 x 0.812 in |
| **all four entries wrapped short (current)** | **1.032 x 0.983 in** |

1.515 in is 1.381 x-units here, so that block reaches x = -0.013 and
there is simply no room between it and the frameshift curve for the
top tier — the top tier's left edge, with `E_COL_C` at its largest
legal value, is -0.273. Shrinking the legend type does not rescue it
either: the old wrap is 1.36 in at 7 pt and 1.27 in at 6.5 pt, both
still far too wide, which is why this rebuild kept 8 pt.

The narrow wrap costs 0.34 in of extra height (10 text lines instead
of 6), and that is what forces the 5 / 2 / 8 tier split above. The
block lands at x -1.393..-0.452, y 0.565..1.099; `E_LEG_TOP` (0.068 in)
is the drop from the guide box's top — which sits on the panel top —
to the first line's ink.

### p-value

`p = 4.2e-07` is right-aligned at the x limit (x = 0.80), exactly as in
figure 4 panel g. **Its row had to come down from 0.78 to 0.40** when
the cell narrowed, and that is measured, not taste:

* the string is 0.685 in = **0.625 x-units** wide here (it was 0.465),
  so parked at x = 0.80 it reaches back to x = 0.175;
* at y = 0.78 the rightmost of the four curves is already at x = 0.302
  (measured on the 400 dpi render), and RBBP8's leader runs out to its
  point at x = 0.20 — the leader printed straight through the "p";
* clearing the curves at that row would need a right edge at x = 0.955,
  which is past the panel.

At y = 0.40 the curves reach only x = 0.047, so the string clears them
by 0.128 x-units (0.140 in, measured 58 px at 400 dpi) and clears every
leader — tier B's longest ends at x = -0.02.

## Gotchas & warnings

- **Never put `limits =` on the scales.** A scale limit drops rows
  *before* `stat_ecdf` runs, so the drawn curve would be the ECDF of a
  truncated sample. The view is set with `coord_cartesian(clip="on")`.
  (Figure 2 panel b was burned by this; see `../../figure2/compose.md`.)
- **Everything positional is derived from a measured panel geometry**
  (`E_IN_PER_X` 1.4718, `E_IN_PER_Y` 1.7600, from the x ticks landing
  294.3 px apart and the y ticks 176.0 px per 0.25 at 400 dpi). If
  `compose.R`'s `ROW_H` or `COL_W` change, re-measure them first; the
  pitch, the bands and the legend block are all in those units.
- **`panel_cached()` round-trips the panel through `saveRDS`.** An
  `aes()` that closes over a local scalar renders BLANK after
  `readRDS`. Everything the layers need travels as a data column
  (`lab_x`, `lab_y`, `seg_x`); the one scalar-ish thing left is the
  per-row `nudge_x` vector for DDB1, which is fine because it is a
  plain numeric stored as a layer parameter, not tidy-eval.
- **ggrepel is blind to other layers.** The four step curves and the
  legend rectangle are fed to the DDB1 layer as zero-text obstacle
  rows. If another repel layer is ever added, feed it the same rows.
- **NHEJ1's italic J descends into the gutter above MAD1L1.** The
  measured band gap there is 3 px rather than 9–11; that is a descender
  doing what descenders do, not a collision. Do not "fix" it by
  opening the pitch — tier C is the tier with the least slack and it
  would be the one to suffer. (NHEJ1 and MAD1L1 are still adjacent
  after the re-split; they are now the two lowest names of tier C
  rather than the top of tier B and the bottom of tier C.)
- **Tier B has only two names and that is fine.** 5 / 2 / 8 looks
  lopsided in the source but reads as one continuous, evenly pitched
  column of names in the render, because all three tiers share the
  same pitch and the inter-tier gutters are small. Do not "balance" it
  without re-checking the legend clearance — tier B's third name has
  nowhere to go.
- **SMAD4 and DEPDC5 are in `HIGHLIGHT` but must not be drawn.** They
  fall in the in-frame set, not the frameshift set, so they have no
  position on the labelled curve. 16 of 18 survive; the build prints
  which two were suppressed.
- **`fread(header = TRUE)` on the cached matrices is load-bearing** —
  fread's header sniffing gets these wrong and returns `V1..V21`.
- The build is byte-reproducible (`panels/e.rds` md5 is stable across
  runs) despite `MatchIt::matchit`, because of `set.seed(E_SEED)`.
- The `Removed 2 rows containing missing values (geom_text)` and
  `cdf_place_labels fig3a: … 2 label-leader collisions remain`
  messages that `compose.R` prints come from **panel a**, not here.

## User preferences

- Quinn: "ensure no labels overlap, everything looks nice and is
  readable, no lines intersect/overlap text"; "check whitespace and
  how the pdf output actually looks, not making guesses".
- **"Stretch the plot, don't shrink the type."** This panel is the
  worked example, twice over. When the row was too short the 7 pt
  version was abandoned and `ROW_H[3]` grown from 2.45 to 2.95 in
  instead. When the cell was later narrowed to 3.23 in, the width was
  bought back from the legend's wrap and the tier split — and the
  measurements say shrinking the type would not even have worked (the
  old legend wrap is still 1.24 x-units wide at 7 pt), so 8 pt stayed.
  `compose.R` is not this panel's to change.
- Audit figures by measurement on the render, not by eye alone — the
  gutters, the legend block and the curve clearances here were all
  fixed from pixel measurements at 400 dpi.

## Input/output

```
cache/panel_e_scaled.tsv  ──┐   (20-col slice, scaled RPM: the curves)
cache/panel_e_raw.tsv     ──┤   (20-col slice, raw RPM: the matching)
cache/panel_e_top_ids.txt ──┼──► panel_e/build.R ──► panel_e (patchwork)
diff_splicing/unified_results/  │                        │
  nkap/20v60/dpsi0.01_fc2.5/    │                        ▼
  df_final.tsv            ──────┘                 panels/e.rds ──► out/figure3.pdf
```

## Verifying a change

```
cd quetzal_code/figures/figure3
/scratch/midway2/qhauck/conda_env/rstudio-server/bin/Rscript --vanilla compose.R e
cp out/figure3.pdf /tmp/mine.pdf        # out/ is shared; crop a COPY
pdftoppm -r 400 -png -singlefile -x 1190 -y 2520 -W 1300 -H 1220 /tmp/mine.pdf /tmp/e
```

Then *look at* `/tmp/e.png`, and re-run the measurements that matter:

* per-tier ink-band gutters — expect band counts of **5 / 2 / 8**, ink
  heights of 32–34 px and gutters of 9–11 px (0.023–0.028 in), with
  the single 3 px exception under NHEJ1's J;
* the clearance table the build prints — expect a curve minimum of
  **0.100** x-units (0.110 in) at MAD1L1 and a legend minimum of
  **0.069** x-units at MAD1L1 / ARID1A (the render is looser than the
  audit: 36 px = 0.090 in of real ink-to-ink air at MAD1L1, because
  `lab_w()` uses the nominal character box);
* the curve-to-`p` gap at the p-value's rows — expect 58 px (0.145 in);
* a cached-vs-fresh render must be **0 differing pixels** (`compose.R e`
  then `compose.R`, crop both, compare) — that is the `panel_cached()`
  `saveRDS` trap.

Calibration for that crop (400 dpi): x = 0 is column 840.5 at
438.67 px/x-unit, y = 0 is row 970.5 at 736.0 px/y-unit; the panel is
columns 223–1237 and rows 134–1041.
