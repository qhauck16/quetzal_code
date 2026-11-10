# figure3/panel_b/build.R

## Purpose & context

Panel b of Figure 3 (NKAP / GSP #13): *"Mutations cluster in specific
regions of NKAP and PRPF8."* Two stacked protein cartoons — **NKAP**
over **PRPF8** — with domain boxes on a **broken** amino-acid axis and
a stack of lollipop circles (no stems) sitting on each body's top edge,
one circle per carrier sample, coloured by TCGA project. Two prose
annotations (a bracket over the NKAP stacks, and a captioned pair of
stubs for the circles hanging *below* the bodies), the WXS/RNA-seq
footnote, and a 16-entry project legend set as a **4 × 4 block in the
dead space above NKAP's label staircase**.

Reference being reproduced:
`tcga/genome_wide/flashier/paper_figures/Figure3-NKAP_improved.pdf`,
row 1 right. That reference was drawn **by hand in Illustrator**.
Nothing here is hand-placed: stack membership comes from the TSV, label
text from a prefix-collapse, x positions from a piecewise aa→inch map
plus two collision sweeps, and every font size, circle diameter and
vertical band from one auto-fit scale factor.

`../compose.R` sources this file and expects a `panel_b` object; it is
wrapped in `panel_cached("b", ...)`.

## Decision rationale

**The panel is drawn in INCHES.** x and y data units are both inches
and `coord_cartesian` pins the view to exactly the drawn size, so a
circle of radius *r* data-units is *r* inches on the page and a text
offset can be derived from a point size without a unit conversion at
every site. This is the whole reason the geometry reads as arithmetic
rather than as tuning.

**`B_CAL_W` / `B_CAL_H` are measured, not computed.** The
cell-minus-margins arithmetic is only an estimate of the drawn panel —
patchwork adds its own margins around a `wrap_elements()` and the title
band is not exactly two lineheights. The two calibration factors
correct it off the 400 dpi render (px per data unit ÷ 400). The
`[panel b] cal:` line printed at build time is the handle:
re-measure if the row height, the column split or the title changes.
Same practice as figure 2's title alignment — see
`../../figure2/compose.md`.

**One auto-fit scalar `S` drives every size** (`geom_for(S)`, bisected
in section 6). `h_req(S)` is the exact sum of the bands the content
needs; the bisection finds the `S` at which that equals the cell. So
"is there dead space?" has a code answer, not an eyeball one — if a
band over-reserves, `S` comes out too small and the type shrinks for
nothing. That is exactly what the caption row used to do (below).

**Text metrics are measured at 100 pt and scaled down.** `cairo_pdf`
reports string widths rounded to WHOLE POINTS, which at 4 pt type
quantises every measurement to ~15% — enough to make the bisection land
on a wrap discontinuity. `str_w()` measures once at 100 pt (rounding is
then 1 part in 1000) and scales; font metrics are linear in size. It is
memoised because the bisection re-measures ~25 times.

**The legend is a 4 × 4 block top-left, not a column on the right.**
This is the 2026-09-23 restructure and it is the reason the panel is no
longer mostly empty. A 1 × 16 column cost 0.67 in — **17% of the cell**
— off the right-hand edge, which (a) starved both cartoons of width,
(b) forced the broken axis to spend 30% of NKAP on a region with no
mutations, and (c) still left the whole upper-left quadrant blank. The
block now lives IN that quadrant, so it costs no height at all and both
proteins get the full `B_W_TOT`. Its type size is solved in section 7
against the width NKAP's label staircase actually leaves free
(`free_x()`), and its columns are justified across that width.

**Broken aa axis.** NKAP's mutations are all in aa 300–415 and PRPF8's
in 1570–1900, but both cartoons must still show the whole protein.
`B_PROTEIN$<gene>$brk` / `$frac` give the piecewise map: NKAP spends
24% of the width on 0–272 and 76% on 272–415; PRPF8 15/24/47/14 on
0–1579 / 1579–1754 / 1754–2000 / 2000–2305. Ticks are placed on GROUP
EDGES, so a tick can never disagree with the box it sits under.

The fractions were re-cut with the legend column gone. PRPF8 gained a
FOURTH segment at 1754 (already a tick and a domain boundary) because
nothing is mutated between 1612 and 1756: at the old 20/60/20 that gap
drew as a 0.9 in empty corridor running the full height of the panel.
Compressing it to 0.24 of the width narrowed the corridor to ~0.5 in
and handed the 1832–1881 cluster the difference. The cost is that
"Endonuclease-like Domain" auto-shrinks a little further inside its
box; don't push `frac[2]` below ~0.22 or it hits the 3.2 pt floor.

**The top-60 rule.** The circles ABOVE a body are the carriers among
the **top 60 samples by GSP #13** (`B_TOPN`), which is what the
bracket's caption says. A sample carrying two changes in the same
window (866663 = R402Q + L313M, 111261 = L1872P + S1887L) contributes
**one** circle, at the more recurrent residue — that is what removes
L313M and S1887L, which the reference does not draw. The four circles
BELOW a body are samples *outside* the top 60 that carry one of the
plotted changes; `stopifnot(nrow(below) == 4L)` is the regression
alarm.

**Stacks are overlap components of residue SPANS, not equal positions.**
`P1880_N1881del`, `N1881_V1883dup` and `N1881insKN` all touch residue
1881 and so share a column (as the reference draws them); N371K and
Q372H are 1 aa apart but disjoint and stay separate. Positions are read
only from `<AA-letter><digits>`, so the "13" in `R402Sfs*13` is not
mistaken for a residue.

**The stack label is a progressive prefix collapse**, again as the
reference does: `R402*` + `R402Q` + `R402Sfs*13` → `R402*/Q/Sfs*13`. A
prefix is only collapsed when it is at least 3 characters and never to
nothing.

**Label wrapping breaks only at `/` and `_`.** Those are the real token
boundaries in an HGVS-ish protein change. `/` separates merged changes
and is the only break a reader can re-join unambiguously, so it is
always tried first; `_` (a residue RANGE inside one change) is used
only when a single `/`-chunk is still too wide on its own.

**ONE x collision sweep, then one clamp.** The sweep takes the max of
(a) circles may not overlap and (b) two rotated labels are parallel
lines and collide unless their PERPENDICULAR separation clears half a
text block each. Then (c): the sweep only pushes right, so the set is
slid back by the overshoot past `B_X_RIGHT` — but never left of the
functional window's own start, and whatever the slide cannot absorb is
returned as `xover` for the auto-fit to shrink away.

**Where the "outside the top 60" caption goes is SOLVED, not chosen.**
It wants to sit under the two circles it annotates so its leaders are
stubs. But it shares the middle band with PRPF8's label staircase, and
if it overlaps PRPF8's tallest labels the `drop` term flips branches
and PRPF8's whole block falls ~0.2 in. Section 5 therefore slides the
caption left from "centred on the circles" in 0.02 in steps until the
labels it overlaps are short enough that `bP[above]` still binds —
i.e. until the position is free. Both it and the footnote are then
CENTRED in their own free slot rather than hung from the band's top.

## Gotchas & warnings

- **`k13` is `GSP #13`.** The Rmds call the loadings `k`; the paper
  calls them GSPs. Everything user-visible here says GSP.
- **`panel_cached()` round-trips through `saveRDS`.** Any `aes()` that
  references a scalar from the enclosing environment renders BLANK
  after `readRDS`. Everything an `aes()` touches here is a data COLUMN,
  including the per-tick `hjust`; keep it that way. See
  `../../figure2/compose.md`.
- **The caption row is NOT a row.** "TCGA samples carrying a plotted
  mutation, outside the top 60" is two lines about 1.2 in wide sitting
  under NKAP's below-body circles. Reserving a full-width band for it
  (the old `note` term in `blk()`) left ~0.2 in of dead white across
  the whole panel between the two cartoons, AND shrank `S` by ~8%
  because the bisection had to pay for it. Now PRPF8's block rises
  until either its own tallest label reaches the caption row's top, or
  a label whose ink actually OVERLAPS the caption's x span reaches the
  caption's bottom (`drop`, section 5).
- **`capx1`/`capx2` are solved in section 5 and drawn in section 7.**
  The caption's x must be the same number in both places or the
  clearance the layout paid for is not the clearance that gets drawn.
  Do not re-derive it at the drawing site.
- **The caption's leaders are stubs off its TOP edge, one per circle**,
  vertical where the circle's x falls inside the caption and a short
  diagonal off the nearest top corner where it does not. They cross
  NKAP's tick row, which is safe *because* the numbers sit on the group
  edges (0 / 272 / 415) and nothing else is in that band. Two other
  shapes were tried and rejected against the render: a horizontal
  spine under both circles drew 0.06 in below the protein body and
  parallel to it, so the body read as double-ruled; a rule-plus-risers
  comb put a 2.4 in horizontal line across the gutter that read as a
  divider between the two proteins.
- **Reserve `max(n*cd + ink)`, never `max(n)*cd + max(ink)`.** The
  tallest stack and the longest label need not be the same stack; the
  second form over-reserves whenever they differ, and it yields no
  per-stack extents — which are what the caption overlap test needs.
- **The end tick numbers are aligned INWARD** (`hj = c(0, …, 1)`). A
  centred "0" hangs half its width outside the cell (the body starts at
  x = 0) and a centred "2305"/"415" hangs half outside the right edge.
  `gname$x` depends on this: it clears the FULL width of the "0", not
  half of it.
- **The legend's size solve must BISECT, not shrink.** `need(pt)` grows
  with the type while `free_x(pt)` *shrinks* with it — a taller block
  reaches further down into the staircase — so a proportional step
  overshoots, and a shrink-only loop then exits at the first size that
  happens to fit. It settled on 3.93 pt where 5.13 pt fits. There is
  exactly one crossing; bisect for it.
- **The legend's ceiling is `G$note_pt`, not `SIZE_LEGEND`.** The free
  width alone allows ~5.6 pt, which sets the key LARGER than the
  mutation labels it explains. Legend text now matches the prose.
- **`B_X_RIGHT` is the drawing area's own right edge less a hairline,
  and the fit enforces it.** Sweep (c) can only slide the set back as
  far as the functional window's start, so it can FAIL; before `xover`
  was fed into the bisection, `R402*/Q/Sfs*13` overran by 0.22 in and
  printed straight through the TCGA-UCEC and TCGA-KIRC legend entries.
  `reach` is the label's true rightmost ink, so this is a real
  clearance and not a fudge.
- **The drawn lineheight `B_LAB_LH` (0.88) is deliberately TIGHTER than
  the `lh` the collision model assumes (1.0).** The lines of one label
  must sit closer together than any two labels do, or a 3-line label
  stops reading as one object. The `0.85 * lh` white-space term in
  sweep (b) enforces the other side of that ordering; below it,
  "insQ" and "P1880_…" ran together as one string.
- **Circles are tangent by construction — and the two sweeps MUST stay
  fused into one pass.** They used to run as two loops, and that
  silently broke the circle constraint: the label pass moves `x[i-1]`
  right *after* the circle pass has spaced `x[i]` off its old position,
  so wherever the perpendicular term came out below `cd` — which
  happens whenever the stack gets SHORTER (`dy < 0`) — the two stacks'
  circles ended up overlapping. N371K's single red circle sat on top of
  R368Q's blue one in every render up to 2026-09-23. One pass taking
  `max(x[i], x[i-1] + cd, x[i-1] + perp)` fixes it. Circles are drawn
  at `r = cd/2 * 0.94`, so neighbours leave a 0.06 `cd` gap.
- **The asterisk marks a residue ALSO called from WXS or RNA-seq**, not
  one called only from the MAF/AlphaMissense column — that is what the
  reference's footnote means. Shape 8's drawn diameter is ~0.0296 in
  per unit of `size`, so it is sized as a fraction of the lollipop it
  sits on rather than in points.
- **The stack order bottom→top is descending GSP #13**, so the
  strongest carrier sits on the body. The reference's own order is not
  reproducible from the current loadings.
- **The bracket's top is clamped under the legend block** (`leg_bot`),
  and `cap1` is set on ONE line whenever the room left of the bracket
  takes it. Both are measured, not eyeballed: cap1 is centred on the
  bracket, so its ink top is `(br_lo + br_hi)/2 + cap1_h/2` and that is
  what the clamp solves for.
- **The whitespace budget, measured.** Largest empty axis-aligned
  rectangle in the drawing area, from a 100 dpi ink scan of the 400 dpi
  render (`blank.R` in the scratchpad — binarise at 0.95, block-max
  downsample, histogram/stack maximal-rectangle):
  **2.92 × 0.22 in = 0.642 in²**, the gutter strip between the two
  prose blocks and PRPF8's labels. Before the restructure it was
  3.09 × 0.40 = 1.24 in² in the upper-left. Widest blank ROW band
  0.05 in, widest blank COLUMN band 0.09 in. The residual gutter is
  `drop - note_h`, i.e. it is set by how tall PRPF8's 3-line labels
  are; narrowing `B_LAB_MAXW0` to 0.52 does shrink it to 0.53 in² but
  breaks "H1875_L1876del/insQ/P" onto four ragged lines with a single
  "P" hanging — not worth it. Re-run the scan after ANY change here.

## User preferences

- Stretch the plot rather than shrink the type: every band is the exact
  ink it holds so the auto-fit can spend the slack on size.
- Labels on the object; percentages/counts never in a legend.
- Only `FONT_FAMILY` and the `SIZE_*` tokens from `../../theme.R`.
- Look at the rendered PDF, never at a prediction of it. The two
  diagnostic `cat()` lines (`fit`, `cal`, `legend`) exist so a render
  can be checked against what the layout thought it was doing.

## Input/output

Reads one file:

```
<FLASHIER>/NKAP_full_plus_variants.tsv
```

one row per sample, with `rail_id`, `k13` (the GSP #13 loading),
`project_id`, and three free-text mutation columns `gene_mutation`,
`wxs_mutation`, `rnaseq_mutation`. Each is a `;`-separated list of
`GENE:change1,change2` blocks; nothing is pre-parsed. Parsing lives in
`../scripts/nkap_variants.R` (shared with panel a) — note that only ONE
source wins per row in `gene_mutation`, so all three columns have to be
searched or a PRPF8 hit can be hidden behind a FAM50A call.

Functional windows: NKAP aa 300–413, PRPF8 aa 1570–1900 (`B_WIN`).

Writes nothing but the `panel_b` object; `compose.R` serialises it to
`figure3/panels/b.rds` and renders `figure3/out/figure3.pdf`.

## File flow schematic

```
NKAP_full_plus_variants.tsv ──► scripts/nkap_variants.R
                                  (gene_changes / extract_all_positions)
                                        │
                                        ▼
                                long  ──┬──► top   (rank <= 60)  ──► stacks ──► circles above
                                        └──► below (rank >  60)  ─────────────► circles below
                                        │
B_PROTEIN (domains, breaks) ──► aa2in ──┴──► x map ──► collision sweeps ──► st$x
                                                                              │
                              geom_for(S) ──bisect──► S ──► every size ────────┤
                                                                              ▼
                                                        panel_b ──► compose.R ──► out/figure3.pdf
```

## Running it

`out/` is shared with the other panels' agents — copy the PDF before
cropping it, or a concurrent `ggsave` will hand `pdftoppm` a truncated
stream.

```
cd /project2/yangili1/qhauck/spliceNMF/quetzal_code/figures/figure3
Rscript --vanilla compose.R b
cp out/figure3.pdf /tmp/mine.pdf
# whole cell, incl. title + tag
pdftoppm -r 400 -png -singlefile -x 1770 -y 0   -W 1630 -H 1190 /tmp/mine.pdf /tmp/b
# DRAWING AREA only (y = B_H_TOT .. 0) -- this is what a blank scan must use
pdftoppm -r 400 -png -singlefile -x 1770 -y 146 -W 1630 -H 1030 /tmp/mine.pdf /tmp/draw
```

Panel b is row 1 right: page x 4.42–8.5 in, y 0–2.95 in. The drawing
area is `B_H_TOT` = 2.576 in tall, i.e. 1030 px at 400 dpi, and starts
146 px down, under the two-line title — **derive that offset from the
cal measurement, don't assume it**: the body top sits at data
`yN_body_t`, so `y0 = row(body top) - (B_H_TOT - yN_body_t) * 400`.
Crop sub-regions at 600–800 dpi for the label / legend gutters.

Two checks that must pass before calling it done:

1. **Cached vs fresh is 0 differing pixels.** `compose.R b` then
   `compose.R`, rasterise both and diff — this is what catches an
   `aes()` that closed over a local scalar.
2. **The blank-rectangle scan.** See the whitespace budget above.
