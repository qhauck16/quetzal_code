# figure3/panel_a/build.R

Figure 3 panel a — "Somatic mutations in NKAP, PRPF8, FAM50A and SLU7
are enriched upon GSP #13". A full CDF on the left with the zoom
region boxed, and on the right that zoom drawn as a dot plot where
every point is a TCGA sample carrying its protein-change label.

**Figure 4 panel a is the same idiom and shares this code verbatim.**
`../../figure4/panel_a/build.R` holds a byte-identical copy of the
block between the two `# ====` rules. Only these two `build.R` files
were in scope when the solver was written, so the code could not live
in `../scripts/`; if you change one, change the other.

## Purpose & context

Panel a is the claim that GSP #13 is a mutation phenotype: four genes'
variants all sit at the top of the loading. The left CDF carries the
statistical claim; the right zoom carries the identity of each sample.
47 labels in a panel 131.6 x 152.5 pt is the whole difficulty — the
panel exists to be read one mutation at a time.

## Decision rationale

### Why ggrepel was removed

ggrepel placed these labels until 2026-09. It cannot work here, and
the reason is structural, not a matter of tuning:

* **ggrepel only avoids its own layer's data.** It cannot see the
  markers, the four `geom_step` staircases or the *other* repel layer.
  It drew the red PRPF8 step line straight through "N1881insKN" and
  knotted six strings (H1875_L1876delinsQ / H1875_L1876del / R289Q /
  X38_splice / E192del / N1881insKN / Y73H / L67P / N1881H,V2883dup)
  at cumulative 88-93%, where the PRPF8, FAM50A and SLU7 classes all
  converge.
* **The band is at capacity.** Cumulative 85-100% is 15% of the y
  range — five text rows — and the labels that belong in it are 560 pt
  of Helvetica at 4.96 pt, i.e. 88% linear occupancy. No force
  parameter solves an 88%-full band; the labels have to leave it.

Things that were tried and helped but were not enough: giving the zoom
its own y range, pinning classes to flanks, dropping `force_pull`,
growing `ROW_H[1]` from 2.55 to 2.95 in. One thing was tried and made
it worse and must not be redone: a deterministic per-flank pre-spread
that assigned each label a y by RANK with a minimum gap. It places
labels by rank rather than near their points, so leaders grow long,
cross into the other class's column, and the NKAP diagonal is
scrambled.

### What replaced it: `cdf_place_labels()`

An explicit 2-D solver. For each label it searches a cost-ordered grid
of candidate slots — slide out in x on either flank, step in y, or sit
centred above/below the point — and takes the cheapest one that is
geometrically clear of *everything*:

* the panel edges,
* every marker in the panel (all 71, not just the labelled ones),
* every segment of all five staircases (`cdf_step_segments()` turns
  `geom_step(direction = "hv")` into explicit axis-aligned segments),
* every label already placed,
* every leader already drawn (a leader may not run through a string,
  and a string may not land on a leader).

Two things are penalties rather than constraints: a leader crossing
another leader (`w_cross`) and a leader passing within `lead_pad` of a
string without entering it (`w_graze`). Both were tried as hard
constraints and both made the result *worse* — they take slots away
from the greedy so early that later labels need the level-2 escape
hatch, and a leader clean through a string is far uglier than either.

The panel has genuinely more label than room in one band, so the
solver's job is to spend the slack that exists elsewhere: the >100%
headroom and the half-empty 68-85% band. That is exactly the space
ggrepel could not reach, because reaching it means moving a label 20 pt
from its point and drawing a leader.

### Greedy, re-ordered, then polished

Greedy placement is order-dependent and "hardest first" is not always
right — a label that squeezes into a near slot can push a neighbour
thirty points away. So:

1. order by local text demand (summed widths of labels within a few
   text rows in y), longest string first;
2. re-run, promoting anything that failed to the front, until nothing
   fails (`rounds = 26`);
3. then shuffle the order `polish = 80` more times from a seeded RNG
   and **keep the cheapest layout seen**, cost being each label's
   displacement from its own point plus a quadratic far term.

Step 3 is what turns "all 47 placed" into "all 47 placed next to their
own points". Without the quadratic term the solver parks one label
50 pt away to save two labels 25 pt each — the same linear cost, one
leader flying across the panel instead of two short ones.

### Relaxation ladder

If a label still has no slot after all that, it gives up constraints
one step at a time (`relax[bad] <- relax[bad] + 1`, never jumping a
level): level 1 lets a leader cross its string, level 2 lets a step
line cross it. Level 2 is the original defect and exists only so that
a label can never be silently dropped. **Current state: level 0, all
47 placed, zero relaxations.** If the solver ever prints a non-zero
"unplaced" or a collision warning, the panel has a visible defect.

### Side preference

Only NKAP asks for one (alternate left/right). Its 20 points are the
long lower diagonal with room on both sides, and alternating is what
keeps it reading as a *diagonal* rather than a column. PRPF8, FAM50A
and SLU7 are free (`pref = 0`). Pinning them to a flank — which the
ggrepel version had to do — is actively harmful with a real solver:
PRPF8's stack sits ~10 pt from the panel's left edge, so "left" is a
flank it does not have, and forcing it there made the solver buy the
flank with 30 pt of vertical travel.

## Gotchas & warnings

* **`A_ZOOM_W_IN` / `A_ZOOM_H_IN` are the DRAWN PANEL size in inches,
  and every pt→data conversion hangs off them.** Measured at 400 dpi:
  731 x 848 px = 1.8275 x 2.12 in. Re-measure after ANY change to
  `ROW_H[1]`, `COL_W$ab`, the `plot_layout(widths = )` inside the
  panel, or the axis furniture. The recipe: render, `pdftoppm -r 400`,
  find the horizontal run of the zoom's x-axis line and the vertical
  run of its left border. Get these wrong and the solver "succeeds"
  while the drawn panel overlaps.
* **Text widths must come from a `cairo_pdf` device.** The figure is
  written with `cairo_pdf`, which resolves "Helvetica" through
  fontconfig; `grDevices::pdf()`'s built-in AFM metrics are ~6%
  narrower, which is enough to let two labels print flush.
  `cdf_text_metrics()` opens a throw-away cairo device and restores
  the previous one.
* **The collision box must be as tall as the FONT, not as the measured
  ink.** `geom_text(vjust = 0.5)` centres each string on its own ink
  box, so "R366Q" (caps only) and "N1881_V1883dup" (underscore below
  the baseline) sit differently inside the same nominal row. Boxes
  sized from cairo's `grobHeight("Rgp")` printed the upper string's
  underscore as a bar through the lower string's digits. `th` is
  therefore `max(measured, lab_pt * 0.92)`.
* **`pad_x` has to be split across the anchor.** `lx` is the TEXT
  anchor and sits `pad_x` INSIDE the box edge (`geom_k`). Anchoring
  the text at the box edge puts all of `pad_x` on the far side and
  none between the string and its own marker — that is what printed
  "Y1506*" flush against its point in figure 4.
* **Leaders are tested against a box inflated by `lead_pad`, and that
  test is a PENALTY, not a constraint.** A hairline running 0.9 pt
  past the final "I" of "R290I" printed as "R290Il". Making it a hard
  constraint cost more than it bought (see above).
* **`cdf_xlim()` fixes the gutters up front.** Solving needs `xlim`,
  and `xlim` would need to know how far labels reach — the circle is
  cut by reserving 40 pt on the left and 16 pt on the right and
  letting long strings overhang the point cloud at a y the solver has
  proved is clear. 40 pt on the left is load-bearing: a sweep over
  (`pad_x`, `mr_extra`, gutter) showed gutters of 34 pt failing at
  several clearance settings and 40/46 succeeding at all of them,
  because PRPF8's whole lower staircase lives ~10 pt from the left
  edge and its labels are 17-40 pt wide.
* **The zoom's x breaks are pinned to the data** (`seq(A_ZOOM_FROM, 0,
  1)`). The gutters widen the x limits past -4, and a tick at -6 in a
  region the zoom does not contain contradicts the box drawn on the
  CDF at -4.
* **`panel_cached()` serialises this panel with `saveRDS`.** An
  `aes()` that closes over a local scalar renders BLANK after
  `readRDS`. Everything the layers need (`lx`, `ly`, `hj`, `sx`…`ey`,
  `text`) is a DATA COLUMN for that reason. Layer *parameters*
  (`size`, `family`, `linewidth`) are captured by value and are safe.
* **The RNG is borrowed and returned.** The polish loop saves
  `.Random.seed`, sets its own (20260923) and restores, so the panel
  is byte-reproducible and nothing else in `compose.R` moves.
* **Panel a's title is set via `plot_annotation()`, not `labs()`**, and
  patchwork spaces the two differently. The 2 pt top margin and the
  zero `plot.margin` are measured corrections against panel b's title
  at 400 dpi; see `../../figure2/compose.md`.
* **The legend is 0.78 x `SIZE_BASE`**, not 0.9. At 0.9 the entry "No
  PRPF8/NKAP/FAM50A/SLU7 variant" is 1.77 in wide and runs under the
  zoom panel's left axis. Figure 4 panel a uses the same size for
  consistency even though it has the room.
* Panel a takes ~90 s to build (26 ordering rounds + 80 polish rounds
  x 47 labels). `Rscript compose.R a` rebuilds only this panel.

## User preferences (Quinn)

* **Stretch the plot, don't shrink the type.** The label size is
  unchanged at `SIZE_BASE * 0.62`; every gain came from packing and
  from `ROW_H[1]` / the panel's internal width split, never from
  smaller text.
* **No label may sit on its own point marker or on a step line; no
  leader may cross another label.** These are hard constraints in the
  solver, and the audit block re-checks the layout that is actually
  drawn rather than trusting the search.
* **No large dead whitespace bands.** Verified by measurement, not by
  eye: `pdftoppm -r 400` then scan for all-white rows/columns inside
  the panel rectangle. The zoom panel has no blank band over 0.1 in in
  either direction; the left column's largest is 0.22 in, between the
  CDF's axis title and the legend.
* **Figure 3 panel a and figure 4 panel a must look like siblings** —
  same label size, same leader style (0.12 pt, grey55, drawn only when
  the label is more than `lead_min_pt` from its marker), same legend
  treatment, same solver.
* Audit figures by measurement on the render; do not reason about how
  the PDF might look.

## Input/output

Reads one TSV and nothing else:

* `$FLASHIER/NKAP_full_plus_variants.tsv` — one row per TCGA sample,
  with `k13`/`logk13` (the GSP #13 loading) and three free-text
  mutation columns. `scripts/nkap_variants.R` parses those into
  `variant_class` (5 levels, priority NKAP > PRPF8 > FAM50A > SLU7 >
  Other) and `proteinChange`; `VARIANT_COL` is the published palette.

Leaves `panel_a` (a `wrap_elements`) for `../compose.R`, which caches
it to `panels/a.rds` and writes `out/figure3.pdf`.

## File flow schematic

```
NKAP_full_plus_variants.tsv ──► classify_variants()  [scripts/nkap_variants.R]
                                        │
                                        ▼
                            a_cum  (logk13, cumulative, variant_class,
                                    proteinChange)
                                        │
                     ┌──────────────────┴──────────────────┐
                     ▼                                     ▼
            panel_a_cdf (full CDF                  a_zoom (logk13 >= -4)
            + zoom rectangle)                              │
                     │                    ┌────────────────┼────────────────┐
                     │                    ▼                ▼                ▼
                     │            cdf_step_segments   markers      a_lab (47 labels)
                     │                    └────────┬───────┘                │
                     │                             ▼                        │
                     │                     cdf_place_labels() ◄─────────────┘
                     │                             │
                     │                   lx, ly, hj, sx..ey
                     │                             ▼
                     │                      panel_a_zoom
                     └──────────┬──────────────────┘
                                ▼
                     panel_a ──► ../compose.R ──► panels/a.rds
                                              └─► out/figure3.pdf
```
