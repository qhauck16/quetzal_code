# figure4/panel_d/build.R

Figure 4 panel d — "LeafCutter results show upregulation of both alt
3'ss and 5'ss". One pie of splice-event types over every differential
junction, percentages ON the wedges, category names in the legend.

## Purpose & context

The panel exists to make one point: for YLPM1 / GSP #30 the effect is
NOT confined to one splice-site class — alt 5'ss (50.1%) and alt 3'ss
(30.3%) both dominate. Compare figure 3 panel d, where NKAP is 88%
alt 3'ss. The two pies are deliberately the same idiom so the contrast
reads at a glance; they share `../../scripts/pie_chart.R`.

Instance: ylpm1 / 20v60 / dpsi0.01_fc2.5, pinned by
`stopifnot(nrow(df) == 2460)` (the published "N = 2460") and by its
type breakdown reproducing the published percentages to the decimal.
NB Quinn's recollection was dPSI > 0.05; the data says 0.01, the same
threshold family as figure 3 panel d.

## Decision rationale

### Percentages are ON the wedges. This is a standing instruction.

Quinn, verbatim:

> "pie charts in all figures should always have the numbers on the
> chart - you can add a rectangle around the ones that are too big for
> their wedge - fix that preference to never make that mistake again"

The legend carries CATEGORY NAMES ONLY. There are also no leader
lines: the previous `pie_with_leaders()` fanned narrow wedges' numbers
out on leaders and that is exactly what was rejected. Now every label
sits on its own wedge's mid-ray and the ones too big for their wedge
get a white box behind the text. See `../../scripts/pie_chart.R`.

What changed visually here: **12.8% moved INSIDE its wedge** (it has
23 degrees of half-angle, which is enough to hold the string once the
fit test is done properly instead of against a `prop >= 0.15`
threshold), and 5.7% / 1.1% became boxed chips straddling the rim on
their own rays instead of leadered callouts.

### Why the title lives on a separate void plot

A title on the pie itself gets pushed around by `coord_equal`, which
reserves whatever space the square needs and then centres the title
over the result. So the pie is an `inset_element` on an empty
`frame` plot that carries the title. Two traps in that frame:

* an EMPTY ggplot has no scale range, so `annotate()` on it centres
  itself instead of honouring its coordinates — `coord_cartesian(xlim
  = c(0,1), ylim = c(0,1))` pins a frame to place "N = 2460" against;
* `labs(x = NULL, y = NULL)` is required or `theme_publication()`
  draws stray "x" / "y" axis titles.

## Gotchas & warnings

- **`D_PIE_R_IN` is MEASURED, never computed.** It is the pie's drawn
  radius in inches, and it converts the label's point size into pie
  radii — which is what decides whether a label fits inside its wedge
  and how far a boxed one has to travel. It cannot be derived from the
  inset box: `coord_equal` makes the pie height-limited and the drawn
  height is only settled by patchwork at layout time (deriving it from
  the box width in figure 3 once gave a 2x error).
- **Measuring it is a FIXED POINT.** A larger constant shrinks the
  label boxes in radii, which shrinks the helper's computed limits,
  which grows the drawn pie. Re-render and re-measure until the number
  reproduces itself. Here it took two passes, 0.685 -> 0.646 ->
  **0.645**, and 0.645 reproduces itself exactly. Method: render, crop
  the cell at 400 dpi, take the extreme pixels matching the five Set2
  fills, radius = width / 2 / 400.
- **The limits are computed from where the labels land**, so changing
  the label text, the point size or the palette moves the pie's drawn
  size and stales `D_PIE_R_IN`. Same for `ROW_H` or the inset box.
- **`panel_cached()` serialises with `saveRDS`.** An `aes()` closing
  over a local scalar renders BLANK after `readRDS`. Every value the
  pie layers need is carried as a data-frame COLUMN. Verified: fresh
  build vs cached re-render is 0 differing pixels at 400 dpi.
- The cell is 3.95 x 2.00 in (row 2 is `c | d` at 1.15 / 1). Nothing
  here was re-fitted for width — unlike figure 3 panel d, this cell
  did not change.

## User preferences

- Percentages on the pie, never in a legend — standing, all figures.
- Look at the rendered PDF at 400+ dpi; never reason about how it
  might look.

## Input/output

```
unified_results/ylpm1/20v60/dpsi0.01_fc2.5/df_final.tsv (2460 rows)
        │
        └── count(type) ──► pie_labels_on_wedges()  ──┐
                       └──► pie_legend_grob()  ───────┤
                                                      ▼
                                   frame (title + "N = ") + 2 insets
                                                      │
                                    panel_d ──► ../compose.R ──► out/figure4.pdf
```

## Running it

```
cd quetzal_code/figures/figure4
Rscript --vanilla compose.R d
cp out/figure4.pdf /tmp/mine4.pdf     # out/ is shared; crop a COPY
pdftoppm -r 400 -png -singlefile -x 1700 -y 1000 -W 1700 -H 830 /tmp/mine4.pdf /tmp/d4
```
