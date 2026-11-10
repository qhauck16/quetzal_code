# figure3/panel_d/build.R

Figure 3 panel d — "GSP #13 phenotype captures consistent displacement
of canonical 3'ss downstream". A 1-nt-bin histogram of where the
cryptic 3'ss falls relative to the annotated one, with an event-type
pie inset in the histogram's empty left half.

## Purpose & context

The histogram is the claim (the cryptic site is displaced *downstream*,
mass piled at +5 to +10 nt); the pie is the qualifier (88% of the
differential junctions are alt 3'ss at all, so the histogram is about
almost all of them).

Both come from ONE file, `df_final.tsv`, for the nkap / 20v60 /
dpsi0.01_fc2.5 instance of the `diff_splicing/unified_results`
Snakemake pipeline. `stopifnot(nrow(df) == 2806)` is what pins the
instance — several parameter combinations exist under
`unified_results/nkap/` and only this one reproduces the published
percentages.

The histogram uses only `type == "alt 3'ss"` rows; **the pie counts
every row**. That asymmetry is deliberate and is what makes 88.06%.

## Decision rationale

### Percentages are ON the wedges. This is a standing instruction.

Quinn, verbatim:

> "pie charts in all figures should always have the numbers on the
> chart - you can add a rectangle around the ones that are too big for
> their wedge - fix that preference to never make that mistake again"

So the legend carries CATEGORY NAMES ONLY, and there are no leader
lines fanning labels out into white space either — the older
`pie_with_leaders()` did that and it was rejected. Every percentage is
drawn on its own wedge's mid-ray; the four that cannot fit inside
their wedge get a white box behind the text. The geometry lives in
`../../scripts/pie_chart.R` (`pie_labels_on_wedges()`), shared with
figure 4 panel d — read that file's header before changing anything
about the labelling.

This panel is the WORST CASE for that helper and is the reason it is
written the way it is: four of five wedges are under 6% (4.95, 0.32,
5.06, 1.60) and all four are adjacent, packed into a 43-degree sector
between 12 and 2 o'clock. Their four labels are separated purely by
RADIUS along their own rays, so they step outward as a fan.

### Why both titles are wrapped the way they are

Row 3 of figure 3 used to be `d | e` and panel d's cell was 4.42 in.
It is now `d | e | f` at 0.35 / 0.38 / 0.27, so the cell is **2.98 in**
and both of the old strings overflowed into panel e:

* the plot title was 2 lines whose longer line drew 2.90 in against
  2.836 in of usable plot width -> now 3 lines, longest 1.94 in;
* the x axis title was one 2.71 in line. Worse than it looks: an axis
  title is centred on the **panel**, not on the plot, and the panel
  starts ~0.63 in in (y-axis text + title), so only ~2.20 in is
  actually usable before the right edge clips. -> now 2 lines.

**Measure these with `cairo_pdf`, not `pdf()`.** The figures are
written with `cairo_pdf`, which resolves "Helvetica" through
fontconfig to a substitute that is 12-18% WIDER on this machine than
the base `pdf()` device's built-in Helvetica AFM metrics. Every "it
fits" computed with `pdf()` metrics is wrong by about a line's worth.

### Why the type in the legend shrank but nothing else did

Quinn's rule is stretch the plot, not shrink the type — but the cell
width is fixed by `compose.R`. The legend is the one element that is
pure furniture, so it went to 6 pt / 2.6 mm keys; at 6.4 pt / 3.0 mm
"skipped exon" pushed the legend block off the panel's right edge. The
axis text, axis titles and title all stay at the figure's tokens.

## Gotchas & warnings

- **`D_PIE_R_IN` is MEASURED, never computed.** It is the pie's drawn
  radius in inches, and it is what converts the label's point size
  into pie radii — i.e. it decides which labels fit inside and how far
  the rest must travel. It cannot be derived from `D_PIE_BOX`:
  `coord_equal` means the drawn size depends on the panel's height,
  which patchwork only settles at layout time. An early attempt to
  derive it from the box width gave 0.73 in against a real 0.31 in.
- **Measuring it is a FIXED POINT, and the map is not smooth.** A
  larger constant makes the label boxes smaller in radii, which
  shrinks the helper's computed limits, which makes the drawn pie
  bigger. Worse, the helper's placement search flips to a different
  discrete arrangement as the constant moves, so the measured value
  wobbles: 0.47 -> 0.483, 0.48 -> 0.474, **0.485 -> 0.486**, 0.49 ->
  0.476, 0.52 -> 0.494. Take the value that REPRODUCES ITSELF. Method:
  set it, `Rscript compose.R d`, crop the cell at 400 dpi, find the
  extreme pixels matching the five Set2 fills, radius = width/2/400.
- **Do not run a sweep of `D_PIE_R_IN` in the background.** It
  `sed`s the file in place; if it is still running when you set the
  final value it will overwrite it with the last value of its list.
- **`D_PIE_BOX` right edge is 0.50 and that is near the limit.** The
  pie's right rim lands at x ~ 649 px (400 dpi, cell-relative) and the
  first histogram bar's column starts at ~640 — they do not collide
  only because the bar at pos = 5 is short enough that the circle has
  already curved away. Widening the box further will put the pie on
  the bars.
- **Row 3's three titles do not line up** (400 dpi, from the row's top
  edge: d 69 px, e 6 px, f 151 px). d's 9.6 pt `plot.title` top margin
  was measured back when row 3 was `d | e`. This is now a three-panel
  argument and e/f are owned elsewhere, so d was left where it is —
  it sits between the other two. `../scripts/title_align.R` is the
  measuring tool.
- **`panel_cached()` serialises with `saveRDS`.** An `aes()` that
  closes over a local scalar renders BLANK after `readRDS`. Everything
  the pie draws is carried as a data-frame COLUMN for that reason.
  Verified: fresh build vs cached re-render is 0 differing pixels at
  400 dpi over the whole cell.
- **The published title says "GSP #20".** Stale Illustrator leftover
  from an older factor key; NKAP is k13 in current code and panel a of
  this figure says #13. Quinn confirmed #13.

## User preferences

- Percentages on the pie, never in a legend — standing, all figures.
- Stretch the plot rather than shrink the type (which is why the
  title took a third line instead of dropping a point size).
- Look at the rendered PDF at 400+ dpi; never reason about how it
  might look.

## Input/output

```
unified_results/nkap/20v60/dpsi0.01_fc2.5/df_final.tsv (2806 rows)
        │
        ├── type == "alt 3'ss", pos = -distance_to_cryptic, clipped to
        │   D_ZOOM = c(-8, 20)  ──► geom_histogram (1 nt bins)
        │
        └── count(type) over ALL rows ──► pie_labels_on_wedges()  ──┐
                                     └──► pie_legend_grob()  ───────┤
                                                                    ▼
                                                    inset_element into the
                                                    histogram's left half
                                                            │
                                            panel_d ──► ../compose.R ──► out/figure3.pdf
```

## Running it

```
cd quetzal_code/figures/figure3
/scratch/midway2/qhauck/conda_env/rstudio-server/bin/Rscript --vanilla compose.R d
cp out/figure3.pdf /tmp/mine3.pdf     # out/ is shared; crop a COPY
pdftoppm -r 400 -png -singlefile -x 0 -y 2480 -W 1190 -H 1250 /tmp/mine3.pdf /tmp/d3
```

`-W 1190` is the cell's exact width (0.35 * 8.5 in * 400 dpi); crop to
exactly that and any overflow into panel e shows up as ink touching
the right edge.
