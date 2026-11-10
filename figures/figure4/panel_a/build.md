# figure4/panel_a/build.R

Figure 4 panel a — "YLPM1 LoF mutants are the top 13 samples on
GSP #30". A full CDF on the left with the zoom region boxed, and on
the right that zoom drawn as a dot plot where every point is a TCGA
sample carrying its protein-change label.

**This is the same idiom as figure 3 panel a and shares its code.**
The block between the two `# ====` rules is a byte-identical copy of
`../../figure3/panel_a/build.R`'s. Only those two `build.R` files were
in scope when the solver was written, so it could not live in
`../scripts/`; if you change one copy, change the other. The long
rationale for the solver — why ggrepel was removed, what the cost
function is, what the relaxation ladder does — lives in
`../../figure3/panel_a/build.md`. This file records only what is
specific to figure 4.

## Purpose & context

Binary claim, not a four-gene one: the panel splits samples on the
table's own `frameshift` flag (115 TRUE of 9800) and shows that the
LoF samples occupy the top of GSP #30. The table also carries a
three-level `plot_color`, which is NOT what the published figure uses
— its legend has two entries.

Twelve labels, not 47, and a wider panel (2.058 x 1.753 in drawn),
so this one is never short of room. Its defects were of a different
kind.

## Decision rationale

### What was wrong

* **ggrepel could not see the markers or the staircases** (it only
  avoids its own layer's data). `D2015Vfs*40` printed on top of its
  own point, the red LoF step line ran through `M1180Cfs*27` and
  `P177Rfs*42`, `W455*` sat on the line, and the three labels lifted
  into the headroom threw long leaders across each other and across
  the navy cloud.
* **Half the panel was empty.** `A_ZLIM` topped out at 112 for labels
  that only ever needed four percentage points of headroom, so the
  upper ~30% of the zoom was a dead band, and the left column had a
  0.4 in blank band between the CDF's axis title and the first legend
  key.

### What was done

* Same `cdf_place_labels()` as figure 3 panel a. With twelve labels it
  finds every one a slot immediately adjacent to its own point:
  **0 leaders, 0 grazes, 0 relaxations.** That is the point of using
  the same machinery — the solver spends nothing it does not have to.
* `A_ZLIM` top is now 101.6 — one clear text row above the 100% row,
  where the unlabelled navy cloud's step line runs the full panel
  width, and nothing more. Bottom is `min(cumulative) - 1.2`.
* Legend cell `heights = c(1, 0.10)` instead of `0.22`.
* x breaks pinned to `seq(-4, 0, 1)` so the label gutters do not get
  ticks.

### Why alternating flanks

The twelve LoF points are ~0.87 percentage points apart, which at this
panel height is a hair more than one text row — so a single flank
would ALMOST work and then fail on the two long strings
(`P1507Afs*6; R1719*; Y1506*` is 70 pt, half the panel's width).
Alternating gives every label its own row and reads as a ladder. It is
a soft preference (`w_flip`); the solver overrides it wherever the
step line or a neighbour is in the way.

## Gotchas & warnings

* **`A_ZOOM_W_IN` / `A_ZOOM_H_IN` are the DRAWN PANEL size in inches**
  (measured at 400 dpi: 823 x 701 px = 2.058 x 1.753 in), and every
  pt→data conversion in the solver hangs off them. Re-measure after
  any change to `ROW_H[1]`, `COL_W$ab`, the internal
  `plot_layout(widths = )` or the axis furniture: render, `pdftoppm
  -r 400`, find the horizontal run of the zoom's x-axis line and the
  vertical run of its left border. Wrong numbers mean the solver
  "succeeds" while the drawn panel overlaps.
* Every trap listed in `../../figure3/panel_a/build.md` applies
  verbatim: cairo text metrics, font-height collision boxes, the
  `pad_x` split across the text anchor, `lead_pad` as a penalty rather
  than a constraint, and `panel_cached()`'s `saveRDS` requiring every
  aes() value to be a DATA COLUMN.
* **The CDF panel keeps a 16 pt left margin.** patchwork puts the "a"
  tag at the cell's top-left and it otherwise prints through the
  y-axis title.
* **The zoom's y range is NOT the CDF's.** Both classes are already
  past 90% at log = -4, so reusing `c(0, 108)` squeezed every point
  into the top fifth. The comment about "shared so both halves'
  gridlines align" was inherited and was never true; the zoom has no
  y-axis text, so nothing is lost.
* **Panel a's title is set via `plot_annotation()`, not `labs()`**;
  patchwork stacks that title on its own plot margin, hence
  `plot.margin = margin(0,0,0,0)`. See `../../figure2/compose.md`.
* This panel builds in ~5 s (twelve labels; the polish loop exits
  early). `Rscript compose.R a` rebuilds only this panel.

## User preferences (Quinn)

* **Stretch the plot, don't shrink the type** — label size is
  `SIZE_BASE * 0.62`, identical to figure 3 panel a.
* **Figure 3 panel a and figure 4 panel a must look like siblings.**
  Same label size, same leader style (0.12 pt grey55, drawn only when
  the label is more than `lead_min_pt` from its marker), same legend
  treatment (`ncol = 1`, bottom, key 2.3 mm, `0.78 * SIZE_BASE` — the
  small size is forced by figure 3's long "No PRPF8/NKAP/FAM50A/SLU7
  variant" entry and is matched here on purpose), same solver, same
  gutter logic, same x-break rule.
* **No large dead whitespace bands**, verified by measurement on the
  render rather than by eye. After the `A_ZLIM` change the zoom
  panel's largest interior blank band is 0.14 in (the bottom axis
  padding) and the left column's is 0.24 in (CDF axis title → legend).
* Audit figures by measurement on the render; do not reason about how
  the PDF might look.

## Input/output

Reads one TSV and nothing else:

* `$FLASHIER/YLPM1_full_plus_variants.tsv` — one row per TCGA sample,
  with `k30`/`logk30` (the GSP #30 loading), a logical `frameshift`
  and a free-text `proteinChange`. Nothing is parsed here; the class
  is `ifelse(frameshift, …)` and the label is `proteinChange` as
  given, with `NA`/`"NA"`/`""` treated as unlabelled.

Leaves `panel_a` (a `wrap_elements`) for `../compose.R`, which caches
it to `panels/a.rds` and writes `out/figure4.pdf`.

## File flow schematic

```
YLPM1_full_plus_variants.tsv
            │  (frameshift flag; proteinChange as-is)
            ▼
   a_cum (logk30, cumulative, variant_class, proteinChange)
            │
   ┌────────┴─────────────────────┐
   ▼                              ▼
panel_a_cdf                a_zoom (logk30 >= -4)
(full CDF +                       │
 zoom rectangle)   ┌──────────────┼──────────────┐
   │               ▼              ▼              ▼
   │       cdf_step_segments   markers    a_lab (12 LoF labels)
   │               └───────┬──────┘              │
   │                       ▼                     │
   │               cdf_place_labels() ◄──────────┘
   │                       │
   │             lx, ly, hj, sx..ey  (all adjacent: no leaders)
   │                       ▼
   │                 panel_a_zoom
   └──────────┬────────────┘
              ▼
     panel_a ──► ../compose.R ──► panels/a.rds
                              └─► out/figure4.pdf
```
