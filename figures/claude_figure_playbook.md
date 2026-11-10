# Publication-figure playbook for Claude (R / ggplot2 / patchwork)

Lessons from building a full paper's main and supplementary figures in code.
Give this to a Claude session at the start of any figure work. It is data-agnostic.
It covers how to build, check and style figures so they come out clean,
centred and consistent. Section 1 is the owner's taste; everything else is
craft that applies to anyone.

---

## 1. Style preferences (defaults, apply without asking)

**Typography**
- One typeface for the whole page. Name a font the PDF device can actually
  find. See §4: "Helvetica" silently became DejaVuSans for weeks.
- Size all text from a few theme tokens (base ~8 pt, axis titles 9 pt,
  panel titles, legend, strip) and never hand-matched numbers. Two labels
  that should match are then identical by construction.
- Never go below ~6 pt. If labels collide, **stretch the plot or page, don't
  shrink the type** (§3).
- Gene and protein symbols are *italic* everywhere: titles, axis labels,
  legends, annotations.
- Panel letters are lowercase **bold** (a, b, c) at top-left. Legends and keys
  get no letter. Only letter panels where the published layout has them. A
  single-panel figure has none.
- Plain weight for track, cohort and strip labels; bold only for titles and
  panel letters.
- Real superscripts for p-values: *p* = 2.7 × 10⁻⁸³ via plotmath, not `2.7e-83`.
  Where a CDF shows exact p-values, don't put a capped "p < 2.22e-16" next to
  them. That is R's display floor, not the value; compute and print the real p.

**Axes and labels**
- Categorical x tick labels at **45°**, never vertical.
- Label the object, not the panel edge. Put a transcript's name on the
  transcript and a mutation's name at its point.
- Scatter and CDF panels are square (`aspect.ratio = 1`).
- Main figures: no gridlines. Supplements: light major gridlines (grey90,
  0.25 lw) read better. Dense small multiples look odd without them.
- Don't clutter with explanatory captions inside the figure. The figure
  caption does that job.

**Colour**
- Use one palette family consistently (tab20 pastels worked well). Take a
  true medium blue over near-navy.
- An entity that appears in several panels gets **one colour token** used
  everywhere, so the uses can't drift.
- Keep the same meaning → same colour across main and supplementary figures.
  A supplement that says "same as Fig 3e" must reuse 3e's colours.
- Diverging heatmaps: RdBu, symmetric limits, `oob = squish`.

**Layout**
- Minimise whitespace aggressively. Dead bands and empty right-hand gutters
  get noticed.
- Centre content rows over full-width elements. Use fixed-width panels (in
  inches) between two flexible `plot_spacer()`s, rather than letting panels
  float left with a void on the right.
- Equivalent elements in different panels (two violins, two pies) are the
  same drawn size.
- Align the tops of panels' *content*, not just their bounding boxes.
- Squat, squarish violins beat tall slivers stretched to the row height.

**Pie charts**
- Percentages go **on the pie** and never in the legend. The legend holds
  category names only.
- Large wedges: plain text inside the wedge.
- Small wedges: a rounded chip **filled with that wedge's own colour**, no
  border and no leader line. Place it on the disc as close as possible to its
  own wedge without overlapping another label, even if it lies across a
  neighbouring wedge. Never fly it off past the rim.

**Genome tracks / sashimi**
- Junction arcs are thin curved lines, not filled shapes. Thickness ∝
  magnitude (factor weight, or √PSI).
- PSI labels sit **inline, mid-arc, on their own arc**, white-backed with no
  border, and off other arcs and coverage bars. Use the same convention in
  every figure.
- Coverage tracks: numeric y-axis with an axis line and just two ticks (0 and
  max), 0 flush with the baseline. Normalise to max coverage and share the
  scale within a comparison.
- Arcs hug the gene body. Use raw genomic coordinates; never shift x to make
  something fit.
- Explain every glyph: add a small key for asterisks and stars (non-coding
  junction, start/stop codon).
- In a grid of repeated panels, put row labels only on the leftmost column.
  Make them transparent elsewhere so every panel keeps the same drawn width.
- Heat-strip plots: no baseline under every row. A constant line reads as
  signal and flattens the apparent enrichment. Break the fill where it is
  zero.

**Data honesty**
- Never let a display filter drop the feature the panel is about. A
  `min_psi` cutoff once hid the exact junction the panel existed to show.
- Keep precision where it matters, e.g. one decimal for sub-10 % values.
- Label exclusions deliberately. If a sample sits high because of a
  *different* driver mutation, don't label it as an example of this one.

---

## 2. Process: how to work

1. **Find the original code and data first.** For figures that already
   exist, e.g. as Illustrator PDFs, locate the Rmd or script that made each
   panel. Reproduce its statistics verbatim, then check your numbers against
   the published ones (p-values, counts, n). Record every deviation and its
   cause.
2. **Split heavy data work from drawing.** Put a cache builder per figure in
   `cache/` that reduces the multi-GB inputs to small TSV/RDS files. Run
   those on the cluster. The figure `build.R` reads only caches and renders
   in seconds, which keeps the edit → render → inspect loop fast.
3. **One script per figure** (`figS<N>/build.R`) sourcing a shared
   `common.R` / `theme.R`. Keep a `build_all.sh` that runs everything,
   continues past failures and prints PASS/FAIL. Keep a README table mapping
   figure → build script → cache → original source.
4. **Document each script** with a companion `script.md` for future sessions:
   purpose, decision rationale, gotchas, owner preferences, inputs/outputs,
   and an ASCII file-flow diagram. Append dated notes when review changes the
   design.
5. **Audit the render yourself, by measurement.** Never hand back a figure
   for the human to list defects.
   - Rasterise at 200–300 dpi (`pdftoppm -r 250`). Crop panels at high zoom
     (`-x -y -W -H`); full-page views hide collisions.
   - Look for text over text, text over lines or points, clipped labels,
     misaligned panels, blank bands, tiny or inconsistent fonts, missing
     italics, bold that isn't bold, and ggplot leftovers ("group" axis
     titles, default legends).
   - Measure with ink masks: per-panel ink extents for alignment, all-blank
     row bands for dead space. Iterate numerically instead of nudging by feel.
   - Read ggplot warnings. "Removed N rows" means something vanished.
   - Check embedded fonts with `pdffonts out/*.pdf`.
6. **Compare against the published original** side by side for content and
   panel letters. Aim to improve the style, not change the content.
7. **Parallelise review with subagents only when authorised.** One agent per
   group of about six figures works well. Keep shared files (theme, helpers)
   off-limits to them, and have them report shared-code issues instead of
   editing. Then spot-check their output yourself.

---

## 3. Layout craft (ggplot2 + patchwork)

- **Derive offsets, never guess them.** Anchor text by the edge that faces
  the thing it clears (`vjust = 0` above an object, then `y = top + gap`).
  Convert text size to data units from the panel's drawn size:
  `clearance = pt * k * y_range / (height_in * 72)`. Constants tuned on one
  render break silently when the panel changes size.
- **Measure text on the device you'll save with.** `cairo_pdf` metrics differ
  from `pdf()`. Wrap string-width measurements in a helper that opens a
  scratch `cairo_pdf`. Text has width, so collision tests must sample the
  whole label footprint, not its centre.
- **Absolute sizes.** Keep row heights in inches with
  `page_h = sum(row_h)`. With relative heights, growing one row rescales
  every other row, and anything that derives spacing from its own height
  goes stale.
- **Keep every top-level panel the same object type** (all
  `wrap_elements(full = ...)`, or all bare ggplots). Mixing them changes tag
  spacing and styling.
- **Panels sharing an axis need identical margins** on the shared side. Drive
  coupled margins from one constant.
- **`plot_layout(heights =)` needs a flat stack.** `(a / b) / c` has 3 rows.
- **Layers drive the panel range.** Off-panel annotation text stretches the
  panel and drags axis labels with it. Pin the view with
  `coord_cartesian(ylim = ..., clip = "off")`.
- **Scale limits delete data before stats.** `scale_x_continuous(limits =)`
  truncates the sample an ECDF is computed on. Use `coord_cartesian()` for
  views.
- **ggrepel only avoids its own layer.** Feed obstacle points into the repel
  layer with `label = ""`, or split labels into layers with disjoint ranges.
  For crowded clusters, a hand-placed stacked label column in empty space
  with non-crossing leaders beats repel. Seed it and set `max.time` high, or
  the layout changes with machine load.
- **Facet order across layers.** Every layer's facet variable must be a
  factor with the same levels, or facets get sorted alphabetically.
- **Don't rasterise native plots** to fit a cell; it distorts the text. Only
  genuinely external art (Illustrator schematics) should be raster.
- Use `patchwork::wrap_elements(x, ignore_tag = TRUE)` for legends and keys so
  they don't consume a panel letter.
- Randomised plotting helpers (e.g. `fastTopics::structure_plot`) need
  `set.seed()`, or every rebuild reshuffles and diffs become meaningless.

---

## 4. Fonts and PDF output (these fail silently)

- **`family = "Helvetica"` under `cairo_pdf` falls back to DejaVuSans** on
  Linux when fontconfig has no Helvetica alias. Use a named metric clone that
  exists, such as `"Nimbus Sans"` (URW). Verify with `pdffonts`.
- **plotmath operators and spaces** (the `=` in `italic(p) == ...`) use the
  device's *symbol* font, which is separate from `family`. Fix it on the
  device:
  ```r
  cairo_pdf_pub <- function(...) grDevices::cairo_pdf(...,
    symbolfamily = grDevices::cairoSymbolFont(FONT_FAMILY, usePUA = FALSE))
  ggsave(f, p, device = cairo_pdf_pub, bg = "white")
  ```
- **Panel letters won't be bold** from `plot_annotation(theme = theme(plot.tag
  = ...))` alone. Patchwork styles each tag from that panel's own theme, so
  wrapped and nested panels ignore it. Set `plot.tag` bold in your base theme
  and globally:
  ```r
  theme_update(text = element_text(family = FONT_FAMILY),
               plot.tag = element_text(size = 10, face = "bold", family = FONT_FAMILY))
  update_geom_defaults("text",  list(family = FONT_FAMILY))
  update_geom_defaults("label", list(family = FONT_FAMILY))
  ```
- **Build where the fonts are.** HPC compute nodes may lack the fonts (and
  tools like `pdftocairo`) that the login node has, so output silently comes
  out in a different face. Render light figure scripts on the login node;
  send only heavy data caches to the scheduler.
- ggtext markdown titles in cached or serialised ggplots need
  `library(ggtext)` attached at render time.
- `fread` can misdetect headers on all-numeric matrices, so pass
  `header = TRUE`.
- Masking: Bioconductor packages mask `dplyr::first/rename/select`. Namespace
  them.

---

## 5. Working with the human

- Take feedback literally and generalise it: "the percentages should be on
  the pie" applies to every pie in every figure, now and later.
- When a request is ambiguous, e.g. "remove the little text", do the
  narrower reading and say what you kept.
- Report deviations from the published version plainly, with their cause
  (data file edited later, rounding, a classification choice).
- Don't guess at renders; look. Don't claim done until you've checked the
  high-dpi output.
