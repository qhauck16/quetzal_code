# figure2/compose.R

## Purpose & context

Regenerates Figure 2 (SF3B1 / genome-wide GSPs) from raw data, as
`figure1/compose.R` does for Figure 1. The reference being reproduced is
`tcga/genome_wide/flashier/paper_figures/Figure2-SF3B1_improved.pdf`,
which was assembled in Illustrator from separately-saved PDFs.

Seven panels and where their code came from:

| Panel | Content | Source Rmd |
|---|---|---|
| a | genome-wide workflow schematic | bespoke Illustrator art |
| b | GSP #1 neuronal-junction ECDF + WDR11 sashimi inset | `flashier/neuronal_gsp.Rmd` |
| c | GSP x covariate correlation heatmap | `flashier/flash_fact_char.Rmd` |
| d | GSP #10 vs SF3B1 hotspot mutations | `flashier/flash_fact_char.Rmd` |
| e | cryptic 3'ss distance histogram + event-type pie | `flashier/spec_flash_fact.Rmd` |
| f | PSI fold-change vs binned GSP #10 | `diff_splicing/sf3b1/top500_alt3_fc_vs_loading.Rmd` |
| g | mutation-enrichment scatter | `flashier/flash_fact_char.Rmd` |

## Decision rationale

**Shared design tokens.** `../theme.R`, `../scripts/raster_utils.R` and
`../scripts/coverage_plot.R` are shared with figure 1 rather than copied.
That is the whole point — typography and palette cannot drift between
figures if there is one source. Figure 1 was updated to source them too;
re-verify figure 1 renders after touching any of them.

**Three panels are cached, four are not.** Panels b, c and g each need
the 3.5 G flashier RDS (`softimpute_flash_300_qc_filtered.RDS`) and, for
c/g, a 536 M mutation table. Those go through
`scripts/build_cache_panel_{b,c,g}.R`, which persist small TSVs into
`cache/`. `compose.R` then reads only the TSVs and rebuilds in seconds.
Panels d/e/f read small cached artifacts directly (834 K – 19 M) and are
fast enough to recompute every build.

**The caches MUST run under SLURM, not on a login node.** Every attempt
on midway3-login1 was SIGKILLed partway through (reproducibly, even with
nothing else running, at only ~7.4 G RSS — so it is the login-node
watchdog, not the OOM killer). See `scripts/submit_caches.sh`; note it
uses the `build` partition because `pi-yangili1`'s midway3 SU balance is
negative and the shared partitions reject the job.

**Things the reference did in Illustrator, redone properly in R:** the
panel e pie inset (was a pasted PDF), the panel d zoom-region rectangle,
panel c's Group legend, panel g's `GSP #N–GENE` label formatting, and
panel a's title (re-set in R so it matches b–g's title styling). Panel
b's sashimi inset is drawn natively from the bigwigs + `.links` with
`coverage_plot.R`'s `sashimi_track_gg` rather than insetting the
pyGenomeTracks render — vector text, consistent with figure 1.

**Panel g's label rule.** The reference labels five points. That set is
reproduced exactly by "top 5 by `-log10(linear_pval)`" — it is *not*
"passes Bonferroni on both tests" (IDH1 is labelled with a near-zero CMH
score; HRAS/RBM10/PRKDC clear both cutoffs but are unlabelled). See
`N_LABEL` in the panel g block.

## File layout

`compose.R` holds only the shared setup (theme sourcing, paths, `ROW_H`,
the panel cache, `missing_panel`, `sf3b1_info`) and the final
`plot_layout`. Each panel lives in its own `panel_<x>/build.R`, sourced
in a fixed order, and must leave a `panel_<x>` object defined. This is
so panels can be worked on independently — including by separate
agents — without contending for one file. Layout (row heights, column
widths, which panel goes where) is decided ONLY in `compose.R`.

## Per-panel build cache

`compose.R` serialises each built panel to `panels/<name>.rds` and
reuses it, so iterating on one panel's layout does not rebuild the
other six:

```
Rscript compose.R          # reuse every cached panel (~18 s)
Rscript compose.R c a      # rebuild only c and a
Rscript compose.R --all    # rebuild everything (~60 s)
```

**Sharp edge:** a panel's `aes()` must not reference variables from the
enclosing environment. Tidy evaluation resolves them from the calling
env, which is gone after `saveRDS`/`readRDS` — the panel then fails at
render with `Problem while computing aesthetics` and the whole figure
comes out blank. Panel c's leader-line layer hit this by using bare
`x_edge`/`x_lab` scalars; the fix is to carry them as COLUMNS of the
layer's data frame. If you add a layer and the cached build suddenly
renders blank, look here first.

Also note the cached object carries its data, so a panel will not pick
up changed upstream data until you name it (or `--all`).

## Gotchas & warnings

- **`k10` is `GSP #10`.** The Rmds call the loadings `k`; the paper calls
  them GSPs. All axis titles here use the GSP naming.
- **Panel c's `max_expl` is computed over tissue + technical *plus* the
  other summary columns.** That looks like a bug but it drives both the
  "Max explained" column and which GSPs get an `[Unknown]` label, so the
  cache builder reproduces it verbatim. Do not "fix" it without checking
  the row labels against the reference.
- **Panel c's cache stores raw covariate names in row labels**
  (`GSP 3 [avg_len]`); `compose.R` maps them through `panel_c_cols.tsv`'s
  `label` column to get the reference's wording (`[Avg read length]`).
- **Panel g's Bonferroni denominator includes a trailing `project_id`
  column**, so it is `50 x (n_genes + 1)`. Kept for parity — changing it
  moves the dotted threshold lines off the published position.
- **Panel f drops the Rmd's `geom_text` annotation layer** (the
  `p=…/FC=…` strings above each box). The reference has no annotations,
  and that layer is what forced the Rmd's y-axis up to ~22.
- **Panel a must keep `preserve_aspect = TRUE`.** Its art is 5.83:1;
  filling the cell stretched the schematic's text ~19% horizontally.
  `ROW_H[1]` is tuned so the CELL aspect matches the art. If you change
  row 1's height, re-derive it (`8.5/5.83` plus the title's height).
  Because `ROW_H` is absolute inches, changing OTHER rows no longer
  disturbs this one.
- **Panel c's row labels sit in ONE column hard against the matrix**,
  with collisions resolved VERTICALLY. An earlier near/far two-column
  scheme had to place the far column past the longest label, pushing
  the label block out to ~45% of the matrix width and squeezing panel
  b's CDF. `MIN_DY` converts the label's point size into row units.
- **Panel c pins its y range with `coord_cartesian(ylim = C_YLIM)`, and
  this is load-bearing.** There are more labels (34) than the matrix is
  tall (50 rows at `MIN_DY` ~ 1.9 = 66 rows needed), so the block spills
  below the matrix. ggplot sizes a panel from EVERY layer's data range,
  so without the pin that spill stretched the panel ~0.6 in, leaving
  dead space under the matrix and dragging the x-axis tick labels down
  with it. With the pin, the spill draws outside the panel
  (`clip = "off"`) into the x-axis label band — harmless, because the
  GSP labels sit at `x > x_lab`, to the right of the last column, while
  the rotated tick labels only occupy the matrix's own x range.
- **Panel titles are aligned BY MEASUREMENT, not by construction.**
  The panels are not the same kind of object — b and d are patchworks
  inside their `wrap_elements` (b carries an `inset_element`, d an
  inner `|` layout), c/e/f/g are plain ggplots — and patchwork adds its
  own plot margin above a nested patchwork's title. Title tops at
  400 dpi from each row's top: b 150 / c 151, d 68 / e 69, f 15 / g 15.
  The corrections are c's `plot.margin` top (12 pt) and d's
  `plot.title` margin top (9 pt). Re-measure after ANY change to a
  panel's wrapper; the helper is in the scratchpad as `titletop.R`
  (crop a 230 px strip at each row boundary with `pdftoppm -r 400`,
  then scan for rows with a wide run of dark ink).
  NB the panel TAGS are already aligned — a/b/d/f have ascenders and
  c/e/g do not, so their ink tops differ by ~7 px at 200 dpi while
  their baselines match. That is type, not misalignment.
- **Panel c needs an explicit left `plot.margin` (22 pt).** It has no
  y-axis text, so nothing else holds the matrix off the cell's left
  edge — at 2 pt the heatmap started to the LEFT of panel c's own "c"
  tag.
- **`C_SPILL_ROWS` caps how far that spill may go** (`205/948` of the
  matrix height, measured at 400 dpi). `wrap_elements()` clips at panel
  c's plot area, so anything past the cap is silently shaved off — at
  the axis band's full `273/948` the last label, GSP 49 [LUSC],
  disappeared entirely. `MIN_DY` is `min(preferred, band/(n-1) * 0.95)`
  so the labels tighten rather than vanish. `compose.R c` prints the
  block's extent and the floor; check that line after any change to
  panel c's height or label count.
- **Panel c's block names (Tissue/Technical/Summary) are text ON the
  strip, not a legend.** The legend block cost ~25% of the panel's
  width for three words. "Summary" is wider than its 4-column block, so
  it is left-anchored at the block start rather than centred.
- **Panel e's percentages are on the wedges** (two narrow ones outside
  with a leader). Quinn: percentages in a legend are never acceptable.
  Slice midpoints are computed explicitly as
  `1 - (cumsum(prop) - prop/2)` (position_stack puts the first level on
  top) so text and leaders share exactly the same y.
- **Panel b's inset has a WDR11 gene body row and a connector** from
  WDR11's point on the TCGA-brain ECDF (its mean PSI in the top-100
  GSP #1 samples, x = 0.0448) to the inset's top-left corner. The inset
  sits lower-right so it does not cover that point; the legend sits in
  the empty pocket below the curves and left of the inset. Arcs use
  `lw_scale = 0.45` -- full-width stroke was far too heavy at inset size.
  2026-10-05 (Quinn): WDR11 is named ONLY under its gene body (no label at
  the ECDF point); the box bottom is dropped by B_BOX_BOT_PAD so that label
  sits inside it; leaders to the box's top corners are dashed; title is the
  single line "GSP #1 represents neuronal junctions".
- **NEVER use `limits =` on panel b's x scale.** A scale limit DROPS
  rows before the stat runs, so `limits = c(0, 0.11)` handed
  `stat_ecdf` a truncated sample (141 junctions cut): the drawn curves
  were the ECDF of that subset, all converging to 1.0 at the right
  edge and sitting above the true ECDF everywhere. `b_wdr11_y` comes
  from `ecdf()` on the FULL vector, so the WDR11 point was plotted off
  its own curve — which is what "the line doesn't connect the proper
  curve" actually was. The view is set with
  `coord_cartesian(xlim = B_XLIM, ylim = B_YLIM, clip = "on")`
  instead. `clip` must stay `"on"`: the curves now genuinely extend
  past 0.11, and with clip off those tails were drawn across the whole
  figure row. The same trap is why `sashimi_track_gg` also switched
  from scale limits to coord.
- **`inset_element()` takes PANEL NPC; the connector is drawn in DATA
  units.** Treating the npc corner as a data coordinate left the line
  stopping short of the inset. `b_npc_to_data()` maps through the
  scale expansion and MUST stay in step with `B_XLIM`/`B_YLIM`/`B_EXP`.
- **The callout is a WEDGE: a light box around the inset plus two
  leaders from the WDR11 point to the box's two TOP corners**
  (`b_wedge`, one data frame so both share a layer). A single leader
  does not work here — with the inset raised there is only ~0.6 in
  between the point and the arcs, so one line either terminates in
  mid-air inside the box or runs nearly horizontally along the
  GTEx-brain curve. The wedge reads as "this point expands into this
  box" and the box edge is what the lines terminate on. The leaders
  necessarily cross the other ECDF curves: the inset sits below all
  four, so anything joining a point on the orange curve to it must.
  The corners come from `b_npc_to_data()` — `inset_element()` takes
  panel npc, the annotations are in data units.
- **`apex_range` must be checked against the track's own peak height.**
  It is a fraction of `y_max`, and `y_max` scales with drawn track
  height — so shortening the control track to 0.58 of a row left its
  `y_max` at only ~1.7x its coverage peak, and the default 0.545 apex
  landed AT the peak. The 98% arc drew as a straight line. It needs
  `c(0.88, 0.82)` to read as an arc at that height.
- **The inset's cohort captions hang BENEATH their own coverage
  track**, left-aligned, in the track's colour (`label_below = TRUE`,
  `label_hjust = 0`, `label_color`, all added to `sashimi_track_gg`).
  A caption placed inside the track floats in the arc headroom and
  ends up nearer the track ABOVE it than its own coverage. `B_GAP_PT`
  is `B_LAB_PT * 1.9` so the gutter actually holds the caption.
- **Matched controls use `PAL$cohort["TCGA-other"]`, not GTEx-brain.**
  They are the non-brain TCGA samples; the old blue named a cohort
  that is not in that track at all.
- **`ROW_H` is in INCHES and `PAGE_H = sum(ROW_H)`.** patchwork treats
  row heights as relative, so the original relative vector meant that
  growing one row silently rescaled every other row -- and several
  panels derive text-collision spacing from their own drawn height
  (`C_CELL_H_IN`, `D_ZOOM_H_IN`), so their constants went stale
  invisibly. With absolute inches, growing row 3 leaves rows 1/2/4
  untouched. The page is 8.5 x 11.7, not 8.5 x 11.
- **Panel d gets the lion's share of row 3** (`COL_W$de = c(1.75, 1)`,
  zoom `D_SPLIT = c(1, 2.1)`, and row 3 is 0.70 in taller than the
  reference's ink band): the mutation labels only become readable with
  that much room for `geom_text_repel` to work in. Quinn's rule here --
  stretch the plot rather than shrink the type.
- **Panel d's labels are split into THREE repel layers at
  `D_TOP_BAND = 97`.** Every class's CDF tops out at exactly 100%, and
  the unlabelled navy "No SF3B1 Mutation" cloud occupies that band with
  a step line running the full panel width. ggrepel only avoids its own
  layer's data, so it cannot see that line -- it drew it straight
  through "H662R". Labels at or above the band are therefore lifted
  into headroom (`d_top_ylim = c(101.5, 108)`) and everything else is
  capped below it (`d_repel_ylim = c(-2, 96)`). Disjoint y ranges also
  guarantee the layers cannot collide with each other.
- **`D_GAP_X`'s 4.5 pt term sets the corridor between the two label
  columns**, which is `2 * D_GAP_X` wide. At 2.5 pt the left and right
  red columns printed flush ("K700E|K700E") where the diagonal is
  steepest.
- **Panel d's title must be set with zero top margin.** It is the only
  panel titled via `plot_annotation()` rather than `labs()`, and
  patchwork stacks that title band on top of its own 5.5 pt plot
  margin -- which put panel d's content 39 px below panel e's in the
  same row. `plot.margin = margin(0,0,0,0)` in the annotation theme
  fixes it; `figure_qa.R` flags the regression if it comes back.
- Panel b's sashimi region is hardcoded (`chr10:120,885,900-120,889,950`)
  to frame the WDR11 junctions; the `.links` span 120,885,938–120,889,894.

## Open questions for Quinn

1. **Panel b's "non-brain" control label is not what the code does.** The
   Rmd filters with
   `!project_id %in% c('TCGA-LGG, TCGA-GBM, TCGA-PCPG')` — one
   comma-joined string, so it never matches and the 400 controls were
   **not** restricted to non-brain samples. The published panel labels
   them "Non-brain TCGA samples (400)". This is a claim/code mismatch,
   not a cosmetic one, so I reproduced the code faithfully and labelled
   the track "Matched controls (400)". Needs a decision.
2. **The published inset's PSI values don't match the current data**
   (published 71.9 / 1.22 / 7.33 / 4.78 / 96.7 vs on-disk 71.2 / 1.27 /
   7.5 / 4.89 / 98.2). The figure appears to predate the current
   outputs. Current values are used here.
3. Panel c's reference legend reads ±0.4 while the code used ±0.8 with
   squish (data truncated to ±0.7 upstream). ±0.4 is used here to match
   the published legend.

## Input/output

```
softimpute_flash_300_qc_filtered.RDS (3.5G) ─┐
tcga_v2_samples.tsv (68M) ───────────────────┤
merged_grch38_with_am.tsv (536M) ────────────┼─► scripts/build_cache_panel_{b,c,g}.R
standardized_covariates.RDS ─────────────────┤    (SLURM, ~minutes-hours)
per-gene res.RDS / dat.tsv / psi_table.tsv ──┘         │
                                                       ▼
                                                  cache/*.tsv  (small)
                                                       │
SF3B1_full_info.tsv (834K) ──────────┐                 │
fact_10_top_junctions.tsv (3.4M) ────┼─────────────────┼─► compose.R ─► out/figure2.pdf
sf3b1_top500_alt3_fc5_psi.tsv (19M) ─┤                 │
panel_a/genome_wide_workflow.pdf ────┘                 │
WDR11 bigwigs + .links ────────────────────────────────┘
```

## Running it

```
bash scripts/submit_caches.sh     # once, on SLURM; writes cache/
Rscript compose.R
```

`compose.R` renders placeholder boxes for any panel whose cache is
missing, so the layout stays reviewable while a cache is still building.

## Later changes (2026-09-23)

Four fixes Quinn asked for after reviewing the rendered PDF.

- **Panel a's schematic is re-set in the figure's own typeface.** The
  Illustrator art's labels are Arial (plus a TeXGyreHeros-Bold title);
  every R panel draws in `FONT_FAMILY` -> Nimbus Sans. `panel_a/build.R`
  now pipes the PDF through `../scripts/pdf_subst_font()`, which
  strips the embedded font subsets with `pdftops -noemb*` and lets
  Ghostscript resolve the names through its Fontmap. Metrically
  compatible, so the artwork does not reflow.
  `genome_wide_workflow_helv.pdf` is GENERATED (regenerated whenever
  the source is newer) — the `.pdf` beside it is still the source.
  The round trip widens hairline strokes (the black arrows carry a
  purple outline) to ~1 px, hence `dpi = 600` at the call site.
  Same treatment in `figure1/compose.R`.
- **Panel b's inset is smaller and its arcs are flatter.** `B_INS` went
  from 0.60 x 0.74 npc to 0.50 x 0.60 npc (~33% less area).
  `apex_range` went from `c(0.84, 0.40)` to `c(0.68, 0.37)` on the
  TCGA track and `c(0.88, 0.82)` to `c(0.78, 0.74)` on the controls.
  Both were back-solved for **~30% less rise above each arc's own
  endpoints**, which is what the eye reads — a flat 30% cut of the raw
  apex fraction would have put the control arc's apex ON its
  endpoints, i.e. a straight chord. `min_arc_rise` stays 0 (off); see
  the trap note in `../scripts/coverage_plot.R`.
- **Panel b's callout box is padded on the left by `B_AXIS_PAD`.**
  `sashimi_track_gg` draws its y-axis at `start - n_pt` in data units
  with clip off, i.e. OUTSIDE the npc region `inset_element()` gives
  it, so the "38.7"/"19.7"/"0" numbers fell outside the callout frame.
  0.055 npc (measured: 65 px overhang against 1400 px per npc unit at
  400 dpi) puts them back inside.
- **Panel c's heatmap grew into its own whitespace.** Two causes, both
  removed: (1) `scale_y_continuous(expand = expansion(add = c(0.6,
  c_top_expand)))` added the headroom a SECOND time on top of
  `coord_cartesian(ylim = C_YLIM)`, which already carries it — ggplot
  applies scale expansion to the coord limits. Now `expand = c(0, 0)`.
  (2) The title's default bottom margin left a 0.13 in blank band
  above the Group strip; `plot.title` now has a zero margin (its TOP
  does not move, so the b/c title alignment is untouched).
  Measured on the 400 dpi render, panel c's tile block went from
  **2.400 in tall to 2.652 in (+10.5%)**, width unchanged at 2.487 in,
  and the largest interior blank row band from **0.300 in to
  0.072 in**. `C_VP_H_IN` was re-measured to `0.672` of the cell
  accordingly (it drives `row_pt`, `MIN_DY` and `C_FLOOR`).
  Rows 3 and 4 of the page are pixel-identical before and after.

## 2026-10-05 layout changes (Quinn)
- `COL_W$bc` 1.12/1.38 -> 1.00/1.50: 2b slightly narrower, 2c starts further left.
- 2c GSP labels sit LEVEL WITH THEIR OWN ROW. Runs of adjacent rows alternate
  near/far columns (phase chosen so long labels go far); far offset = widest
  near label; a 2-D greedy search (boxes + leaders) is the fallback. Label
  widths and the legend are measured on cairo; panel drawn width is derived
  from PAGE_W/COL_W. Prints `[panel c] N labels: near/far, max |dy|` -- expect 0.
- 2d: step lines + all points are blank-label obstacles in the repel layers
  (ggrepel can't see other layers). Near-hotspot (blue) labels are placed
  deterministically left of their point, lifted bottom-up when crowded, capped
  at the local red left-column edge; a blue point inside the navy band goes to
  the headroom layer.

## Zoom connectors (2026-10-05)
The two lines from the CDF zoom box to the zoom panel end at constants in the
CDF's data units (`*_CALL_DX/_TOP/_BOT`). These are MEASURED, not tuned: run
`../scripts/measure_zoom_connectors.py out/<fig>.pdf x0 y0 x1 y1` (crop = page
fractions around the CDF + zoom) after any layout change to that row, and paste
the DX / TOP / BOT it prints. Hand-tuned values had drifted off the corners.
