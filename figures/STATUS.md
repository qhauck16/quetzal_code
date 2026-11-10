# Quetzal paper-figure pipeline — status

Paused 2026-09-22 ~17:10 at a session limit. Resume from here.

    cd quetzal_code/figures/<figureN>
    /scratch/midway2/qhauck/conda_env/rstudio-server/bin/Rscript --vanilla compose.R          # all cached panels
    /scratch/midway2/qhauck/conda_env/rstudio-server/bin/Rscript --vanilla compose.R c d      # rebuild only c and d
    /scratch/midway2/qhauck/conda_env/rstudio-server/bin/Rscript --vanilla compose.R --all

## Where each figure stands

| | a | b | c | d | e | f | g |
|---|---|---|---|---|---|---|---|
| **fig 1** | done | done | done | done (d/e combined) | — | — | — |
| **fig 2** | done | done | done | done | done | done | done |
| **fig 3** | done | done | done | done | done | done | — |
| **fig 4** | done | done | done | done | done | done | done |

Figures 1 and 2 are **finished and signed off**. Do not disturb them except via
the regression check below.

## RESUMED 2026-09-22 19:40 — figure 3 is COMPLETE

All three stopped-agent panels (fig 3 b/c/e) were inspected at 400 dpi and are good.
Figure 3 panel f (the screenshot) is built. Figure 1 re-verified at 0 differing pixels
against its pre-agent baseline, twice, so the shared coverage_plot.R edits are clean.

Fixed on resume:
  * fig 3 panel a — zoom shared the CDF's y range, crushing 47 labels into the top
    third. Gave it its own range, pinned PRPF8+FAM50A left / SLU7 right / NKAP
    alternating, loosened force_pull, and grew ROW_H[1] 2.55 -> 2.95 in.
    RESIDUAL: ~6 labels still knot mid-panel, where all four classes converge at
    ~88-93% cumulative. Four approaches tried: per-class flank pinning, loosened
    force_pull and a taller row all helped; a deterministic per-flank pre-spread
    made it WORSE (it places by rank not proximity, so leaders grow long and cross
    the other column) and was reverted with a comment in build.R so nobody retries
    it. The band is at genuine capacity; the panel is now cleaner than the
    published reference.
  * fig 4 panel c — the stopped agent left two bugs: the tile text was bound to
    `label`, which after its join is c_cols' COLUMN HEADER (so "SF3B1 introns" drew
    once per tile), and scale_size_identity() reads size in MILLIMETRES so 8 pt
    rendered ~2.8x too large. Both fixed; all 24 betas + 3 Ns match the verified table.
  * scripts/coverage_plot.R gained `min_arc_rise` (default 0 = OFF). See the note in
    the file: a caller that back-solves apex_range from chosen midpoint heights --
    which fig 3 panel c does, so it can place PSI labels deterministically -- must not
    have its apexes raised, or its labels get struck through by their own arcs.

## HISTORICAL: unverified agent output on disk (now resolved)

Five subagents were stopped mid-flight. Three had written files that
**compose.R loads without error but that NOBODY HAS LOOKED AT**:

    figure3/panel_b/build.R   595 lines   NKAP/PRPF8 lollipops
    figure3/panel_c/build.R   442 lines   four sashimi quadrants
    figure3/panel_e/build.R   285 lines   expression ECDF

Their last reported states before being stopped:
  * panel b — "font metrics are quantized to whole points at these sizes,
    fixing the measurement"  (mid-fix; treat the label layout as suspect)
  * panel c — "now update the call sites and the axis break rule"  (mid-edit;
    it may also have edited the SHARED `scripts/coverage_plot.R` — see below)
  * panel e — choosing a legend size/placement  (probably closest to done)

**First job on resume: render figure 3 at 400 dpi and LOOK at b, c, e before
trusting any of them.** Each panel's `build.md`, if present, says what its agent
intended.

Two agents (figure 4 panels c and f) wrote nothing. Panel c's agent did report
"beta table reproduces exactly (all 24 cells + all 3 Ns)" before stopping, so
the model spec in its brief is confirmed correct.

## Regression check before trusting anything shared

`scripts/coverage_plot.R` is shared by figures 1, 2 and 3, and figure 3's panel-c
agent had permission to extend it. Figure 1 is now seeded and byte-reproducible,
so any pixel difference is a regression:

    cd figure1 && Rscript --vanilla compose.R      # twice; diff the renders
    cd figure2 && Rscript --vanilla compose.R

Rasterize at 110 dpi and compare; expect **zero** differing pixels.

## Decisions Quinn has already made — do not re-ask

* Published figures carry three factual errors, all to be FIXED in the rebuild:
  - fig 3 panel d title says "GSP #20" → **use #13** (NKAP is k13 everywhere).
  - fig 3 panel c labels SENP2 as "chr2" → **chr3**.
  - fig 3 panel c PRKDC y-axis reads 5.4 → **use 9.4**, the ini's `max_value`.
* fig 4 panel c row label says "log(intron length)" but the model fits
  `scale(length)` on raw bp → label it **"intron length"**. (Raised, not yet
  explicitly confirmed — it is the last open question.)
* x-axis tick labels: **always 45°, never vertical**. Use `angled_x_labels()`.
* No pygenometracks, no ggbio. Own geoms from `scripts/coverage_plot.R`, reading
  raw bigwigs + `.links`.

## Verified data facts (do not re-derive)

* Factor keys, all current: SF3B1 = k10, NKAP = k13, YLPM1 = k30.
  (`k20` / `k65` anywhere is stale.)
* fig 3 d: `unified_results/nkap/20v60/dpsi0.01_fc2.5/df_final.tsv`, 2806 rows.
* fig 4 c: `unified_results/{sf3b1,nkap,ylpm1}/20v60/dpsi0.01_fc2.5/junctions_plus_characteristics.tsv`,
  N = 662 / 1150 / 762. Binomial GLM, all predictors `scale()`d, filtered
  `bp_distance_to_3ss < 100` and `group_by(subclass) |> filter(n() == 2)`.
* fig 4 d: `unified_results/ylpm1/20v60/dpsi0.01_fc2.5/df_final.tsv`, 2460 rows.
  (Quinn recalled dPSI > 0.05; the data says 0.01.)
* fig 3 e: the two 13.4 GB expression matrices are already sliced to the 20
  columns that matter — `figure3/cache/panel_e_{scaled,raw}.tsv`, 24 MB total.
  Never `fread` the originals.
* fig 4 f: the original is a **heatmaply/plotly widget**, printed from Safari.
  It cannot go into patchwork; rebuild natively. Clustering is euclidean +
  complete + seriation OLO, rows only, columns in `tcga_ordered_juncs` order.
  The three annotation bars SHARE one 18-level palette (hardcode it).

## Figure 3 QA — PASSED 2026-09-22 20:15

  figure_qa  2.95,3.35,2.45,2.45  "0.45,0.55;1;0.52,0.48;1"   -> no dead bands
  title_align: row 1 a/b aligned to 2 px (after +2 pt on panel a's
    plot_annotation title -- a is titled via plot_annotation and b via labs(),
    which patchwork spaces differently); row 3 d/e within 6 px.
Figure 3's page is 8.5 x 11.20 in (row 1 grew to 2.95 for panel a's labels).

## Figure 4 panel g — DONE and verified

604 target genes (Leaf2 !Coding & !UTR, dPSI>0.01, FC>2.5, NO alt-3'ss filter),
p = 1.4e-05, 10/10 highlight genes present and labelled. All match the reference
PDF at diff_splicing/expression_CDFs/publication/YLPM1_leaf2_psi0.01_fc2.5_bootstrap.pdf.
Inspected at 400 dpi: no label collides with the curves, ribbon, legend or p-value.
NB its cache goes the OTHER way from figure 3 panel e: the bootstrap draws from the
whole column pool, so the 7.7 GB matrix is sliced by ROW (604 target genes, all
11,173 columns -> cache/panel_g_expr.tsv, 134 MB), and the computed curves/ribbon
are then cached as small TSVs.

## A FOURTH stale label found in the published figures

Figure 4 panel b's track labels read "High GSP #13 (15)" / "Matched Controls (150)".
BOTH parts are wrong — those are the NKAP figure's numbers. The ini
(GSP_pygenome/YLPM1/*/sashimi_tracks.ini) says `High_k30 (13)` and
`Matched_Controls (130)`, and panel a of the same figure says "top 13 samples".
Corrected to "High GSP #30 (13)" / "Matched Controls (130)".
Running tally of published errors corrected: fig3 d "GSP #20"->#13,
fig3 c "chr2"->chr3, fig3 c PRKDC 5.4->9.4, fig4 c "log(intron length)"->intron
length, fig4 b "#13 (15)"->"#30 (13)".

## Figure 4 panel f — DONE and verified

Rebuilt natively (geom_raster + geom_tile strips + dendextend segments); heatmaply
is called ONCE, cache-guarded, only for the 512-junction column order. Row order
reproduces the published 17 exactly with the clean 11-NO / 6-YES split, asserted
with stopifnot.
Two findings from its agent worth keeping:
  * `group_by(cluster_id)` must use the **TCGA** LeafCutter cluster (from the
    psi_files/tcga psi_table.tsv headers), NOT GIST's own cluster from V1. Both
    give a near-identical matrix and the same 11/6 split, but GIST's clusters
    PERMUTE dendrogram rows 2-9. Only the exact published row order catches it.
  * A transpose bug (`z = as.vector(t(mat))` paired with `x = rep(..., each=N_ROW)`)
    renders a plausible but scrambled heatmap. Pinned with corner/interior
    stopifnot spot-checks.

## ALL FOUR FIGURES COMPLETE — 2026-09-22 20:45

figure_qa + title_align pass on both new figures; figure 1 re-verified at 0
differing pixels after every shared-helper change (4 checks total).
  fig 3: 8.5 x 11.20 in, rows 2.95/3.35/2.45/2.45, cols "0.45,0.55;1;0.52,0.48;1"
  fig 4: 8.5 x 11.00 in, rows 2.50/2.00/3.50/3.00, cols "1.15,1;1.15,1;1;1.09,1"
Titles align within 2 px on every multi-panel row of both figures.

## 2026-09-22 21:00 — QUALITY PASS IN PROGRESS

Quinn reviewed figures 3 and 4 and rejected them: "an unbelievable number of errors
in the PSI curves, labels, even positions of the CDF curves themselves, overlapping
text and points". He is right; several panels were signed off from low-dpi renders.
STANDING INSTRUCTION: "Reiterate all night if you need to make these figures
perfect... be sure to check whitespace and how the pdf output actually looks, not
making guesses as to how it might look."

Defects I confirmed myself at 400 dpi before briefing:
  fig4 b  TRIM38's lavender taper renders as a STEPPED/NOTCHED polygon, not a clean
          trapezoid (KDM4C's is correct — compare). Its arcs cross chaotically and
          several terminate mid-intron at the baseline with no coverage. Leftmost
          coverage peak clipped. "1.27"/"2.13" crammed.
  fig3 e  Legend is enormous (4 entries, 2 wrapped) and dominates the plot. Leader
          lines cross other labels AND the curves. Gene labels pair up and read as
          compound words. p-value sits inside the label ladder.
  fig3 b  A below-body circle on the PRPF8 track is CLIPPED by the panel edge; large
          dead whitespace between the two cartoons.
  fig4 e  ASXL1's gene label is drawn ON TOP OF its taper; gene-name side is
          inconsistent across the four loci.
  fig3 a / fig4 a  residual label/point overlaps in the zooms.

Five fix agents launched, one per problem area. All briefed with: LOOK at the
rendered PDF, iterate many times, do not stop at "better".

### Pass 1 results (verified BY ME at 400 dpi, not taken on report)
* **fig3 panel b — DONE.** Legend pitch now fills its column (was capped, leaving a
  0.5 in void); legend type solved against the column width rather than picked
  (4.96 -> 5.46 pt); band reservation made exact (`max(n*cd + ink)` not
  `max(n)*cd + max(ink)`); the caption row replaced by a derived drop that only
  pays for labels whose ink actually overlaps the caption's x-span. Net auto-fit
  scale 0.994 -> 1.078, so everything is ~8.5% larger in the same cell. A
  measured blank-band scan finds no run over 0.085 in.
  Extra defect the agent found and fixed: the larger legend type made the centred
  `2305` tick print on the `*Called via WXS` footnote; end ticks are now aligned
  INWARD (hjust 0 first / 1 last), which also stopped `0` hanging outside the cell.
  NB **my "clipped PRPF8 circle" defect was MY error** — I cropped -H 1080 when
  row 1 is 1190 px. Both below-body circles were always complete.
* **fig4 panel e — DONE.** Gene-name side is now DERIVED, not hardcoded: the name
  is drawn in the gene row that the lavender highlight spans, so it lands on the
  taper whenever the zoom window reaches that end. The rule picks an end the window
  does not reach (tested with the label's real drawn width in bp + 0.07 in) and
  prefers left. TP53BP1/PBRM1/MSH6 -> left, ASXL1 -> right (forced: its window
  starts at tx_lo). Extra fix: consecutive tracks butt flush, so a track's `0` and
  the next track's y_max sat 0.027 in apart and read as one pair; now 0.074 in.

### Pass 2 — NEW defect I found myself at 400 dpi (agent launched)
* **fig3 panel c: the arcs SAG at every locus.** They leave the left exon's coverage
  top, droop to the baseline, and rise into the right exon — drooping lines, not
  arcs — and the PSI labels then float in empty space above them, disconnected.
  Worst at PRKDC (70.8/22.1, 92.4/2.95). I had only verified SENP2 closely before,
  where it is least obvious. THIS IS PART OF WHAT QUINN MEANT BY "errors in the PSI
  curves".
  ROOT CAUSE in `arc_geom()`: `ARC_MID <- c(0.26, 0.08)` is the arc's midpoint
  height as a fraction of y_max, but the endpoints y1/y2 are the exon COVERAGE
  TOPS, which reach ~0.4-0.55 of y_max here. So every chosen midpoint is BELOW its
  own endpoints and the curve is forced to bow downward. The constants were picked
  as if the endpoints sat at the baseline. Fix = make the midpoint relative:
  `tgt <- pmax(y1, y2) + rise`, tiered, then re-solve the label stack above it.
  NB `coverage_plot.R`'s `min_arc_rise` (default 0) exists for this exact problem
  but OVERRIDES a back-solved apex and would break this panel's label placement.

### Pass 2 — fig4 panel b DONE (verified by me at 400 dpi)
* **Orange rule now read straight out of the factorisation; the 3.7 GB RDS is never
  loaded.** `flashier/fact_30_top_junctions.tsv` IS GSP #30's ranked top list (all
  250 rows, so row i = rank i — asserted, not assumed). Each row's argmax
  unannotated junction pins the factor: **KDM4C.k31 rank 40, TRIM38.k21 rank 55**
  (TRIM38's second group is rank 114, outside the top 100). Junctions coloured =
  poisson_NMF_junction_plots.Rmd's own rule, `loadings >= group_top/2` on
  `poisson2multinom(res)$F` restricted to `annotated == 0`.
  Cached to cache/panel_b_gsp30_juncs.tsv.
  RESULT: KDM4C -> 10.7 + 11.2 (as before). TRIM38 -> **6.91 (confirmed by Quinn)
  and 1.27**. The old heuristic's 22.7 is **k24**, a different group, and 8.38 is
  k22 = rank 114 — both correctly black now.
  >> NEEDS QUINN: the rule as stated also makes **1.27** orange (it is k21's
  >> single highest-loading junction, 0.425 vs 6.91's 0.284) and removes 22.7.
  >> The published reference shows 6.91 + 22.7. Confirm which is intended.
* TRIM38's taper was a ~2 degree sliver because its window is 42% of the transcript
  (KDM4C's is 1.2%, which is why only KDM4C looked right). Local `taper_gg()` now
  narrows at NECK_F = 0.55 and drops a vertical neck. That also fixed an unlisted
  defect: the italic gene name is 0.076 in tall in a 0.12 in row, so it always
  spills into the taper row and was sitting ON the lavender at both loci.
* Arcs use `y = apex*sin(pi*t)^0.78` so the feet leave the baseline with diverging
  slope and visibly LAND on the splice site instead of melting into it.
  Arc feet measured against the .links coordinates at 400 dpi: all seven TRIM38
  feet within 1.5 px = 0.004 in.
* PSI labels are now CONTOUR-STYLE — each number sits in a break cut in its own
  arc. A white glyph halo was tried first and renders as a strike THROUGH the
  number (an outline only thickens glyphs), so cutting the arc is the fix. Nothing
  white is painted; the PDF stays vector-clean.
* `saveRDS` round-trip checked: cached-only vs fresh build render 0 differing px.
* `../scripts/coverage_plot.R` was NOT modified (md5 unchanged), so figures 1-3
  cannot have moved and no regression check was owed.

### Pass 2 — fig3 panel e DONE (verified by me at 400 dpi), now re-tuning to 8 pt
* Gene labels rewritten as THREE HAND-PLACED LADDERS, not ggrepel. Every name on
  its own row, right edges pinned per tier, leaders leave rightward from the
  column so a leader can never enter a box in its own tier. p-value moved to
  y = 0.78 right-aligned (panel g's placement). Legend tightened.
  Why not ggrepel: its equilibrium for a packed column is "boxes just touching",
  box.padding is the only lever, and at 0.22 the boxes stop fitting the band and
  repel SILENTLY OVERLAPS them (all of one tier printed as one block).
* It was forced to 7 pt because 16 names need 1.88 in at 8 pt and the panel was
  1.72 in. Per Quinn's "stretch the plot rather than shrink the type" I grew
  **ROW_H[3] 2.45 -> 2.95 in** (figure 3's page is now 8.5 x 11.70) and sent the
  agent back to re-tune. ROW_H is absolute inches so rows 1/2/4 did not move.
  **RESULT: both callouts and legend are back at the shared tokens (8 pt), 16
  names fit with no forcing.** Measured: min gutter 11 px = 0.028 in uniform
  (2 pt leading on 8 pt caps); min curve clearance 0.137 in at MAD2L2 (was
  0.08 in); inter-tier clearance 0.157 in; legend-to-nearest-label 0.081 in.
  Drawn panel went 3.40 x 1.72 -> 3.40 x 2.17 in; E_IN_PER_Y 1.3937 -> 1.7600.
  build.md records the dependency ORDER: if ROW_H[3] is ever shortened,
  re-measure E_IN_PER_Y first, re-derive E_PITCH and the bands, and touch
  E_LAB_PT LAST, not first.
* KNOCK-ON I had to fix: panel d shares row 3, its pie is HEIGHT-limited by
  coord_equal, so the taller row enlarged it and the measured D_PIE_R_IN = 0.31
  went stale — "4.95%" ran off the inset edge. Re-measured to 0.475 (pie is now
  0.95 in across) and re-boxed. This is the third time a measured-geometry
  constant has gone stale after a layout change; they are all commented with
  "re-measure if ROW_H changes".

### Pass 2 — fig3 panel c DONE (all four quadrants verified by me at 400-500 dpi)
* The sag was WORSE than my diagnosis: the cryptic site of each pair sits inside
  the exon plateau, so endpoint heights reach **0.805 of y_max** (PRKDC controls),
  not the 0.4-0.55 I estimated. `arc_geom()` now works in PEAK space:
  `peak = max(y1,y2) + ARC_RISE`, tiered by ARC_SEP, capped at ARC_CEIL. The sine
  amplitude that puts the DRAWN peak there is back-solved with `uniroot` — with a
  steep chord the true maximum is at t~0.35, not the midpoint, so parameterising by
  midpoint height would have overshot y_max by 5%.
  Values: ARC_RISE 0.20, ARC_SEP 0.16, ARC_CEIL 0.93, ARC_RISE_MIN 0.09.
* **Stacking order changed from "shortest junction highest" to "by endpoint
  height".** The longer junction of each pair is anchored deeper in the exon, so
  its chord dominates its partner's at BOTH ends in all eight tracks — the old
  order forced a crossing. `arc_min_gap()` now prints the SIGNED minimum gap every
  build as a regression check; all eight are positive (0.02-0.05).
* PSI labels rewritten: each number sits one clearance from ITS OWN curve (above
  where that fits, below into the empty bowl otherwise) and slides ALONG the arc
  until its box clears everything. Horizontal room is the cheap axis. Consistent
  rule in all eight tracks: upper arc's number above, lower arc's number below.
  Two real bugs fixed en route: (1) `for (t) for (sgn)` with a single `break`
  leaves only the INNER loop in R, so the search returned the LAST acceptable slot
  not the first — every number had drifted to t = 0.78; (2) the ceiling test
  double-counted the pad, costing 0.06 of y_max.
* CHANGE THAT NEEDED SIGN-OFF, accepted: the inline "-11 nt"/"-5 nt" is now drawn
  ONCE per locus, on the High track only. In the control tracks it printed ON TOP
  of the grey exon coverage (verified at 800 dpi) and cannot fit: y_max is the ini
  value and is driven by the controls, whose plateau fills 0.91 of y_max while half
  a line is 0.11. NB the reference is inconsistent here — EZH2 blue only, BRD9
  both — so this makes it a consistent rule. Tell Quinn.

### Known-structural, documented, NOT fixed
* fig3 panel c: PRKDC's two control arcs nearly touch near their SHARED right
  endpoint (signed gap 0.02 of y_max). They share e2 so they must converge; the
  reference converges too. Verified analytically as well as numerically that it is
  not a crossing.
* fig3 panel c: label->arc gaps are ~1.9 pt. Comfortable at 6.24 pt type but not
  generous; widening needs a taller ROW_H[2].

### I SILENTLY DEGRADED FIGURE 2 AND ONLY CAUGHT IT IN THE REGRESSION CHECK
`out/figure2.pdf` as shipped at 19:45 was built while coverage_plot.R still had my
faulty floor `apex_y <- max(apex_rows[i], max(y1,y2) + min_arc_rise*y_max)` with
min_arc_rise = 0 — which still clamps apex up to max(y1,y2), making
`apex - max(y1,y2)` exactly 0 and drawing arcs as STRAIGHT CHORDS. Figure 2's
panel b 71% arc was flattened. The conditional fix at 19:47 restored it, but the
shipped PDF predated the fix.
Caught by rebuilding figure 2 from scratch and diffing against the shipped PDF:
683 differing px, all in panel b's inset arcs. After the rebuild, two consecutive
builds are 0-diff, so figure 2 is correct AND reproducible again.
LESSON: "figure 2 is signed off" is only true of the CODE. After any shared-helper
change, the shipped PDF must be REBUILT, not just left alone. A stale artifact is
not a frozen one.

### FINAL QA — 2026-09-23, ALL PANELS PASS MY OWN 400 dpi INSPECTION
  figure1  0 differing px vs pre-agent baseline
  figure2  0 differing px vs post-fix baseline (and correct again — see above)
  figure3  figure_qa clean; titles row1 a40/b38, row3 d20/e21  (<= 2 px)
  figure4  figure_qa clean; titles row1 30/30, row2 29/29, row4 19/20 (<= 1 px)
  fig3 page 8.5 x 11.70 in; fig4 page 8.5 x 11.00 in
The figure_qa "content tops differ" flags on fig3 row3 / fig4 row2 are TAG and
inter-row bleed, not titles — title_align is the definitive measurement.

### fig3 panel a + fig4 panel a — DONE, ggrepel removed from both
Replaced by `cdf_place_labels()`, a deterministic 2-D solver shared verbatim by
both files. It knows every obstacle ggrepel is blind to: markers, the step
staircase expanded into explicit axis-aligned segments, previously placed labels,
and previously drawn leaders. Per label it searches a cost-ordered grid — slide
out in x on either flank, step in y, OR sit centred above/below the point (that
third anchor was essential: a centred label needs only a row's height, not the
string's width). Greedy is order-dependent, so: order by local text demand ->
promote failures and re-run -> 80 seeded-shuffle polish rounds keeping the
cheapest total displacement. Text widths measured on a throwaway CAIRO device
(pdf()'s AFM metrics are ~6% narrow and let labels print flush).
RESULT fig3 a: 47/47 placed, 26 leaders, 0 collisions. fig4 a: 12/12 placed,
**0 leaders needed**, 0 collisions. Both byte-reproducible.
Two real bugs fixed en route:
  * collision boxes were sized from MEASURED ink, but geom_text(vjust=0.5) centres
    each string on its own ink box, so a caps-only string and one with a
    below-baseline glyph sit differently in the same nominal row — an underscore
    printed as a bar through the digits below. Boxes now max(measured, 0.92 em).
  * pad_x landed entirely on the OUTER side, so there was zero padding between a
    string and its own marker (that is what printed "Y1506*" flush against its
    point). The anchor now sits pad_x inside the edge.
Recorded in both build.md as NOT to retry: hard constraints on leader-crossing and
leader-grazing (they starve the greedy so later labels need the relaxation escape
hatch, and a leader THROUGH a string is far uglier); inflating the leader box in y
as a hard test (stacked labels are ~0.8 pt apart, so it swallows its neighbours);
escalating a failure straight to the relaxation level rather than one step.

### QA STATUS (superseded — see FINAL QA above)
  figure1  0 differing px vs pre-agent baseline (checked 5x across the pass)
  figure2  rebuilt; 2 consecutive builds 0-diff
  figure3  figure_qa clean; row 3 titles aligned to 1 px after nudging panel d's
           top margin to 9.6 pt (panel d had drifted 20 px high when ROW_H[3] grew)
  figure4  figure_qa clean; rows 1, 2 and 4 titles all aligned within 1 px

### 2026-09-23 — FIGURE 3 RESTRUCTURED TO THREE ROWS
Quinn: "3a should be wider" and "arrange 3 so that f is not on its own at the
bottom, make stuff fit together in a rectangle and look nice together".
  ROW_H  c(2.95, 3.35, 3.05)            page 8.5 x 9.35 in  (was 4 rows, 11.70)
  COL_W  ab  = c(0.52, 0.48)            a widened, b takes the loss
         def = c(0.35, 0.38, 0.27)      row 3 is now d | e | f
Panel f no longer sits alone on its own row as a centred square with white space
either side.
KNOCK-ON, as always with this design: d's cell went 4.42 -> 2.98 in (title and
x-axis title both overflowed) and e's went 3.40 -> 3.23 in (its legend started
printing over the CTC1/CHFR/BRD2/BRCA1/ARID1A/RBBP8 ladder). Both are with agents.
VERIFIED BY ME at 400 dpi: panel a is good in the wider cell — 47 labels, no
text-on-text, zoom visibly roomier.

### PIES — DONE, both figures, VERIFIED BY ME at 400 dpi
`pie_with_leaders()` is gone, replaced by `pie_labels_on_wedges()` in
scripts/pie_chart.R. Every percentage now sits on its own wedge's mid-ray; when the
text will not fit inside the wedge it still sits on that ray but gets a filled
rectangle behind it. No percentages in the legend, no outside leader callouts.
Whether a label fits is SOLVED, not guessed: the box centre on the mid-ray must
satisfy an edge constraint (r >= support/sin(halfangle)) and a rim constraint
(quadratic in r); a wedge wider than a half circle is not convex at the apex so only
the rim binds. Boxed labels then slide out along the ray until clear of neighbours.
  fig4 d: 30.3 / 50.1 / 12.8 inside; 5.7 and 1.1 boxed at the rim.
  fig3 d: 88.06 inside; 1.60 / 5.06 / 4.95 / 0.32 all boxed and separated — this is
          the hard case, four wedges under 6% with three of them adjacent.
fig3 d also re-tuned for its narrower cell: title re-wrapped to 3 lines, x-axis
title to 2, D_PIE_R_IN re-measured.
STRING WIDTHS ARE MEASURED ON cairo_pdf, not the base pdf() device. pdf() uses
built-in Helvetica AFM metrics; cairo resolves "Helvetica" through fontconfig to
whatever the machine has, which is 12-18% WIDER here. Using pdf() metrics sized
every box ~15% too small, so text overhung its own white box and the collision test
let neighbours touch. Same class of error as every other measured-geometry bug in
this pipeline.
NB figure 2 panel e keeps its own older pie copy and is FROZEN — pie_chart.R is not
used by it.

### fig3 panel c ZOOM + panel e RE-TUNE — DONE, VERIFIED BY ME (600/400 dpi)
panel c zoom: the decorative geom_curve leader is switched OFF at the call site
(`nt_label = NULL, nt_site = NULL` -- it lives in the shared zoom_callout_gg, which
the agent did not own) and replaced by the FEET OF BOTH JUNCTION ARCS, anchored at
their true splice-site coordinates on the zoomed coverage and stroked by the main
track's own law 0.35 + 0.45*sqrt(psi/100). The nt offset now reads as the horizontal
distance between the two feet, and the thicker foot is the higher-PSI junction. The
label is plain and horizontal with no leader. Non-crossing is ENFORCED: end heights
fanned by ZOOM_ARC_SEP, signed minimum gap measured on the drawn curves each build
with a stop() if negative (0.130 box heights in all four boxes).
Convention documented: the foot is the last quarter of sin(t*pi). A literal mapping
would be a straight diagonal -- over the 56 bp SENP2 shows, the true arc is 0.03 of
y_max from straight -- and the box's y axis is the coverage rescaled to its own max,
so there is no faithful mapping to make.
KNOWN, disclosed: PRKDC's feet are short and steep because its zoom window puts both
donors at 86-91% across the box. Buying room means moving zoom_start, which shrinks
the approved source box on the tracks.

panel e: the 3.40 -> 3.23 in cell broke it HORIZONTALLY, not vertically --
E_IN_PER_X fell 25% (1.4718 -> 1.0967), so everything in x-units got 34% wider in
data units and the legend block reached x = -0.013, leaving no room for the top
tier. Re-measured first, then: E_PITCH is now DERIVED from cap height + gutter
rather than stated (stating it is exactly what goes stale); legend re-wrapped with
all four entries short (measured: the old wrap is 1.515 in, and shrinking the type
does NOT fix it -- 1.36 in at 7 pt, 1.27 at 6.5 -- so 8 pt stayed); tiers re-cut
5/4/6 -> 5/2/8 at the widest natural percentile gaps, because tiers A/B sit in the
legend's x-shadow and tier C's column is right of it.
DEVIATION ACCEPTED: the p-value moved from y = 0.78 to y = 0.40. At this width
"p = 4.2e-07" is 0.625 x-units, so right-aligned it reaches back to x = 0.175, and
at y = 0.78 the rightmost curve is at x = 0.302 with RBBP8's leader running to
x = 0.20 -- the leader printed through the "p". Clearing it at that row would need a
right edge past the panel. At y = 0.40 the curves reach only x = 0.047.

### fig3 panel b RESTRUCTURED — DONE, VERIFIED BY ME at 400 dpi
The 16-entry legend moved out of the right-hand COLUMN into a 4 x 4 BLOCK in the
dead upper-left. That column cost 0.67 in = 17% of the 4.08 in cell and was the root
cause of all three voids: it starved both cartoons of width, forced the broken axis
to spend 30% of NKAP on a mutation-free region, and still left the quadrant blank.
In the quadrant the block costs ZERO height. Its type size is bisected against the
width NKAP's rotated-label staircase actually leaves free inside the block's own y
band (each rotated label is a strip, so its leftmost ink at height y is computable).
Consequences: both bodies 3.23 -> 3.93 in (+22%); NKAP break re-cut 30/70 -> 24/76;
PRPF8 20/60/20 -> 15/24/47/14 with a NEW break at aa 1754 (already both a tick and a
domain boundary, so no tick can disagree with its box) because nothing is mutated
between 1612 and 1756 and that gap drew as a 0.9 in empty corridor; the "outside the
top 60" caption moved under the two circles it annotates with one short stub each,
its x SOLVED by sliding until the placement costs no height.
Largest blank band 3.09 x 0.40 in -> 2.92 x 0.22 in, a 48% reduction.

THREE REAL BUGS IT FOUND THAT I HAD LOOKED STRAIGHT PAST:
 1. `R402*/Q/Sfs*13` was printing THROUGH the TCGA-UCEC and TCGA-KIRC legend
    entries -- it overran B_X_RIGHT by 0.22 in. The leftward slide that was meant to
    fix it is bounded by the functional window's start and can simply FAIL; the
    residual is now returned and the auto-fit bisects on it as well as on height.
 2. Circles from ADJACENT STACKS overlapped -- N371K's red circle sat on top of
    R368Q's blue one in every render to date, including ones I inspected and passed.
    Two collision sweeps ran as separate loops, so the label pass moved x[i-1] right
    AFTER the circle pass had spaced x[i] off its old position; wherever the stack
    got shorter the gap fell below one circle diameter. Now one pass taking
    max(x[i], x[i-1]+cd, x[i-1]+perp).
 3. B_CAL_W / B_CAL_H were stale for the new cell width (+0.65% h, -1.4% v),
    re-measured by scanning the render for the rows carrying the protein bodies'
    >1000 px horizontal rules.

### fig4 panel b ARCS — DONE, VERIFIED BY ME at 400 dpi
Feet now land on the COVERAGE at s1 and e2 (ANCHOR_BP = 0 -- a 10 bp window already
walks up into the exon plateau and floats the foot). The two orange KDM4C arcs meet
exactly on the small cryptic-exon bump between them, which is the visual proof.
ARC_P 0.78 -> 0.50 and that is NOT cosmetic: near a sine an arc leaves its foot
almost linearly, so a long arc is flat near its own feet and a short arc sharing
that foot climbs out through it -- dominating TRIM38's 22.7 needed a peak of 1.09 of
y_max (off the track) at 0.78, and 0.51 at 0.50.
Heights are SOLVED bottom-up in a CONTAINMENT order (Kahn), each amplitude the
smallest meeting rise-above-own-feet, separation from nested neighbours, and a
dominance gap over the shared x; bisection + uniroot. Disjoint pairs get no
constraints, which keeps it affordable. Infeasible dominance constraints are
SKIPPED, not driven to the ceiling (that parks arcs on the ceiling and makes
everything above them cross). Ceiling DERIVED from the label solver's own ceiling.
Signed min gaps printed and asserted every build: +0.0226 / +0.0247 / +0.0300 /
+0.0300 of y_max. 4 forced crossings, all TRIM38 high, all from 1.62 whose interval
PARTIALLY overlaps four others -- that is in the junctions, not the drawing.
Arc-foot accuracy re-verified at 1600 dpi: 0.0002-0.0004 in.
Label fixes that only appear once feet are at coverage height: TWO break sizes (cut
the named arc generously, notch every other, or two neighbours' lines stop either
side of a number and it names neither); protect every arc's LANDING (FOOT_KEEP) or a
wide label on a short arc eats its own arc and sits on its neighbours' landings, so
all three stop in mid-air short of the exon; place most-constrained-first.

### fig4 panel e STILL BASELINE-ANCHORED — agent running
I verified at 600 dpi: MSH6's 85.4 arc leaves y = 0 even though the exon coverage is
tall at that donor. All four loci. I had checked this earlier by grepping for
peak_y/covat, found nothing, and wrongly concluded it used sashimi_track_gg's
anchoring -- it carries its own older duplicate machinery. GREP IS NOT INSPECTION.
Agent is porting panel b's solver across so the two panels in the same figure match.

### PROCESS NOTE — out/figure*.pdf races, and it got worse
With FOUR agents composing concurrently I hit `out/figure3.pdf` at **0 bytes**
(pdftoppm: "Document stream is empty") — an agent was mid-write. Poll for a
plausible file size before copying:
    for i in 1..5; do sz=$(stat -c%s out/figureN.pdf); [ $sz -gt 100000 ] && break; sleep 10; done
then `cp` aside and crop the copy. Never crop out/ directly.
`out/figure3.pdf` is a shared artifact and panels are cached to `panels/*.rds`, so
a compose that STARTED before your rebuild re-emits the PDF with the PREVIOUS
panel. An agent lost a cycle diagnosing "a render that didn't match the build log".
ALWAYS `cp` the PDF aside immediately after rendering and crop the copy.
* fig4 panel b: TRIM38 keeps exactly four arc crossings, all from the 1.62 junction
  (25,974,930-25,983,163) partially overlapping four others. Two arcs cross iff
  their intervals partially overlap, which no tiering can change when both feet sit
  on one baseline — only dropping the junction would remove them.
* fig4 panel e still carries the OLD duplicated arc/label machinery, so it will now
  look different from panel b (which has contour labels and sin^p feet). Lifting
  panel b's solver into coverage_plot.R and adopting it in e is the real fix.
* fig3 panel b's upper-left quadrant is quiet — the broken axis gives 30% of NKAP's
  width to a region with no mutations. Filling it means moving the 16-entry legend
  in as a two-column block and widening the cartoons: a restructure, not a tweak.

## RESOLVED THIS PASS
* "intron length" -> **"log(intron length)"** everywhere (Quinn confirmed). NB the
  fitted model is scale(length) on RAW bp, so label and predictor disagree; kept as
  instructed because the betas must match the published figure. Noted in build.R.
* Taper colour: **LAVENDER #ccccff everywhere** (consistency beats matching the
  figure-4 reference's pink). Quinn: "maintain consistent color".
* fig4 b/e orange junctions: must come from the **Poisson NMF junction groups** —
  the gene's junction group among GSP #30's top 100 groups, colouring the
  high-intensity junctions in it. NOT the geometric smallest-gap heuristic currently
  coded. Quinn confirmed TRIM38's answer is the **6.91** junction.

## 2026-09-23 — REAL DATA BUG IN fig4 PANEL a, FOUND BY QUINN, FIXED
Quinn: "the cdf in 4 is totally misaligned, none of the non-variant score higher
than the top 13 samples on the GSP with YLPM1 mutations".
He was right and it was a DATA bug, not layout. I had classified LoF from the
table's `frameshift` logical (and `plot_color`). Both are INCOMPLETE — a strict
subset that misses 5 samples, THREE OF WHICH ARE IN THE TOP 13 on GSP #30:
    rank 2  X1783_splice              flagged "Other YLPM1 Mutation"
    rank 3  YLPM1:S145fs              flagged "No YLPM1 Mutation"
    rank 9  YLPM1:A252fs,YLPM1:Q437*  flagged "No YLPM1 Mutation"
Drawing those as non-variants put apparent non-variants above the LoF samples and
destroyed the panel's claim.
FIX: classify from `proteinChange` with `fs|\*|splice`. 120 LoF vs the column's
115; no sample the column calls TRUE is excluded. The boundary is then exactly
clean — ranks 1-13 all carry a LoF protein change, rank 14 is the first with none,
so the title is literally true of the data. Zoom now shows 15 labels (was 12).
LESSON: never trust a pre-computed convenience column in these tables without
checking it against the underlying strings. figure 3 panel a already did this
properly via classify_variants(); figure 4 panel a did not, and I did not notice
the inconsistency.

## 2026-09-23 — fig4 panel b arcs: TWO MORE QUINN REQUIREMENTS (agent running)
1. "ensure junctions in the sashimi plots map to coverage height, not to the bottom
   of the figure". Panel b draws `y = apex*sin(pi t)^0.78`, which is ZERO at both
   ends, so every foot sits on the BASELINE. Must anchor at the coverage height at
   each junction's own coordinate, as sashimi_track_gg's peak_y() and figure 3
   panel c's arc_geom() do. CHECKED: fig4 panel e and fig3 panel c already anchor
   correctly; only panel b is wrong.
2. "align junctions ... they can go under other junctions, or flip completely".
   Once feet are at coverage height, a tier whose apex is below max(y1,y2) bows
   DOWNWARD. Follow figure 3 panel c: work in peak space, back-solve the amplitude
   with uniroot, order tiers by ENDPOINT HEIGHT not span, and assert the signed
   minimum gap between every arc pair is positive on every build.

## OPEN QUESTIONS FOR QUINN

(all three previously-open questions were answered on 2026-09-22 21:00 — see
"RESOLVED THIS PASS" above)

## KNOWN DEBT

* fig 4 panels b and e DUPLICATE their arc/label machinery (the agent could not
  touch ../scripts/coverage_plot.R). They must be kept in step by hand, or the
  shared solver should be lifted into coverage_plot.R. Documented in both build.md.
* fig 3 panel a has a residual ~6-label knot (see above).
* `zoom_callout_gg()` / `zoom_source_box()` are still absent from the shared-toolkit
  list in figure2/compose.md.
   fig 4 b/e must be checked for a STALE "High GSP #13" track label — this is a
   GSP #30 figure and the ini says High_k30.
2. Final pass per figure: `scripts/figure_qa.R` for layout, `scripts/title_align.R`
   for per-row title alignment, then a 400 dpi visual read of every panel.
   Figure 3's ROW_H changed (row 1 is now 2.95), so its page is 11.2 in — re-run
   figure_qa with `2.95,3.35,2.45,2.45`.
3. `zoom_callout_gg()` and `zoom_source_box()` were added to scripts/coverage_plot.R
   by figure 3 panel c's agent and are used by figure 4 b/e. They are undocumented
   in figure2/compose.md — add them to the shared-toolkit list when things settle.

## Shared toolkit (read before writing any panel)

    theme.R                  design tokens; theme_publication(), angled_x_labels()
    scripts/coverage_plot.R  bigwig_binned, sashimi_track_gg, gene_body_gg,
                             transcript_zoom_gg, trapezoid_gg,
                             zoom_callout_gg, zoom_source_box
    scripts/pie_chart.R      pie_labels_on_wedges(), pie_legend_grob()
                             (pie_with_leaders() is GONE -- numbers go on the
                              wedge, boxed and border-coloured, never on leaders)
    scripts/raster_utils.R   gg_from_pdf() etc.
    scripts/figure_qa.R      per-cell ink extents, dead bands, content-top drift
    scripts/title_align.R    title ink rows per row-cell, for cross-panel alignment
    figure2/compose.md       THE trap list. Read this first, every time.

## 2026-09-23 (later) — ROW 3 TITLES, PIE CHIPS, FIG 1/2 REGRESSION

Quinn: "remember how I specifically said align titles always and yet figure 3 d e
and f are not aligned? Also, a needs to be widened and the labels fixed, and b
tighter, fix text you can't read, and lines from extra TCGA samples underneath the
gene bodies also need to be a bit shorter and it looks weird if one is perfectly
parallel with page and other isn't, make them both diagonal"

### Row-3 title alignment — FIXED
`title_align.R` measured d 50 / e 50 / f 138 px (88 px apart). Cause: panel f is a
screenshot under `coord_fixed`, so its plot is ~2.5 in tall in a 3.05 in row and
`wrap_elements` CENTRES the 0.55 in of spare height — floating f's title below the
others. Fix is a BIG BOTTOM margin (`margin(4, 6, 32, 6, "pt")`), not a small top
one, which pushes the element back to the top of its cell.
Re-measured: title ink starts on row 1 in all three cells. d, e, f aligned.
Panel d's stale comment (which still said d 69 / e 6 / f 151) is corrected.

### Pie chips now carry their wedge's colour — BOTH FIGURES
Boxed percentages are pushed OUT along their mid-ray, so they no longer touch the
wedge they name. With several sub-6% wedges adjacent — fig 3 panel d stacks four
(5.06 / 1.60 / 4.95 / 0.32) in one corner — the fan was ambiguous and nothing said
which chip belonged to which wedge. `pie_labels_on_wedges()` now defaults
`box_colour = NULL`, meaning **each box's border is drawn in its own wedge's fill**
(linewidth 0.2 -> 0.6); the fill stays white so the number is still black-on-white.
This keeps the numbers ON the chart with NO leader lines, which Quinn has rejected
twice. Verified at 600 dpi on both figures.
Also: fig 3 panel d's `D_PIE_BOX` right edge 0.50 -> 0.465. The pie is
HEIGHT-limited under coord_equal so this does not shrink it — it slides the centre
left by half the trim (~35 px at 600 dpi), which was exactly what was needed to get
the "4.95%" chip off the left edge of the 650-count bar at pos = 5.

### Figure 1 + 2 regression — 0 differing px each
Re-rendered both from scratch and diffed at 200 dpi against the shipped PDFs.
NB **ImageMagick `compare` is not installed on this machine.** Diff in R instead:
    png::readPNG both, apply(abs(x-y), c(1,2), max), count > 1/255.

### Figure 4 — full-page read at 150 dpi + spot crops at 500-600 dpi: PASSES
panel a's zoom places 15 labels with 0 leaders and no collisions; panel b's orange
set is 10.7/11.2 (KDM4C) and 6.91/1.27 (TRIM38); panel e's arcs all leave the exon
tops with their numbers in breaks in their own arcs; panel f's three annotation
leaders are unambiguous. Nothing clipped, no unreadable type found.

### Still with agents
* fig 3 panel a — widen (COL_W$ab beyond 0.52) and re-tune `cdf_place_labels()`
  to cut the 26 leaders; restore the legend entry from 0.78*SIZE_BASE.
* fig 3 panel b — tighten into the narrower cell; fix unreadable type; make BOTH
  "outside the top 60" stubs diagonal and shorter (one is currently axis-parallel
  because the stub logic draws vertical when the circle's x falls inside the
  caption's x-span).

## 2026-09-23 afternoon — Quinn's second batch, DONE (verified by me at 400-600 dpi)
* 3a zoom wider (COL_W$ab 0.58/0.42), gutters 48/20 + lead_max_pt 40 -> 17 leaders;
  callout lines from the CDF box's right corners to the zoom's corners (A_CALL_TOP/BOT).
* 3b PRPF8 frac c(0.22,0.24,0.40,0.14) so NTD/HB/RT/Thumb/X are readable; B_ANGLE 68
  (re-swept at this width; S 1.074 -> 1.064 after title change); leader-to-caption pad.
* Row-3 titles d/e/f 22/21/21 px. MY EARLIER "ALIGNED" WAS FALSE: the title_align strip
  started below the titles. Always crop from the row's TOP EDGE.
* Pies: 2-D chip placement, chips never touch the disc, border in own wedge colour.
  New helper args: ang_tol_deg, x_max, start_deg. fig3 d rotated (start_deg) so the four
  small wedges centre on 12 o'clock; D_PIE_R_IN 0.499 (box-limited, stable).
* Sashimi y-axis (line + 2 ticks, 0 centred on baseline) in coverage_plot.R, fig3 c,
  fig4 b and fig4 e (b/e draw their own; shared one switched off there). fig4 e got a
  4.5 pt high/control gutter (AX_BOT_PT) so the centred "0" no longer prints over the
  control track's y_max.
* 3c: each zoom leads from its own track's source box.
* 4b: TRIM38 orange arcs (1.27, 6.91) bow BELOW the baseline; crossings 4 -> 2.
* 4f heatmap: geom_tile (vector) instead of geom_raster (viewers smoothed the bitmap).
* 2b: inset smaller, no box (Quinn vetoed it earlier), top tier apex 0.60.
* 1a/2a: Illustrator art re-fonted to Nimbus via scripts/pdf_font_subst.R.
* 2c heatmap +10.5% taller.
* GENE ITALICS everywhere: theme.R gene_md() (markdown, for titles/axis/richtext) and
  gene_pm() (plotmath, for LEGENDS -- gridtext legend labels measured in one font and
  drew in another, so words overprinted). GENE_SYMBOLS list lives in theme.R.

## Open / flagged to Quinn
* 2c bottom column labels overlap at any angle (~45 columns in 2.65 in).
* 4c header says "NKAP/PRP8 introns" -- PRP8 is the yeast name; probably PRPF8.
* 1c BDH1 coverage tracks (coverage_track_gg) have no y-axis; not sashimi, left as is.

## 2026-09-23 — PIES NOW MATCH THE PUBLISHED STYLE
Quinn: make every pie look like Figure3-NKAP_improved.pdf panel d. That means
chips FILLED with the wedge colour (rounded, no border) sitting ON the pie; big
wedges plain text. pie_chart.R: on_disc = TRUE (default), r_pref 0.80, rim_cost
penalises covering any sliver's outer rim, gap_pt 2.4, ang_tol 60. fig3 d's
start_deg rotation removed (wedges start at 12 o'clock as published).
D_PIE_R_IN re-converged: fig3 d 0.499, fig4 d 0.580.
