# figure4/panel_b/build.R -- "Top junctions of GSP #30 capture usage of
# deep intronic splice sites".
#
# Two loci stacked vertically (KDM4C above, TRIM38 below). Every locus
# is a patchwork column drawn ENTIRELY from raw bigwigs + MeanPSI
# .links with ../scripts/coverage_plot.R -- no pygenometracks, no
# ggbio, no re-used *_MANE_gene_model.pdf:
#
#   High GSP #30 coverage + arcs      sashimi_track_gg() + own arc layer
#   Matched Controls coverage + arcs  sashimi_track_gg() + own arc layer
#   lavender taper + neck             taper_gg()  (local, see below)
#   full MANE gene model              transcript_zoom_gg()
#   genomic coordinate labels         bespoke row (same x scale)
#
# plus one legend row at the bottom for the orange "unannotated
# junction driving GSP #30" arcs.
#
# See panel_b/build.md for the design notes and traps.

if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
panel_b <- panel_cached("b", local({

PY <- file.path(QROOT, "coverage_tracks/fig4b_tracks")

# Sampled off the reference render: the coverage red is exactly
# #C0392B (== PAL$sf3b1_mut["SF3B1 hotspot mutation"]) and the
# unannotated-junction highlight is #F1A340 (ColorBrewer PuOr).
RED    <- "#C0392B"
GREY   <- "#969696"
ORANGE <- "#F1A340"
LAV    <- "#ccccff"           # same taper colour as figure 3 panel c

# ------------------------------------------------------------------
# Cohort labels. The PUBLISHED panel reads "High GSP #13 (15)" /
# "Matched Controls (150)" -- BOTH halves are stale here. This is the
# YLPM1 figure, the bigwigs are High_k30-*, and GSP_pygenome's ini
# titles are "High_k30 (13)" and "Matched_Controls (130)". (Panel e's
# pygenome run really does have 15/150; panel b's numbers were copied
# from it.) Corrected on both counts.
# ------------------------------------------------------------------
HIGH_LAB <- "High GSP"
CTRL_LAB <- "Matched\nControls"
HIGH_N   <- "#30 (13)"
CTRL_N   <- "(130)"

BW_HIGH <- "High_k30-3-..bw"            # High_k65-* is a stale factor key
BW_CTRL <- "Matched_Controls-3-..bw"
LK_HIGH <- "MeanPSI.High_k30-3-..links"
LK_CTRL <- "MeanPSI.Matched_Controls-3-..links"

# ------------------------------------------------------------------
# Loci. `wide` is the pyGenomeTracks --region used to make these
# bigwigs (flashier/poisson_NMF_junction_plots.Rmd, prep_pygenometrack
# calls at the bottom of the file); the GSP_pygenome inis do not record
# it. y_max is that region's own ini max_value, shared by both tracks
# so the two cohorts are comparable.
#
# The bigwigs contain data for EXACTLY `wide` and nothing outside it
# (checked with rtracklayer::import: chr6:25,973,078-25,984,559 and
# chr9:7,011,427-7,014,212). TRIM38's first exon block therefore
# starts one base inside the region and, drawn hard against the panel
# edge, read as a coverage peak sliced off by the frame. PAD_FRAC adds
# a small empty gutter at both ends of the DRAWN x range -- there is no
# data to un-clip, so the only honest fix is to let the peak sit inside
# the frame instead of on it. Everything that has to stay aligned (the
# taper's bottom corners, the gene model's highlight) uses the same
# padded window.
# ------------------------------------------------------------------
PAD_FRAC <- 0.014

LOCI <- list(
  KDM4C  = list(gene = "KDM4C",  chrom = "chr9",
                wide = c(7011426, 7014212),
                y_max = 6.10458374023438,  label_side = "left"),
  TRIM38 = list(gene = "TRIM38", chrom = "chr6",
                wide = c(25973077, 25984559),
                y_max = 45.8690692901611,  label_side = "left")
)
for (k in names(LOCI)) {
  p <- PAD_FRAC * diff(LOCI[[k]]$wide)
  LOCI[[k]]$view <- LOCI[[k]]$wide + c(-p, p)
}

# ------------------------------------------------------------------
# Row heights. patchwork heights are RELATIVE, so only the ratios here
# are stable -- see ../figure2/compose.md. Everything that needs a
# DRAWN size (the PSI-label solver) is derived from these, never
# restated.
#
# H_TRAP / H_GENE are bigger than they were: the italic gene name is
# ~0.076 in tall and the gene row is only ~0.12 in, so the name ALWAYS
# spills upward out of its own row into the taper row. That spill is
# what put "TRIM38" and "KDM4C" on top of the lavender. The fix is the
# NECK (see taper_gg): below NECK_F the taper is already down to the
# window width, so the bottom of the taper row is empty everywhere
# outside the window -- which is exactly where the gene name sits.
# The budget is checked below and warned about if it stops working.
# ------------------------------------------------------------------
H_TRK  <- 0.40     # one coverage/sashimi track
H_TRAP <- 0.155    # lavender taper
H_GENE <- 0.145    # MANE gene model
H_AXIS <- 0.100    # genomic coordinate labels
# H_LEG was 0.115, which worked out to 0.097 in of drawn row while the
# caption needs 0.088 in of ascender-to-descender -- so with the text
# centred at 0.33 of the row the descenders of "junction driving" fell
# out of panel b's cell and were CLIPPED by the row below. 0.16 gives
# the line 0.132 in and lets it sit centred, at a cost of 1.8% of the
# track height.
H_LEG  <- 0.16     # orange-arc legend line
STACK_H <- 2 * H_TRK + H_TRAP + H_GENE + H_AXIS

# ------------------------------------------------------------------
# DROP_F -- the band BELOW the coverage baseline that TRIM38's high
# track carries so its two orange arcs can bow DOWNWARD.
#
# Quinn, verbatim: "in 4b the yellow junctions in TRIM38 should curve
# underneath the gene body, its too crowded on top - you can lower the
# height of the 94.2 curve on the controls to allow space for the
# yellow junctions". TRIM38's high track was a seven-arc fan with four
# forced crossings; 1.27 and 6.91 are now solved in a mirrored world
# below the baseline (see arc_frame's `dir`), which leaves five black
# arcs up top.
#
# The band is bought as EXTRA ROW HEIGHT, not by squashing the
# coverage: the high row is H_TRK * (1 + DROP_F) tall and its y range
# is (-DROP_F * y_max, y_max), so the coverage part of every track in
# the panel is still exactly TRACK_H_IN. The two loci therefore stay
# directly comparable, which they would not be if TRIM38's coverage
# were drawn 30% shorter than KDM4C's. The whole stack pays for it
# (everything is ~6% smaller); B_SCALE below carries the correction, so
# it cannot go stale.
#
# The arcs CANNOT be drawn into the control track's own row -- these
# are separate patchwork cells and a cell clips its neighbours' ink.
#
# Quinn also offered "you can lower the height of the 94.2 curve on the
# controls to allow space" -- and it CANNOT be lowered, which is worth
# recording so nobody tries again. 94.2 contains 1.38 and SHARES ITS
# RIGHT FOOT with it, so at that end the only thing keeping 1.38 under
# it is how steeply 94.2 leaves the foot, and that steepness is the
# amplitude, i.e. the apex. 0.702 of y_max is already the SMALLEST
# apex that clears 1.38 by ARC_GAP; capping it at 0.50 (tried) just
# turns a solved pair into a crossing, which is the defect Quinn
# rejected the previous version for. The room comes out of the taller
# high row instead, which costs nothing but page.
# ------------------------------------------------------------------
DROP_F    <- 0.52     # of y_max, below the baseline, TRIM38 high only

# ROW_H[["ab"]] is panel b's whole row; the two-line title band eats ~0.385 in
# of it. Scale the nominal stack above onto what is actually left so the
# PSI-label solver's point-size -> data-unit conversion is right.
# MEASURED off the render, not guessed. patchwork heights are relative,
# so this number changes NOTHING about the layout -- its only job is to
# make TRACK_H_IN equal the track's true drawn height, which is what
# the PSI solver converts point sizes through. At 400 dpi panel b's
# cell runs y = 30..1016 px and the first track's top edge is at 170,
# so the content band is 846 px = 2.115 in.
# With the "0" centred on the baseline (see the y-axis block below) it
# hangs half a text height BELOW its own track. Under KDM4C's high
# track that page row belongs to the control track's y_max number, and
# the two printed on top of each other. AX_BOT_PT buys that one
# boundary a gutter; it is taken off CONTENT_H FIRST and handed back as
# extra row height, so the drawn coverage height is identical in all
# four tracks -- the two cohorts of a locus share a y_max and must be
# read against each other. The other three boundaries need nothing:
# TRIM38's high track has the whole drop band under its baseline, and
# a control track's "0" hangs over the top of the lavender taper,
# which is a tint, not ink.
AX_BOT_PT  <- 3.2
B_TITLE_H  <- 0.385
CONTENT_H  <- ROW_H[["ab"]] - B_TITLE_H
B_SCALE    <- (CONTENT_H - AX_BOT_PT / 72) /
              (2 * STACK_H + H_LEG + H_TRK * DROP_F)
H_AXGAP    <- (AX_BOT_PT / 72) / B_SCALE      # that gutter, in row units
TRACK_H_IN <- H_TRK  * B_SCALE
TRAP_H_IN  <- H_TRAP * B_SCALE
GENE_H_IN  <- H_GENE * B_SCALE

# Gutters, in pt. The left gutter holds the y-axis extremes ("45.9"),
# the right one the cohort captions; both are drawn with clip = "off"
# so they live in the margin rather than stealing track width.
PANEL_W <- PAGE_W * COL_W$ab[2] / sum(COL_W$ab)
MARG    <- c(20, 42)                       # c(left, right) in pt
gutter  <- function(bottom = 0) margin(0, MARG[2], bottom, MARG[1], "pt")
TRACK_W_IN <- PANEL_W - sum(MARG) / 72

SZ_AXIS   <- SIZE_BASE * 0.78 / .pt   # "45.9" / "0" and coordinate labels
SZ_COHORT <- SIZE_BASE * 0.88 / .pt
SZ_PSI    <- SIZE_BASE * 0.78 / .pt
SZ_GENE   <- SIZE_BASE * 0.95 / .pt
SZ_LEG    <- SIZE_BASE * 0.85 / .pt

# Cap height + descender of a string, in inches, at a ggplot `size`.
txt_h_in <- function(sz) sz * .pt * 0.72 / 72

# ------------------------------------------------------------------
# The coverage y axis. Quinn, for EVERY sashimi in the paper: "the
# numerical axis should have a line indicating it is an axis with two
# ticks for the numbers, and they are currently misaligned, the 0 sits
# entirely above the bottom of the track, it should sit beneath/even
# with it."
#
# So: a rule spanning 0..y_max (the COVERAGE range -- never the drop
# band below it, which is not part of the scale), a tick at each end,
# and the numbers hung off the ticks. The 0 is CENTRED on its tick, so
# it straddles the baseline instead of floating above it; the top
# number is top-aligned at its tick, because centring it would push
# half the glyphs out of the track and into the cell above, where
# patchwork clips them. Weight matches the gene model's coordinate
# ticks (grey25, 0.3) -- the same thin chrome used everywhere else in
# the figure.
#
# All four numbers live in the left gutter, which is OUTSIDE the drawn
# x window, so the x SCALE has to keep them (a scale limit drops data
# before anything is drawn -- ../figure2/compose.md) and `coord` sets
# the view. AX_GUT is the gutter the scale keeps, in pt.
# ------------------------------------------------------------------
# Chrome matched to `sashimi_track_gg`'s own y axis (figures 1-3):
# grey15, 0.3 pt rule, 2.5 pt ticks, 1.5 pt of air before the number,
# and the number at the SAME size (SZ_AXIS == SZ_PSI here, which is
# that function's default).
AX_COL  <- "grey15"
AX_LW   <- 0.3
AX_OFF  <- 0      # pt, axis rule left of the panel's left edge
AX_TICK <- 2.5    # pt, tick length, leftward from the rule
AX_PAD  <- 1.5    # pt, tick end -> right edge of the number
AX_GUT  <- 20     # pt of left gutter the x scale must keep

# ------------------------------------------------------------------
# .links helpers. sashimi_track_gg filters s1 >= start & e2 <= end;
# reproduce that here so the arcs we reason about are exactly the arcs
# that get drawn. The junction window is `wide`, NOT `view`: the pad
# exists to un-crowd the frame, not to pull a new junction (TRIM38's
# 63.3, a 96 bp intron sitting one base outside) into the figure.
# ------------------------------------------------------------------
read_links <- function(path, start, end, min_psi = 1) {
  lk <- read.table(path, header = FALSE, sep = "\t",
                   col.names = c("c1","s1","e1","c2","s2","e2","psi"))
  lk[lk$s1 >= start & lk$e2 <= end & lk$psi >= min_psi, , drop = FALSE]
}

# PSI exactly as the reference prints it: the .links value at its own
# precision, but never bare-integer ("2" -> "2.0", "99" -> "99.0").
# format() on the whole VECTOR pads to a common width, so each value is
# formatted on its own.
psi_fmt <- function(p) {
  s <- vapply(p, function(x) sprintf("%g", x), character(1))
  ifelse(grepl(".", s, fixed = TRUE), s, paste0(s, ".0"))
}

# ------------------------------------------------------------------
# MANE exons + every annotated splice site of the gene, cached (the
# GENCODE GRanges is 26 MB). `ends`/`starts` are the annotated donor /
# acceptor coordinates over ALL transcripts of the gene, in the
# .links convention: a junction (s1, e2) is annotated iff s1 is an exon
# END and e2 + 1 is an exon START.
# ------------------------------------------------------------------
EXON_RDS <- file.path(CACHE, "panel_b_exons.rds")
if (!file.exists(EXON_RDS)) {
  cat("[panel b] building", EXON_RDS, "(hg38_granges.RDS is 26 MB) ...\n")
  gencode <- readRDS(file.path(QROOT, "inputs/gencode_v46_granges.RDS"))
  mane <- data.table::fread(
    file.path(QROOT, "inputs/MANE.GRCh38.v1.4.summary.txt"),
    data.table = FALSE)
  ex_list <- lapply(names(LOCI), function(g) {
    ex <- gencode[gencode$gene_name == g & gencode$type == "exon"]
    tx <- mane$Ensembl_nuc[mane$symbol == g]
    me <- ex[ex$transcript_id %in% tx]
    stopifnot(length(me) > 0)
    list(mane   = me,
         ends   = sort(unique(GenomicRanges::end(ex))),
         starts = sort(unique(GenomicRanges::start(ex))))
  })
  names(ex_list) <- names(LOCI)
  dir.create(CACHE, showWarnings = FALSE)
  saveRDS(ex_list, EXON_RDS)
}
EXONS <- readRDS(EXON_RDS)

# ------------------------------------------------------------------
# Which junctions are "driving GSP #30" (the orange ones)?
#
# Quinn's rule, verbatim: "TRIM38 orange junctions should be based upon
# the TRIM38 junctions that have high intensities in the TRIM38 poisson
# NMF junction group in the top 100 junction groups that make up GSP
# #30". So the orange set is READ OUT OF THE FACTORISATION, not
# guessed from the geometry of the arcs.
#
# The chain, re-derived here and cached into
# figure4/cache/panel_b_gsp30_juncs.tsv:
#
#  1. GSP #30's top junction groups. `flashier/fact_30_top_junctions.tsv`
#     IS that list -- spec_flash_fact.Rmd's "top 250 GSPs 10 13 30"
#     chunk writes one row per feature of
#     rownames(norm_F)[order(norm_F[, 30], decreasing = TRUE)][1:250],
#     IN RANK ORDER, where norm_F = t(diag(ldf$D) %*% t(ldf$F)) of the
#     3.7 GB softimpute_flash_300_qc_filtered.RDS. That loop can `next`
#     past a feature, which would silently shift every later row's
#     rank -- but the file has all 250 rows, so nothing was skipped and
#     row i is rank i. The check is asserted below, not assumed.
#     KDM4C is rank 40, TRIM38 rank 55. TRIM38 has a SECOND junction
#     group at rank 114, OUTSIDE the top 100 -- that one owns the 8.38
#     junction, which is exactly why 8.38 must not be orange.
#
#  2. Which gene-level factor `kN` each group is. The tsv drops the
#     `.kN` suffix but keeps the group's argmax UNANNOTATED junction,
#     and that identifies the factor uniquely: recompute the
#     per-factor argmax over the gene's unannotated junctions and take
#     the column that reproduces the recorded junction.
#
#  3. Which junctions to colour: poisson_NMF_junction_plots.Rmd -- the
#     Rmd these very tracks came from -- highlights every unannotated
#     junction of the group within 2x of the group's top loading
#     (`loadings[loadings >= top_val / 2]` on poisson2multinom(res)$F).
#     Same rule here, so the panel and the Rmd cannot drift apart.
#
# `res$F` rownames are 1-BASED; .links is BED, so s1 = start - 1 and
# e2 = end. What it yields:
#
#   KDM4C  k31, rank 40:  11.2 (0.1121) and 10.7 (0.1013) -- the
#     group's own top is 0.1217 at chr9:7,172,713-7,174,552, 160 kb
#     outside the drawn window, so both in-window junctions clear
#     top/2 = 0.0609 comfortably.
#   TRIM38 k21, rank 55:  1.27 (0.4254, the group's top) and 6.91
#     (0.2841). NOT 22.7 (k24) and NOT 8.38 (k22, rank 114) -- the two
#     the old smallest-intronic-gap heuristic got wrong.
# ------------------------------------------------------------------
GSP_TSV  <- file.path(CACHE, "panel_b_gsp30_juncs.tsv")
GSP_TOPN <- 100L
GSP_ROWS_EXPECTED <- 250L
TOPJ_TSV <- file.path(QROOT, "gsp_analysis/fact_30_top_junctions.tsv")
PNMF_DIR <- file.path(QROOT, "quetzal_runs/genome_wide")
SNAP_DIR <- file.path(QROOT, "inputs/snaptron/tcga")

build_gsp_tsv <- function() {
  cat("[panel b] building", GSP_TSV, "from the GSP #30 top-group list ...\n")
  top <- data.table::fread(TOPJ_TSV, data.table = FALSE,
                           select = c("gene_name", "chromosome", "start", "end"))
  if (nrow(top) != GSP_ROWS_EXPECTED)
    stop(sprintf(paste("[panel b] %s has %d rows, not %d: features were",
                       "skipped, so a row index is no longer a rank and the",
                       "cut at %d is unsafe."),
                 TOPJ_TSV, nrow(top), GSP_ROWS_EXPECTED, GSP_TOPN))
  top      <- top[seq_len(GSP_TOPN), , drop = FALSE]
  top$rank <- seq_len(GSP_TOPN)
  suppressPackageStartupMessages(library(fastTopics))

  out <- list()
  for (g in names(LOCI)) {
    rows <- top[top$gene_name == g, , drop = FALSE]
    if (nrow(rows) == 0)
      stop("[panel b] ", g, " has no junction group in GSP #30's top ",
           GSP_TOPN)
    chr  <- rows$chromosome[1]
    res  <- fastTopics::poisson2multinom(
              readRDS(file.path(PNMF_DIR, chr, "FastTopics_output", g,
                                "res.RDS")))
    snap <- data.table::fread(
              file.path(SNAP_DIR, chr, "snaptron_output",
                        paste0(g, "_snaptron.tsv")), data.table = FALSE)
    unn  <- paste0(snap$start[snap$annotated == 0], "-",
                   snap$end  [snap$annotated == 0])
    cand <- rownames(res$F) %in% unn
    stopifnot(any(cand))
    Fu   <- res$F[cand, , drop = FALSE]
    amax <- rownames(Fu)[apply(Fu, 2, which.max)]

    for (r in seq_len(nrow(rows))) {
      want <- paste0(rows$start[r], "-", rows$end[r])
      kk   <- which(amax == want)
      if (length(kk) == 0)
        stop("[panel b] ", g, ": no gene factor has ", want,
             " as its top unannotated junction")
      if (length(kk) > 1) kk <- kk[which.max(Fu[want, kk])]
      ld  <- Fu[, kk]
      sel <- which(ld >= max(ld) / 2)
      jc  <- do.call(rbind, strsplit(names(ld)[sel], "-", fixed = TRUE))
      out[[length(out) + 1L]] <- data.frame(
        gene      = g,
        chrom     = chr,
        group     = colnames(res$F)[kk],
        gsp_rank  = rows$rank[r],
        s1        = as.numeric(jc[, 1]) - 1,   # 1-based -> BED/.links
        e2        = as.numeric(jc[, 2]),
        loading   = as.numeric(ld[sel]),
        group_top = max(ld),
        orange    = TRUE,
        stringsAsFactors = FALSE)
    }
  }
  res_df <- do.call(rbind, out)
  res_df <- res_df[order(res_df$gene, -res_df$loading), ]
  dir.create(CACHE, showWarnings = FALSE)
  write.table(res_df, GSP_TSV, sep = "\t", quote = FALSE, row.names = FALSE)
  res_df
}

if (!file.exists(GSP_TSV)) build_gsp_tsv()
GSP30 <- read.table(GSP_TSV, header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE)

orange_keys <- function(gene, hi) {
  g <- GSP30[GSP30$gene == gene & GSP30$orange, , drop = FALSE]
  intersect(paste(hi$s1, hi$e2), paste(g$s1, g$e2))
}

# ------------------------------------------------------------------
# Arc geometry.  (Rewritten; see build.md "Arcs" for the full story.)
#
# EACH FOOT SITS ON THE COVERAGE AT ITS OWN JUNCTION COORDINATE. The
# old code drew y = apex * sin(pi t)^p, which is ZERO at t = 0 and
# t = 1, so every foot landed on the track BASELINE no matter where the
# coverage was -- Quinn's requirement #1, and the root of the "arcs
# dangle off the bottom of the figure" complaint. The curve is now the
# same parametric form `sashimi_track_gg()` and figure 3 panel c use:
#
#     y(t) = (1 - t) y1 + t y2 + sin(t pi) * amp
#
# with y1 / y2 the coverage AT s1 and AT e2. X is untouched: every foot
# is still drawn at the raw .links coordinate, so the arcs land exactly
# on their splice sites (re-verified at 800 dpi, see build.md).
#
# ANCHOR_BP = 0, as in figure 3 panel c: the foot takes the coverage of
# the bin the junction coordinate falls in, which is a bin that IS
# drawn, so the foot sits ON the coverage silhouette and never floats
# above it. Widening the window walks up into the exon plateau (KDM4C's
# 66.2 right foot: 0.13 of y_max at the site, 0.66 fifty bp in) and the
# foot would then hang a tenth of an inch over the fill.
#
# Heights are solved in PEAK space (the drawn maximum of the curve),
# never as an absolute fraction of y_max: with the feet at coverage
# height an absolute apex is routinely BELOW its own feet, the sine
# bump goes negative and the arc is drawn drooping to the baseline and
# back -- Quinn's "flip completely". Every arc here is guaranteed
# peak >= max(y1, y2) + ARC_RISE by construction: the amplitude is
# never solved DOWN from that floor, only up.
#
# The sine AMPLITUDE that puts the DRAWN maximum at a chosen peak is
# back-solved with uniroot: with a steep chord the maximum is a long
# way from t = 0.5 (KDM4C's 10.7 runs 0.53 -> 0.01 and peaks at
# t = 0.34), so assuming the midpoint sags the arc.
#
# ---- Why a plain tier ladder is not enough -------------------------
# Sorting the arcs and handing tier k the peak (k-1) * ARC_SEP above
# tier 1 is what figure 3 panel c does, and it works there because no
# quadrant has more than two arcs. It FAILS here, and not marginally:
# with the feet at coverage height an arc is LOW near its own feet, so
# a short arc nested against the end of a long one pops out through it
# even though its peak is far lower. TRIM38's controls are the clean
# example -- 1.38 and 94.2 share the right foot at 0.20 of y_max, and
# over 1.38's 1266 bp the 9878 bp arc only climbs to 0.27, while 1.38's
# own peak was 0.30. It crossed by 0.074 of y_max.
#
# So the heights are SOLVED, arc by arc, bottom-up in the containment
# order. Arc i's amplitude is the smallest that satisfies all of:
#   * peak_i >= max(y1, y2) + ARC_RISE            (it reads as an arc)
#   * peak_i >= peak_j + ARC_SEP                  for every already
#     placed j that CONTAINS i or is contained by it (the tiers stay
#     visibly apart where they actually run together)
#   * min over their shared x of (y_i - y_j) >= ARC_GAP, same j
#     (dominance: i never passes under j, which is requirement #2)
# Each is monotone in the amplitude, so the binding one is found by
# bisection. ARC_CEIL caps the result; if the cap bites before
# dominance is reached, the pair is reported as a forced crossing
# rather than silently drawn wrong.
#
# ---- Which crossings are genuinely unavoidable ---------------------
# Two arcs whose intervals PARTIALLY overlap (neither nested nor
# disjoint) must cross, whatever their heights. TRIM38's 1.62
# (25,974,930-25,983,163) partially overlaps 1.27, 2.13, 8.38 and 6.91,
# and those four crossings are in the junctions, not in the drawing.
#
# A NESTED pair can also be forced, which is new with coverage-anchored
# feet: if the inner arc's far foot sits below the outer arc's curve
# there, and they share the near foot, the inner arc can only stay
# under the outer one if the outer one is steep enough at the shared
# foot. That is a lower bound on the outer arc's amplitude and it CAN
# exceed the track. When it does the solver leaves the constraint
# unsatisfied and says so, instead of spending the whole budget on an
# amplitude that still does not clear -- doing the latter is what put
# five of KDM4C's and TRIM38's arcs on the ceiling in the first draft
# and made everything above them cross.
#
# Disjoint arcs never interact and are given no constraints at all --
# which is what keeps the budget affordable. KDM4C's 11.2 and 10.7 are
# disjoint, so 11.2 costs the four arcs above it nothing.
#
# ---- Why ARC_P = 0.5 -----------------------------------------------
# The bump is sin(t pi)^ARC_P, not a plain sine. With ARC_P = 1 the
# arc leaves each foot LINEARLY, so a long arc is nearly flat near its
# own feet -- and a short arc that shares that foot climbs straight out
# through it. Making the long arc taller does not help, because the
# sine's tail stays flat: dominating TRIM38's 22.7 over its last
# 1.3 kb needed amplitude 0.96, i.e. a peak of 1.09 of y_max, which is
# off the track. At ARC_P = 0.5 the bump's derivative at the feet
# diverges, the tails are steep, and the same constraint is met at
# amplitude 0.51. Every avoidable crossing in all four tracks then
# solves inside the ceiling.
#
# 0.5 is also the right shape on its own merits: sin(t pi)^0.5 is
# within 0.03 of the half-ellipse sqrt(4 t (1-t)) that pyGenomeTracks
# draws, everywhere. So the arcs read as sashimi arcs, land crisply ON
# their splice sites instead of melting into the coverage, and carry a
# flat top that the contour-style PSI labels can sit on.
# ------------------------------------------------------------------
ANCHOR_BP <- 0
ARC_P     <- 0.50     # bump shape, sin(t pi)^ARC_P (~ a half-ellipse)
ARC_RISE  <- 0.08     # an arc's own bulge above its higher foot
ARC_SEP   <- 0.10     # minimum peak-to-peak separation, nested pairs
ARC_GAP   <- 0.03     # minimum drawn clearance, nested pairs
ARC_TRIM  <- 0.10     # of the shared span ignored at each end (see below)

# ARC_CEIL is DERIVED, not chosen: it is exactly `place_psi()`'s own
# ceiling for a label CENTRE, so the top arc's apex is always a legal
# place for that arc's number. This matters more than it sounds. When
# the ceiling was a round 0.93 the top arc spent most of its length
# ABOVE the label ceiling, its number could only sit near a foot, and
# at KDM4C that is inside the bundle where four arcs leave one splice
# site together -- the break then had two lines entering and two
# leaving and named neither. Both values now come from the same
# expression, so a change to SZ_PSI or the row heights moves them
# together.
lab_ceil_frac <- function(trk_h_in)
  1 - (txt_h_in(SZ_PSI) / 2 + 0.45 / 72) / trk_h_in
ARC_CEIL  <- lab_ceil_frac(TRACK_H_IN)

# ---- the mirrored world below the baseline -------------------------
# A downward arc is solved by MIRRORING: the solver is handed -y1/-y2
# and its own ceiling, so "peak" means depth, "above" means deeper, and
# every constraint (ARC_RISE, ARC_SEP, ARC_GAP, dominance) keeps its
# meaning without a second implementation. `dir` on the returned frame
# maps the solution back at draw time.
#
# DOWN_CEIL is the mirror of ARC_CEIL: the deepest an arc may go and
# still leave its own number a legal seat inside the drop band.
# DOWN_MIN is a FLOOR on the shallowest one -- a down arc has to clear
# the baseline by a whole label half-height before its number can sit
# in a break in it, which ARC_RISE (0.08, measured from the arc's own
# feet, and the right foot of both of these is nearly on the baseline)
# does not on its own guarantee.
DOWN_CEIL <- DROP_F - (1 - ARC_CEIL)
# 2.4 pt, not the 0.45 pt ARC_CEIL uses, because a label has to sit in
# a break in the arc over its whole WIDTH, and with ARC_P = 0.5 the arc
# climbs away from its trough fast. At 1.0 pt TRIM38's 1.27 bottomed
# out 0.014 of y_max shy of its own label's seat and the number was
# pushed off the curve (off = -0.172 in the build log).
DOWN_MIN  <- (txt_h_in(SZ_PSI) / 2 + 2.4 / 72) / TRACK_H_IN

arc_bump <- function(t) sin(t * pi)^ARC_P

# Drawn maximum of (1-t) y1 + t y2 + bump(t) * amp over t in [0, 1].
# The grid is fine because with ARC_P < 1 the maximum can sit close to
# a foot when the chord is steep.
arc_peak_val <- function(y1, y2, amp, tt = seq(0, 1, length.out = 1201))
  max((1 - tt) * y1 + tt * y2 + arc_bump(tt) * amp)

# Amplitude whose drawn maximum is exactly `p` (0 if the chord already
# reaches it, which only happens for a degenerate p <= max(y1, y2)).
amp_for_peak <- function(y1, y2, p) {
  f <- function(a) arc_peak_val(y1, y2, a) - p
  if (f(0) >= 0) 0 else uniroot(f, c(0, 6), tol = 1e-9)$root
}

# 600 points, not 320: with ARC_P = 0.5 the limbs turn hard within a
# few percent of t and a coarse grid visibly chamfers the feet.
# NB y1 / y2 / amp / peak on an arc frame are in the SOLVER's space,
# which for a downward arc is mirrored. `dir` (+1 up, -1 down) maps
# back to the page: every drawn y goes through it.
arc_curve <- function(a, n = 600) {
  t <- seq(0, 1, length.out = n)
  data.frame(x = a$x1 + t * (a$x2 - a$x1),
             y = a$dir *
                 ((1 - t) * a$y1 + t * a$y2 + arc_bump(t) * a$amp) * a$y_max)
}
arc_y_at <- function(a, x) {
  t  <- (x - a$x1) / (a$x2 - a$x1)
  tc <- pmin(pmax(t, 0), 1)
  y  <- a$dir *
        ((1 - tc) * a$y1 + tc * a$y2 + arc_bump(tc) * a$amp) * a$y_max
  ifelse(t < 0 | t > 1, NA_real_, y)
}
# y in FRACTION-of-y_max units, for the solver (which works there).
arc_frac_at <- function(y1, y2, amp, x1, x2, x) {
  t <- (x - x1) / (x2 - x1)
  (1 - t) * y1 + t * y2 + arc_bump(t) * amp
}

# The x at which two arcs are compared. The outer ARC_TRIM of the span
# they share is dropped: these junction sets share endpoints -- four of
# KDM4C's five start at 7,011,879 and three of TRIM38's seven end at
# 25,983,163 -- and at a shared foot the gap is 0 by construction, so
# including it would make every dominance test unsatisfiable.
shared_xs <- function(a1, a2, b1, b2, n = 240) {
  lo <- max(a1, b1); hi <- min(a2, b2)
  if (hi <= lo) return(numeric(0))
  m <- ARC_TRIM * (hi - lo)
  seq(lo + m, hi - m, length.out = n)
}

nests <- function(a1, a2, b1, b2)
  (a1 >= b1 && a2 <= b2) || (b1 >= a1 && b2 <= a2)

# Signed vertical gap, in fractions of y_max, between every pair of
# arcs that is NOT a forced crossing: (upper) - (lower) over the x they
# share. Negative anywhere means an arc passes UNDER one it is tiered
# above -- exactly the defect Quinn rejected. Printed on every build.
arc_min_gap <- function(arcs) {
  n <- nrow(arcs)
  out <- list(gap = Inf, worst = "", forced = character(0))
  if (n < 2) { out$gap <- NA_real_; return(out) }
  for (i in seq_len(n - 1)) for (j in (i + 1):n) {
    xs <- shared_xs(arcs$x1[i], arcs$x2[i], arcs$x1[j], arcs$x2[j])
    if (!length(xs)) next                                  # disjoint
    if (!nests(arcs$x1[i], arcs$x2[i], arcs$x1[j], arcs$x2[j])) {
      out$forced <- c(out$forced,
                      sprintf("%s/%s [partial overlap]", arcs$lab[i],
                              arcs$lab[j]))
      next
    }
    u <- if (arcs$peak[i] > arcs$peak[j]) i else j
    l <- if (arcs$peak[i] > arcs$peak[j]) j else i
    g <- min(arc_frac_at(arcs$y1[u], arcs$y2[u], arcs$amp[u],
                         arcs$x1[u], arcs$x2[u], xs) -
             arc_frac_at(arcs$y1[l], arcs$y2[l], arcs$amp[l],
                         arcs$x1[l], arcs$x2[l], xs))
    if (g <= 0)
      out$forced <- c(out$forced, sprintf("%s under %s [%+.3f, ceiling]",
                                          arcs$lab[u], arcs$lab[l], g))
    else if (g < out$gap) {
      out$gap <- g; out$worst <- sprintf("%s over %s", arcs$lab[u], arcs$lab[l])
    }
  }
  if (is.infinite(out$gap)) out$gap <- NA_real_
  out
}

# `dir = -1` solves the arcs UPSIDE DOWN: the feet keep their true
# coverage heights (ANCHOR_BP = 0, exactly as before -- only the bow
# direction changes) but everything the solver sees is negated, so it
# places them in the band below the baseline with the same
# rise/separation/dominance guarantees. `ceil` is that world's ceiling
# and `min_peak` a floor every arc must reach.
arc_frame <- function(lk, d, y_max, orange_key = character(0),
                      dir = 1, ceil = ARC_CEIL, min_peak = 0) {
  n <- nrow(lk)
  if (n == 0) return(NULL)
  lk <- lk[order(lk$e2 - lk$s1), , drop = FALSE]

  covat <- function(x) {
    w <- d$bin_center >= x - ANCHOR_BP & d$bin_center <= x + ANCHOR_BP
    if (!any(w)) d$coverage[which.min(abs(d$bin_center - x))]
    else         max(d$coverage[w])
  }
  y1   <- dir * vapply(lk$s1, covat, numeric(1)) / y_max
  y2   <- dir * vapply(lk$e2, covat, numeric(1)) / y_max
  E    <- pmax(y1, y2)
  span <- lk$e2 - lk$s1

  # nested_in[a, b]: arc a lies inside arc b, so a is solved (and sits)
  # BELOW it. Kahn's algorithm over that partial order, taking the
  # lowest (E, span) among the arcs whose nested children are all done.
  nested_in <- matrix(FALSE, n, n)
  for (a in seq_len(n)) for (b in seq_len(n))
    nested_in[a, b] <- a != b && lk$s1[a] >= lk$s1[b] &&
                       lk$e2[a] <= lk$e2[b] && span[a] < span[b]
  ord <- integer(0); rem <- rep(TRUE, n)
  for (k in seq_len(n)) {
    ready <- which(rem & vapply(seq_len(n),
                     function(b) !any(rem & nested_in[, b]), logical(1)))
    if (!length(ready)) ready <- which(rem)   # unreachable: strict order
    pick <- ready[order(E[ready], span[ready])][1]
    ord  <- c(ord, pick); rem[pick] <- FALSE
  }

  amp <- rep(NA_real_, n)
  for (k in seq_len(n)) {
    i    <- ord[k]
    cap  <- amp_for_peak(y1[i], y2[i], ceil)
    need <- max(amp_for_peak(y1[i], y2[i], E[i] + ARC_RISE),
                amp_for_peak(y1[i], y2[i], min_peak))
    for (j in ord[seq_len(k - 1)]) {
      if (!nests(lk$s1[i], lk$e2[i], lk$s1[j], lk$e2[j])) next
      xs <- shared_xs(lk$s1[i], lk$e2[i], lk$s1[j], lk$e2[j])
      if (!length(xs)) next
      # (a) stay ARC_SEP above j's peak
      need <- max(need, amp_for_peak(y1[i], y2[i],
                          arc_peak_val(y1[j], y2[j], amp[j]) + ARC_SEP))
      # (b) dominate j by ARC_GAP everywhere they share x. Monotone in
      #     the amplitude (sin(t pi) > 0 in the interior), so bisect.
      yj <- arc_frac_at(y1[j], y2[j], amp[j], lk$s1[j], lk$e2[j], xs)
      gp <- function(a) min(arc_frac_at(y1[i], y2[i], a,
                                        lk$s1[i], lk$e2[i], xs) - yj)
      if (gp(cap) >= ARC_GAP) {
        lo <- need; hi <- cap
        if (gp(lo) < ARC_GAP) {
          for (it in seq_len(60)) {
            mid <- (lo + hi) / 2
            if (gp(mid) >= ARC_GAP) hi <- mid else lo <- mid
          }
          need <- max(need, hi)
        }
      }
      # else: not reachable inside the ceiling. Leave `need` alone --
      # spending the budget on it would not clear the pair and would
      # push every arc above this one through the top of the track.
      # arc_min_gap() reports the pair as a crossing.
    }
    amp[i] <- min(need, cap)
  }
  peak <- vapply(seq_len(n),
                 function(i) arc_peak_val(y1[i], y2[i], amp[i]), numeric(1))

  key <- paste(lk$s1, lk$e2)
  data.frame(x1 = lk$s1, x2 = lk$e2, psi = lk$psi,
             y1 = y1, y2 = y2, amp = amp, peak = peak, y_max = y_max,
             dir = dir, lab = psi_fmt(lk$psi),
             col = ifelse(key %in% orange_key, ORANGE, "black"),
             stringsAsFactors = FALSE)
}

# ------------------------------------------------------------------
# PSI labels: CONTOUR-STYLE, i.e. the number sits in a break in the arc
# it names, exactly the way a contour map labels an isoline.
#
# Why not "just put the number in a clear gap between two arcs": it
# does not fit and cannot be made to. TRIM38's high track carries seven
# tiers across 0.48 * y_max of a 0.33 in track -- 0.026 in apart -- and
# the digits are 0.062 in tall, so a number needs a gap four times
# wider than any that exists. Even KDM4C, with five, is 0.085 of y_max
# between tiers against a label half-height of 0.111. The old solver demanded exactly
# that, failed, and walked its level ladder until it parked "1.27" and
# "2.13" in the empty upper-left corner, nowhere near the arcs they
# name. (A white glyph halo does not rescue it either: an outline only
# thickens the glyphs, so a near-horizontal arc still shows through the
# counters and the gaps between digits -- it renders as a strike
# through the number.)
#
# So the arcs are CUT instead. Every label gets a box; each arc's curve
# is split wherever it enters one, and only the surviving runs are
# drawn. Nothing white is painted, the PDF stays vector-clean, and the
# invariant Quinn asked for holds exactly: no line ever touches a
# glyph. The number lands ON the arc it belongs to, so attribution
# needs no leader and no guesswork.
#
# The placement search is scored, not first-fit. Hard constraints:
# inside the panel, under the ceiling, clear of the coverage area,
# and no overlap with a label already placed. Among what is left it
# minimises
#   CUT_W * (number of OTHER arcs the box would cut)
#     + 3.0 * (label pushed off its own arc)
#     + 3.0 * min(on-page slope of its own arc under the box, 2)
#     + 0.9 * |t - 0.5|
# The SLOPE term is what keeps the number readable as part of its arc:
# a break in a near-horizontal line reads as a labelled line, a break
# in a steep one reads as a line that stopped.
#
# CUT_W stays LOW (1.5) because the neighbour's break is now a hairline
# notch rather than a hole -- see bw2/bh2 in place_psi(). Raising it to
# 3.0 was tried and pushed TRIM38's 2.13 and 8.38 off their own arcs to
# chase positions that cut nothing; with a cheap notch there is nothing
# to chase.
#
# Labels are placed MOST-CONSTRAINED FIRST: in increasing order of
# (arc span / label width), i.e. the number that barely fits on its own
# arc chooses before the ones with a whole track to move along. TRIM38's
# "22.7" is 0.7 kb of label on a 1.3 kb arc that shares its right splice
# site with 52.4 and 1.62; placed last it found no legal slot at all and
# fell through to the fallback, and placed by tier order it landed on
# the other two arcs' landings. First, it takes the one slot it has and
# the roomy arcs work around it. (Tier order -- what this used to do --
# and peak order were both tried.)
# ------------------------------------------------------------------
LAB_TS    <- seq(0.10, 0.90, by = 0.02)
CUT_W     <- 1.5
FOOT_KEEP <- 5.0   # pt of arc that must survive leading into every foot

# Helvetica advance widths: digits are 0.556 em, the period 0.278 em.
# nchar() * one fudged average over-estimates "6.91" by ~25%, and the
# width is what sizes the cut in the arc, so it is worth being exact.
str_w_em <- function(s)
  vapply(strsplit(s, "", fixed = TRUE),
         function(cs) sum(ifelse(cs == ".", 0.278, 0.556)), numeric(1))

# `idx` picks WHICH arcs are placed -- the whole frame is still passed
# in, because the cut count and the protected landings must be scored
# against every arc in the track, not just the ones being placed. The
# up and down fans are placed in separate calls (`below = TRUE` for the
# mirrored ones): their boxes live in disjoint y bands, so neither can
# collide with the other's, and each gets the seat rules that fit it.
place_psi <- function(arcs, d, y_max, trk_h_in, trk_w_in, s, e,
                      idx = seq_len(nrow(arcs)), below = FALSE,
                      floor_y = 0) {
  psi_pt  <- SZ_PSI * .pt
  txt_h   <- txt_h_in(SZ_PSI)
  x_pt    <- (e - s)  / (trk_w_in * 72)     # data-x per pt
  y_pt    <- y_max    / (trk_h_in * 72)     # data-y per pt
  n       <- nrow(arcs)
  hw      <- str_w_em(arcs$lab) / 2 * psi_pt * x_pt
  hh      <- txt_h / 2 * y_max / trk_h_in
  # TWO break sizes, which is the whole trick for reading a crowded
  # track. The arc the number NAMES is cut generously (bw/bh), so the
  # number visibly sits in a gap in that line. Every OTHER arc the box
  # touches is cut as tightly as "no line may touch a glyph" allows
  # (bw2/bh2) -- a hairline notch that the eye reads straight through.
  # That is the contour-map convention: only the labelled isoline
  # breaks. Cutting every arc the same width is what made KDM4C's
  # "1.92" and "4.87" unreadable -- two lines stopped on the left of
  # each number and two started on the right, and nothing said which
  # one the number belonged to.
  bw      <- hw + 1.7 * x_pt
  bh      <- hh + 0.9 * y_pt
  bw2     <- hw + 0.5 * x_pt
  bh2     <- hh + 0.5 * y_pt
  # Extra air demanded BETWEEN two numbers, on top of their boxes:
  # boxes that merely fail to overlap still print as "6.91:8.38".
  air_x   <- 6.0 * x_pt
  air_y   <- 1.0 * y_pt
  # How far a label may be pushed off its own arc when it cannot sit on
  # it (a low arc whose apex is inside the coverage).
  off     <- hh + 1.6 * y_pt
  # Highest label CENTRE whose GLYPHS stay in the track. bh carries the
  # box's air as well as the glyphs, so testing `ly + bh > y_max - bh`
  # would charge the pad twice -- 0.07 of y_max, which with seven tiers
  # is most of a tier.
  ceil_y  <- lab_ceil_frac(trk_h_in) * y_max     # == ARC_CEIL * y_max

  # The protected run leading into each foot, as drawn points.
  prot <- lapply(seq_len(n), function(j) {
    cv  <- arc_curve(arcs[j, ], n = 200)
    dl  <- c(0, cumsum(sqrt((diff(cv$x) / x_pt)^2 + (diff(cv$y) / y_pt)^2)))
    tot <- dl[length(dl)]
    cv[dl < FOOT_KEEP | (tot - dl) < FOOT_KEEP, c("x", "y"), drop = FALSE]
  })

  px <- py <- rep(NA_real_, n)
  # Most-constrained first: increasing (arc span / label width). See the
  # header -- "22.7" has one slot and everything else has a hundred.
  for (i in idx[order(((arcs$x2 - arcs$x1) / (2 * bw))[idx])]) {
    a     <- arcs[i, ]
    done  <- which(!is.na(py))
    best  <- NULL; best_s <- Inf
    # 0 = on its own arc, which is what the contour idiom needs. The
    # half steps are there so a number that CANNOT sit on its arc
    # leaves it by as little as the constraints allow: at +-1 only,
    # TRIM38's 1.38 -- a 25 pt arc carrying a 15 pt number between two
    # 5 pt protected landings -- jumped a whole label-height clear of
    # its curve and read as an unattached number. Half a push still
    # breaks the arc it names, so attribution survives.
    for (sgn in c(0, 0.5, -0.5, 1, -1)) {
      for (tt in LAB_TS) {
        lx <- a$x1 + tt * (a$x2 - a$x1)
        if (lx - bw[i] < s || lx + bw[i] > e) next
        ly <- arc_y_at(a, lx) + sgn * off
        if (below) {
          # In the drop band: the whole box stays under the baseline
          # (a number half in the track reads as a coverage artefact)
          # and inside the band. There is no coverage down here, so
          # the clearance test above does not apply.
          if (ly + bh > -0.4 * y_pt || ly - bh < floor_y + 0.5 * y_pt) next
        } else {
          if (ly > ceil_y || ly - bh < 0) next
          inw  <- d$bin_center >= lx - bw[i] & d$bin_center <= lx + bw[i]
          covb <- if (any(inw)) max(d$coverage[inw]) else 0
          if (ly - bh < covb + 0.6 * y_pt) next
        }
        if (length(done) &&
            any(abs(px[done] - lx) < (bw[done] + bw[i] + air_x) &
                abs(py[done] - ly) < (bh + bh + air_y))) next
        # The break must never eat an arc's LANDING. This panel is ABOUT
        # where the junctions come down, and a number is wider than a
        # short arc: "22.7" spans 700 bp of TRIM38's window against its
        # own arc's 1266, so centred on that arc the box ate the arc
        # whole. Moved off the arc it then sat exactly over the last
        # 7 pt of 52.4 and 1.62, which share 22.7's right splice site --
        # both simply stopped in mid-air short of the exon. So the last
        # FOOT_KEEP points of every arc, at BOTH ends, are protected and
        # no label box may touch them.
        if (any(vapply(seq_len(n), function(j) {
              hwj <- if (j == i) bw[i] else bw2[i]
              hhj <- if (j == i) bh   else bh2
              any(abs(prot[[j]]$x - lx) < hwj & abs(prot[[j]]$y - ly) < hhj)
            }, logical(1)))) next
        # How many OTHER arcs would this box cut?
        probe <- seq(lx - bw[i], lx + bw[i], length.out = 11)
        cuts  <- sum(vapply(setdiff(seq_len(n), i), function(j) {
                   yj <- arc_y_at(arcs[j, ], probe)
                   any(!is.na(yj) & abs(yj - ly) < bh)
                 }, logical(1)))
        # On-page slope of the OWN arc across the label's width.
        rise  <- abs(diff(range(arc_y_at(a, probe), na.rm = TRUE))) /
                 y_max * trk_h_in
        run   <- 2 * bw[i] / (e - s) * trk_w_in
        slope <- if (run > 0) rise / run else 0
        # Leaving the arc at all costs 1.2 on top of the distance, so a
        # number only takes the half step when it has NO seat on its
        # own curve. Without the fixed part TRIM38's 22.7 bought its
        # way off a perfectly good seat for half a push (1.5) to dodge
        # one hairline notch (1.5).
        sc <- CUT_W * cuts + 3.0 * abs(sgn) + 1.2 * (sgn != 0) +
              3.0 * min(slope, 2) + 0.9 * abs(tt - 0.5)
        if (sc < best_s) { best_s <- sc; best <- c(lx, ly) }
      }
    }
    if (is.null(best)) {
      mid  <- (a$x1 + a$x2) / 2
      best <- c(mid, if (below) max(floor_y + bh, arc_y_at(a, mid))
                     else       min(ceil_y, arc_y_at(a, mid) + off))
      warning(sprintf("place_psi: no clear slot for PSI label %s", a$lab),
              call. = FALSE)
    }
    px[i] <- best[1]; py[i] <- best[2]
  }
  # Every column travels as DATA, never as an aes() reference to a
  # scalar in this environment -- a panel restored from
  # panel_cached()'s saveRDS renders blank otherwise
  # (../figure2/compose.md).
  data.frame(x = px[idx], y = py[idx], lab = arcs$lab[idx],
             col = arcs$col[idx], arc = idx,
             bw = bw[idx], bh = bh, bw2 = bw2[idx], bh2 = bh2,
             stringsAsFactors = FALSE)
}

# ------------------------------------------------------------------
# Draw one arc, split wherever it runs into a PSI label's box. NA-ing
# the points would make geom_path drop rows with a warning and join
# across the hole, so the surviving points are grouped into contiguous
# runs and the group travels as a data column.
# ------------------------------------------------------------------
arc_layer <- function(a, boxes, color, linewidth, x_pt, y_pt, min_pt = 3,
                      own = NA_integer_) {
  cv  <- arc_curve(a)
  hit <- rep(FALSE, nrow(cv))
  if (!is.null(boxes) && nrow(boxes) > 0)
    for (j in seq_len(nrow(boxes))) {
      # Wide break for this arc's own number, hairline notch for the rest.
      mine <- !is.na(own) && identical(boxes$arc[j], own)
      hw_j <- if (mine) boxes$bw[j] else boxes$bw2[j]
      hh_j <- if (mine) boxes$bh[j] else boxes$bh2[j]
      hit  <- hit | (abs(cv$x - boxes$x[j]) <= hw_j &
                     abs(cv$y - boxes$y[j]) <= hh_j)
    }
  cv$grp <- cumsum(c(1L, abs(diff(hit))))
  cv     <- cv[!hit, , drop = FALSE]
  if (nrow(cv) == 0) return(NULL)
  # A run only a few points long is the sliver of arc left between two
  # neighbouring numbers. Drawn, it reads as a dash -- "66.2 - 4.87" --
  # not as a line, so drop it and let the two labels share one break.
  len_pt <- vapply(split(cv, cv$grp), function(g)
              sum(sqrt((diff(g$x) / x_pt)^2 + (diff(g$y) / y_pt)^2)),
              numeric(1))
  cv <- cv[as.character(cv$grp) %in% names(len_pt)[len_pt >= min_pt], ,
           drop = FALSE]
  if (nrow(cv) == 0) return(NULL)
  geom_path(data = cv, aes(x = x, y = y, group = grp), inherit.aes = FALSE,
            color = color, linewidth = linewidth, lineend = "round")
}

# ------------------------------------------------------------------
# One coverage / sashimi track.
#
# sashimi_track_gg draws the coverage (and owns the theming, so the
# tracks look identical to figure 3's); its own arc engine is switched
# OFF with min_psi = Inf because it interpolates a two-element
# apex_range over its own short->long ordering, which cannot express
# seven solved per-arc peaks (figure 3 panel c documents the same
# limit). The FORM of the curve is identical to the helper's -- see
# arc_frame above. Its x scale is REPLACED (not extended by coord,
# which would not stop the scale dropping data) so the drawn window is
# `view`, not `wide`.
#
# The arcs are built HERE, not by the caller, because their feet are
# sampled from this track's own coverage.
# ------------------------------------------------------------------
make_track <- function(L, cohort, col, lk_df, orange_key = character(0),
                       caption = NULL, caption_n = NULL,
                       down_key = character(0), drop = 0,
                       ceil = ARC_CEIL, bot = 0) {
  wd    <- file.path(PY, L$gene)
  bw    <- file.path(wd, if (cohort == "high") BW_HIGH else BW_CTRL)
  lk    <- file.path(wd, if (cohort == "high") LK_HIGH else LK_CTRL)
  y_max <- L$y_max
  s <- L$view[1]; e <- L$view[2]
  x_pt  <- (e - s) / (TRACK_W_IN * 72)          # data-x per drawn pt
  y_pt  <- y_max   / (TRACK_H_IN * 72)          # data-y per drawn pt
  y_lo  <- -drop * y_max                        # bottom of the drop band

  # The x SCALE keeps the left gutter (the axis rule, its ticks and the
  # numbers live there) -- a scale limit would DROP them before they
  # were drawn. `coord` sets the view, and carries the drop band.
  # sashimi_track_gg() grew a y axis of its own (for figures 1-3). It
  # draws the rule at `start`, which for this panel is 161 bp INSIDE
  # the drawn window -- PAD_FRAC exists precisely so TRIM38's first
  # exon is not against the frame, and an axis there would stand in the
  # coverage. Panel b draws its own, below, so switch the shared one
  # off -- defensively, in case that argument is not there yet.
  sash <- list(bw_path = bw, chrom = L$chrom,
               start = L$wide[1], end = L$wide[2],
               track_color = col, track_label = "",
               links_path = lk, min_psi = Inf, y_max = y_max,
               label_family = FONT_FAMILY,
               track_height_in = TRACK_H_IN, track_width_in = TRACK_W_IN)
  if ("y_axis" %in% names(formals(sashimi_track_gg))) sash$y_axis <- FALSE
  p <- suppressMessages(
        do.call(sashimi_track_gg, sash) +
        # Right limit keeps the caption gutter too: the cohort captions sit
        # at e + 3 pt, and a scale limit at `e` DROPPED all four (2026-10-06,
        # "Removed 1 row" x4). coord_cartesian still sets the view.
        scale_x_continuous(limits = c(s - AX_GUT * x_pt, e + MARG[2] * x_pt),
                           expand = c(0, 0)) +
        coord_cartesian(xlim = c(s, e), ylim = c(y_lo, y_max),
                        clip = "off"))

  d <- bigwig_binned(bw, L$chrom, L$wide[1], L$wide[2])

  # Two fans: the ordinary one bowing up, and (TRIM38's high track
  # only) the orange pair bowing DOWN into the drop band. They are
  # solved independently -- an up arc lives at y >= its own feet and a
  # down arc at y <= its own, so no pair drawn from different fans can
  # cross except at a foot they share.
  is_dn   <- paste(lk_df$s1, lk_df$e2) %in% down_key
  arcs_up <- arc_frame(lk_df[!is_dn, , drop = FALSE], d, y_max,
                       orange_key, ceil = ceil)
  arcs_dn <- if (any(is_dn))
               arc_frame(lk_df[is_dn, , drop = FALSE], d, y_max, orange_key,
                         dir = -1, ceil = DOWN_CEIL, min_peak = DOWN_MIN)
  arcs  <- rbind(arcs_up, arcs_dn)
  n_up  <- if (is.null(arcs_up)) 0L else nrow(arcs_up)

  if (!is.null(arcs) && nrow(arcs) > 0) {
    # Regression check against "junctions go under other junctions",
    # run per fan (in each fan's own, possibly mirrored, space).
    for (fan in list(list("up", arcs_up), list("down", arcs_dn))) {
      fa <- fan[[2]]
      if (is.null(fa) || nrow(fa) == 0) next
      ag <- arc_min_gap(fa)
      cat(sprintf("[panel b]   %-6s %-4s %-4s feet/peaks (of y_max): %s\n",
                  L$gene, cohort, fan[[1]],
                  paste(sprintf("%s %.3f->%.3f ^%.3f", fa$lab, fa$y1,
                                fa$y2, fa$peak), collapse = " | ")))
      cat(sprintf("[panel b]   %-6s %-4s %-4s min signed gap %s   %s\n",
                  L$gene, cohort, fan[[1]],
                  if (is.na(ag$gap)) "n/a"
                  else sprintf("%+.4f of y_max (%s)", ag$gap, ag$worst),
                  if (length(ag$forced))
                    sprintf("| %d crossing(s): %s", length(ag$forced),
                            paste(ag$forced, collapse = "; "))
                  else "| no crossings"))
      if (!is.na(ag$gap) && ag$gap <= 0)
        stop(sprintf(paste("[panel b] %s %s (%s fan): a SOLVABLE pair still",
                           "crosses (min signed gap %+.4f)."),
                     L$gene, cohort, fan[[1]], ag$gap))
    }

    # Labels are solved FIRST: the arcs are then drawn with a gap at
    # every label box, so nothing overprints anything.
    lab_df <- place_psi(arcs, d, y_max, TRACK_H_IN, TRACK_W_IN, s, e,
                        idx = seq_len(n_up))
    if (nrow(arcs) > n_up)
      lab_df <- rbind(lab_df,
                      place_psi(arcs, d, y_max, TRACK_H_IN, TRACK_W_IN, s, e,
                                idx = (n_up + 1L):nrow(arcs),
                                below = TRUE, floor_y = y_lo))
    # Every number must be sitting IN its own arc. "off" is how far its
    # centre had to leave the curve, in fractions of y_max; anything
    # non-zero means the contour idiom broke down for that label.
    cat(sprintf("[panel b]   %-6s %-4s labels: %s\n", L$gene, cohort,
                paste(sprintf("%s t=%.2f off=%+.3f", lab_df$lab,
                        (lab_df$x - arcs$x1) / (arcs$x2 - arcs$x1),
                        (lab_df$y - vapply(seq_len(nrow(arcs)), function(i)
                           arc_y_at(arcs[i, ], lab_df$x[i]), numeric(1))) /
                          y_max),
                      collapse = " | ")))
    # Black first, orange last, so the junctions the panel is ABOUT are
    # never buried under a neighbouring tier.
    for (i in c(which(arcs$col == "black"), which(arcs$col != "black")))
      p <- p + arc_layer(arcs[i, ], lab_df, arcs$col[i],
                         0.25 + 0.60 * sqrt(arcs$psi[i] / 100),
                         x_pt = (e - s) / (TRACK_W_IN * 72),
                         y_pt = y_max  / (TRACK_H_IN * 72), own = i)
    # The number is drawn in its arc's colour, so an orange arc's PSI
    # reads as part of the legend's claim rather than as a stray number.
    lab_df$col <- ifelse(lab_df$col == "black", "grey15", lab_df$col)
    p <- p +
      geom_text(data = lab_df,
                aes(x = x, y = y, label = lab, color = col),
                hjust = 0.5, vjust = 0.5, size = SZ_PSI,
                family = FONT_FAMILY, inherit.aes = FALSE,
                show.legend = FALSE) +
      scale_color_identity()
  }

  # ---- the y axis: a rule over 0..y_max, a tick at each end, the two
  # numbers hung off the ticks. The 0 is CENTRED on its tick so it sits
  # even with the baseline rather than floating above it; the top
  # number hangs just inside, because centring it would put half its
  # glyphs in the cell above, where patchwork clips. Every coordinate
  # is a DATA COLUMN -- an aes() closing over a local scalar renders
  # blank after panel_cached()'s saveRDS.
  ax_x <- s - AX_OFF * x_pt
  p <- p +
    geom_segment(data = data.frame(x = ax_x, xend = ax_x,
                                   y = 0, yend = y_max),
                 aes(x = x, xend = xend, y = y, yend = yend),
                 inherit.aes = FALSE, color = AX_COL, linewidth = AX_LW,
                 lineend = "square") +
    geom_segment(data = data.frame(x = ax_x, xend = ax_x - AX_TICK * x_pt,
                                   y = c(0, y_max), yend = c(0, y_max)),
                 aes(x = x, xend = xend, y = y, yend = yend),
                 inherit.aes = FALSE, color = AX_COL, linewidth = AX_LW,
                 lineend = "square") +
    geom_text(data = data.frame(
                x   = ax_x - (AX_TICK + AX_PAD) * x_pt,
                y   = c(y_max, 0),
                lab = c(sprintf("%.1f", y_max), "0"),
                vj  = c(1, 0.5)),
              aes(x = x, y = y, label = lab, vjust = vj),
              hjust = 1, size = SZ_AXIS, family = FONT_FAMILY,
              color = "grey15", inherit.aes = FALSE)

  if (!is.null(caption))
    p <- p + geom_text(
      data = data.frame(x = e + 3 * (e - s) / (TRACK_W_IN * 72), y = y_max * 0.5,
                        lab = paste0(caption, "\n", caption_n)),
      # hjust 0 + fixed 3 pt gap (hjust -0.06 shifts each line by its own width)
      aes(x = x, y = y, label = lab), hjust = 0, vjust = 0.5,
      size = SZ_COHORT, lineheight = 1.1, family = FONT_FAMILY,
      color = "grey15", inherit.aes = FALSE)

  p + theme(plot.margin = gutter(bot))
}

# ------------------------------------------------------------------
# Lavender taper WITH A NECK.
#
# ../scripts/coverage_plot.R's trapezoid_gg() runs the taper all the
# way to the bottom of its row, so lavender covers the whole row at
# every x the transcript spans. The italic gene name drawn by
# transcript_zoom_gg() is ~0.076 in tall in a ~0.12 in row, so it
# ALWAYS spills up out of the gene row and into the taper row -- and
# landed on the lavender at both loci.
#
# Here the taper reaches the window width at NECK_F of the row's height
# and then drops straight down as a neck, so below NECK_F the row is
# empty everywhere outside the window -- which is exactly where the
# gene name sits (label_side = "left" at both loci; the window starts
# well right of tx_lo at both). The neck also makes TRIM38's taper read
# as a taper at all: its window is 42% of the transcript, so a
# full-height trapezoid is a shallow sliver, and the two sides looked
# like unexplained risers rather than one funnel.
#
# Kept local rather than added to coverage_plot.R: figures 1-3 depend
# on that file and none of them needs a neck.
# ------------------------------------------------------------------
NECK_F <- 0.55

taper_gg <- function(x_lo, x_hi, win_lo, win_hi, fill = LAV) {
  ggplot() +
    annotate("polygon",
             x = c(x_lo, x_hi, win_hi, win_hi, win_lo, win_lo),
             y = c(1, 1, NECK_F, 0, 0, NECK_F),
             fill = fill, alpha = 1, color = NA) +
    scale_x_continuous(limits = c(x_lo, x_hi), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(plot.margin      = margin(0, 4, 0, 0, "pt"),
          plot.background  = element_rect(fill = NA, color = NA),
          panel.background = element_rect(fill = NA, color = NA))
}

# ------------------------------------------------------------------
# Gene model + coordinate rows.
# ------------------------------------------------------------------
make_gene_rows <- function(L, axis_bottom_pt) {
  ex    <- EXONS[[L$gene]]$mane
  tx_lo <- min(GenomicRanges::start(ex)); tx_hi <- max(GenomicRanges::end(ex))
  span  <- tx_hi - tx_lo
  win   <- c(max(L$view[1], tx_lo), min(L$view[2], tx_hi))
  brks  <- scales::breaks_pretty(4)(c(tx_lo, tx_hi))
  # Drop breaks that would collide with the chromosome label parked at
  # the left edge, or run off the right-hand end of the gene model.
  # Both cutoffs are DERIVED from the labels' real drawn width at this
  # point size rather than being eyeballed offsets.
  lab_w <- 8   * SZ_AXIS * .pt * 0.56 / 72
  chr_w <- 4.5 * SZ_AXIS * .pt * 0.56 / 72
  at    <- (brks - tx_lo) / span * TRACK_W_IN
  keep  <- at - lab_w / 2 > chr_w + 0.05 & at + lab_w / 2 < TRACK_W_IN + 0.22 &
           brks >= tx_lo & brks <= tx_hi
  brks  <- brks[keep]

  # Budget check for the gene name's spill (see H_TRAP above). It sits
  # with its BOTTOM 0.415 row-half-units below the gene row's top edge,
  # so it reaches (txt - 0.415 * GENE_H_IN / 2) above that edge and has
  # NECK_F * TRAP_H_IN of empty taper row to reach into.
  spill <- txt_h_in(SZ_GENE) - 0.415 * GENE_H_IN / 2
  if (spill > NECK_F * TRAP_H_IN)
    warning(sprintf(paste("[panel b] %s: gene name spills %.3f in into the",
                          "taper row but only %.3f in is clear -- raise",
                          "H_TRAP or NECK_F"),
                    L$gene, spill, NECK_F * TRAP_H_IN), call. = FALSE)
  # ... and the name must stay left of the neck.
  name_w <- nchar(L$gene) * 0.52 * SZ_GENE * .pt / 72
  if (name_w > (win[1] - tx_lo) / span * TRACK_W_IN)
    warning(sprintf("[panel b] %s: gene name is wider than the clear strip",
                    L$gene), call. = FALSE)

  trap <- taper_gg(tx_lo, tx_hi, win[1], win[2]) +
    theme(plot.margin = gutter())
  gene <- transcript_zoom_gg(
            ex, tx_lo, tx_hi,
            highlight = win, highlight_fill = LAV, highlight_alpha = 1,
            tx_label = L$gene, exon_h = 0.85, n_arrows = 26,
            label_size = SZ_GENE, label_family = FONT_FAMILY,
            x_ticks = brks, tick_h = 0.55, label_side = L$label_side) +
    theme(plot.margin = gutter())
  axis <- ggplot() +
    geom_text(data = data.frame(x = brks, y = 1,
                                lab = format(brks, scientific = FALSE,
                                             trim = TRUE)),
              aes(x = x, y = y, label = lab), hjust = 0.5, vjust = 1,
              size = SZ_AXIS, family = FONT_FAMILY, color = "grey15") +
    geom_text(data = data.frame(x = tx_lo, y = 1, lab = L$chrom),
              aes(x = x, y = y, label = lab), hjust = 0, vjust = 1,
              size = SZ_AXIS, family = FONT_FAMILY, color = "grey15") +
    scale_x_continuous(limits = c(tx_lo, tx_hi), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, 1), clip = "off") +
    theme_void() +
    theme(plot.margin      = gutter(axis_bottom_pt),
          plot.background  = element_rect(fill = NA, color = NA),
          panel.background = element_rect(fill = NA, color = NA))
  list(trap = trap, gene = gene, axis = axis)
}

# ------------------------------------------------------------------
# Assemble one locus.
# ------------------------------------------------------------------
build_locus <- function(key, axis_bottom_pt = 0) {
  L  <- LOCI[[key]]
  wd <- file.path(PY, L$gene)
  hi <- read_links(file.path(wd, LK_HIGH), L$wide[1], L$wide[2])
  ct <- read_links(file.path(wd, LK_CTRL), L$wide[1], L$wide[2])
  hi <- hi[order(hi$e2 - hi$s1), , drop = FALSE]

  key_or <- orange_keys(L$gene, hi)
  cat(sprintf("[panel b] %-6s %s:%d-%d  y_max %.4f\n",
              key, L$chrom, L$wide[1], L$wide[2], L$y_max))
  cat(sprintf("[panel b]   high PSI %s | ctrl PSI %s\n",
              paste(hi$psi, collapse = ", "), paste(ct$psi, collapse = ", ")))
  cat(sprintf("[panel b]   orange (GSP #30 junction group) -> %s\n",
              if (length(key_or) == 0) "NONE"
              else paste(hi$psi[paste(hi$s1, hi$e2) %in% key_or],
                         collapse = " / ")))

  # Orange marks the driving junctions in the HIGH cohort only -- the
  # controls' copy of a shared endpoint (TRIM38's 1.38) is not what
  # drives the GSP, and the reference draws it black.
  #
  # TRIM38 ONLY: its two orange arcs bow DOWN, into the band below the
  # baseline, and its control track's top arc (94.2) is capped so the
  # two do not close on each other. KDM4C is untouched -- its 10.7 and
  # 11.2 meet on the cryptic-exon bump between them, which is the
  # visual proof of the panel.
  dn <- if (identical(L$gene, "TRIM38")) key_or else character(0)
  trk_hi <- make_track(L, "high", RED,  hi, key_or,        HIGH_LAB, HIGH_N,
                       down_key = dn,
                       drop = if (length(dn)) DROP_F else 0,
                       # only the flush KDM4C boundary needs the gutter
                       bot  = if (length(dn)) 0 else AX_BOT_PT)
  trk_ct <- make_track(L, "ctrl", GREY, ct, character(0),  CTRL_LAB, CTRL_N)
  g      <- make_gene_rows(L, axis_bottom_pt)

  # Returned as a flat LIST, not a patchwork. `/` between two column
  # patchworks merges them into one flat stack and silently loses the
  # outer plot_layout's heights ("Need 7 panels, but ... only provide
  # 3") -- see ../figure2/compose.md. Panel b is a single column, so
  # the honest construction is one flat stack of every row with one
  # heights vector.
  list(trk_hi, trk_ct, g$trap, g$gene, g$axis)
}

# ------------------------------------------------------------------
# Legend row: one orange arc glyph + its caption, in the panel's own
# 0..1 space. Sized from the caption's drawn width so the pair stays
# centred under the tracks whatever the panel width is.
# ------------------------------------------------------------------
LEG_TXT <- "Unannotated junction driving GSP #30"
LEG_H_IN <- H_LEG * B_SCALE
LEG_CAP  <- txt_h_in(SZ_LEG) / LEG_H_IN        # cap height, in row units
if (txt_h_in(SZ_LEG) * 1.30 > LEG_H_IN)        # + descenders
  warning(sprintf(paste("[panel b] legend caption needs %.3f in but the",
                        "legend row is only %.3f in -- its descenders will",
                        "be clipped by the row below. Raise H_LEG."),
                  txt_h_in(SZ_LEG) * 1.30, LEG_H_IN), call. = FALSE)
legend_row <- local({
  txt_w  <- nchar(LEG_TXT) * 0.50 * (SZ_LEG * .pt) / 72 / TRACK_W_IN
  glyph_w <- 0.055
  gap     <- 0.012
  x0      <- 1 - (txt_w + glyph_w + gap)        # right-aligned, as published
  t  <- seq(0, 1, length.out = 100)
  ggplot() +
    # arc_bump() with y1 = y2 IS the arcs' own curve, so the glyph and
    # the thing it explains are the same shape.
    # Glyph baseline and cap height are set from the caption's own
    # drawn size, so the arc sits on the text's baseline and rises to
    # its cap height rather than floating at a guessed height.
    geom_path(data = data.frame(x = x0 + t * glyph_w,
                                y = 0.50 - 0.5 * LEG_CAP +
                                    arc_bump(t) * LEG_CAP),
              aes(x = x, y = y), color = ORANGE, linewidth = 0.8) +
    geom_text(data = data.frame(x = x0 + glyph_w + gap, y = 0.50,
                                lab = LEG_TXT),
              aes(x = x, y = y, label = lab), hjust = 0, vjust = 0.5,
              size = SZ_LEG, family = FONT_FAMILY, color = "grey15") +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(plot.margin      = gutter(),
          plot.background  = element_rect(fill = NA, color = NA),
          panel.background = element_rect(fill = NA, color = NA))
})

rows <- c(build_locus("KDM4C", axis_bottom_pt = 6),
          build_locus("TRIM38"),
          list(legend_row))
# TRIM38's high track is taller by exactly its drop band, so the
# COVERAGE part of all four tracks is still TRACK_H_IN (B_SCALE has the
# matching correction). Rows 1-5 KDM4C, 6-10 TRIM38, 11 the legend.
ROW_HEIGHTS <- c(H_TRK + H_AXGAP,       H_TRK, H_TRAP, H_GENE, H_AXIS,
                 H_TRK * (1 + DROP_F),  H_TRK, H_TRAP, H_GENE, H_AXIS,
                 H_LEG)

wrap_elements(full =
  Reduce(`/`, rows) +
  plot_layout(heights = ROW_HEIGHTS) +
  plot_annotation(
    title = paste0("Top junctions of GSP #30 capture\n",
                   "usage of deep intronic splice sites"),
    theme = theme(plot.title = element_text(
                    size = SIZE_TITLE, face = "bold", hjust = 0.5,
                    family = FONT_FAMILY, margin = margin(0, 0, 3, 0)),
                  plot.margin = margin(0, 0, 0, 0))))
}))
