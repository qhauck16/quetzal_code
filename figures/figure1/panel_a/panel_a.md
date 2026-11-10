# Figure 1a — Quetzal workflow schematic

`fig1a_updated.ai` — Illustrator drawing. No R code produces this panel;
it is an illustrated schematic of a per-gene splice-junction count
matrix flowing into Poisson NMF and out as an `L` × `F` factorisation
(with KRAS as the working example gene).

If regenerating, open the `.ai` file directly in Illustrator.

## Typography (2026-09-23)

The Illustrator art sets its labels in **Arial / Arial-Bold /
Arial-Italic** plus a TeXGyreHeros-Bold title, all embedded as font
subsets. Every R-drawn panel uses `FONT_FAMILY` ("Helvetica"), which
fontconfig resolves to **Nimbus Sans** on this machine — a visibly
different face (R leg, G spur, t/C terminals, digits). Quinn asked for
the schematic to match the rest of the figure.

`compose.R` therefore feeds the art through
`../scripts/pdf_subst_font()` before rasterising it. That re-sets the
text in Nimbus Sans by stripping the embedded subsets
(`pdftops -noembt1 -noembtt -noembcidps -noembcidtt`) and letting
Ghostscript resolve the font names through its own Fontmap. The
families are metrically compatible, so nothing reflows.

**Do not "fix" this by editing the raster.** If the schematic is ever
re-exported from `fig1a_updated.ai`, just re-export it — the
substitution regenerates automatically (it keys off the source's
mtime) and writes `gene_level_workflow_helv.pdf` beside it. That
`_helv.pdf` is a GENERATED file; the `.pdf`/`.ai` remain the source.

Cost of the route: the Ghostscript round trip renders hairline strokes
(the block arrows' outline, and the title's stroke) about a pixel wide
where poppler drew them sub-pixel. `compose.R` rasterises at 600 dpi
rather than 300 so they stay sub-pixel at the drawn size.
