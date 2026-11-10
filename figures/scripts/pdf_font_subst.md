# `pdf_font_subst.R`

## Purpose & context

Panels **1a** and **2a** are Illustrator schematics, brought into the
figure as rasterised vector PDFs by `raster_utils.R::gg_from_pdf()`.
Their labels are real vector text in **Arial / Arial-Bold /
Arial-Italic** (plus a `TeXGyreHeros-Bold` title), embedded as font
subsets. Every R-drawn panel uses `theme.R`'s `FONT_FAMILY`
("Helvetica"), which fontconfig resolves to **Nimbus Sans** on this
account. Arial and Nimbus Sans are metrically compatible but visually
distinct, so the schematics read as a different typeface from the rest
of the page. Quinn asked (2026-09-23) for 1a and 2a to match.

`pdf_subst_font(path)` returns the path of a `<name>_helv.pdf` beside
the source whose text has been re-set in Nimbus Sans. It regenerates
only when the source is newer, so re-exporting the schematic from
Illustrator is all that is needed.

## Decision rationale

The mechanism: `pdftops -noembt1 -noembtt -noembcidps -noembcidtt`
emits PostScript that *references* the fonts by name instead of
embedding the subsets; Ghostscript then resolves those names through
its own Fontmap, which aliases the Arial/Helvetica names onto Nimbus
Sans — exactly the face the R panels draw in. `pdfwrite` re-embeds the
result.

Alternatives considered and why they were rejected:

| Route | Verdict |
| --- | --- |
| Edit the `.ai` in Illustrator | Correct, but not something Claude can do on the cluster; also loses the "regenerate from source" property. |
| Redraw the labels natively in R over the artwork | Needs the original text REMOVED from the PDF first, which needs a PDF content-stream editor. |
| `mutool` / `qpdf` CLI / `pikepdf` / PyMuPDF / Inkscape | **None are installed on midway3.** (The R `qpdf` package is, but it cannot touch fonts.) |
| poppler / cairo font substitution | Neither honours a substitution request when the font is embedded. |

So the `pdftops -noemb*` + Ghostscript round trip is the only route
available here.

## Gotchas & warnings

- **It converts text to outlines.** `pdffonts` on the output is empty
  and `pdftotext` returns nothing. Harmless *here* only because these
  PDFs are rasterised into the figure anyway — do not reuse this
  helper for anything whose text must stay selectable.
- **Hairline strokes get thicker.** The Ghostscript round trip renders
  a sub-pixel stroke at ~1 px, so figure 2a's black block arrows pick
  up a visible purple outline. This is the round trip, NOT the font
  flags — an embedded-font control run shows the same halo. Both call
  sites therefore rasterise at `dpi = 600` instead of 300 so the
  hairline stays sub-pixel relative to the drawn size.
- **Metric compatibility is what makes this safe.** Arial, Helvetica
  and Nimbus Sans share advance widths, so no string changes length
  and the schematic's boxes still fit their labels. A substitution
  into a non-metric-compatible family would reflow the art.
- The `_helv.pdf` files are **generated artifacts**. The `.pdf` (and
  `.ai`) beside them are the source of truth.
- On any failure the function `warning()`s and returns the ORIGINAL
  path, so a figure still builds (with the old font) on a machine
  without `pdftops`/`gs`.

## User preferences

Quinn: typography must be consistent across the whole paper — that is
the entire reason `theme.R` is shared between figures. A schematic in
a different face defeats it.

## Input/output

Input: any vector PDF with embedded text. Output: a sibling
`*_helv.pdf` with the same page geometry, text drawn in the
Ghostscript-resolved (Nimbus Sans) faces.

## File flow schematic

```
panel_a/<name>.pdf ──► pdf_subst_font() ──► panel_a/<name>_helv.pdf
   (Illustrator)        pdftops -noemb*           │
                          + gs pdfwrite           ▼
                                        gg_from_pdf(dpi = 600)
                                                  │
                                                  ▼
                                        raster_cache/*.png ──► compose.R
```
