# Re-set an external vector PDF's text in the figures' own typeface.
#
# WHY: panels 1a and 2a are Illustrator schematics whose labels are vector
# text in embedded Arial subsets. Every R-drawn panel uses `FONT_FAMILY`
# ("Helvetica"), which fontconfig resolves to Nimbus Sans. The two are
# metrically compatible but visually distinct, so the schematics are re-set
# to match the rest of the page.
#
# HOW: `pdftops -noembt1 -noembtt -noembcidps -noembcidtt` writes
# PostScript that REFERENCES the fonts by name instead of embedding
# the subsets; Ghostscript then resolves those names through its own
# Fontmap, which aliases Arial* / Helvetica* onto the Nimbus Sans
# family -- i.e. exactly the face the R panels draw in. The result is
# re-embedded by `pdfwrite`. Because the families are metrically
# compatible, nothing reflows: every string keeps its advance width
# and the artwork's boxes still fit their labels.
#
# (No mutool, qpdf, pikepdf, PyMuPDF or inkscape is available on the build
# host; poppler/cairo cannot substitute a font on their own.)
#
# COSMETIC COST: the PostScript round trip renders the black block arrows'
# hairline stroke about a pixel wide. Rasterise the substituted PDF at a
# higher dpi (the `dpi` argument at the call sites) to keep it sub-pixel.

# Returns the path of the substituted PDF, regenerating it only when
# the source is newer (or `force`). Output lands beside the source as
# `<name>_helv.pdf`; the raster cache keys off the full path, so the
# two variants never collide.
pdf_subst_font <- function(pdf_path, force = FALSE) {
  stopifnot(file.exists(pdf_path))
  out <- sub("\\.pdf$", "_helv.pdf", pdf_path)
  if (!force && file.exists(out) &&
      file.info(out)$mtime >= file.info(pdf_path)$mtime) {
    return(out)
  }
  cat("[font] re-setting ", pdf_path, " -> ", basename(out), "\n", sep = "")
  ps <- tempfile(fileext = ".ps")
  on.exit(unlink(ps), add = TRUE)
  st <- system2("pdftops",
                c("-noembt1", "-noembtt", "-noembcidps", "-noembcidtt",
                  "-level3", shQuote(pdf_path), shQuote(ps)),
                stdout = TRUE, stderr = TRUE)
  if (!file.exists(ps)) {
    warning("pdftops failed for ", pdf_path, "; using the original.\n",
            paste(st, collapse = "\n"))
    return(pdf_path)
  }
  tmp_pdf <- tempfile(fileext = ".pdf")
  on.exit(unlink(tmp_pdf), add = TRUE)
  st <- system2("gs",
                c("-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=pdfwrite",
                  "-dEmbedAllFonts=true",
                  paste0("-sOutputFile=", shQuote(tmp_pdf)), shQuote(ps)),
                stdout = TRUE, stderr = TRUE)
  if (!file.exists(tmp_pdf) || file.info(tmp_pdf)$size < 1000) {
    warning("ghostscript failed for ", pdf_path, "; using the original.\n",
            paste(st, collapse = "\n"))
    return(pdf_path)
  }
  file.copy(tmp_pdf, out, overwrite = TRUE)
  out
}
