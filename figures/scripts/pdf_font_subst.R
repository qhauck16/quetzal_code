# Re-set an external vector PDF's text in the figures' own typeface.
#
# WHY: panels 1a and 2a are Illustrator schematics. Their labels are
# real vector text, but in Arial / Arial-Bold / Arial-Italic (plus a
# TeXGyreHeros-Bold title), all EMBEDDED as font subsets. Every
# R-drawn panel uses `FONT_FAMILY` ("Helvetica"), which fontconfig
# resolves to Nimbus Sans on this machine. Arial and Nimbus Sans are
# metrically compatible but visually distinct (R leg, G spur, t and C
# terminals, digit shapes), so the schematics read as a different
# typeface from the rest of the page. Quinn asked for them to match.
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
# This is the only route available on midway3: there is no mutool,
# qpdf binary, pikepdf, PyMuPDF or inkscape here, and poppler/cairo
# cannot substitute a font on their own.
#
# KNOWN COSMETIC COST: the PostScript round trip makes Ghostscript
# render the black block arrows' hairline stroke (a purple outline in
# the source art) about a pixel wide, where poppler drew it
# sub-pixel. It is a property of the gs round trip, not of the font
# flags -- an embedded-font control run shows the same halo. Rasterise
# the substituted PDF at a higher dpi (see the `dpi` argument at the
# call sites) so the hairline stays sub-pixel relative to the drawn
# size.
.libPaths('/scratch/midway2/qhauck/conda_env/rstudio-server/lib/R/library/')

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
