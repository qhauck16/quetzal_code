# Figure 2, panel a. Sourced by ../compose.R, which supplies the
# shared setup (theme tokens, FLASHIER/CACHE paths, ROW_H, panel_cached,
# missing_panel, sf3b1_info). Must leave a `panel_a` object defined.
#
# Rebuild just this panel:   Rscript compose.R a

# ------------------------------------------------------------------
# a. Workflow schematic (vector art; title re-set here so its
#    typography matches the R-drawn panels)
# ------------------------------------------------------------------
# preserve_aspect = TRUE: the art is 5.83:1 and filling a 4.90:1 cell
# stretched its text ~19% horizontally. ROW_H[1] below is tuned so the
# cell aspect matches the art instead, so nothing is distorted.
#
# The art's labels are Arial; every other panel draws in FONT_FAMILY
# (Helvetica -> Nimbus Sans here). `pdf_subst_font()` re-sets the
# schematic's text in that face -- see ../scripts/pdf_font_subst.R for
# why this is the only workable route and why the dpi is raised.
panel_a <- panel_cached("a", wrap_elements(full =
  (gg_from_pdf(pdf_subst_font("panel_a/genome_wide_workflow.pdf"),
               dpi = 600, preserve_aspect = TRUE) +
     ggtitle("Quetzal's genome-wide workflow") +
     theme(plot.title = element_text(size = SIZE_TITLE, face = "bold",
                                      family = FONT_FAMILY, hjust = 0.5,
                                      margin = margin(0, 0, 3, 0)),
           plot.margin = margin(0, 0, 0, 0, "pt")))))

