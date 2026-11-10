# Figure 2, panel a. Sourced by ../compose.R, which supplies the
# shared setup (theme tokens, FLASHIER/CACHE paths, ROW_H, panel_cached,
# missing_panel, sf3b1_info). Must leave a `panel_a` object defined.
#
# Rebuild just this panel:   Rscript compose.R a

# ------------------------------------------------------------------
# a. Workflow schematic (vector art; title re-set here so its
#    typography matches the R-drawn panels)
# ------------------------------------------------------------------
# preserve_aspect = TRUE; ROW_H[1] in compose.R is tuned so the cell
# aspect matches the art's 5.83:1.
#
# pdf_subst_font() re-sets the schematic's Arial text in FONT_FAMILY
# (see ../scripts/pdf_font_subst.R).
panel_a <- panel_cached("a", wrap_elements(full =
  (gg_from_pdf(pdf_subst_font("panel_a/genome_wide_workflow.pdf"),
               dpi = 600, preserve_aspect = TRUE) +
     ggtitle("Quetzal's genome-wide workflow") +
     theme(plot.title = element_text(size = SIZE_TITLE, face = "bold",
                                      family = FONT_FAMILY, hjust = 0.5,
                                      margin = margin(0, 0, 3, 0)),
           plot.margin = margin(0, 0, 0, 0, "pt")))))

