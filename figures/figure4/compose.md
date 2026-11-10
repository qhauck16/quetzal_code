
## Zoom connectors (2026-10-05)
The two lines from the CDF zoom box to the zoom panel end at constants in the
CDF's data units (`*_CALL_DX/_TOP/_BOT`). These are MEASURED, not tuned: run
`../scripts/measure_zoom_connectors.py out/<fig>.pdf x0 y0 x1 y1` (crop = page
fractions around the CDF + zoom) after any layout change to that row, and paste
the DX / TOP / BOT it prints. Hand-tuned values had drifted off the corners.
