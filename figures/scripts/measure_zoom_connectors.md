# measure_zoom_connectors.py

**Purpose** — Sets the endpoints of the zoom-callout lines in figures 2d, 3a and 4a so
they hit the zoom panel's two left corners exactly. Written 2026-10-05, after Quinn saw
the hand-tuned endpoints missing the corners.

**How** — It renders the PDF with ghostscript at 600 dpi.
- The zoom box's black edges (data x −4…0.25, y 0…102) give the CDF's px→data scale.
- The zoom panel's grey20 left axis line gives the target corners.
- It prints DX, TOP and BOT in CDF data units.

**Gotchas**
- Curves cross the box edges, so vertical runs bridge gaps of up to 40 px.
- Box top and bottom come from the box's horizontal edges, not the vertical runs.
- Keep the crop to one row; a box in another panel can be picked up instead.
- The lines are drawn outside the CDF's pinned axis range, so changing them doesn't move
  the layout, and the measurement stays valid.

```
out/figureN.pdf ─► measure_zoom_connectors.py ─► DX/TOP/BOT ─► <panel>/build.R *_CALL_* ─► compose
```
