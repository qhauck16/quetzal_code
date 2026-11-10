"""Measure where a CDF zoom box's connectors must end so they hit the zoom
panel's left corners exactly (figures 2d, 3a, 4a).

The zoom box is drawn in the CDF's data units: x from ZOOM_FROM to
X_RIGHT (= max log loading + 0.25), y from 0 to 102. Its black edges in the
render therefore give the CDF's px<->data scale. The zoom panel's left axis
line (grey20) gives the two target corners. Output: DX (data-x reach past
the box's right edge), TOP and BOT (data-y of the zoom's top / bottom-left
corners) -- the values for <X>_CALL_DX / _TOP / _BOT.

usage: python measure_zoom_connectors.py figure.pdf x0 y0 x1 y1 [zoom_from x_right]
  x0..y1 = crop of the page holding the CDF + zoom, as page fractions.
"""
import sys, subprocess, tempfile, os
import numpy as np
from PIL import Image

pdf = sys.argv[1]; x0, y0, x1, y1 = map(float, sys.argv[2:6])
zoom_from = float(sys.argv[6]) if len(sys.argv) > 6 else -4.0
x_right   = float(sys.argv[7]) if len(sys.argv) > 7 else 0.25
DPI = 600
png = tempfile.mktemp(suffix=".png")
subprocess.run(["gs", "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=pnggray", f"-r{DPI}",
                "-dFirstPage=1", "-dLastPage=1", f"-sOutputFile={png}", pdf], check=True)
g = np.array(Image.open(png)).astype(int); os.remove(png)
H, W = g.shape
c0, r0 = int(x0 * W), int(y0 * H)
g = g[r0:int(y1 * H), c0:int(x1 * W)]

def vruns(mask, min_len, max_gap=40):
    """Longest vertical run per column, bridging gaps <= max_gap px (a
    curve crossing a box edge breaks it), as (col, top, bottom, length)."""
    out = []
    for c in range(mask.shape[1]):
        rows = np.where(mask[:, c])[0]
        if len(rows) < min_len: continue
        splits = np.where(np.diff(rows) > max_gap)[0]
        segs = np.split(rows, splits + 1)
        seg = max(segs, key=lambda s: s[-1] - s[0])
        if seg[-1] - seg[0] >= min_len:
            out.append((c, int(seg[0]), int(seg[-1]), int(seg[-1] - seg[0])))
    return out

def centre(runs):
    """Merge adjacent columns of one stroke; report its centre column."""
    groups, cur = [], []
    for r in runs:
        if cur and r[0] != cur[-1][0] + 1: groups.append(cur); cur = []
        cur.append(r)
    if cur: groups.append(cur)
    return [(np.mean([r[0] for r in gr]), int(np.median([r[1] for r in gr])),
             int(np.median([r[2] for r in gr])), max(r[3] for r in gr)) for gr in groups]

black = centre(vruns(g < 35, 150))                     # zoom-box edges
if "--debug" in sys.argv: print("black", [(round(a[0]), a[1], a[2]) for a in black])
grey  = centre(vruns((g >= 35) & (g < 90), 150))       # axis lines (grey20)
# Box edges: black strokes whose TOPS match (a curve can extend a stroke
# downward); rightmost such pair. Top/bottom then come from the box's own
# horizontal edges: rows where black spans >= 85% of the box width.
pairs = [(a, b) for i, a in enumerate(black) for b in black[i + 1:]
         if abs(a[1] - b[1]) < 6 and b[0] - a[0] > 20]
L, R = max(pairs, key=lambda p: p[1][0])
span = g[:, int(L[0]):int(R[0]) + 1] < 35
lo, hi = min(L[1], R[1]) - 5, max(L[2], R[2]) + 5
full = np.where(span.mean(1) >= 0.85)[0]
full = full[(full >= lo) & (full <= hi)]
top_px, bot_px = float(full.min()), float(full.max())
sx = (x_right - zoom_from) / (R[0] - L[0])          # data-x per px
sy = 102.0 / (bot_px - top_px)                      # data-y per px
zoom_axis = min((a for a in grey if a[0] > R[0] + 5), key=lambda a: a[0])
DX  = (zoom_axis[0] - R[0]) * sx
TOP = 102.0 + (top_px - zoom_axis[1]) * sy
BOT = (bot_px - zoom_axis[2]) * sy
print(f"box px x {L[0]:.1f}-{R[0]:.1f} y {top_px:.0f}-{bot_px:.0f} | zoom axis x {zoom_axis[0]:.1f} "
      f"y {zoom_axis[1]}-{zoom_axis[2]}")
print(f"DX={DX:.4f} TOP={TOP:.3f} BOT={BOT:.3f}")
