#!/usr/bin/env python
"""Transform a tracks.ini produced by AggregateBigwigsForPlotting.py into a
sashimi-only ini, replicating the post-processing in
updated_single_gene_figures.Rmd::make_sashimi_ini() plus the tweaks used
in tcga/gene_level_poisson/pygenome/*/sashimi_tracks.ini:

  * set every max_value line to the MAXIMUM of all tracks * 1.1, so TCGA
    and GTEx sashimi panels share the same y-scale (the reference
    figure does this so the two are directly comparable)
  * drop the block from "## Loop thru non-shashimi plots" up to (but not
    including) the next "[vlines]" section
  * set the first two `height = ` lines to `height = 10`
  * flip `show_data_range = false` -> true
  * bump `fontsize = 8` -> 14 on the sashimi link labels (bigger PSI %
    text -- reference-style readable)
"""

import re
import sys


def transform(lines):
    lines = list(lines)

    # tracks.ini from AggregateBigwigsForPlotting.py has 4 max_value
    # lines in this order:
    #   [0] TCGA whole-region bigwig     (~38.7 for NF1)
    #   [1] GTEx whole-region bigwig     (~24.4)
    #   [2] TCGA sashimi mock_for_axes   (~7.6)  <-- SASHIMI y-scale
    #   [3] GTEx sashimi mock_for_axes   (~3.3)  <-- SASHIMI y-scale
    #
    # The reference sashimi_tracks.ini at
    # tcga/gene_level_poisson/pygenome/NF1/sashimi_tracks.ini takes
    # max(indices [2:]) * 1.1 and applies it to both surviving bigwig
    # tracks (indices [0] and [1]). That's what gives the "compressed"
    # sashimi view with visible detail down low.
    max_idx = [i for i, l in enumerate(lines)
               if re.match(r'^max_value\s*=\s*[0-9]', l)]
    if len(max_idx) >= 4:
        sashimi_scale_vals = [
            float(re.sub(r'^max_value\s*=\s*', '', lines[i]).strip())
            for i in max_idx[2:]
        ]
        shared = max(sashimi_scale_vals) * 1.1
    elif max_idx:
        vals = [float(re.sub(r'^max_value\s*=\s*', '', lines[i]).strip())
                for i in max_idx]
        shared = max(vals) * 1.1
    else:
        shared = None
    if shared is not None:
        for i in max_idx[:2]:
            lines[i] = f"max_value = {shared}\n"

    start = next((i for i, l in enumerate(lines)
                  if l.startswith("## Loop thru non-shashimi plots")), None)
    end   = next((i for i, l in enumerate(lines)
                  if l.startswith("[vlines]")), None)
    if start is not None and end is not None and start < end:
        del lines[start + 1:end]

    # Track height. 12 gives a wide-enough aspect (~1.4 W/H) that the
    # rasterized sashimi fills a wide panel_e cell with preserve_aspect
    # = TRUE, without leaving huge horizontal whitespace.
    height_idx = [i for i, l in enumerate(lines)
                  if re.match(r'^height\s*=', l)]
    for i in height_idx[:2]:
        lines[i] = "height = 12\n"

    for i, l in enumerate(lines):
        if re.match(r'^show_data_range\s*=\s*false', l):
            lines[i] = "show_data_range = true\n"
        # Larger PSI % labels on the sashimi links
        if re.match(r'^fontsize\s*=\s*8\s*$', l):
            lines[i] = "fontsize = 14\n"

    return lines


def main(in_path, out_path):
    with open(in_path) as fh:
        lines = fh.readlines()
    out = transform(lines)
    with open(out_path, "w") as fh:
        fh.writelines(out)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("usage: make_sashimi_ini.py <tracks.ini> <sashimi_tracks.ini>")
    main(sys.argv[1], sys.argv[2])
