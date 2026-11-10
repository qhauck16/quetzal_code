#!/usr/bin/env python
"""Transform a tracks.ini produced by AggregateBigwigsForPlotting.py into a
sashimi-only ini, replicating the post-processing in
updated_single_gene_figures.Rmd::make_sashimi_ini():

  * scale every max_value line by 1.1 of its OWN value (so TCGA and GTEx
    sashimi panels each get headroom relative to their own coverage, instead
    of being forced to share the larger of the two)
  * drop the block from "## Loop thru non-shashimi plots" up to (but not
    including) the next "[vlines]" section
  * set the first two `height = ` lines to `height = 10`
  * flip `show_data_range = false` -> true
"""

import re
import sys


def transform(lines):
    lines = list(lines)

    max_idx = [i for i, l in enumerate(lines)
               if re.match(r'^max_value\s*=\s*[0-9]', l)]
    for i in max_idx:
        v = float(re.sub(r'^max_value\s*=\s*', '', lines[i]).strip())
        lines[i] = f"max_value = {v * 1.1}\n"

    start = next((i for i, l in enumerate(lines)
                  if l.startswith("## Loop thru non-shashimi plots")), None)
    end   = next((i for i, l in enumerate(lines)
                  if l.startswith("[vlines]")), None)
    if start is not None and end is not None and start < end:
        del lines[start + 1:end]

    height_idx = [i for i, l in enumerate(lines)
                  if re.match(r'^height\s*=', l)]
    for i in height_idx[:2]:
        lines[i] = "height = 10\n"

    for i, l in enumerate(lines):
        if re.match(r'^show_data_range\s*=\s*false', l):
            lines[i] = "show_data_range = true\n"

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
