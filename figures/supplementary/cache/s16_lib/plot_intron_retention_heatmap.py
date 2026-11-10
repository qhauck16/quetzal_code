#!/usr/bin/env python3
"""Intron-retention heatmap-style plot for YLPM1-axis disrupted introns.

For each top gene, draw one horizontal track per column (top 1-13,
top 14-25, matched controls):

  [exon]----[intron region, solid-colored by retention ratio]----[exon]

Exons are small black rectangles. The intron region is filled with a
solid color drawn from a blue (0) -> red (1) colormap. The color value
is

  group_mean_intron_coverage / max(group_mean_intron_coverage across the
                                   three groups)

so the highest-coverage group always lands at 1.0 (red) and the other
two are read off as fractions of it.

No splice-site arcs.
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches

# Reuse the sample-selection / bigwig / classifications logic.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_top25_intron_panels import (
    BW_DIR, INFO_FILE, DF_FINAL, TCGA_META, LEAF2_CLASS,
    FACTOR_COL, SEED, FLANK_BP, CTRL_LOADING_THRESH,
    bw_path, column_mean_bp,
)

# Local cohort definition (overrides the 13/25/N_CTRL split used by the
# coverage-panel script). Heatmap uses one case column (top 1-N_TOP) and
# one matched-control column.
N_TOP_LOCAL  = 15
N_CTRL_LOCAL = 200


def pick_samples_local(info_df, factor_col, seed,
                       n_top=N_TOP_LOCAL, n_ctrl=N_CTRL_LOCAL,
                       ctrl_thresh=CTRL_LOADING_THRESH):
    """Return (col1_ids, ctrl_ids) by factor_col desc.

    col1  = top n_top samples (ranks 1..n_top).
    ctrl  = n_ctrl samples with factor_col < ctrl_thresh, drawn so
            per-project proportions match col1.
    """
    info = info_df.sort_values(factor_col, ascending=False).reset_index(drop=True)
    col1 = info.iloc[:n_top]["rail_id"].tolist()
    col1_proj = info.iloc[:n_top]["project_id"].value_counts(normalize=True)
    pool = info[(info[factor_col] < ctrl_thresh)
                & (~info["rail_id"].isin(col1))]
    print(f"Control eligibility ({factor_col} < {ctrl_thresh}): "
          f"{len(pool)} samples in pool")

    rng = np.random.default_rng(seed)
    # Largest-remainder allocation so the per-cohort counts sum to exactly
    # n_ctrl.
    raw = col1_proj * n_ctrl
    alloc = np.floor(raw).astype(int)
    short = n_ctrl - int(alloc.sum())
    for proj in (raw - alloc).sort_values(ascending=False).index[:short]:
        alloc[proj] += 1
    ctrl = []
    for proj, prop in col1_proj.items():
        n_need = int(alloc[proj])
        avail = pool[pool["project_id"] == proj]
        if len(avail) == 0 or n_need == 0:
            continue
        ctrl.extend(
            avail.sample(n=min(n_need, len(avail)),
                          random_state=rng)["rail_id"].tolist()
        )
    return col1, ctrl

# ---- Visual constants ------------------------------------------------------

ROW_HEIGHT_IN = 0.09
HSPACE        = 0.08   # very little gap between stacked rows
EXON_FRAC     = 0.05   # fraction of panel width per exon box
EXON_HEIGHT   = 0.6    # axis height ("ratio space") at which the intron
                       # trace reaches the top of the exon box -- i.e.
                       # ratio = 1 maps to y0 + EXON_HEIGHT
EXON_BOTTOM   = 0.2    # y position of the exon-box bottom (panel uses
                       # ylim = [0, PANEL_YMAX] in ratio space)
PANEL_YMAX    = 1.4    # axis upper bound -- gives ~1.4x exon height of
                       # headroom for traces with intron > exon
EXON_AVG_BP   = 20     # flanking exon window used as the local "1.0"
                       # ratio reference (mean coverage in the 20 bp
                       # immediately upstream AND downstream of intron)
RATIO_THRESH  = 0.1    # global threshold: bins with ratio < 1/10 of
                       # the local exon coverage are rounded to 0
INTRON_BINS   = 200    # warp the intron coverage trace into this many
                       # bins so all panels render at identical width
# Solid maroon fill for every intron coverage trace -- no color scale.
INTRON_FILL   = "#800000"


def load_filtered_df(df_final_path, leaf2_class_path,
                     apply_coding_utr_filter=True):
    """Load df_final, keep alt 3'ss / alt 5'ss rows, optionally join the
    leafcutter2 classifications and drop Coding=True / UTR=True
    junctions, and return sorted by deltapsi (descending).
    """
    df = pd.read_csv(df_final_path, sep="\t", low_memory=False)
    df = df[df["type"].isin(["alt 3'ss", "alt 5'ss"])].copy()
    df["deltapsi"] = pd.to_numeric(df["deltapsi"], errors="coerce")
    df = df.dropna(subset=["deltapsi", "mane_intron_start", "mane_intron_end"])
    df["mane_intron_start"] = df["mane_intron_start"].astype(int)
    df["mane_intron_end"]   = df["mane_intron_end"].astype(int)
    df["start"] = df["start"].astype(int)
    df["end"]   = df["end"].astype(int)

    if apply_coding_utr_filter:
        cls = pd.read_csv(leaf2_class_path, sep="\t",
                           usecols=["Intron_coord", "Coding", "UTR"],
                           low_memory=False)
        coords = cls["Intron_coord"].str.extract(
            r"^(?P<chrom>[^:]+):(?P<s>\d+)-(?P<e>\d+)$")
        cls["chrom"] = coords["chrom"]
        cls["start"] = pd.to_numeric(coords["s"], errors="coerce") + 1
        cls["end"]   = pd.to_numeric(coords["e"], errors="coerce")
        cls = cls.dropna(subset=["chrom", "start", "end"])
        cls["start"] = cls["start"].astype(int)
        cls["end"]   = cls["end"].astype(int)
        cls["Coding"] = cls["Coding"].astype(str).str.lower().eq("true")
        cls["UTR"]    = cls["UTR"].astype(str).str.lower().eq("true")
        cls = (cls.groupby(["chrom", "start", "end"], as_index=False)
                  .agg({"Coding": "any", "UTR": "any"}))

        df = df.merge(cls, on=["chrom", "start", "end"], how="inner")
        df = df[(~df["Coding"]) & (~df["UTR"])].copy()

    df = df.sort_values("deltapsi", ascending=False).reset_index(drop=True)
    return df


def compute_row_data(intron, columns, bw_lookup, flank_bp):
    """Read BigWig coverage for every column at this intron and return
    the binned coverage trace, per-group mean intron coverage, ratios,
    and the per-column mean coverage in EXON_AVG_BP of flanking exon
    (used as the "1.0" reference when plotting -- intron / exon_mean).
    """
    chrom        = intron["chrom"]
    intron_start = int(intron["mane_intron_start"])
    intron_end   = int(intron["mane_intron_end"])
    win_start    = intron_start - flank_bp
    win_end      = intron_end   + flank_bp
    intron_len   = intron_end - intron_start
    exon_bp      = max(1, min(EXON_AVG_BP, flank_bp))

    group_means   = []
    binned_traces = []
    exon_means    = []
    for _, ids in columns:
        paths = [bw_lookup[rid] for rid in ids if rid in bw_lookup]
        sig = column_mean_bp(paths, chrom, win_start, win_end)
        intron_sig = sig[flank_bp:flank_bp + intron_len]
        group_means.append(float(np.nanmean(intron_sig)))

        # Flanking-exon reference = mean coverage in the last exon_bp bp
        # of the upstream exon AND the first exon_bp bp of the downstream
        # exon.  Coverage values are normalised per sample, so this is a
        # local "what does the exon look like at this gene" baseline.
        up_sig   = sig[flank_bp - exon_bp:flank_bp]
        down_sig = sig[flank_bp + intron_len:flank_bp + intron_len + exon_bp]
        em       = float(np.nanmean(np.concatenate([up_sig, down_sig])))
        exon_means.append(em if em > 0 else 1e-12)

        if len(intron_sig) > 0:
            edges = np.linspace(0, len(intron_sig),
                                 INTRON_BINS + 1).astype(int)
            binned = np.array([
                float(np.mean(intron_sig[edges[i]:edges[i + 1]]))
                if edges[i + 1] > edges[i] else 0.0
                for i in range(INTRON_BINS)
            ])
        else:
            binned = np.zeros(INTRON_BINS)
        binned_traces.append(binned)

    denom  = max(group_means) if max(group_means) > 0 else 1.0
    ratios = [m / denom for m in group_means]
    return {
        "binned_traces": binned_traces,
        "group_means":   group_means,
        "exon_means":    exon_means,
        "ratios":        ratios,
    }


def com_of(binned, strand):
    """Center of mass of a binned coverage trace, expressed in the
    displayed orientation (left = 5', right = 3'). Negative-strand
    introns flip the trace because the panel x-axis is inverted at
    draw time. Returns nan when the trace is all zero.
    """
    sig = np.asarray(binned, dtype=float)
    if strand == "-":
        sig = sig[::-1]
    sig = np.clip(sig, 0.0, None)
    s = float(np.sum(sig))
    if s <= 0:
        return float("nan")
    return float(np.sum(np.arange(len(sig)) * sig) / s)


def plot_retention(introns, row_data, columns, flank_bp, out_pdf,
                   row_height_in=None, hspace=None):
    n_rows = len(introns)
    n_cols = len(columns)
    rh     = ROW_HEIGHT_IN if row_height_in is None else float(row_height_in)
    hsp    = HSPACE        if hspace        is None else float(hspace)
    fig_w  = 13.0
    fig_h  = rh * n_rows + 1.0    # extra height for colorbar + titles
    fig    = plt.figure(figsize=(fig_w, fig_h))
    gs = gridspec.GridSpec(
        n_rows, n_cols, figure=fig,
        wspace=0.04, hspace=hsp,
        left=0.08, right=0.97,
        top=1 - 0.4 / fig_h, bottom=0.7 / fig_h,
    )

    for c in range(n_cols):
        ax_top = fig.add_subplot(gs[0, c])
        ax_top.set_title(columns[c][0], fontsize=10, pad=2)
        ax_top.remove()    # we only wanted the title slot; real axes below

    for r, (intron, rd) in enumerate(zip(introns, row_data)):
        strand        = str(intron.get("strand", "+"))
        gene_name     = intron["gene_name"]
        binned_traces = rd["binned_traces"]
        exon_means    = rd["exon_means"]

        ef   = EXON_FRAC
        ex_h = EXON_HEIGHT
        y0   = EXON_BOTTOM
        x_intron = np.linspace(ef, 1 - ef, INTRON_BINS)

        for c, ((label, _), binned, em) in enumerate(
                zip(columns, binned_traces, exon_means)):
            ax = fig.add_subplot(gs[r, c])
            # Exon boxes: top at ratio = 1 (y0 + ex_h).
            ax.add_patch(mpatches.Rectangle(
                (0, y0), ef, ex_h, color="black"))
            ax.add_patch(mpatches.Rectangle(
                (1 - ef, y0), ef, ex_h, color="black"))
            # Intron coverage as ratio = binned / exon_mean. The trace
            # is plotted in ratio space (1 maps to the exon-box top,
            # > 1 lets the trace exceed the exon box). RATIO_THRESH is
            # a *global* floor: any bin whose ratio < RATIO_THRESH
            # renders at zero -- same threshold across all introns.
            ratio = binned / em
            ratio = np.where(ratio < RATIO_THRESH, 0.0, ratio)
            trace_y = y0 + ex_h * ratio
            # Clip at the axis upper bound so traces don't draw outside.
            trace_y = np.clip(trace_y, y0, PANEL_YMAX)
            ax.fill_between(x_intron, y0, trace_y,
                             color=INTRON_FILL, linewidth=0)

            ax.set_xlim(0, 1)
            ax.set_ylim(0, PANEL_YMAX)
            ax.set_xticks([]); ax.set_yticks([])
            for s in ("top", "right", "bottom", "left"):
                ax.spines[s].set_visible(False)

            if strand == "-":
                ax.invert_xaxis()

            if c == 0:
                ax.text(-0.04, 0.5, gene_name,
                         transform=ax.transAxes, ha="right", va="center",
                         fontsize=7.5, fontstyle="italic")

            if r == 0:
                ax.set_title(label, fontsize=10, pad=4)

    fig.savefig(out_pdf)
    plt.close(fig)
    print(f"Wrote {out_pdf}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--info_file",   default=INFO_FILE)
    ap.add_argument("--df_final",    default=DF_FINAL)
    ap.add_argument("--tcga_meta",   default=TCGA_META)
    ap.add_argument("--bw_dir",      default=BW_DIR)
    ap.add_argument("--factor_col",  default=FACTOR_COL)
    ap.add_argument("--leaf2_class", default=LEAF2_CLASS)
    ap.add_argument("--n_genes",     type=int, default=1,
                    help="number of unique top genes to plot (ignored when "
                         "--all_introns is set)")
    ap.add_argument("--all_introns", action="store_true",
                    help="plot every alt 3'ss / alt 5'ss disrupted intron "
                         "in df_final (skips Coding/UTR filter and gene "
                         "dedup); row height derived to fit --total_height_in")
    ap.add_argument("--total_height_in", type=float, default=10.0,
                    help="target total figure height when --all_introns "
                         "is set (row height = (this - 1) / n_introns)")
    ap.add_argument("--flank_bp",    type=int, default=FLANK_BP)
    ap.add_argument("--seed",        type=int, default=SEED)
    ap.add_argument("--out",         default="intron_retention_heatmap.pdf")
    args = ap.parse_args()

    info = pd.read_csv(args.info_file, sep="\t", low_memory=False)
    if "project_id" not in info.columns:
        meta = pd.read_csv(args.tcga_meta, sep="\t",
                           usecols=["rail_id", "gdc_cases.project.project_id"],
                           low_memory=False)
        meta = meta.rename(columns={"gdc_cases.project.project_id": "project_id"})
        info = info.merge(meta, on="rail_id", how="left")
    info = info.dropna(subset=["rail_id", args.factor_col, "project_id"]).copy()
    info["rail_id"] = info["rail_id"].astype(int)

    col1, ctrl = pick_samples_local(info, args.factor_col, args.seed)
    print(f"col1 (ranks 1-{N_TOP_LOCAL}):      {len(col1)}")
    print(f"ctrl (project-matched, n={N_CTRL_LOCAL}): {len(ctrl)}")

    meta_full = info.set_index("rail_id")["project_id"].to_dict()
    bw_lookup = {}
    for rid in set(col1) | set(ctrl):
        proj = meta_full.get(rid)
        if proj is None:
            continue
        p = bw_path(rid, proj, args.bw_dir)
        if os.path.isfile(p):
            bw_lookup[rid] = p
    print(f"BW files resolved           : {len(bw_lookup)}")

    df = load_filtered_df(args.df_final, args.leaf2_class,
                           apply_coding_utr_filter=not args.all_introns)
    if args.all_introns:
        introns = df.to_dict("records")
        print(f"All disrupted introns        : {len(introns)}  "
              f"(no Coding/UTR filter, no dedup)")
    else:
        top_per_gene = df.drop_duplicates(subset=["gene_name"], keep="first")
        introns = top_per_gene.head(args.n_genes).to_dict("records")
        print(f"Top genes selected          : {len(introns)}")

    columns = [
        (f"top 1-{N_TOP_LOCAL}",                col1),
        (f"matched controls (n={N_CTRL_LOCAL})", ctrl),
    ]

    # Read BigWig coverage once per intron, then sort rows by the
    # display-space center of mass of the top-1-13 (col 1) coverage
    # trace -- rightmost COM at top, leftmost at bottom.
    print(f"Reading BigWig coverage for {len(introns)} introns...")
    row_data = [compute_row_data(it, columns, bw_lookup, args.flank_bp)
                for it in introns]
    coms = [com_of(rd["binned_traces"][0], str(it.get("strand", "+")))
            for it, rd in zip(introns, row_data)]
    order = sorted(range(len(introns)),
                    key=lambda i: (-1e18 if np.isnan(coms[i]) else coms[i]),
                    reverse=True)
    introns  = [introns[i]  for i in order]
    row_data = [row_data[i] for i in order]
    coms     = [coms[i]     for i in order]

    print("Rows ordered by COM(col1):")
    for r, com in zip(introns, coms):
        print(
            f"  {r['gene_name']:>10s}  {r['chrom']}:{int(r['mane_intron_start'])}"
            f"-{int(r['mane_intron_end'])}  {r['type']}  "
            f"Delta_psi = {r['deltapsi']:+.3f}  COM = {com:.1f}"
        )

    row_height_in = None
    hspace        = None
    if args.all_introns and len(introns) > 0:
        row_height_in = max(0.005, (args.total_height_in - 1.0) / len(introns))
        hspace        = 0.0
        print(f"Derived row height          : {row_height_in:.4f} in  "
              f"(target total {args.total_height_in:.1f} in)")

    plot_retention(introns, row_data, columns, args.flank_bp, args.out,
                    row_height_in=row_height_in, hspace=hspace)


if __name__ == "__main__":
    main()
