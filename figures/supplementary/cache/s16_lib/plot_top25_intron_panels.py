# Vendored 2026-10-06 from diff_splicing/ylpm1/ (S16 cache dependency); relative paths made absolute,
# INFO_FILE repointed to the quetzal_code variant table.
#!/usr/bin/env python3
"""Spot-check coverage plot for ONE YLPM1-axis disrupted intron.

import os
from pathlib import Path as _P
QROOT = str(next(p for p in [_P.cwd(), *_P.cwd().parents] if (p / ".quetzal_root").exists()))  # quetzal_code root
Layout: 1 row x 3 columns.
  - col 1: rank 1-13  on k30
  - col 2: rank 14-25 on k30
  - col 3: 200 controls (k30 < 1e-4) drawn from project_id distribution
           proportional to col 1.

The plotted genomic window for every column is the SAME:
  [mane_intron_start - FLANK, mane_intron_end + FLANK]
so the three panels are bp-for-bp comparable. FLANK = 20 nt by default
so the panel includes a sliver of each flanking exon -- the exon
coverage anchors the eye on what "full" coverage looks like for each
column, making the intron-retention difference visible.

Per-sample library-size normalization (divide each sample's per-bp
signal by `bw.header()['sumData']`) is applied before averaging across
the column's samples, so expression-magnitude differences between
samples don't bias the column means.
"""

import os
import re
import sys
import argparse
import numpy as np
import pandas as pd
import pyBigWig
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import PathPatch
from matplotlib.path import Path as MplPath

ALT_COLORS    = {"alt 3'ss": "#2CA02C",   # green
                 "alt 5'ss": "#FF7F0E"}   # orange
EXON_COLOR    = "#1F77B4"                 # blue: exonic flanks
INTRON_COLOR  = "#D62728"                 # red: intronic body


# ---- Defaults --------------------------------------------------------------

BW_DIR     = os.path.join(QROOT, "inputs/bigwigs/tcga")
INFO_FILE  = os.path.join(QROOT, "figures/figure4/cache/YLPM1_full_plus_variants.tsv")  # vt_path("YLPM1_full_plus_variants")
DF_FINAL   = os.path.join(QROOT, "leafcutter_analyses/diff_splicing/unified_results/ylpm1/20v60/dpsi0.01_fc2.5/df_final.tsv")
TCGA_META  = os.path.join(QROOT, "inputs/tcga_v2_samples.tsv")
PSI_DIR    = os.path.join(QROOT, "leafcutter_analyses/psi_files/tcga")
LEAF2_CLASS = os.path.join(QROOT, "leafcutter_analyses/leaf2_classification/clustering/Leaf2_junction_classifications.txt")
FACTOR_COL = "k30"
N_TOP      = 25      # top samples covered by cols 1 + 2
SPLIT_AT   = 13      # col 1: ranks 1..SPLIT_AT; col 2: SPLIT_AT+1..N_TOP
N_CTRL     = 200
SEED       = 123
CTRL_LOADING_THRESH = 1e-4   # control pool restricted to factor_col < this
FLANK_BP   = 50      # +/- bp into each flanking exon
INTRON_BINS = 200    # intron is squashed to this many bins so all
                     # panels render exons + introns at identical widths
ARC_LW_MIN  = 0.4    # arc thickness at PSI = 0
ARC_LW_MAX  = 3.5    # arc thickness at PSI = 1
ROW_HEIGHT_IN = 0.65 # vertical inches per row in the grid


# ---- Sample selection ------------------------------------------------------

def pick_samples(info_df, factor_col, seed):
    """Return (col1_ids, col2_ids, ctrl_ids) by factor_col desc.

    col 1 = ranks 1..SPLIT_AT
    col 2 = ranks SPLIT_AT+1..N_TOP
    col 3 = N_CTRL samples with factor_col < CTRL_LOADING_THRESH,
            per-project proportions matched to col 1.
    """
    info = info_df.sort_values(factor_col, ascending=False).reset_index(drop=True)
    col1 = info.iloc[:SPLIT_AT]["rail_id"].tolist()
    col2 = info.iloc[SPLIT_AT:N_TOP]["rail_id"].tolist()

    col1_proj = info.iloc[:SPLIT_AT]["project_id"].value_counts(normalize=True)
    top_ids   = set(col1) | set(col2)
    pool = info[(info[factor_col] < CTRL_LOADING_THRESH)
                & (~info["rail_id"].isin(top_ids))]
    print(f"Control eligibility ({factor_col} < {CTRL_LOADING_THRESH}): "
          f"{len(pool)} samples in pool")

    rng = np.random.default_rng(seed)
    ctrl = []
    for proj, prop in col1_proj.items():
        n_need = int(round(prop * N_CTRL))
        avail = pool[pool["project_id"] == proj]
        if len(avail) == 0 or n_need == 0:
            continue
        ctrl.extend(
            avail.sample(n=min(n_need, len(avail)),
                         random_state=rng)["rail_id"].tolist()
        )
    return col1, col2, ctrl


# ---- BigWig path lookup ----------------------------------------------------

def cohort_short(project_id):
    return project_id.replace("TCGA-", "")


def bw_path(rail_id, project_id, bw_root):
    return os.path.join(bw_root, cohort_short(project_id), f"{rail_id}.bw")


# ---- Per-bp coverage, per-sample library-size normalized ------------------

# ---- Per-junction column-mean PSI for arc thickness ----------------------

def column_mean_psi(gene_name, chrom, cluster_id, start, end, rail_ids,
                    psi_dir, _cache={}):
    """Mean PSI of one junction across `rail_ids`. Per-gene PSI table is
    cached because we revisit it for the three columns of each row.
    """
    key = (chrom, gene_name)
    if key not in _cache:
        path = os.path.join(psi_dir, chrom, gene_name, "psi_table.tsv")
        if not os.path.isfile(path) or os.path.getsize(path) < 10:
            _cache[key] = None
        else:
            try:
                # First column has no header (just rail_id values) so pandas
                # otherwise pushes them into a phantom column. index_col=0
                # makes rail_id the row index.
                df = pd.read_csv(path, sep="\t", engine="c", index_col=0)
                df.index.name = "rail_id"
                _cache[key] = df
            except Exception:
                _cache[key] = None
    psi = _cache[key]
    if psi is None:
        return float("nan")
    # PSI table columns are keyed by `<numeric_cluster_id>:<start>-<end>`,
    # but the numeric cluster id comes from a different leafcutter run
    # than the one that produced df_final -- the cluster numbering does
    # not line up across runs. Match on the start-end suffix only, since
    # those genomic coordinates uniquely identify the junction.
    suffix = f":{start}-{end}"
    matching = [c for c in psi.columns if c.endswith(suffix)]
    if not matching:
        return float("nan")
    col_key = matching[0]
    wanted = [int(r) for r in rail_ids]
    present = [r for r in wanted if r in psi.index]
    if not present:
        return float("nan")
    return float(pd.to_numeric(psi.loc[present, col_key],
                               errors="coerce").mean(skipna=True))


# ---- x-axis warping (fixed exon-flank fraction, variable intron) ----------

def _exon_frac():
    return FLANK_BP / (2 * FLANK_BP + INTRON_BINS)


def warp_pos(pos, win_start, intron_start, intron_end, win_end):
    """Map a genomic position into the warped [0, 1] x-axis where every
    panel reserves the same fraction of width for the flanking exons.
    """
    ef = _exon_frac()
    if pos < intron_start:
        return (pos - win_start) / max(1, intron_start - win_start) * ef
    if pos > intron_end:
        return (1 - ef) + (pos - intron_end) / max(1, win_end - intron_end) * ef
    return ef + (pos - intron_start) / max(1, intron_end - intron_start) * (1 - 2 * ef)


def warp_signal(sig, win_start, intron_start, intron_end, win_end):
    """Resample a length-(win_end-win_start) per-bp signal onto the warped
    x grid: keep exon flanks at raw resolution, bin-mean the intron to
    INTRON_BINS. Returns three (x, y) tuples so the caller can colour
    exon vs intron regions separately. Each region's first/last points
    overlap the neighbour's so fill_between renders without seams.
    """
    ef = _exon_frac()
    li = intron_start - win_start            # left exon length in indices
    ri = win_end - intron_end                # right exon length in indices
    intron_len = max(0, intron_end - intron_start)

    sig_left   = sig[:li]
    sig_intron = sig[li:li + intron_len]
    sig_right  = sig[li + intron_len:]

    if intron_len > 0:
        edges = np.linspace(0, intron_len, INTRON_BINS + 1).astype(int)
        intron_binned = np.array([
            float(np.mean(sig_intron[edges[i]:edges[i + 1]]))
            if edges[i + 1] > edges[i] else 0.0
            for i in range(INTRON_BINS)
        ])
    else:
        intron_binned = np.zeros(INTRON_BINS)

    # Replace the boundary y-values so the exon fill ends at the same
    # height as the intron fill begins (and vice versa). This is purely
    # cosmetic -- avoids a vertical notch where the three polygons meet.
    sig_left  = sig_left.astype(float).copy()
    sig_right = sig_right.astype(float).copy()
    if intron_len > 0:
        if len(sig_left)  > 0: sig_left[-1] = intron_binned[0]
        if len(sig_right) > 0: sig_right[0] = intron_binned[-1]

    x_left   = np.linspace(0,      ef,        max(li, 1), endpoint=True)
    x_intron = np.linspace(ef,     1 - ef,    INTRON_BINS, endpoint=True)
    x_right  = np.linspace(1 - ef, 1,         max(ri, 1), endpoint=True)

    return (x_left, sig_left), (x_intron, intron_binned), (x_right, sig_right)


def column_mean_bp(bw_paths, chrom, start, end):
    """Per-bp average coverage over [start, end), with per-sample
    library-size normalization (divide by sumData). Returns a length-
    (end - start) float vector.
    """
    L = end - start
    out = np.full((len(bw_paths), L), np.nan, dtype=np.float64)
    for i, p in enumerate(bw_paths):
        if not os.path.isfile(p):
            continue
        try:
            with pyBigWig.open(p) as bw:
                hdr = bw.header()
                denom = hdr.get("sumData", 0.0)
                if denom <= 0:
                    continue
                # bw.values returns a list of per-bp floats (or None for gaps).
                raw = bw.values(chrom, start, end, numpy=True)
                if raw is None:
                    continue
                arr = np.asarray(raw, dtype=np.float64)
                arr[~np.isfinite(arr)] = 0.0
                out[i] = arr / denom
        except Exception:
            continue
    return np.nanmean(out, axis=0)


# ---- Plot ------------------------------------------------------------------

def draw_arc(ax, x_start, y_start, x_end, y_end, peak_y, color, lw=1.5):
    """Bezier arc anchored at (x_start, y_start) and (x_end, y_end), with
    a quadratic control point at the midpoint at height `peak_y`. Lets
    the two endpoints sit at different y-values so the arc joins the
    actual coverage at each splice site rather than a fixed baseline.
    """
    mid_x = (x_start + x_end) / 2.0
    verts = [
        (x_start, y_start),
        (mid_x,   peak_y),
        (x_end,   y_end),
    ]
    codes = [MplPath.MOVETO, MplPath.CURVE3, MplPath.CURVE3]
    patch = PathPatch(MplPath(verts, codes), fill=False,
                      edgecolor=color, lw=lw, clip_on=False)
    ax.add_patch(patch)


def plot_grid(introns, columns, bw_lookup, flank_bp, psi_dir, out_pdf,
              row_height_in=None):
    """One row per intron, 3 columns. Within a row:
       - All panels are plotted on a warped x in [0, 1] so the flanking
         exons take up the same fraction of width across all rows,
         independent of intron length.
       - Each panel's arc end-points sit at the coverage values at the
         splice sites; arc thickness scales with the column-mean PSI of
         the diff-spliced junction.
       - Negative-strand genes are mirrored along x so every row reads
         5' -> 3' left to right.
    """
    n_rows = len(introns)
    n_cols = len(columns)

    rh = ROW_HEIGHT_IN if row_height_in is None else float(row_height_in)
    fig_w = 13.0
    fig_h = rh * n_rows
    fig = plt.figure(figsize=(fig_w, fig_h))
    gs  = gridspec.GridSpec(
        n_rows, n_cols, figure=fig,
        wspace=0.04, hspace=0.4,
        left=0.08, right=0.99,
        top=1 - 0.15 / fig_h, bottom=0.15 / fig_h,
    )

    for r, intron in enumerate(introns):
        chrom         = intron["chrom"]
        intron_start  = int(intron["mane_intron_start"])
        intron_end    = int(intron["mane_intron_end"])
        junc_start    = int(intron["start"])
        junc_end      = int(intron["end"])
        cluster_id    = intron["cluster_id"]
        alt_type      = intron["type"]
        gene_name     = intron["gene_name"]
        strand        = str(intron.get("strand", "+"))
        win_start     = intron_start - flank_bp
        win_end       = intron_end   + flank_bp

        sigs_raw = [column_mean_bp(
                        [bw_lookup[rid] for rid in ids if rid in bw_lookup],
                        chrom, win_start, win_end)
                    for _, ids in columns]
        y_max_row = max(np.nanmax(s) for s in sigs_raw if np.isfinite(np.nanmax(s)))
        y_max_row = float(y_max_row) if y_max_row > 0 else 1e-12

        # Per-column mean PSI for this junction -> arc thickness.
        col_psis = [
            column_mean_psi(gene_name, chrom, cluster_id,
                             junc_start, junc_end, ids, psi_dir)
            for _, ids in columns
        ]

        arc_color  = ALT_COLORS.get(alt_type, "#888888")
        arc_height = 0.04 * y_max_row
        y_upper    = 1.10 * y_max_row

        # Warped junction endpoints (same across columns).
        jx_s = warp_pos(junc_start, win_start, intron_start, intron_end, win_end)
        jx_e = warp_pos(junc_end,   win_start, intron_start, intron_end, win_end)

        for c, (_, ids) in enumerate(columns):
            ax = fig.add_subplot(gs[r, c])
            left, mid, right = warp_signal(sigs_raw[c], win_start,
                                            intron_start, intron_end, win_end)
            # Blue exonic flanks, red intronic body.
            ax.fill_between(left[0],  0, left[1],
                            color=EXON_COLOR,   alpha=0.9, linewidth=0)
            ax.fill_between(mid[0],   0, mid[1],
                            color=INTRON_COLOR, alpha=0.9, linewidth=0)
            ax.fill_between(right[0], 0, right[1],
                            color=EXON_COLOR,   alpha=0.9, linewidth=0)

            # Coverage values AT the splice sites (raw indices into the
            # un-warped per-bp signal); use those as the arc anchor heights.
            idx_s = max(0, min(len(sigs_raw[c]) - 1, junc_start - win_start))
            idx_e = max(0, min(len(sigs_raw[c]) - 1, junc_end   - win_start))
            y_s   = float(np.nan_to_num(sigs_raw[c][idx_s], nan=0.0))
            y_e   = float(np.nan_to_num(sigs_raw[c][idx_e], nan=0.0))
            peak  = max(y_s, y_e) + arc_height

            psi = col_psis[c]
            lw  = ARC_LW_MIN + (ARC_LW_MAX - ARC_LW_MIN) * np.clip(
                0.0 if (psi is None or np.isnan(psi)) else psi, 0.0, 1.0)
            draw_arc(ax, jx_s, y_s, jx_e, y_e,
                     peak_y=peak, color=arc_color, lw=lw)

            ax.set_xlim(0, 1)
            ax.set_ylim(0, y_upper)
            ax.set_xticks([]); ax.set_yticks([])
            for s in ("top", "right", "bottom", "left"):
                ax.spines[s].set_visible(False)

            # Mirror so every gene reads 5'->3' left to right.
            if strand == "-":
                ax.invert_xaxis()

            if c == 0:
                ax.text(
                    -0.04, 0.5, gene_name,
                    transform=ax.transAxes,
                    ha="right", va="center",
                    fontsize=11, fontstyle="italic",
                )

    fig.savefig(out_pdf)
    plt.close(fig)
    print(f"Wrote {out_pdf}")


# ---- Main ------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--info_file",  default=INFO_FILE)
    ap.add_argument("--df_final",   default=DF_FINAL)
    ap.add_argument("--tcga_meta",  default=TCGA_META)
    ap.add_argument("--bw_dir",     default=BW_DIR)
    ap.add_argument("--factor_col", default=FACTOR_COL)
    ap.add_argument("--n_genes",    type=int, default=10,
                    help="number of unique top genes (by best deltapsi) to plot")
    ap.add_argument("--flank_bp",   type=int, default=FLANK_BP)
    ap.add_argument("--seed",       type=int, default=SEED)
    ap.add_argument("--out_dir",    default="ylpm1_introns_meta_plot",
                    help="directory for the three output PDFs")
    ap.add_argument("--leaf2_class", default=LEAF2_CLASS,
                    help="leafcutter2 junction classifications table; "
                         "junctions with Coding=True or UTR=True are dropped")
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

    col1, col2, ctrl = pick_samples(info, args.factor_col, args.seed)
    print(f"col1 (ranks 1-{SPLIT_AT}):   {len(col1)}")
    print(f"col2 (ranks {SPLIT_AT + 1}-{N_TOP}): {len(col2)}")
    print(f"ctrl (project-matched)      : {len(ctrl)}")

    meta_full = info.set_index("rail_id")["project_id"].to_dict()
    bw_lookup = {}
    for rid in set(col1) | set(col2) | set(ctrl):
        proj = meta_full.get(rid)
        if proj is None:
            continue
        p = bw_path(rid, proj, args.bw_dir)
        if os.path.isfile(p):
            bw_lookup[rid] = p
    print(f"BW files resolved           : {len(bw_lookup)}")

    df = pd.read_csv(args.df_final, sep="\t", low_memory=False)
    df = df[df["type"].isin(["alt 3'ss", "alt 5'ss"])].copy()
    df["deltapsi"] = pd.to_numeric(df["deltapsi"], errors="coerce")
    df = df.dropna(subset=["deltapsi", "mane_intron_start", "mane_intron_end"])
    df["mane_intron_start"] = df["mane_intron_start"].astype(int)
    df["mane_intron_end"]   = df["mane_intron_end"].astype(int)
    df["start"] = df["start"].astype(int)
    df["end"]   = df["end"].astype(int)

    # Join leafcutter2 junction classifications and keep only "true intron"
    # junctions -- Coding == False AND UTR == False. The classifications
    # table writes coords as `chr:start-end` with 0-based start, matching
    # end; df_final uses 1-based start, so add 1 to the parsed start to
    # build the join key.
    cls = pd.read_csv(args.leaf2_class, sep="\t",
                       usecols=["Intron_coord", "Coding", "UTR"],
                       low_memory=False)
    coords = cls["Intron_coord"].str.extract(r"^(?P<chrom>[^:]+):(?P<s>\d+)-(?P<e>\d+)$")
    cls["chrom"] = coords["chrom"]
    cls["start"] = pd.to_numeric(coords["s"], errors="coerce") + 1  # 0-based -> 1-based
    cls["end"]   = pd.to_numeric(coords["e"], errors="coerce")
    cls = cls.dropna(subset=["chrom", "start", "end"])
    cls["start"] = cls["start"].astype(int)
    cls["end"]   = cls["end"].astype(int)
    cls["Coding"] = cls["Coding"].astype(str).str.lower().eq("true")
    cls["UTR"]    = cls["UTR"].astype(str).str.lower().eq("true")
    # Some junctions are classified more than once (multiple gene contexts);
    # be conservative -- mark Coding/UTR True if ANY row says so.
    cls = (cls.groupby(["chrom", "start", "end"], as_index=False)
              .agg({"Coding": "any", "UTR": "any"}))

    n_before = len(df)
    g_before = df["gene_name"].nunique()
    df = df.merge(cls, on=["chrom", "start", "end"], how="inner")
    df = df[(~df["Coding"]) & (~df["UTR"])].copy()
    print(f"Leaf2 Coding=F & UTR=F filter: "
          f"{n_before} -> {len(df)} junctions  "
          f"({g_before} -> {df['gene_name'].nunique()} genes)")

    # Top-N unique-gene introns by best (highest deltapsi) junction.
    df = df.sort_values("deltapsi", ascending=False).reset_index(drop=True)
    top_per_gene = df.drop_duplicates(subset=["gene_name"], keep="first")
    top50 = top_per_gene.head(50).to_dict("records")
    top10 = top50[:10]
    print(f"Top-gene introns available  : {len(top_per_gene)}  "
          f"(plotting top 50)")

    columns = [
        (f"top 1-{SPLIT_AT}",          col1),
        (f"top {SPLIT_AT + 1}-{N_TOP}", col2),
        (f"matched controls",          ctrl),
    ]

    os.makedirs(args.out_dir, exist_ok=True)
    out10 = os.path.join(args.out_dir, "top10_intron_panels.pdf")
    out50 = os.path.join(args.out_dir, "top50_intron_panels_squished.pdf")

    # 1. Top 10, default row height.
    plot_grid(top10, columns, bw_lookup, args.flank_bp, PSI_DIR, out10)

    # 2. Top 50 -- same total figure height as top 10, so each row is
    # ROW_HEIGHT_IN / 5.
    plot_grid(top50, columns, bw_lookup, args.flank_bp, PSI_DIR, out50,
              row_height_in=ROW_HEIGHT_IN / 5.0)


if __name__ == "__main__":
    main()
