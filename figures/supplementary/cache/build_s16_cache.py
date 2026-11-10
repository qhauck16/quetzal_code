"""Cache for Supplementary Figure S17 (YLPM1 intron-coverage meta-plot).

import os
from pathlib import Path as _P
QROOT = str(next(p for p in [_P.cwd(), *_P.cwd().parents] if (p / ".quetzal_root").exists()))  # quetzal_code root
Runs the data half of tcga/genome_wide/flashier/diff_splicing/ylpm1/
plot_intron_retention_heatmap.py unchanged (same sample picks, seed,
Coding/UTR filter, per-gene dedup, bigwig normalisation) for the top
N_MAX genes by delta-PSI and writes the binned intron/exon ratio traces
as a long TSV. figS17/build.R picks the top N and orders rows by the
centre of mass of the top-15 trace, as the original does.
Imports the vendored copies in s16_lib/ (plot_intron_retention_heatmap.py,
plot_top25_intron_panels.py). Run with the python3.12 env (has pyBigWig).
"""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "s16_lib"))
import plot_intron_retention_heatmap as H

N_MAX = 100          # published S16 has 100 rows (--n_genes 100)
OUT = os.path.join(QROOT, "figures/supplementary/cache/s16_traces.tsv")

info = pd.read_csv(H.INFO_FILE, sep="\t", low_memory=False)
if "project_id" not in info.columns:
    meta = pd.read_csv(H.TCGA_META, sep="\t",
                       usecols=["rail_id", "gdc_cases.project.project_id"], low_memory=False)
    meta = meta.rename(columns={"gdc_cases.project.project_id": "project_id"})
    info = info.merge(meta, on="rail_id", how="left")
info = info.dropna(subset=["rail_id", H.FACTOR_COL, "project_id"]).copy()
info["rail_id"] = info["rail_id"].astype(int)
col1, ctrl = H.pick_samples_local(info, H.FACTOR_COL, H.SEED)
proj = info.set_index("rail_id")["project_id"].to_dict()
bw = {}
for rid in set(col1) | set(ctrl):
    if rid in proj:
        p = H.bw_path(rid, proj[rid], H.BW_DIR)
        if os.path.isfile(p):
            bw[rid] = p
print(f"col1 {len(col1)}  ctrl {len(ctrl)}  bw {len(bw)}")

df = H.load_filtered_df(H.DF_FINAL, H.LEAF2_CLASS, apply_coding_utr_filter=True)
introns = df.drop_duplicates(subset=["gene_name"], keep="first").head(N_MAX).to_dict("records")
columns = [("top", col1), ("ctrl", ctrl)]
from multiprocessing import Pool
def _one(it):
    return H.compute_row_data(it, columns, bw, H.FLANK_BP)
with Pool(10) as pool:                     # bigwig reads dominate; ~20 s/intron serial
    rds = pool.map(_one, introns)
rows = []
for k, (it, rd) in enumerate(zip(introns, rds)):
    com = H.com_of(rd["binned_traces"][0], str(it.get("strand", "+")))
    for (lab, _), binned, em in zip(columns, rd["binned_traces"], rd["exon_means"]):
        ratio = binned / em
        ratio = np.where(ratio < H.RATIO_THRESH, 0.0, ratio)
        if str(it.get("strand", "+")) == "-":
            ratio = ratio[::-1]          # display orientation: 5' left
        for b, v in enumerate(ratio):
            rows.append((k + 1, it["gene_name"], it["strand"], it["type"],
                         it["deltapsi"], com, lab, b + 1, float(v)))
    print(k + 1, it["gene_name"], f"{com:.1f}", flush=True)
pd.DataFrame(rows, columns=["dpsi_rank", "gene", "strand", "type", "deltapsi",
                            "com", "column", "bin", "ratio"]).to_csv(OUT, sep="\t", index=False)
print("wrote", OUT, "n_col1", len(col1), "n_ctrl", len(ctrl))
