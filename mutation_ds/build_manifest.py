#!/usr/bin/env python3
"""Exact GDC manifest for the mutation_ds samples that lack a .junc file.
import os
from pathlib import Path as _P
QROOT = str(next(p for p in [_P.cwd(), *_P.cwd().parents] if (p / ".quetzal_root").exists()))  # quetzal_code root
One STAR 2-Pass Genome BAM per full aliquot barcode; fails on any barcode with
zero or several matches. Writes download/manifest.tsv (gdc-client format plus a
`dir` column, the gdc_parallel shard layout) and download/rail_map.tsv."""
import csv, json, sys, urllib.request
D = os.path.join(QROOT, "mutation_ds")
# Optional suffix: build_manifest.py null -> samples/download_list_null.tsv ->
# download/{manifest,rail_map}_null.tsv (same bams/ folder).
SUF = f"_{sys.argv[1]}" if len(sys.argv) > 1 else ""
rows = [r for r in csv.DictReader(open(f"{D}/samples/download_list{SUF}.tsv"), delimiter="\t")
        if r["have_junc"] == "FALSE"]
bc2rail = {}
for r in rows: bc2rail.setdefault(r["tcga_barcode"], []).append(r["rail_id"])
filt = {"op": "and", "content": [
    {"op": "in", "content": {"field": "cases.samples.portions.analytes.aliquots.submitter_id", "value": list(bc2rail)}},
    {"op": "=", "content": {"field": "data_type", "value": "Aligned Reads"}},
    {"op": "=", "content": {"field": "analysis.workflow_type", "value": "STAR 2-Pass Genome"}}]}
q = json.dumps({"filters": filt, "size": 5000,
                "fields": "file_id,file_name,md5sum,file_size,state,cases.samples.portions.analytes.aliquots.submitter_id"}).encode()
hits = json.load(urllib.request.urlopen(urllib.request.Request(
    "https://api.gdc.cancer.gov/files", q, {"Content-Type": "application/json"})))["data"]["hits"]
by_bc = {}
for h in hits:
    alq = {a["submitter_id"] for c in h["cases"] for s in c["samples"] for p in s["portions"]
           for an in p["analytes"] for a in an["aliquots"]}
    for bc in alq & set(bc2rail): by_bc.setdefault(bc, []).append(h)
bad = {bc: len(by_bc.get(bc, [])) for bc in bc2rail if len(by_bc.get(bc, [])) != 1}
if bad: sys.exit(f"barcodes without exactly one genome BAM: {bad}")
out = f"{D}/download/bams"
with open(f"{D}/download/manifest{SUF}.tsv", "w") as m, open(f"{D}/download/rail_map{SUF}.tsv", "w") as rm:
    m.write("id\tfilename\tmd5\tsize\tstate\tdir\n")
    rm.write("rail_id\ttcga_barcode\tfile_id\tfilename\tmd5\tsize\n")
    for bc, (h,) in sorted(by_bc.items()):
        m.write(f"{h['file_id']}\t{h['file_name']}\t{h['md5sum']}\t{h['file_size']}\t{h['state']}\t{out}\n")
        for rail in bc2rail[bc]:
            rm.write(f"{rail}\t{bc}\t{h['file_id']}\t{h['file_name']}\t{h['md5sum']}\t{h['file_size']}\n")
print(f"{len(by_bc)} barcodes -> {len(by_bc)} BAMs, {sum(h[0]['file_size'] for h in by_bc.values())/1e9:.0f} GB")
