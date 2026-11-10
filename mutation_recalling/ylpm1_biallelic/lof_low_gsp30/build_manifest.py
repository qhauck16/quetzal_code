# GDC file manifest for the YLPM1 LoF mutants with GSP #30 < 1e-4 (samples.tsv; copied from top13_recall): per case, the tumour RNA-seq BAM of
# the exact aliquot used in TCGA (STAR 2-Pass Genome), the tumour WXS BAM from the SAME sample
# (e.g. -01A; full barcode, never 12-char patient -- the 2026-09 mislabel lesson), and the matched
# normal WXS (blood-derived preferred). Metadata only (no token). Out: manifest.tsv
import csv, json, urllib.request
rows = list(csv.DictReader(open("samples.tsv"), delimiter="\t"))
F = ["file_id", "file_name", "file_size", "experimental_strategy", "analysis.workflow_type",
     "cases.samples.submitter_id", "cases.samples.sample_type",
     "cases.samples.portions.analytes.aliquots.submitter_id"]
out = []
for r in rows:
    case = r["rna_barcode"][:12]; samp = r["rna_barcode"][:16]
    flt = {"op": "and", "content": [
        {"op": "in", "content": {"field": "cases.submitter_id", "value": [case]}},
        {"op": "in", "content": {"field": "data_format", "value": ["BAM"]}},
        {"op": "in", "content": {"field": "experimental_strategy", "value": ["WXS", "RNA-Seq", "WGS"]}}]}
    q = {"filters": flt, "fields": ",".join(F), "size": 200, "format": "json"}
    req = urllib.request.Request("https://api.gdc.cancer.gov/files", data=json.dumps(q).encode(),
                                 headers={"Content-Type": "application/json"})
    hits = json.load(urllib.request.urlopen(req))["data"]["hits"]
    def info(h):
        s = h["cases"][0]["samples"][0]
        al = [a["submitter_id"] for p in s.get("portions", []) for an in p.get("analytes", []) for a in an.get("aliquots", [])]
        return s["submitter_id"], s["sample_type"], al
    pick = {}
    for h in hits:
        sid, st, al = info(h); es = h["experimental_strategy"]; wf = h.get("analysis", {}).get("workflow_type", "")
        if es == "RNA-Seq" and wf == "STAR 2-Pass Genome" and r["rna_barcode"] in al: pick["rna_tumor"] = (h, sid)
        if es in ("WXS", "WGS") and sid == samp and wf in ("BWA with Mark Duplicates and BQSR", "BWA with BQSR"):
            k = es.lower() + "_tumor"
            if k not in pick or h["file_size"] > pick[k][0]["file_size"]: pick[k] = (h, sid)
        if es in ("WXS", "WGS") and "Normal" in st and wf in ("BWA with Mark Duplicates and BQSR", "BWA with BQSR"):
            k = es.lower() + "_normal"; better = "Blood" in st
            if k not in pick or (better and "Blood" not in pick[k][2]): pick[k] = (h, sid, st)
    for k, v in pick.items():
        h = v[0]
        out.append({"rank": r["rank"], "rail_id": r["rail_id"], "role": k, "sample": v[1], "file_id": h["file_id"],
                    "file_name": h["file_name"], "size_gb": round(h["file_size"] / 1e9, 1)})
    print(r["rank"], r["rail_id"], case, sorted(pick))
with open("manifest.tsv", "w") as f:
    w = csv.DictWriter(f, fieldnames=list(out[0]), delimiter="\t"); w.writeheader(); w.writerows(out)
