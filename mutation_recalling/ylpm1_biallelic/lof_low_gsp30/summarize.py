# Summarise the top-13 YLPM1 re-call (call.sbatch): every YLPM1 HIGH/MODERATE consequence (VEP CSQ, first
# YLPM1 RefSeq transcript per variant) per sample x source (wxs / wgs Mutect2 tumour-normal, rna bcftools),
# with FILTER, tumour alt reads / depth / VAF. Out: recall_calls.tsv (+ printed per-sample summary)
import csv, glob, os, re
S = {r["rail_id"]: r for r in csv.DictReader(open("samples.tsv"), delimiter="\t")}
rows = []
for f in sorted(glob.glob("vep/*.vep.vcf")):
    rail, src = os.path.basename(f).split(".")[0].split("_")
    fmt_keys = None; samples = []
    for line in open(f):
        if line.startswith("##INFO=<ID=CSQ"): csq_keys = re.search(r"Format: ([^\"]+)", line).group(1).split("|")
        if line.startswith("#CHROM"): samples = line.rstrip().split("\t")[9:]
        if line.startswith("#"): continue
        c = line.rstrip().split("\t"); info = dict(kv.split("=", 1) if "=" in kv else (kv, True) for kv in c[7].split(";"))
        csq = [dict(zip(csq_keys, x.split("|"))) for x in info.get("CSQ", "").split(",")]
        csq = [x for x in csq if x.get("SYMBOL") == "YLPM1" and x.get("IMPACT") in ("HIGH", "MODERATE")]
        if not csq: continue
        x = csq[0]; keys = c[8].split(":")
        # tumour column: Mutect2 header line ##tumor_sample; else the single RNA sample
        vals = [dict(zip(keys, s.split(":"))) for s in c[9:]]
        tv = vals[0]
        if src in ("wxs", "wgs"):
            tname = [l for l in open(f) if l.startswith("##tumor_sample=")]
            tname = tname[0].strip().split("=", 1)[1] if tname else samples[0]
            tv = vals[samples.index(tname)]
        ad = [int(a) for a in tv.get("AD", "0,0").split(",") if a not in (".", "")] or [0, 0]
        alt = sum(ad[1:]); dp = sum(ad)
        rows.append({"rank": S[rail]["rank"], "rail_id": rail, "project": S[rail]["project_id"], "source": src,
                     "pos": c[1], "ref": c[3], "alt": c[4], "filter": c[6], "consequence": x["Consequence"],
                     "hgvsp": x.get("HGVSp", "").split(":")[-1] or f'{x.get("Amino_acids")}@{x.get("Protein_position")}',
                     "alt_reads": alt, "depth": dp, "vaf": round(alt / dp, 3) if dp else "",
                     "known": S[rail]["proteinChange"]})
rows.sort(key=lambda r: (int(r["rank"]), int(r["pos"]), r["source"]))
with open("recall_calls.tsv", "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
for r in rows:
    print("\t".join(str(r[k]) for k in ("rank", "rail_id", "source", "pos", "filter", "consequence", "hgvsp", "alt_reads", "depth", "vaf")))
