#!/usr/bin/env python3
"""Consolidate the 2026-07 benchmark into clean TSVs for figures + docs.

Reads directly from the result FASTAs (not the buggy compare.sh), plus the QC
stats and the same-individual (v3) results. Writes:
  bench/results_v2.tsv    cross-individual (12 species)
  bench/results_v3.tsv    same-individual (7 species)
  bench/qc.tsv            per-dataset QC (Illumina + HiFi)
  bench/master.tsv        one row/species merging v2 + QC + v3 identity (drives figures)
"""
import subprocess, os
from pathlib import Path

V2 = Path("/path/to/benchmark_v2_data/results")
V2REF = Path("/path/to/benchmark_v2_data/refs")
V3 = Path("/path/to/benchmark_v3/results")
RAW2 = Path("/path/to/benchmark_v2_data/0_rawdata")
BENCH = Path("/path/to/ngs245/bench")

ORDER = dict(Beta_vulgaris="Caryophyllales", Helianthus_annuus="Asterales",
    Solanum_lycopersicum="Solanales", Sesamum_indicum="Lamiales",
    Citrus_sinensis="Sapindales", Glycine_max="Fabales", Fragaria_vesca="Rosales",
    Populus_trichocarpa="Malpighiales", Vitis_vinifera="Vitales",
    Actinidia_chinensis="Ericales", Oryza_sativa="Poales",
    Musa_acuminata="Zingiberales", Lindera_aggregata="Laurales")

V2_SP = ["Beta_vulgaris","Helianthus_annuus","Solanum_lycopersicum","Sesamum_indicum",
    "Citrus_sinensis","Glycine_max","Fragaria_vesca","Populus_trichocarpa",
    "Vitis_vinifera","Actinidia_chinensis","Oryza_sativa","Musa_acuminata"]
V3_SP = ["Oryza_sativa","Beta_vulgaris","Citrus_sinensis","Sesamum_indicum",
    "Musa_acuminata","Populus_trichocarpa","Lindera_aggregata"]

def flen(fa):
    if not fa.exists() or fa.stat().st_size == 0:
        return None
    return sum(len(l.strip()) for l in fa.read_text().splitlines() if not l.startswith(">"))

def blast_pid(q, s):
    if not (q.exists() and s.exists() and q.stat().st_size and s.stat().st_size):
        return None
    r = subprocess.run(["blastn","-query",str(q),"-subject",str(s),
        "-outfmt","6 pident length","-dust","no"], capture_output=True, text=True)
    best = sorted((ln.split("\t") for ln in r.stdout.splitlines() if ln),
                  key=lambda x: int(x[1]), reverse=True)
    return float(best[0][0]) if best else None

def gb_acc(sp):
    f = V2REF / f"{sp}.its.fasta"
    if f.exists() and f.stat().st_size:
        h = f.read_text().splitlines()[0]
        return h[1:].split()[0]
    return "-"

def fmt(x, d="-"):
    return d if x is None else (f"{x:.3f}" if isinstance(x, float) else str(x))

# ---- v2 ----
v2 = {}
for sp in V2_SP:
    nu, ec = V2/sp/"ngs45/nrDNA_45S.fasta", V2/sp/"easy45/consensus.fasta"
    ni, gb = V2/sp/"ngs45/its.fasta", V2REF/f"{sp}.its.fasta"
    nb, eb = flen(nu), flen(ec)
    row = dict(order=ORDER[sp], ngs45_bp=nb, easy45_bp=eb,
        ngs_vs_easy=blast_pid(nu, ec), its_vs_gb=blast_pid(ni, gb), gb=gb_acc(sp),
        outcome=("fail" if nb is None else "partial" if sp=="Musa_acuminata" else "full"))
    v2[sp] = row

# ---- v3 same-individual ----
v3 = {}
for sp in V3_SP:
    nu, ec = V3/sp/"ngs45/nrDNA_45S.fasta", V3/sp/"easy45/consensus.fasta"
    v3[sp] = dict(ngs45_bp=flen(nu), easy45_bp=flen(ec), same_pid=blast_pid(nu, ec))

# ---- QC (parse scratchpad qc.tsv if present, else skip) ----
QCSRC = Path("/tmp/claude-1065/-data06-users-vutrinh-ngs245/b3fc8574-ec91-4490-b6fc-b943d8228dd0/scratchpad/qc.tsv")
qc = {}  # (set,species,type) -> dict
if QCSRC.exists():
    for ln in QCSRC.read_text().splitlines()[1:]:
        p = ln.split("\t")
        if len(p) < 8 or p[0]=="QC" or p[4]=="avg_len":
            continue
        qc[(p[0],p[1],p[2])] = dict(reads=p[3], avglen=p[4], n50=p[5], q30=p[6], avgq=p[7])

def qcget(setname, sp, typ, key):
    d = qc.get((setname, sp, typ))
    return d[key] if d else "-"

# ---- write results_v2.tsv ----
with open(BENCH/"results_v2.tsv","w") as fh:
    fh.write("species\torder\tngs45_bp\teasy45_bp\tngs_vs_easy_pid\tits_vs_GB_pid\tGB_acc\toutcome\n")
    for sp in V2_SP:
        r = v2[sp]
        fh.write("\t".join([sp, r["order"], fmt(r["ngs45_bp"]), fmt(r["easy45_bp"]),
            fmt(r["ngs_vs_easy"]), fmt(r["its_vs_gb"]), r["gb"], r["outcome"]])+"\n")

# ---- write results_v3.tsv ----
with open(BENCH/"results_v3.tsv","w") as fh:
    fh.write("species\torder\tngs45_bp\teasy45_bp\tsame_indiv_pid\n")
    for sp in V3_SP:
        r = v3[sp]
        fh.write("\t".join([sp, ORDER[sp], fmt(r["ngs45_bp"]),
            fmt(r["easy45_bp"]) if r["easy45_bp"] else "no_consensus", fmt(r["same_pid"])])+"\n")

# ---- write master.tsv (drives figures) ----
with open(BENCH/"master.tsv","w") as fh:
    fh.write("species\torder\tngs_vs_easy_pid\toutcome\till_len\till_q30\thifi_n50\thifi_q30\tsame_indiv_pid\n")
    for sp in V2_SP:
        r = v2[sp]
        s3 = v3.get(sp, {})
        fh.write("\t".join([sp, r["order"], fmt(r["ngs_vs_easy"]), r["outcome"],
            qcget("S1_cross",sp,"Illumina","avglen"), qcget("S1_cross",sp,"Illumina","q30"),
            qcget("S1_cross",sp,"HiFi","n50"), qcget("S1_cross",sp,"HiFi","q30"),
            fmt(s3.get("same_pid"))])+"\n")

print("wrote results_v2.tsv, results_v3.tsv, master.tsv")
for f in ["results_v2.tsv","results_v3.tsv","master.tsv"]:
    print("\n==", f, "==")
    print((BENCH/f).read_text())
