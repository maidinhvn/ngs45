#!/usr/bin/env python3
"""Regenerate the benchmark figures from bench/master.tsv + bench/results_v3.tsv.
Outputs PNG (300 dpi) + SVG + PDF into bench/figures/ . Vector formats for the
manuscript; line-art only (no rasterised panels)."""
import csv
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BENCH = Path(__file__).resolve().parent
FIG = BENCH / "figures"; FIG.mkdir(exist_ok=True)
plt.rcParams.update({"font.size": 11, "svg.fonttype": "none", "axes.spines.top": False,
                     "axes.spines.right": False, "font.family": "DejaVu Sans"})
GREEN, ORANGE, GREY = "#2e7d32", "#ef6c00", "#9e9e9e"

def load(f):
    return list(csv.DictReader((BENCH/f).read_text().splitlines(), delimiter="\t"))

def save(fig, name):
    for ext in ("png", "svg", "pdf"):
        fig.savefig(FIG/f"{name}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)

def sci(s):  # Genus_species -> "G. species" italic label
    g, *rest = s.split("_")
    return f"{g[0]}. {' '.join(rest)}"

# ---------- Figure 1: cross-individual concordance ----------
m = load("master.tsv")
m.sort(key=lambda r: (r["outcome"] != "full", r["outcome"] != "partial",
                      -(float(r["ngs_vs_easy_pid"]) if r["ngs_vs_easy_pid"] != "-" else 0)))
easy_bp = {r["species"]: r["easy45_bp"] for r in load("results_v2.tsv")}
labels = [sci(r["species"]) for r in m]
FAILW = 99.12   # short stub width for "no unit" species (visually != 100%)
fig, ax = plt.subplots(figsize=(7.6, 5.0))
for i, r in enumerate(m):
    if r["ngs_vs_easy_pid"] == "-":
        ax.barh(i, FAILW-99.0, left=99.0, color=GREY, height=0.62)
        ax.text(FAILW+0.02, i, f"no unit — easy45 recovers ({easy_bp[r['species']]} bp)",
                va="center", ha="left", fontsize=8.5, color=GREY)
    else:
        v = float(r["ngs_vs_easy_pid"])
        c = GREEN if r["outcome"] == "full" else ORANGE
        ax.barh(i, v-99.0, left=99.0, color=c, height=0.62)
        ax.text(v+0.01, i, f"{v:.2f}"+("  (partial)" if r["outcome"]=="partial" else ""),
                va="center", fontsize=8.5)
ax.set_xlim(99.0, 100.35); ax.set_yticks(range(len(labels)))
ax.set_yticklabels(labels, style="italic", fontsize=10)
ax.invert_yaxis(); ax.set_xlabel("ngs45 unit identity to HiFi/easy45 consensus (%)")
ax.axvline(100.0, color="k", lw=0.6, ls=":")
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color=GREEN, label="full unit (8/12)"),
                   Patch(color=ORANGE, label="partial — hybrid (1/12)"),
                   Patch(color=GREY, label="ngs45 no unit → easy45 recovers (3/12)")],
          loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=3, fontsize=8.5, frameon=False)
ax.set_title("Cross-individual benchmark (12 species / 12 orders)", fontsize=11)
save(fig, "Figure1_concordance")

# ---------- Figure 2: same-individual identity ----------
v3 = [r for r in load("results_v3.tsv")]
labels2 = [sci(r["species"]) for r in v3]
vals2, cols2, txt = [], [], []
for r in v3:
    if r["same_indiv_pid"] == "-":
        vals2.append(100.0); cols2.append(GREY); txt.append("easy45 no consensus")
    else:
        v = float(r["same_indiv_pid"]); vals2.append(v); cols2.append(GREEN)
        txt.append("100.00 (0 mismatch)" if v == 100.0 else f"{v:.3f}")
fig, ax = plt.subplots(figsize=(7.2, 3.6))
y = range(len(labels2))
ax.barh(list(y), [v-99.5 for v in vals2], left=99.5, color=cols2, height=0.6)
ax.set_xlim(99.5, 100.1); ax.set_yticks(list(y)); ax.set_yticklabels(labels2, style="italic", fontsize=10)
ax.invert_yaxis(); ax.set_xlabel("ngs45 unit identity to easy45 consensus, SAME individual (%)")
ax.axvline(100.0, color="k", lw=0.6, ls=":")
for i, (v, t, c) in enumerate(zip(vals2, txt, cols2)):
    ax.text((99.52 if c == GREY else v+0.005), i, t, va="center", ha="left",
            fontsize=8.5, color=("white" if c == GREY else "k"),
            fontweight=("bold" if c == GREY else "normal"))
ax.set_title("Same-individual validation — no intraspecific confound", fontsize=11)
save(fig, "Figure2_same_individual")

# ---------- Figure 3: read length vs outcome ----------
fig, ax = plt.subplots(figsize=(7.6, 4.6))
oc_c = {"full": GREEN, "partial": ORANGE, "fail": "#c62828"}
# per-species label offset (dx,dy in points) + alignment, to avoid collisions
OFF = {"Actinidia_chinensis": (7, 9, "left"), "Fragaria_vesca": (6, -13, "left"),
       "Vitis_vinifera": (8, 2, "left"), "Helianthus_annuus": (8, -2, "left"),
       "Solanum_lycopersicum": (8, 4, "left"), "Musa_acuminata": (8, 2, "left")}
for r in m:
    x = float(r["ill_len"]); yv = float(r["ill_q30"])
    ax.scatter(x, yv, s=75, color=oc_c[r["outcome"]], edgecolor="k", lw=0.4, zorder=3)
    dx, dy, ha = OFF.get(r["species"], (5, 4, "left"))
    ax.annotate(sci(r["species"]), (x, yv), fontsize=7.5, style="italic", ha=ha,
                xytext=(dx, dy), textcoords="offset points")
ax.axvline(150, color=GREY, ls="--", lw=0.8)
ax.text(150, 100.6, "150 bp", fontsize=8, color=GREY, ha="center")
ax.set_ylim(77, 101.5); ax.set_xlim(72, 262)
ax.set_xlabel("Illumina read length (bp)"); ax.set_ylabel("Illumina bases ≥ Q30 (%)")
ax.legend(handles=[Patch(color=GREEN, label="full"), Patch(color=ORANGE, label="partial"),
                   Patch(color="#c62828", label="fail (short reads can't span)")],
          loc="lower right", fontsize=8.5, frameon=False)
ax.set_title("Recovery vs read length & quality", fontsize=11)
save(fig, "Figure3_readlength")

print("figures written to", FIG)
for p in sorted(FIG.glob("*.png")):
    print(" ", p.name)
