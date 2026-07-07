#!/usr/bin/env python3
"""Dotplot: cross-individual Illumina input size vs ngs45 processing time."""
import csv
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

BENCH = Path(__file__).resolve().parent
FIG = BENCH / "figures"; FIG.mkdir(exist_ok=True)
plt.rcParams.update({"font.size": 11, "svg.fonttype": "none", "axes.spines.top": False,
                     "axes.spines.right": False, "font.family": "DejaVu Sans"})
GREEN, ORANGE, RED = "#2e7d32", "#ef6c00", "#c62828"

# ngs45 wall-clock seconds (definitive successful/final run; Beta has no v2 timing)
TIME = {"Helianthus_annuus":323,"Solanum_lycopersicum":313,"Sesamum_indicum":553,
        "Citrus_sinensis":984,"Glycine_max":351,"Fragaria_vesca":266,
        "Populus_trichocarpa":1073,"Vitis_vinifera":373,"Actinidia_chinensis":169,
        "Oryza_sativa":853,"Musa_acuminata":243}
OUTCOME = {r[0]: r[7] for r in csv.reader((BENCH/"results_v2.tsv").read_text().splitlines(), delimiter="\t")}

# Illumina input as on-disk file size (R1+R2 gzipped), in GB
RAW = Path("/path/to/benchmark_v2_data/0_rawdata")
qc = {}
for sp in TIME:
    b = 0
    for m in ("R1", "R2"):
        f = RAW / sp / f"{sp}_{m}.fastq.gz"
        if f.exists():
            b += f.stat().st_size
    if b:
        qc[sp] = b / 1e9   # GB (R1+R2 gzipped)

col = {"full":GREEN,"partial":ORANGE,"fail":RED}
def sci(s):
    p = s.split("_"); return f"{p[0][0]}. {' '.join(p[1:])}"

# per-species label offset (dx,dy pt, ha) to declutter
OFF = {"Solanum_lycopersicum": (0, 10, "center"), "Helianthus_annuus": (8, 0, "left"),
       "Fragaria_vesca": (8, -5, "left"), "Glycine_max": (0, -13, "center"),
       "Vitis_vinifera": (-6, 8, "right"), "Musa_acuminata": (8, 0, "left")}
fig, ax = plt.subplots(figsize=(9.6, 5.2))
xs, ys = [], []
for sp, sec in TIME.items():
    if sp not in qc:
        continue
    x, y = qc[sp], sec/60.0
    xs.append(x); ys.append(y)
    oc = OUTCOME.get(sp, "full")
    ax.scatter(x, y, s=90, color=col[oc], edgecolor="k", lw=0.5, zorder=3)
    dx, dy, ha = OFF.get(sp, (6, 3, "left"))
    ax.annotate(sci(sp), (x, y), fontsize=8, style="italic", ha=ha,
                xytext=(dx, dy), textcoords="offset points")
ax.set_xlabel("Illumina input file size (GB, R1+R2 gzipped)")
ax.set_ylabel("ngs45 processing time (min)")
ax.set_title("Cross-individual: data size vs ngs45 runtime", fontsize=11)
ax.legend(handles=[Patch(color=GREEN,label="full unit"),
                   Patch(color=ORANGE,label="partial"),
                   Patch(color=RED,label="fail (fast, at S3)")],
          loc="upper left", fontsize=9, frameon=False)
ax.margins(0.12)
for ext in ("png","svg","pdf"):
    fig.savefig(FIG/f"Figure4_size_vs_time.{ext}", dpi=300, bbox_inches="tight")
plt.close(fig)
print("wrote Figure4_size_vs_time; points:", len(xs))
print("Gbp range:", round(min(xs),2), "-", round(max(xs),2), "| min range:", round(min(ys),1), "-", round(max(ys),1))
