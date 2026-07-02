#!/usr/bin/env python3
"""Generate manuscript figures for the ngs45 benchmark.

Reads bench/results_summary.tsv (+ benchmark_modern.tsv) and writes PNG+SVG to
bench/figures/.
"""
import csv, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures"); os.makedirs(FIG, exist_ok=True)

rows = list(csv.DictReader(open(os.path.join(HERE, "results_summary.tsv")), delimiter="\t"))
def f(x):
    try: return float(x)
    except: return None
for r in rows:
    r["short"] = r["species"].split("_")[0]
    r["rl_old"] = f(r["readlen_old"]); r["rl_mod"] = f(r["readlen_modern"])
    r["pid"] = f(r["unitPID"]); r["ribo"] = f(r["ribo_sites"])
    r["ok_old"] = r["ngs45_old"] == "Y"; r["ok_mod"] = r["ngs45_modern"] == "Y"

GREEN, RED, BLUE, GREY = "#2ca02c", "#d62728", "#1f77b4", "#888888"

def save(fig, name):
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("wrote", name)

# ---- Fig 1: read length determines short-read success --------------------
fig, ax = plt.subplots(figsize=(7, 4.5))
FAILY = 99.30  # dedicated row for "not assembled"
for r in rows:
    for rl, ok, mk, lab in [(r["rl_old"], r["ok_old"], "o", "old"),
                            (r["rl_mod"], r["ok_mod"], "s", "modern")]:
        if rl is None: continue
        if ok and r["pid"] is not None:
            ax.scatter(rl, r["pid"], c=GREEN, marker=mk, s=70, edgecolor="k", lw=.5, zorder=3)
        else:
            ax.scatter(rl, FAILY, c=RED, marker="x", s=70, zorder=3)
# arrows old->modern for rescued species
for r in rows:
    if r["rl_old"] and r["rl_mod"] and not r["ok_old"] and r["ok_mod"] and r["pid"]:
        ax.annotate("", xy=(r["rl_mod"], r["pid"]), xytext=(r["rl_old"], FAILY),
                    arrowprops=dict(arrowstyle="->", color=GREY, lw=1, alpha=.6))
ax.axvline(150, ls="--", c=BLUE, lw=1.2)
ax.text(152, 99.55, "150 bp\nthreshold", color=BLUE, fontsize=8, va="center")
# annotate Musa (fails even at 309 bp) — short label above its point, no long arrow
ax.annotate("Musa: fails even at 309 bp\n(high rDNA heterozygosity)",
            xy=(309, FAILY), xytext=(305, 99.40), fontsize=7.5, color=RED, ha="right",
            arrowprops=dict(arrowstyle="->", color=RED, lw=.8))
ax.set_xlabel("Illumina read length (bp)"); ax.set_ylabel("ngs45 unit identity to HiFi (%)")
ax.set_ylim(99.25, 100.1); ax.set_xlim(60, 330); ax.set_yticks([99.3, 99.5, 99.75, 100.0])
ax.set_yticklabels(["not\nassembled", "99.5", "99.75", "100"])
ax.set_title("Read length ≥150 bp enables short-read 45S assembly")
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([],[],marker="o",color="w",mfc=GREEN,mec="k",label="old run"),
                   Line2D([],[],marker="s",color="w",mfc=GREEN,mec="k",label="modern run"),
                   Line2D([],[],marker="x",color=RED,ls="",label="not assembled")],
          fontsize=8, loc="center right", bbox_to_anchor=(1.0, 0.62))
save(fig, "Figure1_readlength")

# ---- Fig 2: success rate old vs modern vs HiFi ---------------------------
fig, ax = plt.subplots(figsize=(4.6, 4))
n_old = sum(r["ok_old"] for r in rows); n_mod = sum(r["ok_mod"] for r in rows)
bars = ax.bar(["ngs45\n(old <150 bp)", "ngs45\n(modern ≥150 bp)", "easy45\n(HiFi)"],
              [n_old, n_mod, 10], color=[RED, GREEN, BLUE], edgecolor="k")
for b, v in zip(bars, [n_old, n_mod, 10]):
    ax.text(b.get_x()+b.get_width()/2, v+.15, f"{v}/10", ha="center", fontweight="bold")
ax.set_ylabel("species with full 45S unit recovered"); ax.set_ylim(0, 11)
ax.set_title("Recovery rate (10 species, 9 orders)")
save(fig, "Figure2_success_rate")

# ---- Fig 3: ngs45 (short-read) vs easy45 (HiFi) concordance --------------
ok = [r for r in rows if r["ok_mod"] and r["pid"] is not None]
ok.sort(key=lambda r: r["pid"])
fig, ax = plt.subplots(figsize=(6.5, 4))
cols = [GREEN if r["pid"] >= 99.9 else "#66bd63" for r in ok]
ax.barh([r["short"] for r in ok], [r["pid"] for r in ok], color=cols, edgecolor="k")
for i, r in enumerate(ok):
    ax.text(r["pid"]-.02, i, f'{r["pid"]:.2f}%', va="center", ha="right", fontsize=7.5, color="w")
ax.set_xlim(99.5, 100.02); ax.axvline(100, ls=":", c=GREY)
ax.set_xlabel("ngs45 unit identity to HiFi consensus (%)")
ax.set_title("Short-read ↔ HiFi concordance (modern input)")
save(fig, "Figure3_concordance")

# ---- Fig 4: ribotype heterozygosity vs consensus identity ----------------
pts = [r for r in rows if r["ok_mod"] and r["ribo"] is not None and r["pid"] is not None]
fig, ax = plt.subplots(figsize=(6, 4.3))
xs = [r["ribo"] for r in pts]; ys = [r["pid"] for r in pts]
ax.scatter(xs, ys, c=BLUE, s=70, edgecolor="k", zorder=3)
# stagger labels so near-identical points (e.g. Beta/Glycine at 100%) don't overlap
_off = {"Beta": (5, 4), "Glycine": (5, -12), "Solanum": (-8, 6), "Sesamum": (6, -12)}
for r in pts:
    dx, dy = _off.get(r["short"], (6, 3))
    ax.annotate(r["short"], (r["ribo"], r["pid"]), textcoords="offset points",
                xytext=(dx, dy), fontsize=7.5)
# trend line
if len(xs) > 2:
    import numpy as np
    b, a = np.polyfit(xs, ys, 1)
    xr = np.array([min(xs), max(xs)])
    ax.plot(xr, a + b*xr, ls="--", c=RED, lw=1, label=f"trend (slope {b:.3f}%/site)")
    ax.legend(fontsize=8)
ax.set_xlabel("ribotype_sites (intragenomic heterozygosity)")
ax.set_ylabel("ngs45 unit identity to HiFi (%)")
ax.set_title("More ribotype heterogeneity → blended consensus, lower identity")
save(fig, "Figure4_ribotype")

print("all figures ->", FIG)
