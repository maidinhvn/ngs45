#!/usr/bin/env python3
"""Publication-quality pipeline diagram for ngs45 (PNG 600 dpi + SVG)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import os

plt.rcParams.update({"font.family": "DejaVu Sans", "svg.fonttype": "none"})
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "docs", "figures")

# colours
CORE  = "#dbeafe"; CORE_E = "#1d4ed8"      # blue fill / edge
OPT   = "#f1f5f9"; OPT_E  = "#64748b"      # grey (optional)
IO    = "#1e293b"                          # dark slate (input/output)
TXT   = "#0f172a"
CHIP  = "#1e3a8a"

# (id, title, desc, tool, kind)  kind: io / core / opt
rows = [
    ("", "Illumina paired-end reads", "+ bundled 45S seed (Arabidopsis) · Rfam SSU/LSU CMs", "", "io"),
    ("S0", "QC", "quality & adapter trimming", "cutadapt", "opt"),
    ("S1", "Bait", "iterative recruitment of rDNA reads against the 45S seed", "bowtie2", "core"),
    ("S2", "Assemble", "multi-k de Bruijn assembly  (k < read length, no --careful)", "SPAdes", "core"),
    ("S3", "Resolve", "locate the rDNA contig/scaffold → cut one repeat unit", "BLAST", "core"),
    ("S4", "Boundary", "orient + trim to the mature transcribed unit (18S 5′ – 26S 3′)", "infernal · Rfam CM", "core"),
    ("S5", "Annotate", "delimit 18S · ITS1 · 5.8S · ITS2 · 26S  +  ITS barcode", "ITSx", "core"),
    ("S6", "Variants", "map reads back → ribotype-heterozygosity sites", "bwa · samtools · bcftools", "opt"),
    ("S7", "Report", "summary table + human-readable report", "", "core"),
    ("", "Outputs", "nrDNA_45S.fasta · annotation.gff3 · its.fasta · ribotype_variants.tsv · summary/report", "", "io"),
]

n = len(rows)
BH, GAP = 1.0, 0.62                 # box height, gap
W, X = 8.6, 0.7                     # box width, left x
top = n * (BH + GAP)
fig_h = top + 0.4
fig, ax = plt.subplots(figsize=(7.4, fig_h))
ax.set_xlim(0, 10); ax.set_ylim(0, top + 0.2); ax.axis("off")

centers = []
for i, (sid, title, desc, tool, kind) in enumerate(rows):
    y = top - (i + 1) * (BH + GAP) + GAP
    cx, cy = X + W / 2, y + BH / 2
    centers.append((cx, y, y + BH))
    if kind == "io":
        box = FancyBboxPatch((X, y), W, BH, boxstyle="round,pad=0.02,rounding_size=0.12",
                             fc=IO, ec="none", zorder=2)
        ax.add_patch(box)
        ax.text(cx, cy + 0.17, title, color="white", ha="center", va="center",
                fontsize=13, fontweight="bold")
        ax.text(cx, cy - 0.20, desc, color="#cbd5e1", ha="center", va="center", fontsize=9)
    else:
        fc, ec = (CORE, CORE_E) if kind == "core" else (OPT, OPT_E)
        ls = "-" if kind == "core" else (0, (4, 2))
        box = FancyBboxPatch((X, y), W, BH, boxstyle="round,pad=0.02,rounding_size=0.10",
                             fc=fc, ec=ec, lw=1.6, ls=ls, zorder=2)
        ax.add_patch(box)
        # stage badge
        ax.add_patch(FancyBboxPatch((X + 0.18, y + 0.20), 0.72, BH - 0.40,
                     boxstyle="round,pad=0.01,rounding_size=0.08",
                     fc=ec, ec="none", zorder=3))
        ax.text(X + 0.54, cy, sid, color="white", ha="center", va="center",
                fontsize=12.5, fontweight="bold", zorder=4)
        # title + description
        ax.text(X + 1.15, cy + 0.19, title, color=TXT, ha="left", va="center",
                fontsize=13, fontweight="bold")
        ax.text(X + 1.15, cy - 0.20, desc, color="#334155", ha="left", va="center", fontsize=9.3)
        # tool chip (right)
        if tool:
            ax.text(X + W - 0.22, cy, tool, color=CHIP, ha="right", va="center",
                    fontsize=8.8, style="italic")
        if kind == "opt":
            ax.text(X + W - 0.22, y + 0.16, "optional", color=OPT_E, ha="right",
                    va="center", fontsize=7.2, style="italic")

# arrows between consecutive boxes
for i in range(n - 1):
    _, ytop_i, _ = centers[i]              # unused
    cx = X + W / 2
    y_from = centers[i][1]                 # bottom of box i
    y_to = centers[i + 1][2]               # top of box i+1
    ax.add_patch(FancyArrowPatch((cx, y_from), (cx, y_to),
                 arrowstyle="-|>", mutation_scale=16, lw=1.8, color="#475569", zorder=1))

ax.text(X, top + 0.02, "ngs45 pipeline", ha="left", va="bottom",
        fontsize=15, fontweight="bold", color=TXT)
ax.text(X + W, top + 0.05, "45S nrDNA from Illumina short reads", ha="right",
        va="bottom", fontsize=9.5, color="#475569", style="italic")

fig.tight_layout(pad=0.4)
for ext, dpi in (("png", 600), ("svg", 600)):
    fig.savefig(os.path.join(OUT, f"pipeline.{ext}"), dpi=dpi, bbox_inches="tight")
print("wrote docs/figures/pipeline.png (600 dpi) + pipeline.svg")
