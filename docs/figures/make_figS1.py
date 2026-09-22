"""Supplementary Figure S1 — the easy45 workflow (PacBio HiFi) in full.

Simplified relative to the earlier version: single linear flow (no crossing or
dashed connectors), no internal stage codes, no script or output file names, no
command-line flags. Terminology follows the manuscript ("core sequence",
"sequence variants"); the earlier hybrid / allopolyploidy wording is dropped.
Dashed outlines mark optional steps.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fig_style import (new_canvas, box, plain, arrow, save,
                       INPUT_FC, INPUT_EC, KEY_FC, KEY_EC, OUT_FC, OUT_EC)

fig, ax = new_canvas(175, 184)

CX, BW, BH = 0.50, 0.66, 0.058
FS = 7.3

# ---------------- inputs ----------------
ax.text(0.035, 0.988, "Input", ha="left", va="center", fontsize=7.0,
        color="#555555", fontweight="bold")
for x, t in ((0.20, "PacBio HiFi reads"),
             (0.50, "Bundled 45S nrDNA anchor\n($\\it{Arabidopsis}$, GenBank OR453402)"),
             (0.80, "Rfam SSU and LSU\ncovariance models")):
    box(ax, x, 0.950, 0.275, 0.044, t, fc=INPUT_FC, ec=INPUT_EC, fs=7.0)

STEPS = [
    ("Deplete organelle reads  (minimap2)", "dashed"),
    ("Recruit reads carrying 45S nrDNA by mapping to the anchor  (minimap2)", "solid"),
    ("Locate 18S, 5.8S and 26S in each read and keep the reads that span\n"
     "the whole region  (barrnap)", "solid"),
    ("Excise one copy per spanning read and orient it 18S → 26S", "solid"),
    ("Group the excised copies by sequence similarity  (VSEARCH, 97 % identity)", "solid"),
    ("Build a consensus for each group and trim it to the 18S and 26S\n"
     "boundaries  (abPOA, Infernal)", "key"),
    ("Screen candidate sequence variants: keep read-supported substitutions,\n"
     "discard homopolymer and short-tandem-repeat indels, exclude sequences\n"
     "less than 90 % identical to the primary consensus", "solid"),
    ("Recover the intergenic spacer from reads spanning 26S to the next 18S", "dashed"),
    ("Annotate the component regions and extract the ITS barcode  (ITSx)", "solid"),
]

Y0, DY = 0.878, 0.0855
ys = [Y0 - i * DY for i in range(len(STEPS))]

arrow(ax, CX, 0.928, ys[0] + BH / 2)
for i, ((txt, kind), y) in enumerate(zip(STEPS, ys)):
    h = BH if txt.count("\n") < 2 else BH + 0.016
    if kind == "key":
        box(ax, CX, y, BW, h, txt, fc=KEY_FC, ec=KEY_EC, fs=FS)
    else:
        box(ax, CX, y, BW, h, txt, fs=FS, ls=kind)
    if i < len(STEPS) - 1:
        arrow(ax, CX, y - h / 2, ys[i + 1] + BH / 2)

ax.text(CX + BW / 2 + 0.012, ys[0], "optional", ha="left", va="center",
        fontsize=7.0, style="italic", color="#666666")
ax.text(CX + BW / 2 + 0.012, ys[7], "optional", ha="left", va="center",
        fontsize=7.0, style="italic", color="#666666")

# ---------------- outputs ----------------
arrow(ax, CX, ys[-1] - BH / 2, 0.122)
ax.text(0.035, 0.108, "Output", ha="left", va="center", fontsize=7.0,
        color="#555555", fontweight="bold")
for x, t in ((0.155, "45S nrDNA core\nsequence"),
             (0.385, "Candidate sequence\nvariants"),
             (0.615, "ITS barcode and\nper-region FASTA"),
             (0.845, "Intergenic spacer,\nannotation, run report")):
    box(ax, x, 0.068, 0.215, 0.062, t, fc=OUT_FC, ec=OUT_EC, fs=7.0)

save(fig, os.path.join(os.path.dirname(os.path.abspath(__file__)), "FigureS1"))
