"""Figure 1 — overview of easy45 and ngs45 (main text).

Deliberately minimal: four steps per workflow, no tool flags, no internal stage
codes or output file names (those live in Supplementary Figs S1-S2).
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fig_style import (new_canvas, box, plain, arrow, panel_frame, save,
                       INPUT_FC, INPUT_EC, KEY_FC, KEY_EC, OUT_FC, OUT_EC)

fig, ax = new_canvas(175, 104)          # double-column width

PA, PB = (0.020, 0.487), (0.513, 0.980)
CA, CB = sum(PA) / 2, sum(PB) / 2
PANEL_TOP, PANEL_BOT = 0.975, 0.352
BW, BH = 0.415, 0.070                    # step box width / height

panel_frame(ax, PA[0], PANEL_BOT, PA[1], PANEL_TOP, "A", "easy45 – PacBio HiFi")
panel_frame(ax, PB[0], PANEL_BOT, PB[1], PANEL_TOP, "B", "ngs45 – Illumina paired-end")

# one-line rationale for why the two strategies differ
ax.text(CA, 0.902, "a single read spans the whole region → no assembly",
        ha="center", va="center", fontsize=7.0, style="italic", color="#444444")
ax.text(CB, 0.902, "reads are too short to span it → assemble, then resolve",
        ha="center", va="center", fontsize=7.0, style="italic", color="#444444")

Y_IN = 0.838
Y = [0.726, 0.618, 0.510, 0.402]

box(ax, CA, Y_IN, BW, 0.060, "PacBio HiFi reads", fc=INPUT_FC, ec=INPUT_EC, fs=7.6, weight="bold")
box(ax, CB, Y_IN, BW, 0.060, "Illumina paired-end reads", fc=INPUT_FC, ec=INPUT_EC, fs=7.6, weight="bold")

steps_a = [
    "Recruit reads carrying 45S nrDNA\n(minimap2)",
    "Locate the rRNA genes and keep reads\nspanning 18S to 26S (barrnap)",
    "Cluster the spanning units and build\na consensus (VSEARCH, abPOA)",
]
steps_b = [
    "Recruit rDNA reads by iterative\nmapping (Bowtie 2)",
    "Assemble the recruited reads\n(SPAdes, coverage-capped)",
    "Resolve a single repeat unit\n(BLAST+)",
]
shared_step = "Trim to the 18S and 26S boundaries;\nannotate ITS (Infernal, ITSx)"

for i, (ta, tb) in enumerate(zip(steps_a, steps_b)):
    box(ax, CA, Y[i], BW, BH, ta)
    box(ax, CB, Y[i], BW, BH, tb)
# the final step is identical in both workflows -> highlighted
box(ax, CA, Y[3], BW, BH, shared_step, fc=KEY_FC, ec=KEY_EC)
box(ax, CB, Y[3], BW, BH, shared_step, fc=KEY_FC, ec=KEY_EC)

for c in (CA, CB):
    arrow(ax, c, Y_IN - 0.030, Y[0] + BH / 2)
    for i in range(3):
        arrow(ax, c, Y[i] - BH / 2, Y[i + 1] + BH / 2)
    arrow(ax, c, Y[3] - BH / 2, 0.273)

# shared foundation + shared output
plain(ax, 0.5, 0.232, 0.960, 0.074,
      "Same bundled $\\it{Arabidopsis}$ 45S anchor and Rfam SSU/LSU covariance models\n"
      "→ both workflows define the 18S and 26S boundaries identically", fs=7.3)
arrow(ax, 0.5, 0.193, 0.170)

box(ax, 0.5, 0.108, 0.960, 0.108,
    "Standardized 45S nrDNA core sequence  (18S–ITS1–5.8S–ITS2–26S)\n"
    "ITS barcode, FASTA for each component region, run report\n"
    "optional: variable sites within the rDNA copies",
    fc=OUT_FC, ec=OUT_EC, fs=7.5)

save(fig, os.path.join(os.path.dirname(os.path.abspath(__file__)), "Figure1"))
