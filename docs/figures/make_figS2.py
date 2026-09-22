"""Supplementary Figure S2 — the ngs45 workflow (Illumina) in full.

Same simplification as Figure S1: one linear flow, no internal stage codes, no
script or output file names, no command-line flags, and manuscript terminology
("core sequence", "variable sites") in place of the earlier ribotype wording.
Dashed outlines mark optional steps.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fig_style import (new_canvas, box, arrow, save,
                       INPUT_FC, INPUT_EC, KEY_FC, KEY_EC, OUT_FC, OUT_EC)

fig, ax = new_canvas(175, 150)

CX, BW, BH = 0.50, 0.66, 0.058
FS = 7.3

# ---------------- inputs ----------------
ax.text(0.035, 0.982, "Input", ha="left", va="center", fontsize=7.0,
        color="#555555", fontweight="bold")
for x, t in ((0.20, "Illumina paired-end reads"),
             (0.50, "Bundled 45S nrDNA anchor\n($\\it{Arabidopsis}$, GenBank OR453402)"),
             (0.80, "Rfam SSU and LSU\ncovariance models")):
    box(ax, x, 0.941, 0.275, 0.044, t, fc=INPUT_FC, ec=INPUT_EC, fs=7.0)

STEPS = [
    ("Trim adapters and low-quality bases  (cutadapt)", "dashed"),
    ("Recruit rDNA reads by mapping to the anchor, repeating until the\n"
     "recruited set stops growing  (Bowtie 2)", "solid"),
    ("Assemble the recruited reads over several k-mer sizes, with the\n"
     "coverage passed to the assembler capped  (SPAdes)", "solid"),
    ("Identify the contig or scaffold carrying a complete repeat and\n"
     "excise a single copy  (BLAST+)", "solid"),
    ("Orient the copy, trim it to the 18S and 26S boundaries and remove\n"
     "tandem-duplication artefacts  (Infernal)", "key"),
    ("Annotate the component regions and extract the ITS barcode  (ITSx)", "solid"),
    ("Map the recruited reads back to the core sequence and flag variable\n"
     "sites (minor allele ≥ 10 % at depth ≥ 20)  (BWA-MEM)", "dashed"),
]

Y0, DY = 0.868, 0.108
ys = [Y0 - i * DY for i in range(len(STEPS))]

arrow(ax, CX, 0.919, ys[0] + BH / 2)
for i, ((txt, kind), y) in enumerate(zip(STEPS, ys)):
    h = BH if "\n" not in txt else BH + 0.010
    if kind == "key":
        box(ax, CX, y, BW, h, txt, fc=KEY_FC, ec=KEY_EC, fs=FS)
    else:
        box(ax, CX, y, BW, h, txt, fs=FS, ls=kind)
    if i < len(STEPS) - 1:
        arrow(ax, CX, y - h / 2, ys[i + 1] + BH / 2 + 0.005)

for idx in (0, 6):
    ax.text(CX + BW / 2 + 0.012, ys[idx], "optional", ha="left", va="center",
            fontsize=7.0, style="italic", color="#666666")

# ---------------- outputs ----------------
arrow(ax, CX, ys[-1] - BH / 2 - 0.005, 0.125)
ax.text(0.035, 0.112, "Output", ha="left", va="center", fontsize=7.0,
        color="#555555", fontweight="bold")
for x, t in ((0.155, "45S nrDNA core\nsequence"),
             (0.385, "ITS barcode and\nper-region FASTA"),
             (0.615, "Variable sites\n(optional)"),
             (0.845, "Annotation and\nrun report")):
    box(ax, x, 0.070, 0.215, 0.062, t, fc=OUT_FC, ec=OUT_EC, fs=7.0)

save(fig, os.path.join(os.path.dirname(os.path.abspath(__file__)), "FigureS2"))
