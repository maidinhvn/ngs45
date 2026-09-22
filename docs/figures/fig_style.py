"""Shared drawing helpers for the easy45 / ngs45 manuscript figures.

Journal-quality output: vector PDF + 600 dpi PNG, sans-serif, print- and
grayscale-safe fills, no internal stage codes / file names / CLI flags.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle

# --- typography -------------------------------------------------------------
# Liberation Sans is metrically identical to Arial (the usual journal house
# font), so the layout matches an Arial-set figure exactly; Nimbus Sans is the
# Helvetica equivalent. Real Arial/Helvetica are used if the system has them.
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans",
                        "Nimbus Sans", "DejaVu Sans"],
    "pdf.fonttype": 42,      # embed TrueType (text stays selectable/editable)
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    # mathtext (used for the italic species names) otherwise silently falls back
    # to DejaVu, which would mix two typefaces in one figure.
    "mathtext.fontset": "custom",
    "mathtext.rm": "Liberation Sans",
    "mathtext.it": "Liberation Sans:italic",
    "mathtext.bf": "Liberation Sans:bold",
    "mathtext.cal": "Liberation Sans:italic",
    "mathtext.sf": "Liberation Sans",
    "mathtext.tt": "Liberation Sans",
})

MIN_PT = 7.0                 # journal minimum font size at final printed size

# --- palette (light fills keep text black and legible in grayscale) ---------
INPUT_FC, INPUT_EC = "#dce8f6", "#5b7fa6"   # blue-grey  : data in
STEP_FC,  STEP_EC  = "#fdf4dd", "#c08a2e"   # warm cream : processing
KEY_FC,   KEY_EC   = "#f6d79a", "#b9761a"   # amber      : the step that differs
OUT_FC,   OUT_EC   = "#e2f0e2", "#4f8a53"   # green      : results
BAND_FC,  BAND_EC  = "#eef2f7", "#8fa4bd"   # shared foundation band
ARROW_C            = "#8a6a2f"


def box(ax, x, y, w, h, text, fc=STEP_FC, ec=STEP_EC, fs=7.4, weight="normal",
        radius=0.012, lw=0.9, ha="center", ls="solid"):
    """Rounded box with centred, wrapped text. (x, y) = centre.

    ls="dashed" marks an optional step.
    """
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                                boxstyle=f"round,pad=0,rounding_size={radius}",
                                linewidth=lw, facecolor=fc, edgecolor=ec,
                                linestyle=ls, zorder=2))
    ax.text(x, y, text, ha=ha, va="center", fontsize=fs, zorder=3,
            fontweight=weight, linespacing=1.35, color="black")


def plain(ax, x, y, w, h, text, fc=BAND_FC, ec=BAND_EC, fs=7.2, weight="normal", lw=0.9):
    """Square-cornered band (shared foundation / notes)."""
    ax.add_patch(Rectangle((x - w / 2, y - h / 2), w, h, linewidth=lw,
                           facecolor=fc, edgecolor=ec, zorder=2))
    ax.text(x, y, text, ha="center", va="center", fontsize=fs, zorder=3,
            fontweight=weight, linespacing=1.4, color="black")


def arrow(ax, x, y0, y1, c=ARROW_C, lw=1.15):
    """Vertical arrow from y0 down to y1."""
    ax.add_patch(FancyArrowPatch((x, y0), (x, y1), arrowstyle="-|>",
                                 mutation_scale=9, linewidth=lw,
                                 color=c, zorder=1, shrinkA=0, shrinkB=0))


def elbow(ax, x0, y0, x1, y1, c=ARROW_C, lw=1.15):
    """Right-angle connector: down from (x0,y0), across, then into (x1,y1)."""
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>",
                                 mutation_scale=9, linewidth=lw, color=c, zorder=1,
                                 connectionstyle="angle,angleA=-90,angleB=180,rad=0",
                                 shrinkA=0, shrinkB=0))


def panel_frame(ax, x0, y0, x1, y1, label, title, fs_label=9.5, fs_title=8.6):
    """Panel outline with an 'A  title' header."""
    ax.add_patch(FancyBboxPatch((x0, y0), x1 - x0, y1 - y0,
                                boxstyle="round,pad=0,rounding_size=0.012",
                                linewidth=0.9, facecolor="white",
                                edgecolor="#9aa5b1", zorder=0))
    ax.text(x0 + 0.018, y1 - 0.035, label, ha="left", va="center",
            fontsize=fs_label, fontweight="bold", color="black", zorder=3)
    ax.text(x0 + 0.055, y1 - 0.035, title, ha="left", va="center",
            fontsize=fs_title, fontweight="bold", color="black", zorder=3)


def new_canvas(w_mm, h_mm):
    """Figure sized in millimetres (Bioinformatics: 175 mm double column)."""
    fig = plt.figure(figsize=(w_mm / 25.4, h_mm / 25.4))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    return fig, ax


def save(fig, stem):
    """Write the full publication set.

    PDF / EPS : vector — the preferred submission format for line art
    TIFF      : 1200 dpi, LZW — for journals that require raster line art
    PNG       : 600 dpi — for placing in the Word manuscript / review
    """
    fig.savefig(f"{stem}.pdf", facecolor="white")
    fig.savefig(f"{stem}.eps", facecolor="white")
    fig.savefig(f"{stem}.tif", dpi=1200, facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(f"{stem}.png", dpi=600, facecolor="white")
    plt.close(fig)

    # Flatten the TIFF to RGB on white: an alpha channel is not accepted by
    # most production systems and can print as an unexpected background.
    from PIL import Image
    Image.MAX_IMAGE_PIXELS = None
    with Image.open(f"{stem}.tif") as im:
        if im.mode in ("RGBA", "LA", "P"):
            im = im.convert("RGBA")
            flat = Image.new("RGB", im.size, "white")
            flat.paste(im, mask=im.split()[-1])
            flat.save(f"{stem}.tif", compression="tiff_lzw", dpi=(1200, 1200))
    import os
    sizes = {e: os.path.getsize(f"{stem}.{e}") / 1e6 for e in ("pdf", "eps", "tif", "png")}
    print(f"  {os.path.basename(stem)}: "
          + ", ".join(f"{e} {s:.1f} MB" for e, s in sizes.items()))
