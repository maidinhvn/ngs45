# Pipeline figure sources

Scripts that draw the workflow figures, kept here so the figures can be
regenerated and edited rather than only existing as images.

    python3 make_fig1.py     # both workflows, side by side (overview)
    python3 make_figS1.py    # easy45 (PacBio HiFi), full
    python3 make_figS2.py    # ngs45 (Illumina), full

`fig_style.py` holds the shared palette, box/arrow helpers and the output
settings. Each script writes four files next to itself:

| Output | Use |
|---|---|
| `.pdf` / `.eps` | vector, for submission or printing |
| `.tif` | 1200 dpi, LZW, RGB — journals that require raster line art |
| `.png` | 600 dpi, for previewing or dropping into a document |

Sized for a 175 mm double-column page, Arial-metric font (Liberation Sans),
minimum 7 pt type, and legible in greyscale.

The generated image files are not tracked; run a script to recreate them.
