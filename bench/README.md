# ngs45 benchmark (2026-07) — data & scripts

Re-evaluation of **ngs45 v0.2.0** against **easy45 / PacBio HiFi** and GenBank. The v1
artefacts were removed; this folder holds the current source data, scripts and
figures. Narrative + tables: [`../docs/BENCHMARK.md`](../docs/BENCHMARK.md); figures:
[`../docs/FIGURES.md`](../docs/FIGURES.md); accessions: [`../docs/DATA_ACCESSIONS.md`](../docs/DATA_ACCESSIONS.md).

## Panels

- **Cross-individual** (12 species / 12 angiosperm orders): Illumina and HiFi from
  *different* individuals — [`manifest_v2.tsv`](manifest_v2.tsv),
  results in [`results_v2.tsv`](results_v2.tsv).
- **Same-individual** (7 species): HiFi + Illumina from the *same* BioSample —
  results in [`results_v3.tsv`](results_v3.tsv). Removes the intraspecific confound.
- **Application** (8 *Polyscias* samples, in-house): run separately via `ngs45 batch`;
  validated against GenBank ITS. Exercises the v0.2.0 tandem-dup QC (*P. cumingiana*
  26S 5944 → 5817 bp).

## Files

| File | Contents |
|---|---|
| `manifest_v2.tsv` | cross-individual accessions (Illumina + HiFi per species) |
| `results_v2.tsv` | cross-individual: ngs45 vs easy45 vs GenBank ITS + outcome |
| `results_v3.tsv` | same-individual: ngs45 vs easy45 identity |
| `qc.tsv` | per-dataset QC (read length, Q30, N50) from `seqkit stats -a` |
| `master.tsv` | merged per-species table that drives the figures |
| `consolidate.py` | rebuilds results_v2/v3, qc, master from the run outputs |
| `make_figures.py` | renders `figures/Figure{1,2,3}.*` (PNG 300 dpi + SVG + PDF) |

## Reproduce the tables & figures

```bash
python bench/consolidate.py     # -> results_v2.tsv, results_v3.tsv, master.tsv
python bench/make_figures.py    # -> bench/figures/
```

## Headline

Same-individual: 6/7 concordant, 5 base-identical (0 mismatch), mean ≈ 99.98 %.
Cross-individual: 8/12 full recovery at 99.76–100 % vs HiFi (9/12 with ≥150 bp reads),
1 partial (*Musa*), 2 short-read limits (*Helianthus*, *Vitis*) that easy45 recovers.
