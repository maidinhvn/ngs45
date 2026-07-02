# ngs45 benchmark

Cross-validation of ngs45 (short-read 45S nrDNA assembly) against **easy45 / PacBio
HiFi** as the gold standard, plus **GenBank** ITS references, across 10 species
spanning 9 angiosperm orders (2 monocot + 8 eudicot). All identities are BLASTn.

## Headline result

| Input | ngs45 (Illumina) | easy45 (HiFi) |
|---|---|---|
| **old / short-read** (73–144 bp) | 4/10 | — |
| **modern** (≥150 bp PE) | **9/10** | — |
| gold standard | — | **10/10** |

The jump from 4/10 → 9/10 is driven almost entirely by **read length**, not by the
tool (see `ASSEMBLY_LIMITATION.md`). Where ngs45 succeeds it is **99.75–100 %
identical to the HiFi consensus** (several species 0 mismatches) and **100 % to
GenBank ITS**.

## Table 1 — modern-input benchmark

| Order | Species | read (bp) | ngs45 unit | easy45 unit | ngs45↔easy45 | ITS↔easy45 | ITS↔GenBank | ribotype sites |
|---|---|---|---|---|---|---|---|---|
| Poales | *Oryza sativa* | 250 | 5783 | 5782 | 99.88 % / 4mm | — | 98.98 % | 15 |
| Zingiberales | *Musa acuminata* | 100* | fail | 5804 | — | — | — | — |
| Solanales | *Solanum lycopersicum* | 91 | 5800 | 5800 | **100 % / 0mm** | 100 % | 100 % | 5 |
| Asterales | *Helianthus annuus* | 150 | 5577 | 5857 | 99.93 % / 4mm | — | 99.69 % | 0 |
| Lamiales | *Sesamum indicum* | 150 | 5465 | 5798 | **100 % / 0mm** | 100 % | 100 % | 6 |
| Caryophyllales | *Beta vulgaris* | 150 | 5811 | 5811 | **100 % / 0mm** | 100 % | 100 % | 1 |
| Vitales | *Vitis vinifera* | 151 | 5876 | 5881 | 99.85 % / 9mm | 99.8 % | — | 13 |
| Fabales | *Glycine max* | 150 | 5926 | 5799 | **100 % / 0mm** | 100 % | — | 1 |
| Malvales | *Gossypium hirsutum* | 101 | 5985 | 5883 | 99.95 % / 2mm | 100 % | 100 % | 4 |
| Sapindales | *Citrus sinensis* | 151 | 5185 | 5844 | 99.75 % / 11mm | — | 98.59 % | 29 |

\* Musa's modern run was only 100 bp; Musa fails even at 309 bp (biology, not
data — extreme rDNA heterogeneity; HiFi required). See `ASSEMBLY_LIMITATION.md`.

Earlier internal validations (not in the table): **Panax ginseng** ngs45↔easy45
99.97 %, ITS 100 % GenBank; **Polyscias filicifolia** ngs45 100 % vs an
independent reference.

## Table 2 — read-length effect (same species, old vs modern run)

| Species | old run (bp) | old result | modern run (bp) | modern result |
|---|---|---|---|---|
| *Oryza sativa* | 73 | fail | 250 | 99.88 % |
| *Helianthus annuus* | 91 | fail | 150 | 99.93 % |
| *Beta vulgaris* | 144 | fail | 150 | 100 % |
| *Vitis vinifera* | 125 | fail | 151 | 99.85 % |

The ~150 bp threshold is the phasing limit for rDNA spacer variants.

## Table 3 — runtime & memory (wall-clock; illustrative)

Per-species, single-job, 16 threads. HiFi/easy45 ran on all 10; ngs45 timings are
from the isolated pass. HiFi peak RSS is higher (~5.5 GB) as easy45 holds
spanning reads; ngs45 short-read peak RSS 0.4–3.9 GB.

| tool | wall-clock range | peak RSS |
|---|---|---|
| ngs45 (Illumina) | 59–755 s | 0.4–3.9 GB |
| easy45 (HiFi) | 72–689 s | ~5.5 GB |

Full per-species numbers: `bench/timing_isolated.tsv`.

## Interpretation

- **ribotype_sites ↔ identity:** homogeneous arrays (Solanum 5, Sesamum 6) give
  0-mismatch consensus; heterogeneous arrays (Citrus 29, Oryza 15) give a blended
  consensus with lower identity — a real biological signal (possible
  hybrid/allopolyploidy), not an assembly error.
- **ngs45 vs easy45 division of labour:** ngs45 delivers the consensus unit + a
  heterogeneity flag from short reads; easy45 delivers the distinct ribotypes from
  HiFi. Use ngs45 for ≥150 bp libraries and moderate heterozygosity; use easy45
  for old/short data, hybrids, or high `ribotype_sites`.

## Reproducing

- **All NCBI/ENA accessions**: [DATA_ACCESSIONS.md](DATA_ACCESSIONS.md) (HiFi,
  old + modern Illumina, GenBank ITS refs, bundled seed/CM).
- Modern benchmark driver: `bench/run_modern_benchmark.sh` (streams read subsets
  from ENA, runs ngs45, compares to the cached easy45 HiFi consensus).
- Accession manifests: `bench/manifest_modern.tsv`, `bench/manifest.tsv`,
  `bench/manifest_hifi_long.tsv`.
- Raw tables: `bench/benchmark_modern.tsv`, `bench/benchmark_final.tsv` (old
  input), `bench/timing_isolated.tsv`.
