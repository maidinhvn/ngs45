# Data availability — NCBI/ENA/DDBJ accessions

All public sequencing data used to benchmark **ngs45** (Illumina) and **easy45**
(PacBio HiFi), plus the GenBank references used for validation. Runs were fetched
from INSDC as rDNA-sufficient subsets (Illumina 5 M read pairs; HiFi 120 k reads —
rDNA is high-copy, so a random subset retains ample coverage; full runs are at the
accessions below). Prefixes: `SRR` NCBI, `ERR` ENA, `DRR` DDBJ; `SAMN/SAMEA/SAMD` =
BioSample.

## Bundled reference data (shipped in the package)

| Item | Source | Accession |
|---|---|---|
| 45S recruitment/orientation seed | *Arabidopsis thaliana* T2T NOR2, one 45S repeat | GenBank **OR453402.1** (36450–46501) |
| Mature-boundary covariance models | Rfam eukaryotic SSU + LSU rRNA | Rfam **RF01960** (SSU), **RF02543** (LSU) |

## Table S1 — Cross-individual benchmark (12 species, 12 angiosperm orders)

Illumina and HiFi here are from **different individuals** (different BioSamples) of
each species; tests recovery of the species' 45S unit across broad taxonomy.

QC below is for the subset actually used (Illumina 5 M pairs, HiFi 120 k reads):
Illumina read length (bp) and % bases ≥ Q30; HiFi read-N50 (kb) and % bases ≥ Q30.

| Species | Order | Illumina (ngs45) | Ill. len / Q30 | HiFi (easy45) | HiFi N50 / Q30 | GenBank ITS ref |
|---|---|---|---|---|---|---|
| *Beta vulgaris* | Caryophyllales | SRR29552944 | 150 bp / 93 % | SRR37382116 | 17.1 kb / 95 % | PZ361546.1 |
| *Helianthus annuus* | Asterales | SRR8888860 | 150 bp / 91 % | SRR14782853 | 21.8 kb / 100 % | PP232226.1 |
| *Solanum lycopersicum* | Solanales | DRR040154 | 91 bp / 98 % | SRR15243717 | 19.1 kb / 96 % | PZ387679.1 |
| *Sesamum indicum* | Lamiales | DRR806862 | 150 bp / 95 % | SRR21601246 | 17.1 kb / 96 % | OM746326.1 |
| *Citrus sinensis* | Sapindales | SRR33854505 | 151 bp / 89 % | SRR27236983 | 19.9 kb / 96 % | PV589351.1 |
| *Glycine max* | Fabales | SRR11929797 | 150 bp / 94 % | SRR38619744 | 15.1 kb / 95 % | PP060534.1 |
| *Fragaria vesca* | Rosales | ERR14008837 | 100 bp / 94 % | SRR38138057 | 11.1 kb / 96 % | PZ361553.1 |
| *Populus trichocarpa* | Malpighiales | SRR37518787 | 151 bp / 92 % | SRR37518457 | 14.5 kb / 97 % | JQ898636.1 |
| *Vitis vinifera* | Vitales | ERR15002908 | 151 bp / 79 % | ERR17353932 | 7.6 kb / 95 % | OQ519992.1 |
| *Actinidia chinensis* | Ericales | SRR29894970 | **85 bp** / 94 % | SRR24236049 | 9.3 kb / 98 % | OR516785.1 |
| *Oryza sativa* | Poales | DRR160520 | 250 bp / 95 % | SRR13280199 | 19.2 kb / 100 % | PZ572868.1 |
| *Musa acuminata* | Zingiberales | SRR25581090 | 100 bp / 100 % | SRR23425448 | 13.2 kb / 98 % | PZ437810.1 † |

All ITS references verified on NCBI (organism + rDNA content confirmed). Notes:
PZ387679.1 (*Solanum lycopersicum* var. *cerasiforme*) and † PZ437810.1 (*Musa
acuminata* var. *rutilipes*) are infraspecific varieties of the benchmark species;
JQ898636.1 (Populus, 18S-partial) and PP060534.1 (Glycine, 5.8S-partial) span only
part of the ITS barcode.

## Table S2 — Same-individual validation (7 species)

HiFi and Illumina from the **same individual** (same BioSample) — removes the
intraspecific-variation confound, so any ngs45-vs-easy45 difference is method-only.

| Species | Order | HiFi (easy45) | HiFi N50 / Q30 | Illumina (ngs45) | Ill. len / Q30 | BioSample |
|---|---|---|---|---|---|---|
| *Beta vulgaris* | Caryophyllales | SRR21631919 | 18.2 kb / 95 % | SRR21631918 | 150 bp / 91 % | SAMN30877666 |
| *Citrus sinensis* | Sapindales | DRR762555 | 13.8 kb / 95 % | DRR762558 | 150 bp / 97 % | SAMD01669405 |
| *Sesamum indicum* | Lamiales | SRR21601246 | 17.1 kb / 96 % | SRR21384248 | 150 bp / 94 % | SAMN30610312 |
| *Populus trichocarpa* | Malpighiales | SRR37518065 | 15.7 kb / 97 % | SRR37518062 | 151 bp / 89 % | SAMN56350631 |
| *Oryza sativa* | Poales | ERR11838731 | 19.6 kb / 95 % | ERR11838732 | 100 bp / 79 % | SAMEA112646987 |
| *Musa acuminata* | Zingiberales | SRR24235289 | 16.1 kb / 97 % | SRR24235287 | 151 bp / 100 % | SAMN34246735 |
| *Lindera aggregata* | Laurales | SRR31617073 | 17.9 kb / 93 % | SRR31617071 | 150 bp / 96 % | SAMN43038876 |

## Table S3 — Read-length control (*Actinidia chinensis*)

Same species, two Illumina libraries of different read length, both compared to the
same HiFi/easy45 consensus (SRR24236049). Isolates the effect of read length.

| Library | Instrument | Read length | Year | ngs45 outcome |
|---|---|---|---|---|
| SRR29894970 | HiSeq 2000 | 85 bp | 2024 | fail (best contig 3603 bp) |
| SRR28002090 | NovaSeq 6000 | 150 bp | 2024 | full 5834 bp unit, 99.93 % vs HiFi |

## Table S4 — Application datasets (in-house sequencing)

Generated in-house (Illumina WGS / PacBio HiFi), not from NCBI. Internal sample IDs
below; GenBank accessions to be assigned on submission.

**Polyscias — 8 samples (ngs45, Illumina):**

| Species | Internal ID |
|---|---|
| *Polyscias balfouriana* | IM241218-8 |
| *Polyscias cumingiana* | IM250605-1 |
| *Polyscias filicifolia* | IM240312-1 |
| *Polyscias guilfoylei* (green) | IM241218-7 |
| *Polyscias guilfoylei* var. *quinquefolia* | IM241218-6 |
| *Polyscias guilfoylei* (white) | IM241218-5 |
| *Polyscias scutellaria* | IM241218-4 |
| *Polyscias serrata* | IM241218-9 |

**Centella asiatica — 5 accessions (easy45, HiFi):** SNU-CA-4, -10, -11, -14, -18.
Validated against an independent in-house 45S reconstruction.
