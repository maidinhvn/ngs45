# Data availability — NCBI/ENA accessions

All sequencing data are public on the NCBI SRA / ENA. Read length (bp) is the
per-read length actually used (streamed subset).

## Bundled reference data (shipped in the package)

| Item | Source | Accession |
|---|---|---|
| 45S recruitment/orientation seed | *Arabidopsis thaliana* T2T NOR2, one 45S repeat | GenBank **OR453402.1** (36450–46501) |
| Mature-boundary covariance models | Rfam eukaryotic SSU + LSU rRNA | Rfam **RF01960** (SSU), **RF02543** (LSU) |

## Main benchmark — 10 species / 9 orders

Per species: the HiFi run (gold standard, easy45); the *old/short* Illumina run
(baseline, 4/10); the *modern ≥150 bp* Illumina run (final, 9/10); a
representative GenBank ITS reference.

| Order | Species | HiFi (easy45) | Illumina — old (bp) | Illumina — modern (bp) | GenBank ITS ref |
|---|---|---|---|---|---|
| Poales | *Oryza sativa* | SRR13280199 | ERR009651 (73) | DRR160520 (250) | via query |
| Zingiberales | *Musa acuminata* | SRR23425448 | SRR6489399 (309) | SRR25581090 (100)† | EU433925.1 |
| Solanales | *Solanum lycopersicum* | SRR15243717 | SRR1572286 (100) | DRR040154 (91) | AF244747.1 |
| Asterales | *Helianthus annuus* | SRR14782853 | SRR5160794 (91) | SRR8888860 (150) | via query |
| Lamiales | *Sesamum indicum* | SRR21601246 | ERR710504 (96) | DRR806862 (150) | AF169853.1 |
| Caryophyllales | *Beta vulgaris* | SRR37382116 | SRR6315553 (144) | SRR29552944 (150) | AY858597.1 |
| Vitales | *Vitis vinifera* | ERR17353932 | SRR4210153 (125) | ERR15002908 / ERR13382344 (151) | JX290091.1 |
| Fabales | *Glycine max* | SRR38619744 | SRR1463414 (100) | SRR11929797 (150) | via query |
| Malvales | *Gossypium hirsutum* | SRR38842326 | SRR1580594 (101) | — (used old 101) | U12719.1 |
| Sapindales | *Citrus sinensis* | SRR27236983 | ERR466624 | SRR33854505 (151) | via query |

† Musa's modern run was only 100 bp; Musa also fails at 309 bp (old run) — an
intrinsic-heterozygosity case that needs HiFi (see `ASSEMBLY_LIMITATION.md`).

Notes:
- **HiFi runs** are PacBio HiFi/CCS; for easy45 a subset was streamed. For
  Solanum and Sesamum the original short HiFi run was replaced with a long-read
  HiFi run (SRR15243717, SRR21601246) so reads span the unit.
- **GenBank ITS ref** = a representative accession; the actual comparison used the
  full set of records returned by
  `"<species>"[Organism] AND internal transcribed spacer[Title] AND 400:900[SLEN]`.
  "via query" = records were retrieved by that query (no single representative
  pinned; ITS-vs-GenBank not scored for that species).

## Additional validation species (early cross-checks)

| Species | HiFi | Illumina | Reference |
|---|---|---|---|
| *Panax ginseng* (Araliaceae) | SRR35147962 | SRR5196586 | GenBank ITS **U41680.1**, **U41682.1** |
| *Polyscias filicifolia* (Araliaceae) | — | IM240312-1 (in-house WGS) | in-house *P. filicifolia* 45S nrDNA reference |

## Derived / in-silico datasets (no new download)

| Analysis | Source accession | Derivation |
|---|---|---|
| Read-length titration (controlled test that read length — not dataset identity — sets the ~150 bp threshold) | **DRR160520** (*Oryza sativa*, 250 bp) | reads truncated in-silico to 250 / 200 / 150 / 130 / 110 / 90 bp; everything else (individual, coverage, chemistry) held constant. Output: `bench/trunc_titration.tsv` |

## Machine-readable manifests (in `bench/`)

| File | Contents |
|---|---|
| `manifest.tsv` | old benchmark: HiFi + short Illumina accessions + sizes |
| `manifest_hifi_long.tsv` | long-read HiFi accession per species (used by easy45) |
| `manifest_modern.tsv` | modern ≥150 bp Illumina accession per species |
| `benchmark_modern.tsv` | final results (accession, read length, identities) |
