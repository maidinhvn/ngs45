# ngs45

**Recover the 45S nrDNA transcribed unit from Illumina short reads.**

![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
![Python ≥3.9](https://img.shields.io/badge/python-%E2%89%A53.9-blue.svg)
<!-- after pushing, add:
![CI](https://github.com/maidinhvn/ngs45/actions/workflows/ci.yml/badge.svg) -->

The 45S nrDNA (ETS–18S–ITS1–5.8S–ITS2–26S, followed by the IGS) is a high-copy
tandem array, so genome-skimming and WGS libraries carry it at very deep
coverage. Short reads cannot span a repeat, so ngs45 *assembles* the unit: it
baits rDNA reads against a conserved 45S seed, assembles them with SPAdes,
resolves a single repeat unit out of the tandem-array graph, orients and
annotates it, and (optionally) maps reads back to quantify ribotype
heterogeneity — a signal of hybridisation / allopolyploidy.

## Pipeline

```
S0 QC            cutadapt quality/adapter trim            (optional, --trim)
S1 Bait          iterative bowtie2 recruitment: seed -> +recruited-reads each
                 round, walking from the conserved genes out into ETS/ITS/IGS
S2 Assemble      SPAdes multi-k de Bruijn assembly of the recruited reads
S3 Resolve       BLAST seed -> find rDNA contig; self-BLAST -> tandem period;
                 cut exactly one monomer
S4 Boundary      orient to + strand; Rfam SSU/LSU CMs (cmsearch) -> trim to the
                 mature transcribed unit (18S 5' -> 26S 3'), dropping ETS/IGS
S5 Annotate      ITSx -> SSU/ITS1/5.8S/ITS2/LSU tiling -> GFF + ITS barcode
S6 Variants      map reads back -> minor-allele ribotype sites  (optional)
S7 Report        summary table + human-readable report
```

The S4 mature-boundary trim reuses the **same bundled Rfam CM as
[easy45](https://github.com/maidinhvn/easy45)**, so both tools define the 18S/26S termini
identically and their units are directly comparable.

## Install

```bash
conda env create -f environment.yml
conda activate ngs45
ngs45 check-deps
```

The Python package stays pure-Python; all heavy tools (bowtie2, SPAdes, seqkit,
BLAST, barrnap, ITSx, bwa/samtools/bcftools) are conda dependencies.

## Usage

```bash
ngs45 run -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -o out/ -t 16
```

| Flag | Required? | Meaning |
|------|-----------|---------|
| `-1, --reads1` | **yes** | Illumina R1 (FASTQ[.gz]) |
| `-2, --reads2` | no | Illumina R2 (paired-end; omit for single-end) |
| `-o, --outdir` | no | output folder (default `ngs45_out/`) |
| `-s, --seed-ref` | no | custom 45S seed (default: bundled Arabidopsis T2T unit) |
| `-r, --organelle-ref` | no | plastid+mito genomes to deplete before baiting |
| `--trim` | no | quality/adapter-trim reads first (cutadapt) |
| `--bait-rounds` | no | iterative-baiting rounds (default 3) |
| `--call-variants` | no | report intragenomic ribotype heterogeneity (S6) |
| `-t, --threads` | no | threads (default 4) |

Run `ngs45 run --help` for all parameters.

## Outputs (`outdir/`)

| File | Contents |
|------|----------|
| `nrDNA_45S.fasta` | the mature transcribed unit (18S 5' → 26S 3') |
| `annotation.gff3` | 18S / ITS1 / 5.8S / ITS2 / 26S coordinates (ITSx tiling) |
| `its.fasta` | ITS1–5.8S–ITS2 barcode (for BLAST-based species ID) |
| `its_parts.fasta` | ITS1, 5.8S, ITS2 as separate records |
| `regions/` | per-region FASTAs (all five 45S regions) |
| `summary.tsv` | lengths / coordinates / GC |
| `report.txt` | human-readable run report |
| `ribotype_variants.tsv` | ribotype sites (only with `--call-variants`) |

## Relationship to easy45

ngs45 is the **short-read** counterpart of
[easy45](https://github.com/maidinhvn/easy45) (HiFi long reads). easy45 recovers ribotypes
assembly-free because one HiFi read spans a whole unit; ngs45 must assemble the
unit because short reads cannot. They share the same 45S seed and annotation
philosophy.

## Benchmark & scope

Cross-validated against HiFi/easy45 and GenBank on 10 species across 9 angiosperm
orders. With **modern ≥150 bp PE reads ngs45 recovers the unit for 9/10 species**,
**99.75–100 % identical to the HiFi consensus** (several at 0 mismatches) and
100 % to GenBank ITS. See [docs/BENCHMARK.md](docs/BENCHMARK.md) and [figures](docs/FIGURES.md); all data IDs in [docs/DATA_ACCESSIONS.md](docs/DATA_ACCESSIONS.md).

Scope ([docs/ASSEMBLY_LIMITATION.md](docs/ASSEMBLY_LIMITATION.md), [docs/QC.md](docs/QC.md)):

- **Not read-length limited.** A controlled titration recovers the same unit down
  to 60 bp on clean data; the apparent "≥150 bp" pattern across public runs is a
  cross-dataset confound (QC metrics do not predict success).
- **The intrinsic limit is rDNA heterozygosity** — reported as `ribotype_sites`
  (`--call-variants`). When it is high (hybrids/allopolyploids, e.g. *Musa*), the
  short-read consensus blends ribotypes; prefer HiFi + easy45.

## Install from source

```bash
git clone https://github.com/maidinhvn/ngs45.git
cd ngs45
conda env create -f environment.yml   # brings in SPAdes, bowtie2, ITSx, ...
conda activate ngs45
ngs45 check-deps
```

## Documentation

- [docs/BENCHMARK.md](docs/BENCHMARK.md) — cross-validation vs HiFi/easy45 + GenBank
- [docs/FIGURES.md](docs/FIGURES.md) — benchmark figures
- [docs/QC.md](docs/QC.md) — dataset QC + read-length titration
- [docs/ASSEMBLY_LIMITATION.md](docs/ASSEMBLY_LIMITATION.md) — scope & limits
- [docs/DATA_ACCESSIONS.md](docs/DATA_ACCESSIONS.md) — all NCBI/ENA accessions

## Citation

If you use ngs45, please cite it (see [CITATION.cff](CITATION.cff)) and the
accompanying manuscript (in preparation).

## License

MIT. See [LICENSE](LICENSE).
