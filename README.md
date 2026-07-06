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
conda install -c conda-forge -c bioconda ngs45
ngs45 check-deps
```

The Python package stays pure-Python; all heavy tools (bowtie2, SPAdes, seqkit,
BLAST, barrnap, ITSx, bwa/samtools/bcftools) are pulled in automatically as
conda dependencies. To build the exact development environment from source
instead, see [Install from source](#install-from-source) below.

## Usage

ngs45 works in **two stages**, both inside the single `run` command — a flag
decides whether the second one happens:

**Stage 1 — assemble the 45S unit** *(always runs)*. Baits rDNA reads, assembles
them (SPAdes), resolves one repeat unit, trims it to the mature 18S→26S
transcribed unit with the Rfam CMs, and annotates the ITS barcode. This is all
most users need.

```bash
ngs45 run -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -o out/ -t 16
```

**Stage 2 — ribotype heterogeneity** *(optional — add `--call-variants`)*. Maps
the recruited reads back onto the assembled unit and flags positions where the
thousands of rDNA copies disagree — the signal of intragenomic ribotype variation
(hybridisation / allopolyploidy), written to `ribotype_variants.tsv`. It reuses
the stage-1 products (no re-assembly) and maps the **full-depth** read set.

```bash
ngs45 run -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -o out/ -t 16 --call-variants
```

Both stages run in this one command. On a whole-genome library, stage 1 is
~15 min and stage 2 adds ~1 min; omit `--call-variants` if you only need the
unit / ITS barcode.

| Flag | Required? | Meaning |
|------|-----------|---------|
| `-1, --reads1` | **yes** | Illumina R1 (FASTQ[.gz]) |
| `-2, --reads2` | no | Illumina R2 (paired-end; omit for single-end) |
| `-o, --outdir` | no | output folder (default `ngs45_out/`) |
| `-s, --seed-ref` | no | custom 45S seed (default: bundled Arabidopsis T2T unit) |
| `--call-variants` | no | **stage 2**: report intragenomic ribotype heterogeneity |
| `--trim` | no | quality/adapter-trim reads first (cutadapt) |
| `--bait-rounds` | no | iterative-baiting rounds (default 3) |
| `--max-cov` | no | cap assembly depth (default `2000`; `0` disables) — see note |
| `-t, --threads` | no | threads (default 4) |

Run `ngs45 run --help` for all parameters.

> **Speed / coverage note.** The rDNA array is 10⁴–10⁵× deep — useless for
> assembly and the dominant runtime cost. `--max-cov` downsamples the **assembly
> input only** (to ~2000× of a repeat by default), cutting SPAdes from ~70 min to
> ~2 min with an essentially identical unit (99.98 % id). Stage 2 variant calling always uses
> the **full** read depth, so ribotype sensitivity is unaffected by this cap.

### Batch mode

Process a whole folder of samples in one command — ideal for an overnight run:

```bash
ngs45 batch -i rawdata/ -o out/ -t 8 --call-variants
```

`batch` auto-detects the layout: paired FASTQs in one flat folder (mates named
`*_R1/_R2` **or** `*_1/_2`) and/or one subfolder per sample. Each sample is written
to `out/<sample>/`; a failed sample is logged and skipped, samples already done are
skipped (so an interrupted run just re-launches), and a `batch_summary.tsv`
aggregates every sample's unit length, ITS/region lengths, ribotype sites and
timing. All `run` options (`-t`, `--call-variants`, `--max-cov`, …) apply.

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

Cross-validated against **easy45 / PacBio HiFi** and GenBank across 12 species / 12
angiosperm orders, on two tiers ([docs/BENCHMARK.md](docs/BENCHMARK.md), [figures](docs/FIGURES.md);
data IDs in [docs/DATA_ACCESSIONS.md](docs/DATA_ACCESSIONS.md)):

- **Same individual** (HiFi + Illumina from the same plant, 7 species): **6/7
  concordant, 5 base-identical (0 mismatch), mean ≈ 99.98 %** — the short-read unit
  *is* the HiFi molecule, with no intraspecific confound.
- **Cross-individual** (12 species): ngs45 recovers the full unit for **8/12 at
  99.76–100 % identity to HiFi** (9/12 given ≥150 bp reads), 1 partial (*Musa*,
  hybrid), and 2 short-read limits (*Helianthus*, *Vitis*) that easy45 recovers.

Scope ([docs/ASSEMBLY_LIMITATION.md](docs/ASSEMBLY_LIMITATION.md), [docs/QC.md](docs/QC.md)):

- **Read length ≥150 bp is needed but not sufficient.** A same-species control
  (*Actinidia*: 85 bp fails, 150 bp gives the full unit at 99.93 %) shows short reads
  must be long enough; but *Helianthus* and *Vitis* fail even at 150 bp because their
  transcribed-unit spacers are too divergent for short reads to span — use HiFi/easy45.
- **rDNA heterozygosity** — reported as `ribotype_sites` (`--call-variants`). When it
  is high (hybrids/allopolyploids), the short-read consensus blends ribotypes;
  prefer HiFi + easy45 for full ribotype phasing.

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
