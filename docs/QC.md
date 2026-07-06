# Quality control of the benchmark datasets

Every sequencing subset used in the benchmark was QC'd, both to document the inputs
and because it clarifies *why* some short-read runs fail — and shows that read length
and base quality alone do not predict the outcome.

## Method

`seqkit stats -a` on each subset (Illumina 5 M read pairs; HiFi 120 k reads): read
length, read count, %≥Q20 / %≥Q30, mean Phred, N50, GC%. Full per-dataset table:
[`../bench/qc.tsv`](../bench/qc.tsv) (Illumina + HiFi for the cross-individual and
same-individual panels).

## Key finding — QC does not separate success from failure

Illumina inputs (cross-individual panel), with ngs45 outcome:

| species | read (bp) | %≥Q30 | ngs45 |
|---|---|---|---|
| *Solanum lycopersicum* | 91 | 98 | ✅ full |
| *Fragaria vesca* | 100 | 94 | ✅ full |
| *Beta vulgaris* | 150 | 93 | ✅ full |
| *Actinidia chinensis* | **85** | 94 | ❌ fail |
| *Helianthus annuus* | 150 | 91 | ❌ fail |
| *Vitis vinifera* | 151 | **79** | ❌ fail |

Short-but-clean libraries succeed (*Solanum* 91 bp, *Fragaria* 100 bp) while some
150 bp libraries fail. So neither read length nor %Q30 alone predicts recovery — the
two axes intermix (see [FIGURES.md](FIGURES.md) Fig 3).

## What actually sets the requirement

- **Read length interacts with spacer divergence.** A same-species control
  (*Actinidia*: 85 bp fails, 150 bp gives the full unit at 99.93 %) shows reads must
  be long enough to bridge that species' ETS/ITS/IGS — but the length needed depends
  on how divergent those spacers are. *Helianthus* and *Vitis* fail even at 150 bp:
  their spacers are too divergent for any short read to span → HiFi/easy45.
- **rDNA heterozygosity** — high `ribotype_sites` (`--call-variants`) means the array
  carries multiple ribotypes; the short-read consensus blends them (a real
  hybrid/allopolyploidy signal, not an assembly error). See
  [ASSEMBLY_LIMITATION.md](ASSEMBLY_LIMITATION.md).

## Practical QC guidance for ngs45 users

- Run standard QC, but a passing report does **not** guarantee assembly and a poor
  one does not preclude it.
- The analysis-relevant checks are in the ngs45 log: S1 `bait pairs` (were enough
  rDNA reads recruited?) and `ribotype_sites` (high → switch to HiFi + easy45).
- Prefer **≥150 bp** paired reads; for divergent-spacer taxa or hybrids, use HiFi.
