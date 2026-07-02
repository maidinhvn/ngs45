# Operating envelope of ngs45 (short-read 45S nrDNA assembly)

> This document supersedes an earlier draft that wrongly called the failures a
> "fundamental limitation of short-read assembly". A cross-tool benchmark showed
> the failures were driven by **input read length**, not by the tool — see below.

ngs45 recovers the 45S nrDNA transcribed unit from Illumina reads by de-novo
assembly (SPAdes) of recruited rDNA reads. Its success is governed by **two
independent factors**, both understood and quantified on a 10-species benchmark
that used HiFi/easy45 as the gold standard.

## Factor 1 — read length (dominant, technology-driven)

rDNA spacers (ITS1/ITS2/ETS) carry many closely-spaced sequence variants (see
Factor 2). Each variant is a *bubble* in the de-Bruijn graph. To lay down one
consistent path the assembler must **phase** across those bubbles, which requires
reads long enough to span them. Short reads cannot, so the graph fragments into
per-gene contigs.

Empirically the threshold is **~150 bp**. The same species flip from fail→success
purely by swapping an old short-read run for a modern ≥150 bp run:

| Species | OLD run | result | MODERN run | result (vs HiFi) |
|---|---|---|---|---|
| Oryza sativa | 73 bp | fail | **250 bp** | **99.88 %** |
| Helianthus annuus | 91 bp | fail | **150 bp** | **99.93 %** |
| Beta vulgaris | 144 bp | fail | **150 bp** | **100 %** |
| Vitis vinifera | 125 bp | fail | **151 bp** | **99.85 %** |

Benchmark-wide: **4/10 with old short-read input → 9/10 with modern ≥150 bp
input.** Read length correlates with sequencing generation (old Illumina GA/HiSeq
= 36–100 bp; modern MiSeq/HiSeq-v4/NovaSeq = 150–300 bp), which is why published
studies using modern libraries assemble 45S nrDNA routinely.

**Guidance:** use paired-end reads **≥150 bp** (ideally 250 bp). Very short/old
runs (<125 bp) are expected to fail on all but the most homogeneous rDNA.

## Factor 2 — intragenomic ribotype diversity (residual, biology-driven)

The rDNA array holds hundreds–thousands of copies; incomplete concerted evolution
leaves them non-identical (*ribotypes*), especially in the spacers and in
hybrids/allopolyploids. ngs45 reports this as `ribotype_sites` (map-back minor-
allele sites, `--call-variants`). Higher diversity → the majority-vote consensus
blends ribotypes → lower identity to any single reference:

| Species | ribotype_sites | ITS vs GenBank |
|---|---|---|
| Solanum | 5 | 100 % |
| Oryza | 15 | 98.98 % |
| Citrus | 29 | 98.59 % |

At the extreme this defeats short reads regardless of length. **Musa acuminata**
(banana) fails even with **309 bp** reads — the longest run in the benchmark — and
also defeats GetOrganelle; only HiFi/easy45 recovers it (5804 bp). Its rDNA is so
heterogeneous that no short read spans/phases the divergent copies.

**Guidance:** when `ribotype_sites` is high (roughly >20), treat the consensus as
approximate and prefer **HiFi + easy45**, which reads each ribotype as one intact
molecule (no phasing needed) and can enumerate the distinct variants.

## What short reads can and cannot deliver

- **Can:** the consensus 45S unit (base-identical to HiFi where the array is
  homogeneous) + a *site-level* heterogeneity profile (`ribotype_sites`,
  per-site allele fractions).
- **Cannot:** the distinct, full-length ribotype haplotypes — no short read (or
  read pair) spans the ~6 kb unit to phase them. That requires long reads.

## Summary
ngs45 is reliable for taxa sequenced with **≥150 bp reads and moderate rDNA
heterogeneity** — there it matches HiFi and GenBank (99.75–100 %). Use HiFi/easy45
for old/short-read data, for hybrids/allopolyploids, or whenever `ribotype_sites`
is high and the individual ribotypes matter.
