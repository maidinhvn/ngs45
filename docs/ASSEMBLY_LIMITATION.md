# Scope & limits of ngs45 (short-read 45S nrDNA assembly)

> Revision history — two earlier conclusions in this file were wrong and were
> corrected by controlled experiments:
> 1. First draft: "a fundamental limitation of short-read assembly." Wrong —
>    modern-read runs recover 9/10.
> 2. Second draft: "requires ≥150 bp reads." Also wrong — a controlled
>    read-length titration shows ngs45 works down to **60 bp** on clean data
>    (`QC.md`, `bench/trunc_titration.tsv`).
> The account below is what the controlled experiments actually support.

## What is (and isn't) the limit

**Read length is NOT a limit** (on adequate data). One clean run
(*Oryza sativa* DRR160520) truncated in-silico to 250→60 bp — everything else held
constant — gave the identical 5783 bp unit at 99.88 % to HiFi at **every** length,
down to 60 bp. So the "~150 bp threshold" seen when comparing old vs modern *runs*
was a **cross-dataset confound**, not causation.

**Standard QC does not predict success either.** Across the old Illumina runs,
assembled and failed cases overlap completely in (read length × base quality)
space — e.g. *Sesamum* (mean Q 7.6) assembles while *Beta*/*Vitis* (mean Q ~25)
fail. See `QC.md` / Figure 5. Public runs differ simultaneously in length, error,
coverage, individual, contamination and rDNA heterozygosity, so a post-hoc
attribution to any single factor is not defensible.

**The one intrinsic limit — extreme intragenomic rDNA heterozygosity.**
*Musa acuminata* fails even with 309 bp reads (the longest run tested) and also
defeats GetOrganelle; only HiFi/easy45 recovers it (5804 bp). Its rDNA array is so
heterogeneous that no short read spans/phases the divergent copies. ngs45 flags
this as a high `ribotype_sites` count (`--call-variants`); as heterozygosity rises
the majority-vote consensus blends ribotypes and its identity to any single
reference drops (Figure 4). This is biology, not a fixable assembly parameter.

## What short reads can and cannot deliver
- **Can:** the consensus 45S unit (base-identical to HiFi where the array is
  homogeneous) + a *site-level* heterozygosity profile (`ribotype_sites`,
  per-site allele fractions).
- **Cannot:** the distinct, full-length ribotype haplotypes — no short read (or
  pair) spans the ~6 kb unit to phase them. That needs long reads.

## Practical guidance
- ngs45 recovers the unit from adequate short-read data across diverse taxa,
  matching HiFi and GenBank at 99.75–100 %. It is not read-length limited in any
  range you are likely to encounter (≥60 bp works on clean data).
- Some individual public runs fail for dataset-specific reasons that QC does not
  reveal; if a run fails, try another library for that species.
- For hybrids/allopolyploids or whenever `ribotype_sites` is high, prefer
  **HiFi + easy45**, which reads each ribotype as one intact molecule.
