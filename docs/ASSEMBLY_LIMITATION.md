# Scope & limits of ngs45 (short-read 45S nrDNA assembly)

> Revision history — earlier drafts of this file over-generalised from single-species
> tests and have been corrected against the 2026-07 cross-individual + same-individual
> benchmark:
> 1. "a fundamental limitation of short-read assembly" — too strong; ngs45 recovers
>    8/12 species (9/12 with ≥150 bp reads), base-identical to HiFi on the same plant.
> 2. "works down to 60 bp" — that came from truncating one easy species (*Oryza*) and
>    does not generalise: *Actinidia* fails at 85 bp and succeeds at 150 bp.
> 3. "*Musa* fails — intrinsic heterozygosity" — wrong: on same-individual data ngs45
>    recovers *Musa*'s unit at 100 % (0 mismatch); the earlier partial was a data
>    artefact. The consensus is recoverable; only ribotype *phasing* needs HiFi.

## What sets the limit — divergent spacers, via read length

Short reads recover the unit by assembling through the transcribed spacers
(ETS/ITS1/ITS2) and, for baiting, the IGS. Whether they can depends on **how
divergent those spacers are** relative to the read length:

- **Read length is necessary but not sufficient.** Same-species control on
  *Actinidia chinensis* (same HiFi reference): an **85 bp** library fails (best contig
  3603 bp), a **150 bp** library gives the full 5834 bp unit at 99.93 %. So reads must
  be long enough — but 150 bp is not always enough: *Helianthus* and *Vitis* fail at
  150–151 bp because their spacers are too divergent for any short read to span.
- **QC alone does not predict it** ([QC.md](QC.md)): short-but-clean libraries succeed
  (*Solanum* 91 bp, *Fragaria* 100 bp) while some 150 bp libraries fail. Read length
  and %Q30 intermix across success/failure.

When short reads cannot span the divergent spacers, ngs45 stops cleanly at S3 (best
rDNA contig < 4000 bp) rather than emit a wrong unit — and **easy45/HiFi recovers
those species in full** (each HiFi read spans the whole repeat).

## What short reads can and cannot deliver

- **Can:** the consensus 45S unit — base-identical to the HiFi consensus where the
  array is homogeneous (5/7 same-individual species at 0 mismatch) — plus a
  *site-level* heterozygosity profile (`ribotype_sites`, per-site allele fractions
  from `--call-variants`).
- **Cannot:** the distinct, full-length ribotype haplotypes. No short read (or pair)
  spans the ~6 kb unit to phase divergent copies into separate molecules — that needs
  long reads (easy45). This holds even when the consensus itself is recovered (e.g.
  *Musa*: consensus recovered, but per-ribotype phasing still needs HiFi).

## Practical guidance

- Use **≥150 bp** paired reads. ngs45 then recovers the unit across diverse taxa,
  matching HiFi/GenBank at 99.76–100 %.
- If a run fails at S3 (contig too short), the species' spacers are likely too
  divergent for short reads — use **HiFi + easy45** for that species.
- For hybrids/allopolyploids or whenever `ribotype_sites` is high, prefer HiFi +
  easy45 for full ribotype phasing.
