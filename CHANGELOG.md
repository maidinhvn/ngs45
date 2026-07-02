# Changelog

## v0.1.0 (2026-07-02)

First public release.

- Seven-stage pipeline (QC → bait → assemble → resolve → CM mature-boundary trim
  → annotate → ribotype variants → report), self-contained (orchestrates SPAdes;
  no external-assembler dependency).
- Short-read-tuned assembly: k-mers capped below read length, `--careful` off by
  default, unit resolved from contigs **and** scaffolds.
- Mature 18S 5′/26S 3′ boundaries via the bundled Rfam SSU/LSU covariance models
  (identical to easy45, so units are directly comparable).
- Outputs: `nrDNA_45S.fasta`, `annotation.gff3`, `its.fasta` / `its_parts.fasta`,
  per-region FASTAs, `summary.tsv`, `report.txt`, `ribotype_variants.tsv`.
- Benchmark vs HiFi/easy45 + GenBank across 10 species / 9 orders; controlled
  read-length titration and QC of all datasets (see `docs/`).
