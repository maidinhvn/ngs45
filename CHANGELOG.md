# Changelog

## v0.2.0 (2026-07-06)

- **`batch` subcommand** — run a folder of samples in one command (auto-detects
  `*_R1/_R2` or `*_1/_2`, flat or subfolder; resume; `batch_summary.tsv`).
- **Coverage cap** (`--max-cov`, default 2000) downsamples the *assembly input only*
  — SPAdes ~70 min → ~2 min for an essentially identical unit; variant calling keeps
  full depth, so ribotype sensitivity is unaffected.
- **Tandem-duplication QC** (S4) — self-BLAST collapses a spurious adjacent internal
  repeat from a unit-junction mis-assembly; reported as `qc_tandem_dup_bp`.
- **Robustness (surfaced by real-data benchmarking):** SPAdes `--phred-offset 33`
  (avoids BayesHammer offset-detection aborts); bait runaway guard (an extension
  round recruiting >30 % of the library falls back to the prior round); clear error
  on empty contigs.
- **Benchmark refreshed:** cross-individual (12 species / 12 orders) + same-individual
  (6 species, same BioSample; 5/6 base-identical to the HiFi consensus). See
  `docs/BENCHMARK.md`.

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
