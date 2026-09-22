# Changelog

## v0.3.1 (2026-09-22)

- **A partial unit is now reported instead of emitted silently (S4).** When the
  assembly does not span a flanking gene and no adjacent contig carries the missing
  end, ngs45 still emits the best unit it has — but until now it did so without
  saying anything, so a user could take a unit whose 26S was half-length for a
  complete one (*Musa acuminata*: 3864 bp, 26S 1466 bp of the expected ~3390).
  S4 now re-probes the delivered sequence with the covariance models and, if a
  flanking gene is still incomplete, logs a warning, prints it in `report.txt`, and
  records which end in the new `unit_truncated_end` column of `summary.tsv` (empty
  when both genes are full length). Verified on 18 recovered units: the one partial
  case is flagged and none of the 17 complete ones is. Detection only — the emitted
  sequence is unchanged.

- **Baiting no longer truncates a unit, and stays bounded on deep libraries (S1).** The
  v0.3.0 depth-based early stop could end baiting before a gene end's (species-specific,
  seed-invisible) reads arrived in a later round, silently truncating a unit on some
  libraries (e.g. *Populus* same-individual: 5453 bp with a truncated 18S 5′ vs the correct
  5794 bp). The early stop is removed: baiting now always runs to new-read convergence, so
  no round is dropped. To keep that affordable on deep rDNA (which recruits millions of
  reads and blew up bowtie2 index time/memory — e.g. wheat: hours / >14 GB), the reads used
  to *build* each round's reference are capped (`--bait-ref-cap`, default 200000); a capped
  subsample still spans the unit densely, so recruitment is unchanged. Restores correct,
  bounded-memory recovery; the truncation was found by the cross-tool comparison against easy45.

- **SPAdes crash resilience (S2).** SPAdes can abort non-deterministically
  (non-zero exit) on high-TE / low-complexity recruited reads — the same input
  then assembles cleanly on a fresh re-run (seen on *Helianthus*: exit 255 once,
  then a 5857 bp unit). S2 now retries (`--spades-retries`, default 2) before
  failing, and raises a clear error instead of a raw traceback. Removes a
  reproducibility hazard.

- **Docs — troubleshooting a unit that looks wrong** (`docs/TROUBLESHOOTING.md`,
  summarised in the README). ngs45 halts only when the assembly clearly fails to
  span a unit; it can otherwise splice divergent rDNA paralogs or IGS into the
  output and still exit `OK`. Documents how to spot this (unit length vs the
  ~5.8 kb benchmark range, a discontinuous or low-identity BLAST to a related
  45S, ribotype sites clustered at anomalously low depth) and the workaround —
  re-assemble with a tighter `--max-cov` (200; titrate 2000 → 500 → 200), which
  caps the assembly input only, leaving variant-calling depth untouched.
  Documentation only; no code change.

## v0.3.0 (2026-08-23)

- **Faster baiting (S1) — depth early-stop.** Iterative baiting now stops once the
  recruited set already saturates the assembly coverage cap (default `bait_stop_cov`
  6000× = 3× `--max-cov`), after ≥1 extension round. On deep libraries this skips a
  final round that the runaway guard would discard anyway, so the output is
  identical while S1 (95–99 % of runtime) drops ~10× (e.g. wheat ~5 h → ~40 min).
- **Gene-end graft (S4).** When SPAdes fragments the array at a coverage transition
  and the best contig is missing the 18S 5′ or 26S 3′ end (cmsearch reports a
  one-sided `trunc`), the missing end is grafted from an adjacent rDNA contig by
  their overlap. Rescues species like *Sesamum* (5344 → 5798 bp, 100 % to the HiFi
  unit); a no-op when the monomer already spans both genes.
- **`--reference-guided` fallback (opt-in).** When de novo assembly cannot span a
  unit, map the baited reads to a (close) `--seed-ref` and call a consensus — a
  usable 45S for NGS-only libraries with no HiFi option. Reference-biased and
  clearly labelled (FASTA header, report banner, `summary.tsv`); default off.
- **Guiding failure message.** A too-short assembly now points the user at a closer
  `--seed-ref` (a congener 45S rescued a *Vitis* run 2904 → 6739 bp) and, as a last
  resort, `--reference-guided`.
- **Contaminant warning (S3).** Flags a picked contig that looks unlike a plant
  nuclear 45S (GC/identity to the seed) — a fungal endophyte or plastid rRNA
  recruited via the conserved genes.
- **Low-depth warning (S3).** Warns when the picked contig's coverage is low.
- **barrnap pinned `>=0.9,<1.0`** + a `check-deps` probe that fails fast if barrnap
  lacks the eukaryote DB (bioconda's `barrnap 1.10.6` ships only bac/arc/fun).
- **Pipeline figure updated** (`docs/pipeline.{dot,png,svg}`).

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
