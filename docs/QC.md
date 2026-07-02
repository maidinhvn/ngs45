# Quality control of all datasets

Every sequencing dataset used to develop and benchmark ngs45 was QC'd. This step
matters because it settled *why* some public short-read runs failed — and showed
that the intuitive explanations (read length, base quality) do **not** hold up.

## Method
`bench/run_qc.sh` samples up to 200,000 reads per file and computes, with
`seqkit stats -a`: read length, read count, mean base quality (Phred), Q20/Q30
fractions, GC%, N50. Long-read (HiFi) files are handled the same way.
Full table: [`../bench/qc_all_datasets.tsv`](../bench/qc_all_datasets.tsv)
(29 datasets: 10 old-Illumina, 9 modern-Illumina, 10 HiFi).

## Key finding — standard QC does NOT predict assembly success

![Figure 5](figures/Figure5_QC.png)

For the old Illumina runs, assembled (green) and not-assembled (red) points are
**intermixed** in (read length × base quality) space:

| example | read (bp) | mean Q | ngs45 |
|---|---|---|---|
| *Sesamum* | 96 | **7.6** (worst) | ✅ assembled |
| *Gossypium* | 101 | 11.9 | ✅ |
| *Beta* | 144 | 24.3 (good) | ❌ not assembled |
| *Vitis* | 125 | 26.0 (good) | ❌ |
| *Musa* | 309 | 30.0 (excellent) | ❌ |

Neither read length nor base quality separates the two groups. Cross-dataset
comparison of public runs is **confounded** — each run differs simultaneously in
length, error profile, coverage, individual/cultivar, contamination and rDNA
heterozygosity — so no single QC metric explains, or predicts, the outcome.

## The one causal result — controlled read-length titration

![Figure 6](figures/Figure6_titration.png)

To isolate read length from every other variable, one clean modern run
(*Oryza sativa* DRR160520, 250 bp) was truncated in-silico to
250/200/150/130/110/90/73/60 bp — same reads, same individual, same
coverage-in-reads. ngs45 recovered the **identical 5783 bp unit at 99.879 % to
HiFi at every length, down to 60 bp**
([`../bench/trunc_titration.tsv`](../bench/trunc_titration.tsv)).

**Conclusion:** ngs45 is **not read-length limited** on adequate data; the
apparent "~150 bp threshold" seen when comparing old vs modern *runs* was a
cross-dataset confound, not causation. Where specific old runs fail, the cause is
dataset-specific and multi-factorial, and QC metrics alone do not pinpoint it. The
genuinely intrinsic limit is *Musa*-type extreme rDNA heterozygosity, which
requires HiFi regardless of read length (see `ASSEMBLY_LIMITATION.md`).

## Practical QC guidance for ngs45 users
- Standard QC (FastQC/`seqkit stats -a`) is still worth running, but a passing
  report does **not** guarantee assembly, and a poor one does not preclude it.
- The analysis-relevant check is whether enough rDNA reads are recruited (S1
  `bait pairs` in the log) and whether `ribotype_sites` (`--call-variants`) is
  high — high heterozygosity is the real reason to switch to HiFi.
