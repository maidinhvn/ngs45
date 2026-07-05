# ngs45 pipeline — step-by-step

Textual companion to [`pipeline.png`](pipeline.png) (source
[`pipeline.dot`](pipeline.dot)). Describes the workflow as implemented in
`src/ngs45/` (v0.2.0). Stages are chained by `src/ngs45/pipeline.py`; each writes
to the work directory and records progress so `--resume` can skip finished stages.

## Inputs

1. **Illumina paired-end reads** (`FASTQ[.gz]`) — the user's short-read library
   (R1/R2; single-end is also accepted).
2. **Bundled 45S seed** — the *Arabidopsis* T2T 45S unit (GenBank OR453402), used
   to bait rDNA reads and to orient the assembly.
3. **Bundled Rfam SSU+LSU CMs** — covariance models RF01960 (18S / SSU) and
   RF02543 (26S / LSU), used to trim the unit to its mature gene boundaries.

## Stages (S0 → S7)

**S0 — QC / trim reads** *(optional, `--trim`)*
Adapter/quality trimming with **cutadapt**. Skipped unless `--trim` is given.
`reads → cleaned reads`

**S1 — Bait rDNA reads**
Iterative recruitment with **bowtie2**: map reads to the 45S seed, collect the
hits, extend the reference, and repeat over several rounds to pull in the full
rDNA read set.
`reads + 45S seed → recruited reads (rDNA)`

**S2 — Assemble recruited reads** *(core step)*
De-novo assembly with **SPAdes** using a multi-k ladder auto-filtered to the read
length (k < read length). *Optional* coverage cap (`--max-cov`, on by default,
disable with `--max-cov 0`): the rDNA array is 10^4–10^5× deep, so the assembly
input is downsampled to a sane depth to keep SPAdes fast — the full read set is
untouched and still used by S6.
`recruited reads → contigs + scaffolds`

**S3 — Resolve one repeat unit**
Identify the rDNA-bearing sequence by **BLAST** of the 45S seed against the
contigs + scaffolds, then cut a single repeat unit (monomer).
`contigs + scaffolds → one repeat unit`

**S4 — Orient & trim to the mature unit (18S 5′ → 26S 3′)**
Orient and trim to the mature transcribed unit with **cmsearch** against the Rfam
SSU+LSU CMs. Then **collapse tandem-dup (QC, always on)**: the unit is self-BLASTed
and any spurious adjacent internal tandem duplication (an assembly artifact across
a unit junction) is collapsed. Reports `qc_tandem_dup_bp`.
`one repeat unit + CMs → mature 18S–26S unit`

**S5 — Annotate + ITS barcode**
Delimit the five regions (18S / ITS1 / 5.8S / ITS2 / 26S) with **ITSx**; if ITSx
delimits nothing, fall back to **barrnap** (rRNA genes only).
`mature unit → consensus unit + annotation + ITS barcode`

**S6 — Ribotype heterozygosity** *(optional, `--call-variants`)*
Map the **full-depth** recruited reads back onto the unit with **bwa | samtools |
bcftools** and count positions where the many rDNA copies disagree.
`consensus unit → ribotype variants`

> **ribotype_sites (`--call-variants`):** map reads back to the unit; a position
> with **minor allele ≥ 10 % at depth ≥ 20** is a ribotype-variant site. A high
> count means a **heterogeneous array** — a hybridisation / allopolyploidy signal.

**S7 — Summary + report**
Write the summary table and the human-readable report.

## Outputs (`outdir/`)

| File | Contents |
|---|---|
| `nrDNA_45S.fasta` | mature transcribed unit (18S–26S) |
| `its.fasta` + `its_parts.fasta` | ITS barcode (whole + ITS1/5.8S/ITS2) |
| `regions/` | the five 45S regions as separate FASTAs |
| `annotation.gff3` | region coordinates |
| `ribotype_variants.tsv` | ribotype sites (only with `--call-variants`) |
| `summary.tsv` + `report.txt` | summary table + run report |

## Batch mode

`ngs45 batch -i <folder> -o <out>` runs this same pipeline over every sample in a
folder — auto-detecting paired FASTQs (`*_R1/_R2` or `*_1/_2`) in a flat folder
and/or one subfolder per sample. Each sample goes to `out/<sample>/`; failures are
logged and skipped, finished samples are skipped on re-launch (resume), and
`out/batch_summary.tsv` aggregates every sample's summary + timing + status. All
`run` options apply per sample.
