# Troubleshooting — a recovered unit that looks wrong

ngs45 always returns a sequence. It stops with an error only when the assembly
clearly fails to span a unit (see [ASSEMBLY_LIMITATION.md](ASSEMBLY_LIMITATION.md));
it does **not** fail loudly when the assembly instead wanders into divergent rDNA
paralogs or the IGS and splices them into the unit. Such a run exits `OK`, so the
output is worth a look before you use it.

## How to spot a bad unit

**1. Length.** `unit_len` in `summary.tsv` should be ~5.8 kb — 5,780–5,940 bp
across the angiosperms in our benchmark (see [BENCHMARK.md](BENCHMARK.md)).
Running several samples of the same species makes an outlier obvious.

**2. Alignment to a reference.** BLAST the unit against a conspecific or closely
related 45S sequence. A good unit gives **one continuous alignment at high
identity**. Warning signs:

| Pattern | Likely cause |
|---|---|
| alignment stops early | truncated unit (missing one end) |
| tail at markedly lower identity | chimeric join to a divergent paralog |
| internal segment well below the surrounding identity | spurious insert (paralog / IGS) |

**3. Depth of ribotype sites.** In `ribotype_variants.tsv`, check the `depth`
column. Genuine ribotype variation sits at the unit's normal depth; variant sites
clustered in a region of much lower depth usually mark an assembly artifact
rather than biology.

## How to fix it

Re-assemble with a tighter coverage cap:

```bash
ngs45 run -1 R1.fq.gz -2 R2.fq.gz -o out_lowcov --max-cov 200
```

The rDNA array is 10⁴–10⁵× deep and its copies are not identical. At the default
cap (2000) rare divergent paralogs can still carry enough reads for SPAdes to
build them into the unit; a tighter cap pushes the assembly graph toward the
consensus repeat. If 200 does not help, titrate (2000 → 500 → 200).

`--max-cov` limits the **assembly input only** — S6 variant calling still uses
the full read depth, so ribotype sensitivity is unaffected.

## Caveat

A tighter cap biases the assembly toward the dominant ribotype. For a sample
that genuinely carries two co-dominant arrays (a hybrid or an allopolyploid),
inspect `ribotype_variants.tsv` rather than relying on the consensus alone — the
consensus of two divergent parental arrays is not a sequence that exists in the
organism.

This is an empirical workaround, not a guaranteed fix for every failure mode.
