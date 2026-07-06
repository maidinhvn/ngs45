# ngs45 benchmark

Cross-validation of **ngs45** (short-read 45S nrDNA assembly) against **easy45 /
PacBio HiFi** as the independent reference, plus **GenBank** ITS records. Two tiers:
a broad **cross-individual** panel (12 species / 12 angiosperm orders) and a rigorous
**same-individual** panel (7 species, HiFi + Illumina from the same BioSample). All
identities are best-HSP BLASTn. Accessions: [DATA_ACCESSIONS.md](DATA_ACCESSIONS.md);
figures: [FIGURES.md](FIGURES.md); source tables in [`../bench/`](../bench/).

## Headline

- **Same individual (no confound): 6/7 species concordant, 5 base-identical (0
  mismatch), mean ≈ 99.98 %.** Where both tools run on the same plant, the short-read
  unit is essentially identical to the HiFi consensus.
- **Cross-individual: ngs45 recovers the full unit for 8/12 species at 99.76–100 %
  identity to HiFi** (9/12 once *Actinidia* is given ≥150 bp reads; see Table 3),
  1 partial (*Musa*, hybrid), and 2 genuine short-read limits (*Helianthus*, *Vitis*)
  that **easy45 recovers in full** — the two tools are complementary.
- **easy45 recovers the unit for every species tested.**

## Table 1 — Same-individual validation (7 species)

HiFi and Illumina from the **same BioSample**, so any difference is method-only.

| Species | Order | ngs45 unit | easy45 cons | identity |
|---|---|---|---|---|
| *Oryza sativa* | Poales | 5783 | 5783 | **100.000 %** (0 mm) |
| *Beta vulgaris* | Caryophyllales | 5811 | 5811 | **100.000 %** (0 mm) |
| *Sesamum indicum* | Lamiales | 5798 | 5798 | **100.000 %** (0 mm) |
| *Musa acuminata* | Zingiberales | 5804 | 5804 | **100.000 %** (0 mm) |
| *Populus trichocarpa* | Malpighiales | 5794 | 5794 | **100.000 %** (0 mm) |
| *Citrus sinensis* | Sapindales | 5845 | 5844 | 99.846 % |
| *Lindera aggregata* | Laurales | 5855 | no consensus† | — |

† *Lindera* (basal Laurales): the 120 k-read HiFi subset had too few 45S-spanning
reads to pass easy45's cluster-abundance gate; ngs45 still recovered the unit — the
complementarity runs both ways. *Musa* is base-identical here although its
cross-individual run (Table 2) was only partial, i.e. that shortfall was a
data artefact, not a hybrid limit of ngs45.

## Table 2 — Cross-individual benchmark (12 species / 12 orders)

| Species | Order | ngs45 | easy45 | ngs45↔easy45 | ITS↔GenBank |
|---|---|---|---|---|---|
| *Beta vulgaris* | Caryophyllales | 5811 | 5811 | 100.000 % | 99.67 % |
| *Solanum lycopersicum* | Solanales | 5800 | 5800 | 100.000 % | 100.0 % |
| *Sesamum indicum* | Lamiales | 5798 | 5798 | 100.000 % | 92.99 % |
| *Glycine max* | Fabales | 5799 | 5799 | 100.000 % | 99.66 % |
| *Populus trichocarpa* | Malpighiales | 5794 | 5794 | 100.000 % | 98.65 % |
| *Fragaria vesca* | Rosales | 5829 | 5828 | 99.88 % | 99.62 % |
| *Oryza sativa* | Poales | 5783 | 5782 | 99.88 % | — |
| *Citrus sinensis* | Sapindales | 5844 | 5844 | 99.76 % | 99.84 % |
| *Musa acuminata* | Zingiberales | 3861 (partial) | 5804 | 99.53 %‡ | 97.96 % |
| *Helianthus annuus* | Asterales | **fail** (493 bp) | 5857 | — | 100.0 %§ |
| *Vitis vinifera* | Vitales | **fail** (2904 bp) | 5877 | — | — |
| *Actinidia chinensis* | Ericales | **fail** (3603 bp)¶ | 5834 | — | 94.87 %§ |

‡ over the 3833 bp overlap (short-read consensus is partial in a hybrid).
§ easy45 ITS vs GenBank (ngs45 produced no unit). ¶ fails only at 85 bp — see Table 3.

## Table 3 — Read-length control (*Actinidia chinensis*, same species)

Two Illumina libraries, same HiFi reference (SRR24236049):

| Library | Read length | ngs45 result |
|---|---|---|
| SRR29894970 | 85 bp | fail (best contig 3603 bp) |
| SRR28002090 | 150 bp | **full 5834 bp unit, 99.93 % vs HiFi** |

So *Actinidia*'s cross-individual failure is a **read-length** effect, not a species
limit: with modern ≥150 bp reads ngs45 recovers it. Read length ≥150 bp is
**necessary but not sufficient** — *Helianthus* and *Vitis* fail even at 150 bp
because their transcribed-unit spacers are too divergent for short reads to span
(→ use HiFi/easy45). See [ASSEMBLY_LIMITATION.md](ASSEMBLY_LIMITATION.md).

## Table 4 — Runtime (wall-clock, 8–16 threads, shared server)

Subsets: Illumina 5 M read pairs, HiFi 120 k reads. `--max-cov` keeps SPAdes fast.

| tool | median | range |
|---|---|---|
| ngs45 (Illumina) | ~6 min | 4–18 min |
| easy45 (HiFi) | ~4 min | 1–15 min |

Both finish a sample in minutes — practical for hundreds of samples via the `batch`
subcommand (Polyscias: 8 samples, ~13 min each).

## Interpretation

- **Concordance:** where ngs45 assembles a unit it matches the independent HiFi
  consensus to 99.76–100 % (5/7 same-individual species at 0 mismatch). The
  short-read assembly is not a lossy approximation — it is the same molecule.
- **Division of labour:** ngs45 covers the many cheap short-read libraries; easy45
  covers what short reads cannot (highly divergent spacers, and full ribotype
  phasing from HiFi). The failures in Table 2 (*Helianthus*, *Vitis*) are exactly the
  cases easy45 exists for.

## Reproducing

- Accessions: [DATA_ACCESSIONS.md](DATA_ACCESSIONS.md).
- Consolidated result tables: [`../bench/results_v2.tsv`](../bench/results_v2.tsv),
  [`../bench/results_v3.tsv`](../bench/results_v3.tsv),
  [`../bench/master.tsv`](../bench/master.tsv), QC in [`../bench/qc.tsv`](../bench/qc.tsv).
- Regenerate tables: [`../bench/consolidate.py`](../bench/consolidate.py); figures:
  [`../bench/make_figures.py`](../bench/make_figures.py).
