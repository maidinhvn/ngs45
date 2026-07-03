# Assembly QC — internal tandem-duplication check (S4)

Short-read de-novo assembly of a high-copy tandem array occasionally produces a
unit with a **spurious internal tandem duplication**: the mature-boundary trim
(S4, `cmsearch` LSU model) can over-extend across a unit junction when the
neighbouring copy contributes an identical stretch, duplicating a short segment
inside the 26S. Because the rRNA genes are single-copy, any such perfect,
*adjacent* internal repeat is a mis-assembly, not biology.

## Detection & correction
After the S4 trim, ngs45 self-aligns the unit (`blastn` unit-vs-unit) and, while
an off-diagonal hit is found that is

- at least `dup_min_len` bp (default 40),
- at least `dup_min_ident` % identical (default 97), and
- **adjacent** — the second copy starts immediately after the first
  (self-alignment offset == copy length),

it drops one copy. The removed length(s) are logged (`S4: collapsed … tandem
duplication`), written to `report.txt`, and summed into the `qc_tandem_dup_bp`
column of `summary.tsv`. Clean units are never modified.

This is the read-remap / self-validation idea used by mature organelle
assemblers (e.g. NOVOPlasty breaks contigs at mis-assembled repeat boundaries);
here the conserved single-copy nature of rRNA lets a self-alignment alone flag
the artifact deterministically.

## Worked example — *Polyscias cumingiana* (IM250605-1)
The S3 monomer for this sample was a 12,386 bp scaffold spanning ~one repeat plus
the 3′ tail of the neighbouring 26S. The S4 LSU trim over-extended by exactly one
127 bp motif, giving a 5,944 bp unit whose 26S carried a perfect 127 bp tandem
duplication (`GTGACGCGCATG…GACTCTAGT` × 2, 100 % identity). Read-depth across the
duplicated segment dropped ~50 % (single-copy reads spread over two identical
copies) — the mis-assembly signature.

The QC step collapses the duplicate to the correct **5,817 bp** unit, which is
99.8 % identical (11 mismatches, 0 gaps) to the congeneric *P. balfouriana* unit
— exactly the within-genus conservation expected. The congener assembled cleanly
and is left untouched.

## Limitation / future work
This catches duplications that survive into the final unit. The gold-standard
resolution is graph-based (one cycle of the assembly graph = one repeat unit, as
in Ribotin / GetOrganelle), which needs the assembly graph and, for phasing
divergent ribotype *morphs*, long reads (HiFi/ONT). ngs45 targets short reads, so
it validates the assembled unit rather than phasing morphs.
