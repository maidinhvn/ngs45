"""S3 — resolve a single 45S repeat unit (monomer) from the assembly.

rDNA is a tandem array, so SPAdes typically collapses it into one (or few)
high-coverage contigs that are either ~one unit long or several units of the
same sequence concatenated. We recover exactly one period:

  1. BLAST the 45S seed against the contigs to find the rDNA-bearing contig(s).
  2. Self-BLAST that contig: a tandem array self-aligns off the main diagonal at
     an offset equal to the repeat period. The smallest consistent offset in the
     plausible unit-length window is the monomer length.
  3. Cut one period. If the contig is already ~one unit (no internal repeat), keep
     it whole.

Input keys:  {"contigs"}
Output keys: {"monomer_raw", "unit_len", "rdna_contig"}
"""

from __future__ import annotations

import logging
from collections import Counter

from ..config import Config
from ..external import run as sh
from ..io import read_fasta, write_fasta

log = logging.getLogger("ngs45")

_FMT = "6 qseqid sseqid pident length qstart qend sstart send bitscore"


def _blast(query, subject, workdir, tag, perc_identity=None, extra=None) -> list[dict]:
    """Run blastn query-vs-subject, return parsed tabular rows."""
    out = workdir / f"s3_{tag}.tsv"
    cmd = ["blastn", "-query", str(query), "-subject", str(subject),
           "-outfmt", _FMT, "-dust", "no"]
    if perc_identity is not None:
        cmd += ["-perc_identity", str(perc_identity)]
    if extra:
        cmd += extra
    cp = sh(cmd)
    out.write_text(cp.stdout)
    rows = []
    for line in cp.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 9:
            continue
        rows.append({
            "q": f[0], "s": f[1], "pid": float(f[2]), "len": int(f[3]),
            "qs": int(f[4]), "qe": int(f[5]), "ss": int(f[6]), "se": int(f[7]),
            "bits": float(f[8]),
        })
    return rows


def _pick_rdna_contig(config: Config, sources) -> tuple[str, str]:
    """Pick the best rDNA sequence across several FASTAs (contigs + scaffolds).

    Among sequences carrying a substantial seed hit (>=800 bp aligned, i.e. they
    hold most of the conserved genes) we take the *longest* — that is the one most
    likely to span a whole unit (a scaffold often joins gene contigs across the
    spacer gaps that the raw contigs leave broken). If nothing clears that bar we
    fall back to the single most seed-aligned contig.
    """
    best = None           # (seqlen, aligned_bp, id, seq, label)
    fallback = None       # (aligned_bp, id, seq, label)
    n_hit = 0
    for label, fasta in sources:
        if fasta is None:
            continue
        rows = _blast(config.seed_ref, fasta, config.workdir, f"seed_vs_{label}",
                      perc_identity=config.min_gene_ident)
        if not rows:
            continue
        aligned = Counter()
        for r in rows:
            aligned[r["s"]] += r["len"]
        n_hit += len(aligned)
        seqs = dict(_seq_by_id(fasta))
        for cid, bp in aligned.items():
            seq = seqs.get(cid, "")
            if fallback is None or bp > fallback[0]:
                fallback = (bp, cid, seq, label)
            if bp >= 800 and (best is None or len(seq) > best[0]):
                best = (len(seq), bp, cid, seq, label)

    if best is not None:
        slen, bp, cid, seq, label = best
        log.info("S3: rDNA sequence = %s [%s] %d bp (%d bp seed-aligned; %d hits total)",
                 cid, label, slen, bp, n_hit)
        return cid, seq
    if fallback is not None:
        bp, cid, seq, label = fallback
        log.info("S3: rDNA sequence = %s [%s] %d bp (weak: %d bp seed-aligned)",
                 cid, label, len(seq), bp)
        return cid, seq
    raise RuntimeError("S3: no contig/scaffold aligned to the 45S seed — assembly "
                       "may have failed or reads were off-target.")


def _seq_by_id(fasta):
    for name, seq in read_fasta(fasta):
        yield name.split()[0], seq


def _tandem_period(config: Config, contig_id, seq) -> int | None:
    """Smallest self-alignment offset in the unit-length window, or None."""
    tmp = config.workdir / "s3_contig.fasta"
    write_fasta([(contig_id, seq)], tmp)
    rows = _blast(tmp, tmp, config.workdir, "self", perc_identity=95)
    offsets = []
    for r in rows:
        off = r["ss"] - r["qs"]
        if off <= 0:                                   # keep upper triangle only
            continue
        if r["len"] < config.unit_min_len // 2:        # ignore short incidental hits
            continue
        if config.unit_min_len <= off <= config.unit_max_len:
            offsets.append(off)
    if not offsets:
        return None
    # The fundamental period is the smallest well-supported offset; round to
    # tolerate a few indels between copies before taking the mode.
    binned = Counter(round(o / 10) * 10 for o in offsets)
    period = min(o for o, _ in binned.most_common() if binned[o] == max(binned.values()))
    return period


def run(config: Config, state: dict) -> dict:
    sources = [("contigs", state["contigs"]),
               ("scaffolds", state.get("scaffolds"))]
    cid, seq = _pick_rdna_contig(config, sources)

    if len(seq) < config.unit_min_len:
        raise RuntimeError(
            f"S3: best rDNA sequence {cid} is only {len(seq)} bp (< {config.unit_min_len}); "
            "assembly did not span a unit (short reads through divergent spacers). "
            "See docs/ASSEMBLY_LIMITATION.md.")

    period = _tandem_period(config, cid, seq)
    if period:
        monomer = seq[:period]
        log.info("S3: tandem period = %d bp -> cut one monomer", period)
    elif len(seq) <= config.unit_max_len:
        monomer = seq
        log.info("S3: contig is ~one unit (%d bp), no internal repeat -> keep whole", len(seq))
    else:
        monomer = seq[:config.unit_max_len]
        log.warning("S3: no clean tandem period on a %d bp contig; truncated to "
                    "%d bp (inspect manually).", len(seq), config.unit_max_len)

    out = config.workdir / "s3_monomer_raw.fasta"
    write_fasta([(f"{cid}_monomer len={len(monomer)}", monomer)], out)
    return {"monomer_raw": out, "unit_len": len(monomer), "rdna_contig": cid}
