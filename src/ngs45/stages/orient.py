"""S4 — orient and trim to the mature 45S transcribed unit (18S 5' -> 26S 3').

The raw monomer is one repeat period cut at an arbitrary phase and possibly on
the minus strand. We reduce it to the mature transcribed region — the conserved,
alignable 18S-ITS1-5.8S-ITS2-26S that downstream comparison / barcoding wants —
dropping the ETS/IGS.

Boundaries come from Rfam SSU+LSU covariance models via ``cmsearch`` (the same
bundled CM easy45 uses), so both tools define the mature 18S 5' / 26S 3' termini
*identically* and their units are directly comparable. To be robust to the phase
of the S3 cut (which may fall inside a gene), we search a doubled copy of the
period and take the first full SSU->LSU span.

The full repeat (with ETS/IGS) is written under work/ as a best-effort
intermediate, not a headline output.

Input keys:  {"monomer_raw"}
Output keys: {"monomer", "full_repeat"}
"""

from __future__ import annotations

import logging
import subprocess

from ..config import Config
from ..io import read_fasta, write_fasta, revcomp

log = logging.getLogger("ngs45")


def _cmsearch_trunc(config: Config, seq: str, tag: str):
    """Like ``_cmsearch`` but also captures the truncation flag (tblout col 11).

    Returns [(kind, sfrom, sto, strand, score, trunc)] where ``trunc`` is
    cmsearch's own flag: ``"no"``, ``"5'"``, ``"3'"`` or ``"5'&3'"``. A gene end
    left on an adjacent contig makes the SSU/LSU hit report ``5'``/``3'``.
    """
    fa = config.workdir / f"s4_{tag}.fasta"
    tbl = config.workdir / f"s4_{tag}.tbl"
    write_fasta([("s", seq)], fa)
    subprocess.run(["cmsearch", "--noali", "--cpu", str(config.threads),
                    "--tblout", str(tbl), str(config.cm_ref), str(fa)],
                   check=True, capture_output=True, text=True)
    hits = []
    with open(tbl) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) < 15:
                continue
            model, sfrom, sto = f[2], int(f[7]), int(f[8])
            strand, trunc, score = f[9], f[10], float(f[14])
            kind = "SSU" if "SSU" in model else ("LSU" if "LSU" in model else None)
            if kind:
                hits.append((kind, min(sfrom, sto), max(sfrom, sto),
                             strand, score, trunc))
    return hits


def _rdna_donor_contigs(config: Config, state, exclude_id):
    """rDNA-bearing contigs/scaffolds (>=300 bp seed-aligned) other than the
    monomer's own source contig — candidate donors of a missing gene end."""
    from .resolve import _blast, _seq_by_id  # reuse S3 blast/parse helpers
    donors = []
    for label in ("contigs", "scaffolds"):
        fasta = state.get(label)
        if fasta is None:
            continue
        rows = _blast(config.seed_ref, fasta, config.workdir,
                      f"repair_seed_vs_{label}", perc_identity=config.min_gene_ident)
        if not rows:
            continue
        aligned = {}
        for r in rows:
            aligned[r["s"]] = aligned.get(r["s"], 0) + r["len"]
        seqs = dict(_seq_by_id(fasta))
        for cid, bp in aligned.items():
            if cid == exclude_id or bp < 300:
                continue
            donors.append((cid, seqs.get(cid, "")))
    return donors


def _graft(config: Config, seq: str, donors, end: str):
    """Extend ``seq`` at its 5' or 3' end using the donor contig that overhangs
    it, joined on their overlap. ``end`` is "5" or "3". Returns the extended
    sequence, or the original if no clean overhang is found.
    """
    q = config.workdir / "s4_repair_query.fasta"
    write_fasta([("m", seq)], q)
    best_ext = ""
    for cid, dseq in donors:
        if not dseq:
            continue
        for oriented in (dseq, revcomp(dseq)):
            d = config.workdir / "s4_repair_donor.fasta"
            write_fasta([("d", oriented)], d)
            cp = subprocess.run(
                ["blastn", "-query", str(q), "-subject", str(d),
                 "-outfmt", "6 qstart qend sstart send length pident", "-dust", "no"],
                check=True, capture_output=True, text=True)
            for line in cp.stdout.splitlines():
                f = line.split("\t")
                if len(f) < 6:
                    continue
                qs, qe, ss, se, ln = (int(x) for x in f[:5])
                pid = float(f[5])
                # adjacent SPAdes contigs overlap by ~k (up to 127 bp); require a
                # clean, high-identity + overlap rather than a long one
                if ss > se or ln < 80 or pid < 98.0:
                    continue
                if end == "5" and qs <= 5:
                    overhang = ss - qs                    # donor bases before seq start
                    cand = oriented[:overhang]
                    if 0 < overhang and len(cand) > len(best_ext) \
                            and len(cand) + len(seq) <= config.unit_max_len:
                        best_ext = cand
                elif end == "3" and qe >= len(seq) - 5:
                    overhang = len(oriented) - se         # donor bases after seq end
                    cand = oriented[se:]
                    if 0 < overhang and len(cand) > len(best_ext) \
                            and len(cand) + len(seq) <= config.unit_max_len:
                        best_ext = cand
    if not best_ext:
        return seq
    return best_ext + seq if end == "5" else seq + best_ext


def _repair_truncated_monomer(config: Config, seq: str, state):
    """Graft a conserved gene end (SSU 5' / LSU 3') that SPAdes left on an
    adjacent contig back onto the oriented monomer.

    SPAdes fragments the rDNA array at coverage transitions (ETS/18S), so the
    best contig can miss the 18S 5' or 26S 3' end; S3 then yields a truncated
    unit. cmsearch flags this (trunc 5'/3'); when it does *and* no full-length
    gene hit exists on the monomer, we stitch the missing end from a donor
    contig. A no-op for the usual case (monomer already spans both genes).
    Returns (seq, repaired?).
    """
    # Probe the SINGLE monomer (not doubled): a good assembly cut mid-gene shows
    # a 5'- AND a 3'-truncated hit for that gene (the two ends reconstitute it),
    # so it is NOT flagged; only a genuinely absent gene end leaves a lone
    # one-sided truncation. (Doubling would forge a spurious full hit across the
    # ETS/IGS junction.)
    hits = _cmsearch_trunc(config, seq, "repair_probe")
    ssu = [h for h in hits if h[0] == "SSU" and h[3] == "+"]
    lsu = [h for h in hits if h[0] == "LSU" and h[3] == "+"]
    need5 = bool(ssu) and all("5'" in h[5] for h in ssu)   # only 5'-truncated SSU
    need3 = bool(lsu) and all("3'" in h[5] for h in lsu)   # only 3'-truncated LSU
    if not (need5 or need3):
        return seq, False
    donors = _rdna_donor_contigs(config, state, state.get("rdna_contig"))
    if not donors:
        return seq, False
    repaired = False
    if need5:
        new = _graft(config, seq, donors, "5")
        if len(new) > len(seq):
            log.info("S4: SSU 5' end was on an adjacent contig; grafted %d bp "
                     "(SPAdes split the array) -> monomer %d bp",
                     len(new) - len(seq), len(new))
            seq, repaired = new, True
    if need3:
        new = _graft(config, seq, donors, "3")
        if len(new) > len(seq):
            log.info("S4: LSU 3' end was on an adjacent contig; grafted %d bp "
                     "-> monomer %d bp", len(new) - len(seq), len(new))
            seq, repaired = new, True
    return seq, repaired


def _cmsearch(config: Config, seq: str, tag: str):
    """Run cmsearch(SSU+LSU) on `seq`; return [(kind, sfrom, sto, strand, score)]."""
    fa = config.workdir / f"s4_{tag}.fasta"
    tbl = config.workdir / f"s4_{tag}.tbl"
    write_fasta([("s", seq)], fa)
    subprocess.run(["cmsearch", "--noali", "--cpu", str(config.threads),
                    "--tblout", str(tbl), str(config.cm_ref), str(fa)],
                   check=True, capture_output=True, text=True)
    hits = []
    with open(tbl) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) < 15:
                continue
            model, sfrom, sto, strand, score = f[2], int(f[7]), int(f[8]), f[9], float(f[14])
            kind = "SSU" if "SSU" in model else ("LSU" if "LSU" in model else None)
            if kind:
                hits.append((kind, min(sfrom, sto), max(sfrom, sto), strand, score))
    return hits


def _best(hits, kind, strand="+"):
    cand = [h for h in hits if h[0] == kind and h[3] == strand]
    return max(cand, key=lambda h: h[4]) if cand else None


def _collapse_tandem_dup(config: Config, seq: str):
    """Collapse spurious internal tandem duplications in the unit.

    rRNA genes are single-copy, so a near-identical segment that sits *immediately*
    after its source (diagonal self-alignment offset == copy length) is a
    mis-assembly: SPAdes/repeat-resolution occasionally over-extends the LSU trim
    across a unit junction, duplicating a short stretch (see docs/ASSEMBLY QC). We
    self-BLAST the unit and, while such an adjacent internal repeat >= dup_min_len
    at >= dup_min_ident%% remains, drop one copy. Returns (clean_seq, [removed_len]).
    """
    removed = []
    fa = config.workdir / "s4_selfdup.fasta"
    for _ in range(20):                       # bounded; one copy removed per pass
        write_fasta([("u", seq)], fa)
        cp = subprocess.run(
            ["blastn", "-query", str(fa), "-subject", str(fa),
             "-outfmt", "6 qstart qend sstart send length pident", "-dust", "no"],
            check=True, capture_output=True, text=True)
        best = None                            # (offset, qstart) of shortest tandem
        for line in cp.stdout.splitlines():
            f = line.split("\t")
            if len(f) < 6:
                continue
            qs, qe, ss, se, ln = (int(x) for x in f[:5])
            pid = float(f[5])
            off = ss - qs                      # diagonal offset (upper triangle)
            if (off > 0 and ln >= config.dup_min_len and pid >= config.dup_min_ident
                    and abs(ss - (qe + 1)) <= max(3, off // 5)):  # 2nd copy adjacent
                if best is None or off < best[0]:
                    best = (off, qs)
        if best is None:
            break
        off, qs = best
        seq = seq[:qs - 1] + seq[qs - 1 + off:]
        removed.append(off)
    return seq, removed


def run(config: Config, state: dict) -> dict:
    _name, seq = read_fasta(state["monomer_raw"])[0]
    out = config.workdir / "s4_monomer.fasta"
    full_out = config.workdir / "full_repeat.fasta"

    # 1. orient: find the strand of the best SSU on the raw period
    hits = _cmsearch(config, seq, "orient")
    ssu_any = max([h for h in hits if h[0] == "SSU"], key=lambda h: h[4], default=None)
    if ssu_any is None:
        log.warning("S4: cmsearch found no SSU on the monomer; leaving as-is")
        write_fasta([("nrDNA_45S_unit", seq)], out)
        write_fasta([("nrDNA_full_repeat", seq)], full_out)
        return {"monomer": out, "full_repeat": full_out}
    if ssu_any[3] == "-":
        seq = revcomp(seq)
        log.info("S4: unit was on minus strand -> reverse-complemented")

    # 1b. if SPAdes split the array and a conserved gene end (18S 5' / 26S 3')
    #     landed on an adjacent contig, graft it back before trimming. No-op when
    #     the monomer already spans both genes (the usual case).
    seq, _repaired = _repair_truncated_monomer(config, seq, state)

    # 2. search a doubled period so a full SSU->LSU span exists regardless of
    #    where S3 cut the monomer
    period = len(seq)
    dbl = seq + seq
    hits = _cmsearch(config, dbl, "double")
    ssu = _best(hits, "SSU")
    if ssu is None:
        log.warning("S4: no plus-strand SSU after orientation; keeping full period")
        write_fasta([("nrDNA_45S_unit", seq)], out)
        write_fasta([("nrDNA_full_repeat", seq)], full_out)
        return {"monomer": out, "full_repeat": full_out}
    ssu_from = ssu[1]

    # 3. the mature 26S 3' end = best LSU that ends downstream of the 18S start,
    #    within one period
    lsu_cands = [h for h in hits
                 if h[0] == "LSU" and h[3] == "+" and h[2] > ssu_from
                 and (h[2] - ssu_from) <= period * 1.2]
    lsu = max(lsu_cands, key=lambda h: h[4]) if lsu_cands else None

    # full repeat = the period rotated to the mature 18S start
    rot = (ssu_from - 1) % period
    write_fasta([("nrDNA_full_repeat", seq[rot:] + seq[:rot])], full_out)

    if lsu is None:
        log.warning("S4: cmsearch found no downstream LSU; emitting full period "
                    "(inspect manually)")
        write_fasta([("nrDNA_45S_unit", seq[rot:] + seq[:rot])], out)
        return {"monomer": out, "full_repeat": full_out}

    unit = dbl[ssu_from - 1:lsu[2]]
    unit, dup_removed = _collapse_tandem_dup(config, unit)
    if dup_removed:
        log.warning("S4: collapsed %d internal tandem duplication(s) %s bp "
                    "(assembly artifact across a unit junction) -> unit %d bp",
                    len(dup_removed), dup_removed, len(unit))
    write_fasta([("nrDNA_45S_unit", unit)], out)
    log.info("S4: CM mature-boundary trim -> transcribed unit %d bp "
             "(18S 5' .. 26S 3'); full repeat %d bp", len(unit), period)
    return {"monomer": out, "full_repeat": full_out, "qc_dup_removed": dup_removed}
