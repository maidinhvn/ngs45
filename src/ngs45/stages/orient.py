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
