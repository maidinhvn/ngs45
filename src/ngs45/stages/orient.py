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
    write_fasta([("nrDNA_45S_unit", unit)], out)
    log.info("S4: CM mature-boundary trim -> transcribed unit %d bp "
             "(18S 5' .. 26S 3'); full repeat %d bp", len(unit), period)
    return {"monomer": out, "full_repeat": full_out}
