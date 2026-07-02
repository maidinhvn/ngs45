"""S2 — assemble the recruited reads with SPAdes.

Two lessons adopted from mature short-read organelle/nr assemblers (e.g.
GetOrganelle) so ngs45 stays self-contained yet competitive:

  1. **k-mers must be < read length.** A fixed 21..127 ladder is wrong for short
     libraries (k=127 is meaningless for 73 bp reads). We detect the read length
     and keep only k values that fit, so SPAdes builds a usable graph.
  2. **No `--careful` by default.** Its mismatch corrector collapses the many
     divergent rDNA ribotypes and fragments the unit; plain multi-k assembly keeps
     the graph intact. (`--careful` stays available via config for clean isolates.)

S3 later picks the unit from whichever of contigs.fasta / scaffolds.fasta carries
the longest rDNA-spanning sequence (scaffolds join contigs across spacer gaps with
paired-end links, so they often span a unit the raw contigs do not).

Input keys:  {"bait_r1", "bait_r2"}
Output keys: {"contigs", "scaffolds", "graph", "spades_dir"}
"""

from __future__ import annotations

import logging

from ..config import Config
from ..external import run as sh

log = logging.getLogger("ngs45")

_K_LADDER = [21, 33, 55, 77, 99, 127]


def _read_length(path) -> int:
    """Median-ish read length via seqkit stats (avg_len, col 7)."""
    cp = sh(["seqkit", "stats", "-T", str(path)])
    lines = cp.stdout.strip().splitlines()
    if len(lines) < 2:
        return 0
    try:
        return int(round(float(lines[1].split("\t")[6])))
    except (ValueError, IndexError):
        return 0


def _klist_for(readlen: int, configured: str) -> str:
    """k values that fit the reads (k <= readlen-6, odd). Falls back to configured."""
    if readlen <= 0:
        return configured
    ks = [k for k in _K_LADDER if k <= readlen - 6]
    if not ks:
        ks = [21]
    return ",".join(str(k) for k in ks)


def run(config: Config, state: dict) -> dict:
    r1 = state["bait_r1"]
    r2 = state.get("bait_r2")
    spades_dir = config.workdir / "s2_spades"

    readlen = _read_length(r1)
    klist = _klist_for(readlen, config.spades_klist)

    cmd = ["spades.py", "-o", str(spades_dir), "-k", klist, "-t", str(config.threads)]
    if config.spades_careful:
        cmd += ["--careful"]
    if r2 is not None:
        cmd += ["-1", str(r1), "-2", str(r2)]
    else:
        cmd += ["-s", str(r1)]

    log.info("S2: SPAdes assembly (readlen~%d -> k=%s, careful=%s)",
             readlen, klist, config.spades_careful)
    sh(cmd)

    contigs = spades_dir / "contigs.fasta"
    scaffolds = spades_dir / "scaffolds.fasta"
    graph = spades_dir / "assembly_graph_with_scaffolds.gfa"
    if not graph.exists():
        graph = spades_dir / "assembly_graph.gfa"
    if not contigs.exists():
        raise RuntimeError(f"SPAdes produced no contigs at {contigs}")
    log.info("S2: contigs=%s scaffolds=%s", contigs.name,
             scaffolds.name if scaffolds.exists() else "none")
    return {"contigs": contigs,
            "scaffolds": scaffolds if scaffolds.exists() else None,
            "graph": graph if graph.exists() else None,
            "spades_dir": spades_dir}
