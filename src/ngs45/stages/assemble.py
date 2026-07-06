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


def _n_seqs(path) -> int:
    """Read count via seqkit stats (num_seqs, col 4)."""
    cp = sh(["seqkit", "stats", "-T", str(path)])
    lines = cp.stdout.strip().splitlines()
    try:
        return int(lines[1].split("\t")[3])
    except (ValueError, IndexError):
        return 0


def _cap_coverage(config: Config, r1, r2, readlen):
    """Downsample the baited reads to ~assemble_max_cov x of a nominal rDNA repeat.

    The rDNA array is a very high-copy tandem repeat, so baiting recruits 10^4-10^5 x
    coverage — useless for assembly and the dominant time cost (SPAdes ~35x slower).
    We cap the *assembly input* only; the full bait set is untouched for S6 depth.
    Returns (r1, r2) possibly pointing at downsampled copies.
    """
    if config.assemble_max_cov <= 0 or readlen <= 0:
        return r1, r2
    n_pairs = _n_seqs(r1)
    if n_pairs <= 0:
        return r1, r2
    ref = config.unit_min_len * 2.5          # nominal full repeat (~10 kb) for coverage
    per_pair = readlen * (2 if r2 is not None else 1)
    cov = n_pairs * per_pair / ref
    if cov <= config.assemble_max_cov:
        return r1, r2
    frac = round(config.assemble_max_cov / cov, 6)
    ds1 = config.workdir / "s2_ds_R1.fastq.gz"
    sh(["seqkit", "sample", "-p", str(frac), "-s", "100", "-o", str(ds1), str(r1)])
    ds2 = None
    if r2 is not None:
        ds2 = config.workdir / "s2_ds_R2.fastq.gz"
        sh(["seqkit", "sample", "-p", str(frac), "-s", "100", "-o", str(ds2), str(r2)])
    log.info("S2: baited depth ~%.0fx of a ~%.0f kb repeat -> downsampled %d -> ~%d "
             "read pairs (cap %dx) for assembly", cov, ref / 1000, n_pairs,
             int(n_pairs * frac), config.assemble_max_cov)
    return ds1, ds2


def run(config: Config, state: dict) -> dict:
    r1 = state["bait_r1"]
    r2 = state.get("bait_r2")
    spades_dir = config.workdir / "s2_spades"

    readlen = _read_length(r1)
    klist = _klist_for(readlen, config.spades_klist)
    r1, r2 = _cap_coverage(config, r1, r2, readlen)

    # Force Phred+33: modern Illumina is always Phred+33, but BayesHammer's
    # auto-detection ("Failed to determine offset!") aborts on libraries with a
    # narrow quality range. Pinning the offset makes assembly robust across runs.
    cmd = ["spades.py", "-o", str(spades_dir), "-k", klist, "-t", str(config.threads),
           "--phred-offset", "33"]
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
    if not contigs.exists() or contigs.stat().st_size == 0:
        raise RuntimeError(
            "S2: SPAdes produced no contigs — the recruited read set was too "
            "sparse or malformed to assemble. Check that the input reads are "
            "intact and R1/R2 are in sync (equal read counts).")
    log.info("S2: contigs=%s scaffolds=%s", contigs.name,
             scaffolds.name if scaffolds.exists() else "none")
    return {"contigs": contigs,
            "scaffolds": scaffolds if scaffolds.exists() else None,
            "graph": graph if graph.exists() else None,
            "spades_dir": spades_dir}
