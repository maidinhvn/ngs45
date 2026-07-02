"""S0 — optional quality / adapter trimming with cutadapt.

Only runs when ``--trim`` is set (see pipeline gating). rDNA is high-copy, so
even aggressive trimming leaves ample coverage. Illumina TruSeq adapters are
trimmed from both mates; low-quality 3' ends and too-short reads are removed.

Input keys:  none (uses config.reads1/reads2)
Output keys: {"reads1": <trimmed R1>, "reads2": <trimmed R2 or None>}
"""

from __future__ import annotations

import logging

from ..config import Config
from ..external import run as sh

log = logging.getLogger("ngs45")

# Illumina TruSeq universal adapters (read-through into adapter on short inserts).
_ADAPT_FWD = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
_ADAPT_REV = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"


def run(config: Config, state: dict) -> dict:
    r1_in = state.get("reads1", config.reads1)
    r2_in = state.get("reads2", config.reads2)
    out1 = config.workdir / "s0_trim_R1.fastq.gz"

    cmd = [
        "cutadapt", "-j", str(config.threads),
        "-q", str(config.trim_quality),
        "-m", str(config.min_read_len),
        "-a", _ADAPT_FWD,
    ]
    if r2_in is not None:
        out2 = config.workdir / "s0_trim_R2.fastq.gz"
        cmd += ["-A", _ADAPT_REV, "-o", str(out1), "-p", str(out2),
                str(r1_in), str(r2_in)]
    else:
        out2 = None
        cmd += ["-o", str(out1), str(r1_in)]

    log.info("S0: cutadapt trim (q>=%d, min_len=%d)", config.trim_quality, config.min_read_len)
    sh(cmd)
    return {"reads1": out1, "reads2": out2}
