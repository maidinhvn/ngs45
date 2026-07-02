"""I/O helpers: gzip-aware open, format sniffing, small FASTA utilities.

Sequence parsing uses Biopython (the only non-trivial pip dependency). Kept in
one module so format handling lives in one place.
"""

from __future__ import annotations

import gzip
from pathlib import Path


def open_maybe_gzip(path, mode: str = "rt"):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, mode)
    return open(path, mode)


def detect_format(path) -> str:
    """Return 'fasta' or 'fastq' by sniffing the first record character."""
    with open_maybe_gzip(path) as fh:
        first = fh.read(1)
    if first == ">":
        return "fasta"
    if first == "@":
        return "fastq"
    raise ValueError(f"Unrecognised sequence format: {path}")


def read_fasta(path) -> list[tuple[str, str]]:
    """Return [(header_without_gt, sequence), ...]; header keeps the full line."""
    records: list[tuple[str, str]] = []
    name = None
    chunks: list[str] = []
    with open_maybe_gzip(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks)))
                name = line[1:]
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        records.append((name, "".join(chunks)))
    return records


def write_fasta(records, path, width: int = 70) -> None:
    """Write [(header, seq), ...] as wrapped FASTA."""
    with open(path, "w") as out:
        for name, seq in records:
            out.write(f">{name}\n")
            for i in range(0, len(seq), width):
                out.write(seq[i:i + width] + "\n")


_COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]
