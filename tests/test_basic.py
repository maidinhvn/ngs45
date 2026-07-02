"""Lightweight unit tests for ngs45 pure-Python helpers (no external tools)."""

from __future__ import annotations

from pathlib import Path

from ngs45.config import Config, DEFAULT_SEED
from ngs45.io import read_fasta, write_fasta, revcomp, detect_format


def test_revcomp():
    assert revcomp("ACGT") == "ACGT"
    assert revcomp("AAAC") == "GTTT"
    assert revcomp("NnACgt") == "acGTnN"


def test_fasta_roundtrip(tmp_path: Path):
    recs = [("seq1 desc", "ACGT" * 20), ("seq2", "GGGCCC")]
    p = tmp_path / "x.fasta"
    write_fasta(recs, p)
    back = read_fasta(p)
    assert back[0][0] == "seq1 desc"
    assert back[0][1] == "ACGT" * 20
    assert back[1][1] == "GGGCCC"
    assert detect_format(p) == "fasta"


def test_seed_present_and_parseable():
    assert DEFAULT_SEED.exists()
    recs = read_fasta(DEFAULT_SEED)
    assert len(recs) == 1
    name, seq = recs[0]
    assert "45S" in name
    assert len(seq) > 8000        # one Arabidopsis 45S unit incl. IGS


def test_config_defaults_and_paths():
    c = Config(reads1="r1.fq.gz", reads2="r2.fq.gz")
    assert c.workdir == c.outdir / "work"
    assert c.spades_klist == "21,33,55,77,99,127"
    assert Path(c.seed_ref) == DEFAULT_SEED
    c2 = Config(reads1="r1.fq", spades_k="33,55")
    assert c2.spades_klist == "33,55"
    assert c2.reads2 is None
