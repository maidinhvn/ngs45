"""S1 — iterative baiting of rDNA reads (GetOrganelle/MITObim-style).

Short reads can't span a repeat, so we recover the unit by *growing* a read set:

  round 0: map the whole library against the conserved 45S seed. The 18S/5.8S/26S
           genes are conserved across angiosperms and recruit the coding-region
           reads; keeping both mates of every aligned pair pulls in the spacer-side
           mates that overhang into ETS/ITS/IGS.
  round n: rebuild the bait from seed + reads recruited so far, remap. Each round
           walks a read-length further into the divergent spacers, so ITS/ETS/IGS
           reads that never matched the Arabidopsis seed get recruited via overlap.

Stops when a round adds < ``bait_converge`` new reads, or after ``bait_rounds``.

Input keys:  {"reads1", "reads2"} (optional; else config.reads*)
Output keys: {"bait_r1", "bait_r2", "bait_n_pairs"}
"""

from __future__ import annotations

import logging
import subprocess

from ..config import Config
from ..external import run as sh
from ..io import read_fasta, write_fasta

log = logging.getLogger("ngs45")


def _aligned_names(sam_path) -> set[str]:
    """Collect QNAMEs of mapped records from a --no-unal SAM (skip @headers)."""
    names: set[str] = set()
    with open(sam_path) as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            q, flag = line.split("\t", 2)[:2]
            if int(flag) & 4:            # unmapped (belt-and-suspenders; --no-unal already drops these)
                continue
            names.add(q)
    return names


def _recruit_once(config: Config, ref_fasta, r1, r2, round_dir) -> set[str]:
    """Build a bowtie2 index for ``ref_fasta`` and return names of aligned reads."""
    round_dir.mkdir(parents=True, exist_ok=True)
    idx = round_dir / "bait_idx"
    sh(["bowtie2-build", "--threads", str(config.threads), "-q",
        str(ref_fasta), str(idx)])

    sam = round_dir / "aln.sam"
    # --local + a constant min-score = 2 * bait_min_len (local match bonus is 2),
    # i.e. recruit any read with >= bait_min_len bp aligned to the current bait.
    cmd = ["bowtie2", "-x", str(idx), "--local", "--no-unal",
           "-k", "1", "-p", str(config.threads),
           "--score-min", f"C,{2 * config.bait_min_len},0"]
    if r2 is not None:
        cmd += ["-1", str(r1), "-2", str(r2)]
    else:
        cmd += ["-U", str(r1)]
    with open(sam, "w") as out:
        subprocess.run(cmd, check=True, stdout=out, stderr=subprocess.PIPE, text=True)
    return _aligned_names(sam)


def _extract_pairs(config: Config, names_file, r1, r2, out1, out2) -> None:
    """Pull both mates of every recruited read name from the original library."""
    sh(["seqkit", "grep", "-j", str(config.threads), "-f", str(names_file),
        str(r1), "-o", str(out1)])
    if r2 is not None:
        sh(["seqkit", "grep", "-j", str(config.threads), "-f", str(names_file),
            str(r2), "-o", str(out2)])


def _fq_to_fasta_records(config: Config, *fastqs) -> list:
    """Concatenate FASTQ files and return [(name, seq)] via seqkit fq2fa."""
    records: list = []
    for i, fq in enumerate(fastqs):
        if fq is None:
            continue
        cp = sh(["seqkit", "fq2fa", str(fq)])
        # parse the FASTA text seqkit emitted on stdout
        name = None
        seq: list[str] = []
        for line in cp.stdout.splitlines():
            if line.startswith(">"):
                if name is not None:
                    records.append((f"m{i}_{name}", "".join(seq)))
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
        if name is not None:
            records.append((f"m{i}_{name}", "".join(seq)))
    return records


def run(config: Config, state: dict) -> dict:
    r1 = state.get("reads1", config.reads1)
    r2 = state.get("reads2", config.reads2)

    ref = config.seed_ref
    prev_n = 0
    bait_r1 = config.workdir / "s1_bait_R1.fastq.gz"
    bait_r2 = config.workdir / "s1_bait_R2.fastq.gz" if r2 is not None else None
    names: set[str] = set()

    for rnd in range(config.bait_rounds):
        round_dir = config.workdir / f"s1_round{rnd}"
        names = _recruit_once(config, ref, r1, r2, round_dir)
        n = len(names)
        grew = (n - prev_n) / n if n else 0.0
        log.info("S1 round %d: %d reads recruited (+%.1f%% new)", rnd, n, grew * 100)
        if n == 0:
            raise RuntimeError(
                "S1 recruited 0 reads — check --seed-ref matches your reads' "
                "kingdom, or lower --bait-min-len.")

        names_file = round_dir / "names.txt"
        names_file.write_text("\n".join(sorted(names)) + "\n")
        _extract_pairs(config, names_file, r1, r2,
                       round_dir / "rec_R1.fastq.gz",
                       round_dir / "rec_R2.fastq.gz" if r2 else None)

        # converged, or last round: this round's reads are the final bait set
        if rnd == config.bait_rounds - 1 or (prev_n and grew < config.bait_converge):
            round_dir.joinpath("rec_R1.fastq.gz").replace(bait_r1)
            if r2 is not None:
                round_dir.joinpath("rec_R2.fastq.gz").replace(bait_r2)
            break

        # else: grow the bait = seed + recruited reads, and go again
        ext_ref = round_dir / "bait_ref.fasta"
        seed_records = read_fasta(config.seed_ref)
        rec = _fq_to_fasta_records(config, round_dir / "rec_R1.fastq.gz",
                                   round_dir / "rec_R2.fastq.gz" if r2 else None)
        write_fasta(seed_records + rec, ext_ref)
        ref = ext_ref
        prev_n = n

    # optional subsample to cap depth (rDNA can be extremely deep)
    if config.subsample and _count_seqs(config, bait_r1) > config.subsample:
        log.info("S1: subsampling to %d pairs", config.subsample)
        _subsample_pairs(config, bait_r1, bait_r2)

    n_pairs = _count_seqs(config, bait_r1)
    log.info("S1: final bait set = %d read pairs", n_pairs)
    return {"bait_r1": bait_r1, "bait_r2": bait_r2, "bait_n_pairs": n_pairs}


def _count_seqs(config: Config, fq) -> int:
    if fq is None:
        return 0
    cp = sh(["seqkit", "stats", "-T", str(fq)])
    # header line + one data line; num_seqs is column 4 (1-based)
    lines = cp.stdout.strip().splitlines()
    if len(lines) < 2:
        return 0
    return int(lines[1].split("\t")[3].replace(",", ""))


def _subsample_pairs(config: Config, r1, r2) -> None:
    """Deterministically subsample both mates to config.subsample, keeping pairs."""
    for fq in (r1, r2):
        if fq is None:
            continue
        tmp = fq.with_suffix(".sub.fastq.gz")
        sh(["seqkit", "sample", "-s", "11", "-n", str(config.subsample),
            str(fq), "-o", str(tmp)])
        tmp.replace(fq)
