"""Batch mode — run the ngs45 pipeline over every sample in a folder.

Auto-detects two input layouts (both may be mixed):
  (a) flat            — paired FASTQ files in one folder, mates named
                        ``*_R1/_R2`` or ``*_1/_2`` (``.fastq``/``.fq``, optionally
                        gzipped); the sample name is the shared prefix.
  (b) subfolder       — one directory per sample holding its pair; the sample
                        name is the directory name.

Each sample is written to ``outdir/<sample>/``. One sample failing does not stop
the batch (the error is logged and saved). Samples whose unit already exists are
skipped (resume), so an interrupted overnight run can simply be re-launched.
A ``batch_summary.tsv`` aggregates every sample's summary + timing + status.
"""

from __future__ import annotations

import logging
import re
import time
import traceback
from pathlib import Path

from .config import Config, DEFAULT_SEED
from .pipeline import run_pipeline

log = logging.getLogger("ngs45")

_FQ = re.compile(r"\.(fastq|fq)(\.gz)?$", re.I)
_MATE = re.compile(r"^(.*?)[._]R?([12])$")   # sample_R1 / sample_1 / sample.R2 ...

# columns emitted by stages/report.py summary.tsv (kept explicit for a stable table)
_SUMCOLS = ["unit_len", "GC_percent", "bait_pairs", "annot_source", "18S_len",
            "ITS1_len", "5.8S_len", "ITS2_len", "26S_len", "ITS_barcode_len",
            "ribotype_sites", "qc_tandem_dup_bp"]


def _mate_of(fname: str):
    base = _FQ.sub("", fname)
    m = _MATE.match(base)
    return (m.group(1), m.group(2)) if m else (None, None)


def _pairs_in(folder: Path, name: str | None = None):
    """R1/R2 pairs among FASTQ files directly inside `folder`."""
    by_sample: dict = {}
    for f in sorted(p for p in folder.iterdir() if p.is_file() and _FQ.search(p.name)):
        s, mate = _mate_of(f.name)
        if s is not None:
            by_sample.setdefault(s, {})[mate] = f
    out = []
    for s, d in by_sample.items():
        if "1" in d and "2" in d:
            out.append((name or s, d["1"], d["2"]))
    return out


def discover_samples(indir: Path):
    """Return [(sample_name, r1, r2)] auto-detecting flat + subfolder layouts."""
    samples = list(_pairs_in(indir))                        # flat pairs in indir
    for sub in sorted(p for p in indir.iterdir() if p.is_dir()):
        samples += _pairs_in(sub, name=sub.name)            # subfolder-per-sample
    seen, uniq = set(), []
    for s in samples:
        if s[0] not in seen:
            seen.add(s[0])
            uniq.append(s)
    return uniq


def run_batch(args) -> int:
    indir, outroot = Path(args.indir), Path(args.outdir)
    outroot.mkdir(parents=True, exist_ok=True)
    samples = discover_samples(indir)
    if not samples:
        log.error("batch: no paired FASTQ samples found under %s", indir)
        return 1
    resume = not args.no_resume
    log.info("batch: %d sample(s) found under %s", len(samples), indir)

    rows = []
    for i, (name, r1, r2) in enumerate(samples, 1):
        out = outroot / name
        if resume and (out / "nrDNA_45S.fasta").exists():
            log.info("[%d/%d] %s — already done, skipping", i, len(samples), name)
            rows.append((name, out, "skipped(done)", 0))
            continue
        log.info("[%d/%d] %s  (R1=%s)", i, len(samples), name, Path(r1).name)
        t0 = time.time()
        try:
            cfg = Config(
                reads1=r1, reads2=r2, seed_ref=args.seed_ref or DEFAULT_SEED,
                outdir=out, threads=args.threads, trim=args.trim,
                bait_rounds=args.bait_rounds, subsample=args.subsample,
                spades_k=args.spades_k, assemble_max_cov=args.max_cov,
                call_variants=args.call_variants, resume=resume,
            )
            run_pipeline(cfg)
            status = "OK"
        except Exception as e:  # one sample must not sink the batch
            status = f"FAIL:{type(e).__name__}"
            out.mkdir(parents=True, exist_ok=True)
            (out / "batch_error.txt").write_text(traceback.format_exc())
            log.error("[%d/%d] %s FAILED: %s (see batch_error.txt)", i, len(samples), name, e)
        rows.append((name, out, status, int(time.time() - t0)))

    _write_summary(outroot, rows)
    ok = sum(1 for _, _, s, _ in rows if s == "OK")
    log.info("batch: done — %d OK, %d skipped, %d failed of %d",
             ok, sum(1 for r in rows if r[2].startswith("skipped")),
             sum(1 for r in rows if r[2].startswith("FAIL")), len(rows))
    return 0


def _write_summary(outroot: Path, rows):
    out = outroot / "batch_summary.tsv"
    with open(out, "w") as fh:
        fh.write("\t".join(["sample", *_SUMCOLS, "seconds", "status"]) + "\n")
        for name, sdir, status, secs in rows:
            cols = {}
            stsv = Path(sdir) / "summary.tsv"
            if stsv.exists():
                lines = stsv.read_text().splitlines()
                if len(lines) >= 2:
                    cols = dict(zip(lines[0].split("\t"), lines[1].split("\t")))
            fh.write("\t".join([name, *[cols.get(k, "") for k in _SUMCOLS],
                                str(secs), status]) + "\n")
    log.info("batch: wrote %s", out)
