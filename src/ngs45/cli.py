"""Command-line entry point for ngs45."""

from __future__ import annotations

import argparse
import logging
import sys

from . import __version__
from .config import Config, DEFAULT_SEED
from .external import DependencyError, check_dependencies
from .pipeline import run_pipeline


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="ngs45",
        description="Recover the 45S nrDNA transcribed unit from Illumina short reads.",
    )
    p.add_argument("--version", action="version", version=f"ngs45 {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    # --- run ----------------------------------------------------------------
    r = sub.add_parser("run", help="run the full pipeline")
    r.add_argument("-1", "--reads1", required=True, help="Illumina R1 (FASTQ[.gz])")
    r.add_argument("-2", "--reads2", default=None, help="Illumina R2 (FASTQ[.gz])")
    r.add_argument("-s", "--seed-ref", default=None,
                   help="45S seed for baiting/orientation "
                        "(default: bundled Arabidopsis T2T 45S unit)")
    r.add_argument("-o", "--outdir", default="ngs45_out")
    r.add_argument("-t", "--threads", type=int, default=4)
    r.add_argument("--trim", action="store_true", help="quality/adapter trim with cutadapt (S0)")
    r.add_argument("--bait-rounds", type=int, default=3)
    r.add_argument("--subsample", type=int, default=0,
                   help="cap recruited pairs (0 = keep all)")
    r.add_argument("--spades-k", default="auto", help="SPAdes k-mer list (default auto)")
    r.add_argument("--max-cov", type=int, default=2000,
                   help="cap baited depth for assembly to ~this x of an rDNA repeat "
                        "(0 = no cap). The array is 10^4-10^5x deep; capping keeps "
                        "S2 fast without losing the unit (full reads still used for S6)")
    r.add_argument("--call-variants", action="store_true",
                   help="map reads back and report ribotype heterogeneity (S6)")
    r.add_argument("--no-resume", action="store_true", help="ignore previous run state")
    r.add_argument("-v", "--verbose", action="store_true")

    # --- check-deps ---------------------------------------------------------
    c = sub.add_parser("check-deps", help="verify external tools are installed")
    c.add_argument("--optional", action="store_true",
                   help="also require optional tools (cutadapt, bwa, samtools, bcftools)")
    return p


def _setup_logging(verbose: bool) -> None:
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
    )


def main(argv: list | None = None) -> int:
    args = _build_parser().parse_args(argv)

    if args.command == "check-deps":
        try:
            found = check_dependencies(include_optional=args.optional)
        except DependencyError as e:
            print(e, file=sys.stderr)
            return 1
        for tool, path in found.items():
            print(f"  {'OK ' if path else 'MISSING'}  {tool:14s} {path or ''}")
        print("All required tools found.")
        return 0

    if args.command == "run":
        _setup_logging(args.verbose)
        try:
            check_dependencies(include_optional=args.call_variants or args.trim)
        except DependencyError as e:
            print(e, file=sys.stderr)
            return 1
        config = Config(
            reads1=args.reads1,
            reads2=args.reads2,
            seed_ref=args.seed_ref or DEFAULT_SEED,
            outdir=args.outdir,
            threads=args.threads,
            trim=args.trim,
            bait_rounds=args.bait_rounds,
            subsample=args.subsample,
            spades_k=args.spades_k,
            assemble_max_cov=args.max_cov,
            call_variants=args.call_variants,
            resume=not args.no_resume,
        )
        run_pipeline(config)
        return 0

    return 1


if __name__ == "__main__":
    raise SystemExit(main())
