"""Run configuration and default parameters for the ngs45 pipeline.

All user-tunable thresholds live here so they can be surfaced on the CLI and
overridden in one place. The :class:`Config` object is built once in ``cli.py``
and threaded through every stage.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

# Recruitment/orientation seed shipped with the package: one intact 45S repeat
# unit (ETS-18S-ITS1-5.8S-ITS2-26S-IGS) from the Arabidopsis thaliana T2T NOR2
# assembly (GenBank OR453402). The conserved 18S/5.8S/26S genes drive baiting
# across angiosperms; the spacers are recovered by iterative extension.
DEFAULT_SEED = Path(__file__).parent / "data" / "seed_arabidopsis_45S.fasta"

# Rfam covariance models (SSU RF01960 + LSU RF02543). cmsearch against these
# defines the mature 18S 5' / 26S 3' termini structure-first, independent of any
# single reference. This is the SAME bundled CM easy45 uses, so both tools trim
# the transcribed unit to identical mature boundaries (directly comparable).
DEFAULT_CM = Path(__file__).parent / "data" / "rrna_ssu_lsu.cm"


@dataclass
class Config:
    # --- inputs -----------------------------------------------------------
    reads1: Path                       # R1 (paired-end Illumina)
    reads2: Path | None = None         # R2; None => single-end / interleaved handled upstream
    seed_ref: Path = DEFAULT_SEED      # 45S seed for baiting + orientation
    cm_ref: Path = DEFAULT_CM          # Rfam SSU+LSU CMs for mature-boundary trim (S4)

    # --- outputs ----------------------------------------------------------
    outdir: Path = Path("ngs45_out")
    workdir: Path | None = None        # intermediate files; defaults to <outdir>/work
    threads: int = 4

    # --- Stage 0: QC ------------------------------------------------------
    trim: bool = False                 # run cutadapt quality/adapter trim
    trim_quality: int = 20
    min_read_len: int = 50

    # --- Stage 1: iterative baiting ---------------------------------------
    bait_rounds: int = 3               # max extension rounds
    bait_min_len: int = 40             # bowtie2 --local: report only if this many bp align
    bait_converge: float = 0.02        # stop when new-read fraction < this
    subsample: int = 0                 # cap recruited pairs (0 = no cap); rDNA is deep

    # --- Stage 2: assembly ------------------------------------------------
    spades_k: str = "auto"             # "auto" => 21,33,55,77,99,127 (filtered to < read length in S2)
    spades_careful: bool = False       # off by default: --careful's mismatch corrector
                                       # fragments variant-rich rDNA (see stages/assemble.py)
    min_cov: float = 0.0               # drop assembly-graph nodes below this k-mer coverage

    # --- Stage 3: monomer resolution --------------------------------------
    unit_min_len: int = 4000           # a 45S transcribed unit is typically 5-8 kb
    unit_max_len: int = 20000          # ceiling incl. long IGS
    min_gene_ident: float = 70.0       # blast %id for a seed-gene hit on a contig

    # --- Stage 4: mature-boundary trim + assembly QC ----------------------
    dup_min_len: int = 40              # collapse an internal tandem duplication in the
    dup_min_ident: float = 97.0        #   unit that is >= this long at >= this % id and
                                       #   sits immediately next to its source. rRNA
                                       #   genes are single-copy, so such a perfect
                                       #   adjacent repeat is a mis-assembly across a
                                       #   unit junction, not biology (see docs).

    # --- Stage 6: ribotype variants (optional) ----------------------------
    call_variants: bool = False
    var_min_freq: float = 0.10         # minor-allele fraction to report a ribotype site
    var_min_depth: int = 20

    # --- behaviour --------------------------------------------------------
    resume: bool = True                # skip stages whose outputs already exist
    keep_intermediate: bool = True

    def __post_init__(self) -> None:
        self.reads1 = Path(self.reads1)
        if self.reads2 is not None:
            self.reads2 = Path(self.reads2)
        self.seed_ref = Path(self.seed_ref)
        self.cm_ref = Path(self.cm_ref)
        self.outdir = Path(self.outdir)
        if self.workdir is None:
            self.workdir = self.outdir / "work"
        self.workdir = Path(self.workdir)

    @property
    def spades_klist(self) -> str:
        return "21,33,55,77,99,127" if self.spades_k == "auto" else self.spades_k
