"""Thin wrappers around the external bioinformatics tools ngs45 depends on.

Design contract (driven by the conda-packaging goal):
    The Python package contains *no* bundled binaries. Every heavy tool is a
    conda run-dependency declared in environment.yml. This module is the only
    place that shells out, and it fails *fast and clearly* at startup if a
    required tool is missing from PATH — never mid-pipeline.
"""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass


@dataclass(frozen=True)
class Tool:
    name: str            # executable name on PATH
    version_args: tuple  # args that print version info (for the startup check)
    optional: bool = False


# Tools used by the pipeline. `optional` ones are only needed for certain modes.
REQUIRED_TOOLS = (
    Tool("cutadapt", ("--version",), optional=True),   # S0, only with --trim
    Tool("bowtie2", ("--version",)),                   # S1 baiting
    Tool("bowtie2-build", ("--version",)),
    Tool("spades.py", ("--version",)),                 # S2 assembly
    Tool("seqkit", ("version",)),
    Tool("blastn", ("-version",)),                     # S3/S4 locate unit & genes
    Tool("makeblastdb", ("-version",)),
    Tool("cmsearch", ("-h",)),                         # S4 mature-boundary trim (infernal)
    Tool("barrnap", ("--version",)),                   # S4 orientation (strand of 18S)
    Tool("ITSx", ("-h",)),                             # S5 ITS delimitation
    Tool("bwa", ("",), optional=True),                 # S6 map-back (variant mode)
    Tool("samtools", ("--version",), optional=True),   # S6
    Tool("bcftools", ("--version",), optional=True),   # S6
)


class DependencyError(RuntimeError):
    """Raised when a required external tool is missing from PATH."""


def check_dependencies(include_optional: bool = False) -> dict[str, str | None]:
    """Verify required tools are on PATH. Returns {tool: resolved_path_or_None}.

    Raises DependencyError listing every missing tool at once. Optional tools
    are only enforced when ``include_optional`` is set.
    """
    found: dict[str, str | None] = {}
    missing: list[str] = []
    for tool in REQUIRED_TOOLS:
        path = shutil.which(tool.name)
        found[tool.name] = path
        enforce = (not tool.optional) or include_optional
        if path is None and enforce:
            missing.append(tool.name)
    if missing:
        raise DependencyError(
            "Missing required tool(s): "
            + ", ".join(missing)
            + "\nInstall everything with: conda env create -f environment.yml"
        )
    return found


def run(cmd: list, **kwargs) -> subprocess.CompletedProcess:
    """Run an external command, raising CalledProcessError on non-zero exit.

    Output is captured by default so callers can inspect / log it.
    """
    kwargs.setdefault("check", True)
    kwargs.setdefault("text", True)
    kwargs.setdefault("capture_output", True)
    return subprocess.run([str(c) for c in cmd], **kwargs)
