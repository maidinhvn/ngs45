"""Pipeline orchestration: chains the stages and handles resume.

Each stage is a callable ``fn(config, state) -> outputs`` living in
``ngs45.stages``. The pipeline records completed stages in a manifest file
inside the work directory so a re-run can skip finished stages (``--resume``).
"""

from __future__ import annotations

import json
import logging

from .config import Config
from .stages import (
    annotate,
    assemble,
    bait,
    orient,
    qc,
    report,
    resolve,
    variant,
)

log = logging.getLogger("ngs45")

# Ordered (key, human label, callable). Each stage reads/writes files under
# config.workdir and returns a dict of named output paths merged into `state`.
STAGES = [
    ("qc",       "S0 QC / trim reads",           qc.run),
    ("bait",     "S1 bait rDNA reads",           bait.run),
    ("assemble", "S2 assemble recruited reads",  assemble.run),
    ("resolve",  "S3 resolve 45S monomer",       resolve.run),
    ("orient",   "S4 orient & linearize unit",   orient.run),
    ("annotate", "S5 annotate & ITS barcode",    annotate.run),
    ("variant",  "S6 ribotype variants",         variant.run),
    ("report",   "S7 summary & report",          report.run),
]


def _manifest_path(config: Config):
    return config.workdir / "manifest.json"


def _load_manifest(config: Config) -> dict:
    p = _manifest_path(config)
    if p.exists():
        return json.loads(p.read_text())
    return {}


def _save_manifest(config: Config, manifest: dict) -> None:
    _manifest_path(config).write_text(json.dumps(manifest, indent=2, default=str))


def run_pipeline(config: Config) -> dict:
    """Execute all stages in order, returning the accumulated outputs."""
    config.outdir.mkdir(parents=True, exist_ok=True)
    config.workdir.mkdir(parents=True, exist_ok=True)

    manifest = _load_manifest(config) if config.resume else {}
    state: dict = manifest.get("outputs", {})

    for key, label, fn in STAGES:
        if key == "qc" and not config.trim:
            log.info("[skip] %s (no --trim)", label)
            continue
        if key == "variant" and not config.call_variants:
            log.info("[skip] %s (no --call-variants)", label)
            continue
        if config.resume and key in manifest.get("done", []):
            log.info("[resume] %s already complete", label)
            continue
        log.info("[run] %s", label)
        outputs = fn(config, state)
        state.update(outputs or {})
        manifest.setdefault("done", []).append(key)
        manifest["outputs"] = state
        _save_manifest(config, manifest)

    log.info("Pipeline complete. Results in %s", config.outdir)
    return state
