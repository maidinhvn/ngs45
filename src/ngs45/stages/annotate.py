"""S5 — annotate the transcribed unit and write the deliverables.

Mirrors easy45's annotation contract so the two tools' outputs are
interchangeable: ITSx (``-t Tracheophyta``) delimits SSU/ITS1/5.8S/ITS2/LSU as
one contiguous tiling, parsed into a GFF3 and the standard ITS barcode. On top
of easy45 we also drop each of the five regions to its own FASTA (many users
just want the five 45S regions to compare) and a one-row summary table.

Final files in config.outdir:
  nrDNA_45S.fasta     the mature transcribed unit (18S 5' -> 26S 3')
  annotation.gff3     SSU/ITS1/5.8S/ITS2/LSU features (ITSx positions)
  its.fasta           ITS1-5.8S-ITS2 barcode (BLAST/GenBank)
  its_parts.fasta     ITS1, 5.8S, ITS2 as separate records
  regions/            18S / ITS1 / 5.8S / ITS2 / 26S FASTAs
  summary.tsv         lengths / GC / coordinates
  report.txt          human-readable run report

Input keys:  {"monomer", "full_repeat"(opt), "bait_n_pairs"(opt), "n_ribotype_sites"(opt)}
Output keys: {"final_fasta", "gff", "its", "its_parts", "summary", "report"}
"""

from __future__ import annotations

import logging
import re
import subprocess

from ..config import Config
from ..external import run as sh
from ..io import read_fasta, write_fasta

log = logging.getLogger("ngs45")

# ITSx positions-file tag -> (GFF feature name, region-fasta basename).
# Feature names match easy45 exactly (LSU labelled 28S_rRNA) so GFFs are mergeable.
_REGION = {
    "SSU":  ("18S_rRNA", "18S"),
    "ITS1": ("ITS1", "ITS1"),
    "5.8S": ("5_8S_rRNA", "5.8S"),
    "ITS2": ("ITS2", "ITS2"),
    "LSU":  ("28S_rRNA", "26S"),
}
_POS_RE = re.compile(r"(SSU|ITS1|5\.8S|ITS2|LSU): (\d+)-(\d+)")


def _itsx(config: Config, fasta):
    """Run ITSx; return (positions[list of (tag,start,end)], {region_base: seq})."""
    prefix = config.workdir / "s5_itsx"
    sh(["ITSx", "-i", str(fasta), "-o", str(prefix),
        "-t", "Tracheophyta", "--cpu", str(config.threads),
        "--preserve", "T", "--save_regions", "all",
        "--graphical", "F", "--complement", "F"])

    pos_file = prefix.with_name(prefix.name + ".positions.txt")
    positions = []
    if pos_file.exists():
        for line in pos_file.read_text().splitlines():
            for tag, s, e in _POS_RE.findall(line):
                if int(e) >= int(s):
                    positions.append((tag, int(s), int(e)))

    regions = {}
    for tag, (_feat, base) in _REGION.items():
        suffix = {"SSU": "SSU", "ITS1": "ITS1", "5.8S": "5_8S",
                  "ITS2": "ITS2", "LSU": "LSU"}[tag]
        f = prefix.with_name(prefix.name + f".{suffix}.fasta")
        if f.exists() and f.stat().st_size:
            rec = read_fasta(f)
            if rec:
                regions[base] = rec[0][1]
    return positions, regions


def _barrnap_positions(config: Config, fasta):
    """Fallback annotation if ITSx yields nothing: barrnap rRNA genes only."""
    cp = sh(["barrnap", "--kingdom", "euk", "--threads", str(config.threads),
             "--quiet", str(fasta)])
    name_map = {"18S": "SSU", "5_8S": "5.8S", "5.8S": "5.8S", "28S": "LSU", "26S": "LSU"}
    positions = []
    for line in cp.stdout.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split("\t")
        if len(f) < 9:
            continue
        attr = f[8]
        tag = next((v for k, v in name_map.items() if k in attr), None)
        if tag:
            positions.append((tag, int(f[3]), int(f[4])))
    return positions


def _gc(seq: str) -> float:
    seq = seq.upper()
    gc = seq.count("G") + seq.count("C")
    atgc = gc + seq.count("A") + seq.count("T")
    return round(100 * gc / atgc, 2) if atgc else 0.0


def run(config: Config, state: dict) -> dict:
    _n, unit = read_fasta(state["monomer"])[0]

    final_fasta = config.outdir / "nrDNA_45S.fasta"
    write_fasta([("nrDNA_45S_unit", unit)], final_fasta)

    positions, regions = _itsx(config, final_fasta)
    source = "ngs45"
    if not positions:
        positions = _barrnap_positions(config, final_fasta)
        source = "barrnap"
        log.warning("S5: ITSx delimited nothing; annotating rRNA genes with barrnap")

    # --- GFF3 (ITSx tiling; easy45-compatible feature names) ----------------
    gff = config.outdir / "annotation.gff3"
    coords = {}
    with open(gff, "w") as g:
        g.write("##gff-version 3\n")
        for tag, s, e in positions:
            feat = _REGION.get(tag, (tag, tag))[0]
            coords[tag] = (s, e)
            ftype = "rRNA" if feat.endswith("_rRNA") else "misc_feature"
            g.write(f"nrDNA_45S_unit\t{source}\t{ftype}\t{s}\t{e}\t.\t+\t.\tName={feat}\n")

    # --- per-region FASTAs (all five) --------------------------------------
    regdir = config.outdir / "regions"
    regdir.mkdir(exist_ok=True)
    for base, seq in regions.items():
        write_fasta([(base, seq)], regdir / f"{base}.fasta")

    # --- ITS barcode outputs (easy45 names) --------------------------------
    its = config.outdir / "its.fasta"
    its_parts = config.outdir / "its_parts.fasta"
    i1, c58, i2 = coords.get("ITS1"), coords.get("5.8S"), coords.get("ITS2")
    if i1 and i2 and i1[0] <= i2[1]:
        write_fasta([("nrDNA_45S_unit_ITS ITS1-5.8S-ITS2", unit[i1[0] - 1:i2[1]])], its)
        parts = []
        for label, span in (("ITS1", i1), ("5.8S", c58), ("ITS2", i2)):
            if span:
                parts.append((f"nrDNA_45S_unit_{label}", unit[span[0] - 1:span[1]]))
        write_fasta(parts, its_parts)
    else:
        log.warning("S5: incomplete ITS (ITS1/ITS2); ITS barcode not written")
        its = its_parts = None

    # --- stash annotation facts for the final report stage (S7) -------------
    reg_len = lambda base: len(regions.get(base, "")) or _span_len(coords, base)
    annot = {
        "unit_len": len(unit), "GC_percent": _gc(unit), "annot_source": source,
        "coords": {tag: [s, e] for tag, (s, e) in coords.items()},
        "region_lens": {b: reg_len(b) for b in ("18S", "ITS1", "5.8S", "ITS2", "26S")},
        "ITS_barcode_len": (i2[1] - i1[0] + 1) if (i1 and i2) else "",
    }

    log.info("S5: wrote %s (unit %d bp, %d features via %s, %d region FASTAs)",
             final_fasta.name, len(unit), len(positions), source, len(regions))
    return {"final_fasta": final_fasta, "gff": gff, "its": its,
            "its_parts": its_parts, "annot": annot}


def _span_len(coords, base):
    tag = {"18S": "SSU", "ITS1": "ITS1", "5.8S": "5.8S", "ITS2": "ITS2", "26S": "LSU"}[base]
    span = coords.get(tag)
    return (span[1] - span[0] + 1) if span else ""
