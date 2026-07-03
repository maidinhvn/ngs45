"""S7 — write the final summary table and human-readable report.

Runs last so it can see everything, including the optional ribotype count from
S6 (which is produced *after* annotation). Emitting the report here — rather than
inside S5 — avoids baking a stale variant count into the summary.

Input keys:  {"final_fasta", "annot", "bait_n_pairs"(opt), "n_ribotype_sites"(opt)}
Output keys: {"summary", "report"}
"""

from __future__ import annotations

import logging

from ..config import Config

log = logging.getLogger("ngs45")

# GFF feature label per region column, for the report block
_LABEL = {"18S": "18S_rRNA", "ITS1": "ITS1", "5.8S": "5_8S_rRNA",
          "ITS2": "ITS2", "26S": "28S_rRNA"}
_TAG = {"18S": "SSU", "ITS1": "ITS1", "5.8S": "5.8S", "ITS2": "ITS2", "26S": "LSU"}


def run(config: Config, state: dict) -> dict:
    annot = state.get("annot") or {}
    rlens = annot.get("region_lens", {})
    coords = annot.get("coords", {})
    n_sites = state.get("n_ribotype_sites")

    # --- summary.tsv --------------------------------------------------------
    cols = {
        "unit_len": annot.get("unit_len", ""),
        "GC_percent": annot.get("GC_percent", ""),
        "bait_pairs": state.get("bait_n_pairs", ""),
        "annot_source": annot.get("annot_source", ""),
        "18S_len": rlens.get("18S", ""), "ITS1_len": rlens.get("ITS1", ""),
        "5.8S_len": rlens.get("5.8S", ""), "ITS2_len": rlens.get("ITS2", ""),
        "26S_len": rlens.get("26S", ""),
        "ITS_barcode_len": annot.get("ITS_barcode_len", ""),
        "ribotype_sites": n_sites if n_sites is not None else "",
        "qc_tandem_dup_bp": sum(state.get("qc_dup_removed") or []),
    }
    summary = config.outdir / "summary.tsv"
    with open(summary, "w") as s:
        s.write("\t".join(cols.keys()) + "\n")
        s.write("\t".join(str(v) for v in cols.values()) + "\n")

    # --- report.txt ---------------------------------------------------------
    reads = f"{config.reads1}" + (f" , {config.reads2}" if config.reads2 else "")
    L = [
        "ngs45 run report", "=" * 40,
        f"reads:      {reads}",
        f"seed:       {config.seed_ref.name}",
        f"unit:       nrDNA_45S transcribed unit, {annot.get('unit_len', '?')} bp, "
        f"GC {annot.get('GC_percent', '?')}%",
        f"annotation: {annot.get('annot_source', '?')}",
        f"bait pairs: {state.get('bait_n_pairs', '?')}",
        "",
        "features (1-based on nrDNA_45S.fasta):",
    ]
    for base in ("18S", "ITS1", "5.8S", "ITS2", "26S"):
        span = coords.get(_TAG[base])
        if span:
            s, e = span
            L.append(f"  {_LABEL[base]:10s} {s:>6}-{e:<6}  ({e - s + 1} bp)")
    if n_sites is not None:
        tail = "  (heterogeneous array — possible hybrid/allopolyploid)" if n_sites else ""
        L += ["", f"ribotype-variant sites: {n_sites}{tail}"]
    dup = state.get("qc_dup_removed") or []
    if dup:
        L += ["", f"QC: collapsed {len(dup)} internal tandem duplication(s) "
              f"{dup} bp — assembly artifact removed from the unit."]
    report = config.outdir / "report.txt"
    report.write_text("\n".join(L) + "\n")

    log.info("S7: wrote %s and %s", summary.name, report.name)
    return {"summary": summary, "report": report}
