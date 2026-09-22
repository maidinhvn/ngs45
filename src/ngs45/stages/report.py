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
        # "" when both flanking genes are full length; otherwise names the
        # incomplete end, so a partial unit is reported rather than silent.
        "unit_truncated_end": state.get("unit_truncated_end", ""),
        "reference_guided": bool(state.get("reference_guided")),
        "contaminant_suspect": bool(state.get("contaminant_suspect")),
        "seed_identity_pct": state.get("contaminant_ident", ""),
        "GC_vs_seed": state.get("contaminant_gc", ""),
    }
    summary = config.outdir / "summary.tsv"
    with open(summary, "w") as s:
        s.write("\t".join(cols.keys()) + "\n")
        s.write("\t".join(str(v) for v in cols.values()) + "\n")

    # --- report.txt ---------------------------------------------------------
    reads = f"{config.reads1}" + (f" , {config.reads2}" if config.reads2 else "")
    L = [
        "ngs45 run report", "=" * 40,
    ]
    if state.get("reference_guided"):
        L += [
            "*** REFERENCE-GUIDED RESULT (de novo assembly did not span a unit) ***",
            "    Baited reads were mapped to the seed and a consensus called. This is",
            "    reference-biased: uncovered / divergent positions follow the seed, so",
            "    the spacers (ITS/ETS) may reflect the reference, not this sample.",
            "    Trust the conserved genes; treat spacers/barcode with caution unless",
            "    the seed is a close relative. Not an unbiased de novo reconstruction.",
            "",
        ]
    if state.get("contaminant_suspect"):
        L += [
            "*** WARNING: picked contig may be a CONTAMINANT, not your target ***",
            f"    GC {state.get('contaminant_gc', '?')}% and "
            f"{state.get('contaminant_ident', '?')}% identity to the plant seed "
            "look unlike a plant nuclear 45S (fungal endophyte / plastid rRNA "
            "recruited via the conserved genes). Verify the organism.",
            "",
        ]
    L += [
        f"reads:      {reads}",
        f"seed:       {config.seed_ref.name}",
        f"unit:       nrDNA_45S transcribed unit, {annot.get('unit_len', '?')} bp, "
        f"GC {annot.get('GC_percent', '?')}%",
        *( [f"WARNING:    INCOMPLETE unit - {state['unit_truncated_end']} is truncated; "
            f"treat this sequence as partial"]
           if state.get("unit_truncated_end") else [] ),
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
