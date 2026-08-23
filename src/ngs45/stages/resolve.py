"""S3 — resolve a single 45S repeat unit (monomer) from the assembly.

rDNA is a tandem array, so SPAdes typically collapses it into one (or few)
high-coverage contigs that are either ~one unit long or several units of the
same sequence concatenated. We recover exactly one period:

  1. BLAST the 45S seed against the contigs to find the rDNA-bearing contig(s).
  2. Self-BLAST that contig: a tandem array self-aligns off the main diagonal at
     an offset equal to the repeat period. The smallest consistent offset in the
     plausible unit-length window is the monomer length.
  3. Cut one period. If the contig is already ~one unit (no internal repeat), keep
     it whole.

Input keys:  {"contigs"}
Output keys: {"monomer_raw", "unit_len", "rdna_contig"}
"""

from __future__ import annotations

import logging
import subprocess
from collections import Counter

from ..config import Config
from ..external import run as sh
from ..io import read_fasta, write_fasta

log = logging.getLogger("ngs45")


def _reference_guided_monomer(config: Config, state: dict) -> str:
    """Opt-in fallback (``--reference-guided``) when de novo assembly can't span a
    unit: map the baited reads to the (ideally congeneric) ``--seed-ref`` and call
    a consensus with bcftools.

    REFERENCE-BIASED: positions with no read support stay at the reference and any
    divergent intragenomic array collapses to the major allele — a *usable* 45S
    for NGS-only libraries where de novo fails, not an unbiased reconstruction.
    Reuses S6's map-back (bwa mem | samtools sort). Needs bwa/samtools/bcftools.
    """
    from .variant import _map_back           # reuse the S6 map-back
    ref = config.seed_ref
    r1, r2 = state.get("bait_r1"), state.get("bait_r2")
    if r1 is None:
        raise RuntimeError("reference-guided fallback needs S1's baited reads.")
    bam = config.workdir / "s3_refguided.bam"
    _map_back(config, ref, r1, r2, bam)
    vcf = config.workdir / "s3_refguided.vcf.gz"
    p_pile = subprocess.Popen(
        ["bcftools", "mpileup", "-f", str(ref), "-Ou", str(bam)],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    with open(vcf, "wb") as out:
        p_call = subprocess.run(["bcftools", "call", "-mv", "-Oz"],
                                stdin=p_pile.stdout, stdout=out,
                                stderr=subprocess.DEVNULL)
    p_pile.stdout.close()
    p_pile.wait()
    if p_call.returncode != 0:
        raise RuntimeError("S3 reference-guided: bcftools mpileup|call failed "
                           "(is bcftools installed?).")
    sh(["bcftools", "index", "-f", str(vcf)])
    cons = config.workdir / "s3_refguided.fasta"
    with open(cons, "w") as out:
        subprocess.run(["bcftools", "consensus", "-f", str(ref), str(vcf)],
                       stdout=out, stderr=subprocess.DEVNULL, check=True)
    _name, seq = read_fasta(cons)[0]
    return seq

_FMT = "6 qseqid sseqid pident length qstart qend sstart send bitscore"


def _blast(query, subject, workdir, tag, perc_identity=None, extra=None) -> list[dict]:
    """Run blastn query-vs-subject, return parsed tabular rows."""
    out = workdir / f"s3_{tag}.tsv"
    cmd = ["blastn", "-query", str(query), "-subject", str(subject),
           "-outfmt", _FMT, "-dust", "no"]
    if perc_identity is not None:
        cmd += ["-perc_identity", str(perc_identity)]
    if extra:
        cmd += extra
    cp = sh(cmd)
    out.write_text(cp.stdout)
    rows = []
    for line in cp.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 9:
            continue
        rows.append({
            "q": f[0], "s": f[1], "pid": float(f[2]), "len": int(f[3]),
            "qs": int(f[4]), "qe": int(f[5]), "ss": int(f[6]), "se": int(f[7]),
            "bits": float(f[8]),
        })
    return rows


def _pick_rdna_contig(config: Config, sources) -> tuple[str, str]:
    """Pick the best rDNA sequence across several FASTAs (contigs + scaffolds).

    Among sequences carrying a substantial seed hit (>=800 bp aligned, i.e. they
    hold most of the conserved genes) we take the *longest* — that is the one most
    likely to span a whole unit (a scaffold often joins gene contigs across the
    spacer gaps that the raw contigs leave broken). If nothing clears that bar we
    fall back to the single most seed-aligned contig.
    """
    best = None           # (seqlen, aligned_bp, id, seq, label)
    fallback = None       # (aligned_bp, id, seq, label)
    n_hit = 0
    for label, fasta in sources:
        if fasta is None:
            continue
        rows = _blast(config.seed_ref, fasta, config.workdir, f"seed_vs_{label}",
                      perc_identity=config.min_gene_ident)
        if not rows:
            continue
        aligned = Counter()
        for r in rows:
            aligned[r["s"]] += r["len"]
        n_hit += len(aligned)
        seqs = dict(_seq_by_id(fasta))
        for cid, bp in aligned.items():
            seq = seqs.get(cid, "")
            if fallback is None or bp > fallback[0]:
                fallback = (bp, cid, seq, label)
            if bp >= 800 and (best is None or len(seq) > best[0]):
                best = (len(seq), bp, cid, seq, label)

    if best is not None:
        slen, bp, cid, seq, label = best
        log.info("S3: rDNA sequence = %s [%s] %d bp (%d bp seed-aligned; %d hits total)",
                 cid, label, slen, bp, n_hit)
        return cid, seq
    if fallback is not None:
        bp, cid, seq, label = fallback
        log.info("S3: rDNA sequence = %s [%s] %d bp (weak: %d bp seed-aligned)",
                 cid, label, len(seq), bp)
        return cid, seq
    raise RuntimeError("S3: no contig/scaffold aligned to the 45S seed — assembly "
                       "may have failed or reads were off-target.")


def _seq_by_id(fasta):
    for name, seq in read_fasta(fasta):
        yield name.split()[0], seq


def _tandem_period(config: Config, contig_id, seq) -> int | None:
    """Smallest self-alignment offset in the unit-length window, or None."""
    tmp = config.workdir / "s3_contig.fasta"
    write_fasta([(contig_id, seq)], tmp)
    rows = _blast(tmp, tmp, config.workdir, "self", perc_identity=95)
    offsets = []
    for r in rows:
        off = r["ss"] - r["qs"]
        if off <= 0:                                   # keep upper triangle only
            continue
        if r["len"] < config.unit_min_len // 2:        # ignore short incidental hits
            continue
        if config.unit_min_len <= off <= config.unit_max_len:
            offsets.append(off)
    if not offsets:
        return None
    # The fundamental period is the smallest well-supported offset; round to
    # tolerate a few indels between copies before taking the mode.
    binned = Counter(round(o / 10) * 10 for o in offsets)
    period = min(o for o, _ in binned.most_common() if binned[o] == max(binned.values()))
    return period


def _looks_foreign(config: Config, cid, seq) -> tuple[bool, float, float]:
    """Flag a picked rDNA contig that looks unlike a plant *nuclear* 45S — i.e. a
    recruited contaminant pulled in via the conserved 18S/26S (a fungal endophyte,
    or bacterial-type plastid rRNA). Plant nuclear rDNA is GC-rich (~55-62 %) and
    stays >~92 % identical to the plant seed's conserved genes; fungal rDNA is
    AT-richer (~45-50 %) and only ~80-85 % identical (this is exactly how easy45's
    'foreign' filter caught the Vitis endophyte). Returns (suspect, gc, identity).
    """
    gc = 100.0 * (seq.count("G") + seq.count("C")) / max(1, len(seq))
    tmp = config.workdir / "s3_pick_vs_seed.fasta"
    write_fasta([(cid, seq)], tmp)
    rows = _blast(config.seed_ref, tmp, config.workdir, "pick_vs_seed",
                  perc_identity=config.min_gene_ident)
    tot = sum(r["len"] for r in rows)
    ident = sum(r["pid"] * r["len"] for r in rows) / tot if tot else 0.0
    suspect = gc < 50.0 or gc > 66.0 or ident < 90.0
    return suspect, round(gc, 1), round(ident, 1)


def run(config: Config, state: dict) -> dict:
    sources = [("contigs", state["contigs"]),
               ("scaffolds", state.get("scaffolds"))]
    cid, seq = _pick_rdna_contig(config, sources)

    # low-depth warning: SPAdes node names carry k-mer coverage (..._cov_<X>)
    import re as _re
    _m = _re.search(r"cov_([\d.]+)", cid)
    if _m and float(_m.group(1)) < config.min_cov_warn:
        log.warning(
            "S3: picked rDNA contig %s has low coverage (%.1fx < %.0fx) — few rDNA "
            "reads (small/noisy library). The unit may be incomplete or unreliable; "
            "consider more Illumina, or --seed-ref a closer 45S.",
            cid, float(_m.group(1)), config.min_cov_warn)

    suspect, gc, ident = _looks_foreign(config, cid, seq)
    if suspect:
        log.warning(
            "S3: picked contig %s looks unlike a plant nuclear 45S (GC %.1f%%, "
            "%.1f%% identity to the seed) — it may be a recruited CONTAMINANT "
            "(fungal endophyte / plastid rRNA), not your target. Verify the "
            "organism; if wrong, deplete organelles (S0) or use a cleaner library.",
            cid, gc, ident)

    if len(seq) < config.unit_min_len:
        if config.reference_guided:
            log.warning(
                "S3: de novo assembly spanned only %d bp; --reference-guided is "
                "set -> mapping baited reads to %s for a REFERENCE-GUIDED "
                "consensus. This is reference-biased and collapses intragenomic "
                "heterogeneity; use only when de novo fails and HiFi is not an "
                "option, and pass a phylogenetically close --seed-ref.",
                len(seq), config.seed_ref.name)
            seq = _reference_guided_monomer(config, state)
            cid = "reference_guided"
            out = config.workdir / "s3_monomer_raw.fasta"
            write_fasta([(f"{cid} len={len(seq)}", seq)], out)
            return {"monomer_raw": out, "unit_len": len(seq),
                    "rdna_contig": cid, "reference_guided": True}
        raise RuntimeError(
            f"S3: best rDNA sequence {cid} is only {len(seq)} bp "
            f"(< {config.unit_min_len}); de novo assembly did not span a unit. "
            "This happens on small or hybrid-divergent libraries whose short reads "
            "don't bridge the species-specific spacers, so the bundled (distant) "
            "seed recruits mainly the conserved genes.\n"
            "  Fix without new sequencing: pass a phylogenetically closer 45S via "
            "--seed-ref — a congener's GenBank 45S, or an easy45 HiFi unit if you "
            "have one. A closer seed recruits the species' own ITS/ETS/IGS reads in "
            "round 0 and usually rescues the assembly (e.g. a Vitis test went from "
            "2904 bp -> 6739 bp with a matched seed).\n"
            "  Last resort (NGS-only, no HiFi): re-run with --reference-guided plus "
            "a close --seed-ref to map the reads onto that reference and emit a "
            "reference-guided consensus (biased; heterogeneity collapsed).\n"
            "  See docs/ASSEMBLY_LIMITATION.md.")

    period = _tandem_period(config, cid, seq)
    if period:
        monomer = seq[:period]
        log.info("S3: tandem period = %d bp -> cut one monomer", period)
    elif len(seq) <= config.unit_max_len:
        monomer = seq
        log.info("S3: contig is ~one unit (%d bp), no internal repeat -> keep whole", len(seq))
    else:
        monomer = seq[:config.unit_max_len]
        log.warning("S3: no clean tandem period on a %d bp contig; truncated to "
                    "%d bp (inspect manually).", len(seq), config.unit_max_len)

    out = config.workdir / "s3_monomer_raw.fasta"
    write_fasta([(f"{cid}_monomer len={len(monomer)}", monomer)], out)
    return {"monomer_raw": out, "unit_len": len(monomer), "rdna_contig": cid,
            "contaminant_suspect": suspect, "contaminant_gc": gc,
            "contaminant_ident": ident}
