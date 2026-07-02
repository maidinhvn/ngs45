"""S6 — ribotype heterogeneity (optional, ``--call-variants``).

Because the rDNA array holds thousands of copies, positions where recruited
reads disagree beyond sequencing error mark *intragenomic ribotype variants* —
a hallmark of hybridisation / allopolyploidy where divergent parental arrays
coexist. We map the bait reads back to the single-unit consensus and flag every
position whose minor allele exceeds ``var_min_freq`` at sufficient depth.

Input keys:  {"final_fasta", "bait_r1", "bait_r2"}
Output keys: {"variants": <TSV>, "n_ribotype_sites": int}
"""

from __future__ import annotations

import logging
import subprocess

from ..config import Config
from ..external import run as sh

log = logging.getLogger("ngs45")


def _map_back(config: Config, ref, r1, r2, bam):
    # regenerate the .fai every run: the consensus can change between runs and a
    # stale faidx makes bcftools read the reference as Ns (fai_retrieve errors)
    sh(["samtools", "faidx", str(ref)])
    sh(["bwa", "index", str(ref)])
    mem = ["bwa", "mem", "-t", str(config.threads), str(ref), str(r1)]
    if r2 is not None:
        mem.append(str(r2))
    p_mem = subprocess.Popen(mem, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    p_sort = subprocess.Popen(
        ["samtools", "sort", "-@", str(config.threads), "-o", str(bam)],
        stdin=p_mem.stdout, stderr=subprocess.DEVNULL)
    p_mem.stdout.close()
    if p_sort.wait() != 0 or p_mem.wait() != 0:
        raise RuntimeError("S6: bwa mem | samtools sort failed")
    sh(["samtools", "index", str(bam)])


def _scan_minor_alleles(config: Config, ref, bam, vcf):
    """bcftools mpileup with per-allele depth; split multiallelics; write VCF."""
    p_pile = subprocess.Popen(
        ["bcftools", "mpileup", "-f", str(ref), "-a", "FORMAT/AD",
         "-d", "1000000", "-Ou", str(bam)],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    with open(vcf, "w") as out:
        p_norm = subprocess.run(
            ["bcftools", "norm", "-m-", "-f", str(ref), "-Ov"],
            stdin=p_pile.stdout, stdout=out, stderr=subprocess.DEVNULL)
    p_pile.stdout.close()
    p_pile.wait()
    if p_norm.returncode != 0:
        raise RuntimeError("S6: bcftools mpileup|norm failed")


def _parse_sites(vcf, min_depth, min_freq):
    sites = []
    with open(vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            chrom, pos, _id, ref, alt = f[0], f[1], f[2], f[3], f[4]
            if alt in (".", "<*>", ""):
                continue
            fmt, sample = f[8].split(":"), f[9].split(":")
            if "AD" not in fmt:
                continue
            ad = sample[fmt.index("AD")].split(",")
            try:
                depths = [int(x) for x in ad if x not in (".", "")]
            except ValueError:
                continue
            total = sum(depths)
            if total < min_depth or len(depths) < 2:
                continue
            minor = sorted(depths)[-2] / total      # 2nd-largest allele fraction
            if minor >= min_freq:
                sites.append((chrom, pos, ref, alt, total, round(minor, 3)))
    return sites


def run(config: Config, state: dict) -> dict:
    ref = state["final_fasta"]
    r1 = state.get("bait_r1")
    r2 = state.get("bait_r2")
    if r1 is None:
        log.warning("S6: no bait reads in state; skipping variant calling")
        return {"variants": None, "n_ribotype_sites": 0}

    bam = config.workdir / "s6_mapback.bam"
    vcf = config.workdir / "s6_sites.vcf"
    _map_back(config, ref, r1, r2, bam)
    _scan_minor_alleles(config, ref, bam, vcf)
    sites = _parse_sites(vcf, config.var_min_depth, config.var_min_freq)

    out = config.outdir / "ribotype_variants.tsv"
    with open(out, "w") as o:
        o.write("contig\tpos\tref\talt\tdepth\tminor_freq\n")
        for s in sites:
            o.write("\t".join(str(x) for x in s) + "\n")

    log.info("S6: %d ribotype-variant sites (minor>=%.2f, depth>=%d)%s",
             len(sites), config.var_min_freq, config.var_min_depth,
             " — heterogeneous array (possible hybrid/allopolyploid)" if sites else "")
    return {"variants": out, "n_ribotype_sites": len(sites)}
