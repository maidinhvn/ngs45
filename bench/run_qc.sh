#!/bin/bash
# QC of every dataset used to develop ngs45. Samples up to N reads per file and
# reports read length, depth, quality (AvgQual, Q20/Q30 %), GC, N50 via seqkit.
# Writes bench/qc_all_datasets.tsv.
set -u
BENCH=/path/to/ngs245/bench
E=/path/to/.conda/envs; export PATH="$E/easy45/bin:/usr/bin:$PATH"
N=200000            # reads sampled per file for the QC table (fast, representative)
OUT=$BENCH/qc_all_datasets.tsv

# accession lookups
declare -A ILLU_OLD ILLU_MOD HIFI OUT_OLD OUT_MOD
while IFS=$'\t' read -r order sp hifi hm illu im; do
  [ "$order" = DONE_SELECT ] && continue
  s=${sp// /_}; ILLU_OLD[$s]=$illu; HIFI[$s]=$hifi
done < "$BENCH/manifest.tsv"
while IFS=$'\t' read -r sp acc rl; do [ "$sp" = DONE ] && continue; ILLU_MOD[$sp]=$acc; done < "$BENCH/manifest_modern.tsv"
# use the long-HiFi accession actually used for easy45
while IFS=$'\t' read -r order sp srr avg mb; do [ "$order" = DONE_RESELECT ] && continue; HIFI[${sp// /_}]=$srr; done < "$BENCH/manifest_hifi_long.tsv"
# outcomes from figure1_points.tsv
while IFS=$'\t' read -r sp inp rl acc oc; do
  [ "$sp" = species ] && continue
  [ "$inp" = old ] && OUT_OLD[$sp]=$oc; [ "$inp" = modern ] && OUT_MOD[$sp]=$oc
done < "$BENCH/figure1_points.tsv"

printf 'role\tspecies\taccession\tn_reads_sampled\tavg_readlen\tQ20pct\tQ30pct\tAvgQual\tGCpct\tN50\tngs45_outcome\n' > "$OUT"

qc_row(){  # role species accession fastq outcome
  local role="$1" sp="$2" acc="$3" f="$4" oc="$5"
  [ -s "$f" ] || { echo "  [skip] $sp $role: no file"; return; }
  local line
  line=$(seqkit head -n "$N" "$f" 2>/dev/null | seqkit stats -aT 2>/dev/null | awk -F'\t' '
    NR==1{for(i=1;i<=NF;i++)h[i]=$i}
    NR==2{for(i=1;i<=NF;i++)v[h[i]]=$i}
    END{printf "%s\t%s\t%s\t%s\t%s\t%s\t%s",v["num_seqs"],v["avg_len"],v["Q20(%)"],v["Q30(%)"],v["AvgQual"],v["GC(%)"],v["N50"]}')
  printf '%s\t%s\t%s\t%s\t%s\n' "$role" "$sp" "$acc" "$line" "$oc" >> "$OUT"
  echo "  [ok] $role $sp"
}

echo "=== QC: Illumina (old) ==="
for f in "$BENCH"/data/*/illu_1.fastq; do sp=$(basename "$(dirname "$f")"); qc_row "Illumina_old" "$sp" "${ILLU_OLD[$sp]:-NA}" "$f" "${OUT_OLD[$sp]:-NA}"; done
echo "=== QC: Illumina (modern) ==="
for f in "$BENCH"/modern_data/*/illu_1.fastq; do sp=$(basename "$(dirname "$f")"); qc_row "Illumina_modern" "$sp" "${ILLU_MOD[$sp]:-NA}" "$f" "${OUT_MOD[$sp]:-NA}"; done
echo "=== QC: HiFi ==="
for f in "$BENCH"/data/*/hifi.fastq; do sp=$(basename "$(dirname "$f")"); qc_row "HiFi" "$sp" "${HIFI[$sp]:-NA}" "$f" "easy45_ok"; done

echo "DONE" ; echo "-> $OUT"
