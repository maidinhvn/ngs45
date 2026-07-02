#!/bin/bash
# HiFi rescue: for species whose easy45 (HiFi) run FAILED in the concurrent round
# — typically because the chosen HiFi run had reads too short to span a full 45S
# unit — swap in a LONG-READ HiFi run (from manifest_hifi_long.tsv) so the
# subsequent isolated pass can build the HiFi consensus.
set -u
BASE=/path/to/ngs245
BENCH=$BASE/bench
DATA=$BENCH/data
RES=$BENCH/results
LONG=$BENCH/manifest_hifi_long.tsv
LOG=$BENCH/rescue.log
TMP=$BENCH/tmp; mkdir -p "$TMP"

E=/path/to/.conda/envs
export PATH="$E/sra/bin:$E/easy45/bin:/usr/bin:$PATH"
MINAVG=6500     # if current HiFi avg read length is below this, it can't span a unit

now(){ date '+%F %T'; }
log(){ echo "[$(now)] $*" | tee -a "$LOG"; }
slugof(){ echo "$1" | tr ' ' '_'; }
avglen(){ [ -s "$1" ] && seqkit stats -T "$1" 2>/dev/null | tail -1 | cut -f7 || echo 0; }

log "=== HiFi rescue START ==="
[ -s "$LONG" ] || { log "no manifest_hifi_long.tsv; abort"; exit 0; }

grep -vE '^DONE_RESELECT|^$' "$LONG" | while IFS=$'\t' read -r order sp srr avg mb; do
  slug=$(slugof "$sp"); d=$DATA/$slug
  # already have a good HiFi consensus from the concurrent round? skip.
  if [ -s "$RES/$slug/easy45_out/consensus.fasta" ]; then
    log "[$sp] easy45 already succeeded — no rescue"; continue
  fi
  cur=$(avglen "$d/hifi.fastq")
  if awk "BEGIN{exit !($cur+0 >= $MINAVG)}"; then
    log "[$sp] current HiFi avg=${cur}bp is long enough; easy45 failed for another reason — keeping it"
    continue
  fi
  if [ "$srr" = NA ]; then
    log "[$sp] no long-read HiFi run (>=7kb) available on SRA — HiFi comparison not possible"; continue
  fi
  log "[$sp] rescue: current avg=${cur}bp -> downloading long HiFi $srr (avg=${avg}bp, ${mb}MB)"
  rm -rf "$d/$srr"
  timeout 6h prefetch "$srr" -O "$d" --max-size 30g > "$d/pf_hifi_long.log" 2>&1
  timeout 6h fasterq-dump "$d/$srr/$srr.sra" -e 16 -O "$d" -t "$TMP" > "$d/fq_hifi_long.log" 2>&1
  if [ -s "$d/${srr}.fastq" ]; then
    mv -f "$d/${srr}.fastq" "$d/hifi.fastq"
    rm -rf "$d/$srr"
    log "[$sp] rescued HiFi -> avg now $(avglen "$d/hifi.fastq")bp"
  else
    log "[$sp] rescue download FAILED for $srr"
  fi
done
touch "$BENCH/RESCUE_DONE"
log "=== HiFi rescue DONE ==="
