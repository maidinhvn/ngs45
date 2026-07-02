#!/bin/bash
# Fast HiFi acquisition via ENA streaming: pull only the first N reads of each
# species' long-read HiFi run with `curl | zcat | head` (head closes the pipe, so
# curl stops early and only ~a few hundred MB is transferred instead of 10-25 GB).
# rDNA is high-copy, so N reads still hold ~1000 unit-spanning reads for easy45.
# Skips species whose already-downloaded HiFi is long enough (avg >= MINAVG).
set -u
BASE=/path/to/ngs245
BENCH=$BASE/bench
DATA=$BENCH/data
LONG=$BENCH/manifest_hifi_long.tsv
LOG=$BENCH/hifi_subset.log
E=/path/to/.conda/envs
export PATH="$E/easy45/bin:/usr/bin:$PATH"
NREADS=50000
MINAVG=6500
MAXPAR=5

now(){ date '+%F %T'; }
log(){ echo "[$(now)] $*" | tee -a "$LOG"; }
slugof(){ echo "$1" | tr ' ' '_'; }
avglen(){ [ -s "$1" ] && seqkit stats -T "$1" 2>/dev/null | tail -1 | cut -f7 || echo 0; }
gate(){ while [ "$(jobs -rp | wc -l)" -ge "$MAXPAR" ]; do sleep 5; done; }

ena_url(){  # accession -> first fastq_ftp url (https), or empty
  curl -s --max-time 60 "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$1&result=read_run&fields=fastq_ftp&format=tsv" 2>/dev/null \
    | awk 'NR>1{print $NF}' | tr ';' '\n' | grep -E 'fastq.gz$' | head -1
}

one(){  # $1 species  $2 srr
  local sp="$1" srr="$2" d; d=$DATA/$(slugof "$sp"); mkdir -p "$d"
  local cur; cur=$(avglen "$d/hifi.fastq")
  if awk "BEGIN{exit !($cur+0 >= $MINAVG)}"; then log "[$sp] existing HiFi avg=${cur}bp OK — skip"; return; fi
  if [ "$srr" = NA ]; then log "[$sp] no long HiFi accession — skip"; return; fi
  local url; url=$(ena_url "$srr")
  if [ -z "$url" ]; then log "[$sp] no ENA fastq url for $srr — skip"; return; fi
  log "[$sp] streaming ${NREADS} reads of $srr from ENA ..."
  curl -s --max-time 3600 "https://$url" 2>>"$d/ena.log" | zcat 2>/dev/null | head -n $((NREADS*4)) > "$d/hifi.fastq"
  if [ -s "$d/hifi.fastq" ]; then
    log "[$sp] done -> avg=$(avglen "$d/hifi.fastq")bp n=$(( $(wc -l < "$d/hifi.fastq")/4 )) size=$(du -h "$d/hifi.fastq"|cut -f1)"
  else
    log "[$sp] ENA stream FAILED for $srr"
  fi
}

log "=== HiFi subset (ENA) START (N=$NREADS, max $MAXPAR parallel) ==="
grep -vE '^DONE_RESELECT|^$' "$LONG" | while IFS=$'\t' read -r order sp srr avg mb; do
  gate; one "$sp" "$srr" &
done
wait
touch "$BENCH/HIFI_SUBSET_DONE"
log "=== HiFi subset DONE ==="
