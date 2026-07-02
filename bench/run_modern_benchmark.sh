#!/bin/bash
# Re-benchmark ngs45 (with the self-contained assembly fixes) on MODERN Illumina
# runs (>=150 bp PE), streamed as subsets from ENA. Reuses the already-computed
# easy45/HiFi consensus per species as the gold standard. Writes benchmark_modern.tsv.
set -u
BASE=/path/to/ngs245
BENCH=$BASE/bench
DATA=$BENCH/modern_data
RES=$BENCH/modern_results
MAN=$BENCH/manifest_modern.tsv
OUT=$BENCH/benchmark_modern.tsv
LOG=$BENCH/modern_bench.log
mkdir -p "$DATA" "$RES"

E=/path/to/.conda/envs
export PATH="$E/easy45/bin:$E/getorganelle_env/bin:$E/bcftools_env/bin:/usr/bin:$PATH"
export PYTHONPATH="$BASE/src"
PY=$E/easy45/bin/python
NREADS=1000000     # read pairs to stream per species
WORKERS=3

now(){ date '+%F %T'; }
log(){ echo "[$(now)] $*" | tee -a "$LOG"; }
best_hit(){ local q="$1" s="$2" r; [ -s "$q" ] && [ -s "$s" ] || { echo "NA NA NA"; return; }
  r=$(blastn -query "$q" -subject "$s" -outfmt "6 pident length mismatch bitscore" 2>/dev/null | sort -k4 -nr | head -1)
  [ -n "$r" ] && echo "$r" | awk '{print $1,$2,$3}' || echo "NA NA NA"; }
seqlen(){ [ -s "$1" ] && seqkit stats -T "$1" 2>/dev/null | tail -1 | cut -f5 || echo NA; }
ena_urls(){ curl -s --max-time 60 "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$1&result=read_run&fields=fastq_ftp&format=tsv" 2>/dev/null | awk 'NR>1{print $NF}'; }

process(){  # $1 species  $2 acc
  local sp="$1" acc="$2" d="$DATA/$1" r="$RES/$1"
  mkdir -p "$d" "$r"
  if [ ! -s "$d/illu_1.fastq" ]; then
    local urls; urls=$(ena_urls "$acc")
    local u1 u2; u1=$(echo "$urls"|tr ';' '\n'|grep '_1.fastq'|head -1); u2=$(echo "$urls"|tr ';' '\n'|grep '_2.fastq'|head -1)
    [ -z "$u1" ] && { log "[$sp] no ENA _1 url for $acc"; return; }
    log "[$sp] streaming $NREADS pairs of $acc"
    curl -s --max-time 1800 "https://$u1" 2>/dev/null | zcat 2>/dev/null | head -n $((NREADS*4)) > "$d/illu_1.fastq" &
    curl -s --max-time 1800 "https://$u2" 2>/dev/null | zcat 2>/dev/null | head -n $((NREADS*4)) > "$d/illu_2.fastq" &
    wait
  fi
  local rl; rl=$(seqkit stats -T "$d/illu_1.fastq" 2>/dev/null|tail -1|cut -f7)
  log "[$sp] ngs45 (modern reads, readlen~$rl)"
  rm -rf "$r/ngs45_out"
  $PY -m ngs45.cli run -1 "$d/illu_1.fastq" -2 "$d/illu_2.fastq" -o "$r/ngs45_out" \
      -t 12 --bait-rounds 3 --call-variants -v > "$r/ngs45.log" 2>&1 || true
  local u=$r/ngs45_out/nrDNA_45S.fasta  its=$r/ngs45_out/its.fasta
  local ez=$BENCH/results_isolated/$sp/easy45_out/consensus.fasta
  local ezits=$BENCH/results_isolated/$sp/easy45_out/its.fasta
  local ref=$BENCH/results/$sp/its_ref.fasta
  local A B C ribo
  A=$(best_hit "$u" "$ez"); B=$(best_hit "$its" "$ez"); C=$(best_hit "$its" "$ref")
  ribo=$( [ -s "$r/ngs45_out/summary.tsv" ] && tail -1 "$r/ngs45_out/summary.tsv"|awk -F'\t' '{print $NF}' || echo NA )
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$sp" "$acc" "$rl" "$(seqlen "$u")" "$(seqlen "$ez")" "$ribo" \
    "$(echo $A|cut -d' ' -f1)/$(echo $A|cut -d' ' -f3)mm" "$(echo $B|cut -d' ' -f1)" "$(echo $C|cut -d' ' -f1)" >> "$OUT"
  log "[$sp] done: ngs45=$(seqlen "$u")bp vs easy45 pid=$(echo $A|cut -d' ' -f1)"
}

log "=== MODERN benchmark START ==="
printf 'species\tacc\treadlen\tngs45_unit\teasy45_unit\tribo_sites\tunitPID/mm\tITS_ngsVSez\tITS_ngsVSgb\n' > "$OUT"
# read the manifest into an array first: a `... | while` pipe puts the loop in a
# subshell whose background jobs the final `wait` cannot see (the earlier version
# marked DONE before stragglers finished). mapfile + main-shell loop fixes it.
mapfile -t ROWS < <(grep -vE '^DONE|^$' "$MAN")
for row in "${ROWS[@]}"; do
  IFS=$'\t' read -r sp acc rl <<< "$row"
  [ "$acc" = NA ] && { log "[$sp] no modern run — skipped"; continue; }
  while [ "$(jobs -rp | wc -l)" -ge "$WORKERS" ]; do sleep 5; done
  process "$sp" "$acc" &
done
wait   # in the main shell -> reliably waits for every process() before DONE
{ echo "ngs45 MODERN-input benchmark — $(now)"; echo; column -t -s$'\t' "$OUT"; } > "$BENCH/BENCHMARK_MODERN.txt"
touch "$BENCH/MODERN_DONE"
log "=== MODERN benchmark DONE -> $OUT ==="
