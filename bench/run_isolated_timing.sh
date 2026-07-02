#!/bin/bash
# ============================================================================
# Isolated timing pass — clean, reproducible per-species runtimes for the paper.
# Runs ONE assembly at a time (no concurrency, no download competing), fixed
# 16 threads, identical parameters to the concurrent run. Records wall-clock,
# CPU-time (user+sys) and peak RSS via /usr/bin/time -v.
# Prereq: bench/data/<sp>/{illu_1,illu_2,hifi}.fastq already downloaded.
# ============================================================================
set -u
BASE=/path/to/ngs245
BENCH=$BASE/bench
DATA=$BENCH/data
RES=$BENCH/results_isolated
LOG=$BENCH/isolated.log
ISO=$BENCH/timing_isolated.tsv
MAN=$BENCH/manifest.tsv
mkdir -p "$RES"

E=/path/to/.conda/envs
export PATH="$E/easy45/bin:$E/getorganelle_env/bin:$E/bcftools_env/bin:$E/edirect_env/bin:/usr/bin:$PATH"
export PYTHONPATH="$BASE/src"
PY="$E/easy45/bin/python"
EASY45="$E/easy45/bin/easy45"
T=16
TIME=/usr/bin/time

now(){ date '+%Y-%m-%d %H:%M:%S'; }
log(){ echo "[$(now)] $*" | tee -a "$LOG"; }
slugof(){ echo "$1" | tr ' ' '_'; }
seqlen(){ [ -s "$1" ] && seqkit stats -T "$1" 2>/dev/null | tail -1 | cut -f5 || echo NA; }
best_hit(){ local q="$1" s="$2" r; [ -s "$q" ] && [ -s "$s" ] || { echo "NA NA NA"; return; }
  r=$(blastn -query "$q" -subject "$s" -outfmt "6 pident length mismatch bitscore" 2>/dev/null | sort -k4 -nr | head -1)
  [ -n "$r" ] && echo "$r" | awk '{print $1,$2,$3}' || echo "NA NA NA"; }
FINAL=$BENCH/benchmark_final.tsv

# run a command under /usr/bin/time -v; append one timing row.
timed(){  # $1 species  $2 data_type  $3 accession  $4 timefile ; rest = command
  local sp="$1" dt="$2" acc="$3" tf="$4"; shift 4
  local t0 wall usr sys rss cpu
  t0=$(date +%s)
  "$TIME" -v -o "$tf" "$@" >> "${tf%.time}.log" 2>&1
  wall=$(( $(date +%s) - t0 ))
  usr=$(awk -F': ' '/User time/{print $2}' "$tf" 2>/dev/null)
  sys=$(awk -F': ' '/System time/{print $2}' "$tf" 2>/dev/null)
  rss=$(awk -F': ' '/Maximum resident set size/{print $2}' "$tf" 2>/dev/null)
  cpu=$(awk "BEGIN{printf \"%.0f\", ${usr:-0}+${sys:-0}}")
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$sp" "$dt" "$acc" "$wall" "$cpu" "${rss:-NA}" | tee -a "$ISO"
  log "  [$sp] $dt: wall=${wall}s cpu=${cpu}s rss=${rss}kb"
  LAST_WALL=$wall
}

log "=== ISOLATED timing pass START (sequential, T=$T) ==="
printf 'species\tdata_type\taccession\twall_sec\tcpu_sec\tmax_rss_kb\n' > "$ISO"
printf 'order\tspecies\tngs45_unit\tngs45_ITS\tribo_sites\teasy45_unit\teasy45_ITS\tunitPID/mm\tITS_ngsVSez\tITS_ngsVSgb\tITS_ezVSgb\tngs45_wall_s\teasy45_wall_s\n' > "$FINAL"

grep -vE '^DONE_SELECT|^$' "$MAN" | while IFS=$'\t' read -r order sp hifi hm illu im; do
  slug=$(slugof "$sp"); d=$DATA/$slug; r=$RES/$slug; mkdir -p "$r"
  log "--- $sp ---"
  NG_WALL=NA; EZ_WALL=NA; LAST_WALL=NA
  # ngs45 on Illumina
  if [ -s "$d/illu_1.fastq" ]; then
    rm -rf "$r/ngs45_out"
    timed "$sp" "Illumina/ngs45" "$illu" "$r/ngs45.time" \
      $PY -m ngs45.cli run -1 "$d/illu_1.fastq" -2 "$d/illu_2.fastq" \
        -o "$r/ngs45_out" -t "$T" --bait-rounds 3 --call-variants
    NG_WALL=$LAST_WALL
  else
    log "  [$sp] no Illumina data, skip ngs45"
  fi
  # easy45 on HiFi (uses rescued long-read HiFi if the rescue step swapped it in)
  if [ -s "$d/hifi.fastq" ]; then
    rm -rf "$r/easy45_out"
    timed "$sp" "HiFi/easy45" "$hifi" "$r/easy45.time" \
      "$EASY45" run -i "$d/hifi.fastq" -o "$r/easy45_out" -t "$T"
    EZ_WALL=$LAST_WALL
  else
    log "  [$sp] no HiFi data, skip easy45"
  fi
  # ---- final comparison row (biology + isolated timing) ----
  nu=$(seqlen "$r/ngs45_out/nrDNA_45S.fasta"); ni=$(seqlen "$r/ngs45_out/its.fasta")
  ez=$(seqlen "$r/easy45_out/consensus.fasta"); ei=$(seqlen "$r/easy45_out/its.fasta")
  ribo=$( [ -s "$r/ngs45_out/summary.tsv" ] && tail -1 "$r/ngs45_out/summary.tsv" | awk -F'\t' '{print $NF}' || echo NA )
  ref=$BENCH/results/$slug/its_ref.fasta
  if [ ! -s "$ref" ]; then
    esearch -db nuccore -query "\"$sp\"[Organism] AND internal transcribed spacer[Title] AND 400:900[SLEN]" </dev/null 2>/dev/null \
      | efetch -format fasta 2>/dev/null | head -n 800 > "$ref" 2>/dev/null
  fi
  A=$(best_hit "$r/ngs45_out/nrDNA_45S.fasta" "$r/easy45_out/consensus.fasta")
  B=$(best_hit "$r/ngs45_out/its.fasta" "$r/easy45_out/its.fasta")
  C=$(best_hit "$r/ngs45_out/its.fasta" "$ref")
  D=$(best_hit "$r/easy45_out/its.fasta" "$ref")
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$order" "$sp" "$nu" "$ni" "$ribo" "$ez" "$ei" \
    "$(echo $A|cut -d' ' -f1)/$(echo $A|cut -d' ' -f3)mm" \
    "$(echo $B|cut -d' ' -f1)" "$(echo $C|cut -d' ' -f1)" "$(echo $D|cut -d' ' -f1)" \
    "$NG_WALL" "$EZ_WALL" >> "$FINAL"
  log "  [$sp] final: ngs45=${nu}bp easy45=${ez}bp unitPID=$(echo $A|cut -d' ' -f1) ITS_vs_GB=$(echo $C|cut -d' ' -f1)"
done

{
  echo "=== ngs45 benchmark — FINAL (isolated, sequential, T=$T) — $(now) ==="
  echo; echo "## Biology + isolated timing:"; column -t -s$'\t' "$FINAL"
  echo; echo "## Per-tool timing (wall/cpu/RSS):"; column -t -s$'\t' "$ISO"
} > "$BENCH/BENCHMARK_FINAL.txt"
touch "$BENCH/ISOLATED_DONE"
log "=== ISOLATED timing pass DONE -> $FINAL , $ISO ==="
