#!/bin/bash
# ============================================================================
# ngs45 overnight benchmark driver
#   Phase 1: download HiFi + Illumina for every species in manifest.tsv (parallel)
#   Phase 2: run ngs45 (Illumina) + easy45 (HiFi) per species (bounded parallel)
#   Phase 3: fetch GenBank ITS ref + cross-compare -> benchmark_results.tsv
# Designed to run UNATTENDED overnight: never aborts on a single failure,
# logs everything with timestamps, resumable (skips finished steps).
# ============================================================================
set -u
BASE=/path/to/ngs245
BENCH=$BASE/bench
DATA=$BENCH/data
RES=$BENCH/results
TMP=$BENCH/tmp; mkdir -p "$DATA" "$RES" "$TMP"
MAN=$BENCH/manifest.tsv
TABLE=$BENCH/benchmark_results.tsv
LOG=$BENCH/overnight.log
STATUS=$BENCH/STATUS.txt

# ---- composite PATH across the conda envs that hold each tool --------------
E=/path/to/.conda/envs
export PATH="$E/easy45/bin:$E/getorganelle_env/bin:$E/bcftools_env/bin:$E/sra/bin:$E/edirect_env/bin:/usr/bin:$PATH"
export PYTHONPATH="$BASE/src"
PY="$E/easy45/bin/python"
NGS45="$PY -m ngs45.cli"
EASY45="$E/easy45/bin/easy45"

# Low-footprint settings: another user has priority, so run everything at the
# lowest CPU/IO scheduling priority and stay nearly sequential. Slower, but only
# consumes otherwise-idle cycles.
T=16           # threads per tool
WORKERS=6      # assemble up to 6 species concurrently (128 cores, machine idle)
LP=""          # MAX SPEED: no nice/ionice throttle (priority window is over)

now(){ date '+%Y-%m-%d %H:%M:%S'; }
log(){ echo "[$(now)] $*" | tee -a "$LOG"; }
slugof(){ echo "$1" | tr ' ' '_'; }
# count how many of the given PIDs are still alive (assembly-slot accounting)
alive_count(){ local c=0 p; for p in "$@"; do kill -0 "$p" 2>/dev/null && c=$((c+1)); done; echo "$c"; }
# block until a species' inputs are on disk (Illumina always; HiFi unless NA)
wait_for_data(){
  local sp="$1" hifi="$2" d; d=$DATA/$(slugof "$sp"); local n=0
  while [ ! -s "$d/illu_1.fastq" ] && [ $n -lt 2400 ]; do sleep 15; n=$((n+1)); done
  if [ "$hifi" != NA ]; then
    n=0; while [ ! -s "$d/hifi.fastq" ] && [ $n -lt 2400 ]; do sleep 15; n=$((n+1)); done
  fi
}

# best blastn hit: prints "pident length mismatch" (highest bitscore) or "NA NA NA"
best_hit(){
  local q="$1" s="$2"
  [ -s "$q" ] && [ -s "$s" ] || { echo "NA NA NA"; return; }
  local r; r=$(blastn -query "$q" -subject "$s" -outfmt "6 pident length mismatch bitscore" 2>/dev/null | sort -k4 -nr | head -1)
  [ -n "$r" ] && echo "$r" | awk '{print $1,$2,$3}' || echo "NA NA NA"
}
seqlen(){ [ -s "$1" ] && seqkit stats -T "$1" 2>/dev/null | tail -1 | cut -f5 || echo NA; }

# ---------------------------------------------------------------------------
download_one(){
  local sp="$1" hifi="$2" illu="$3" d; d=$DATA/$(slugof "$sp"); mkdir -p "$d"
  echo "DOWNLOADING $(now)" > "$d/STATE"
  # Illumina (paired)
  if [ "$illu" != NA ] && [ ! -s "$d/illu_1.fastq" ]; then
    log "  [$sp] prefetch Illumina $illu"
    $LP timeout 4h prefetch "$illu" -O "$d" --max-size 25g > "$d/pf_illu.log" 2>&1
    $LP timeout 4h fasterq-dump "$d/$illu/$illu.sra" --split-3 -e "$T" -O "$d" -t "$TMP" > "$d/fq_illu.log" 2>&1
    [ -s "$d/${illu}_1.fastq" ] && mv -f "$d/${illu}_1.fastq" "$d/illu_1.fastq"
    [ -s "$d/${illu}_2.fastq" ] && mv -f "$d/${illu}_2.fastq" "$d/illu_2.fastq"
    rm -rf "$d/$illu" "$d/${illu}.fastq"   # drop .sra + singletons
  fi
  # HiFi (single)
  if [ "$hifi" != NA ] && [ ! -s "$d/hifi.fastq" ]; then
    log "  [$sp] prefetch HiFi $hifi"
    $LP timeout 6h prefetch "$hifi" -O "$d" --max-size 25g > "$d/pf_hifi.log" 2>&1
    $LP timeout 6h fasterq-dump "$d/$hifi/$hifi.sra" -e "$T" -O "$d" -t "$TMP" > "$d/fq_hifi.log" 2>&1
    [ -s "$d/${hifi}.fastq" ] && mv -f "$d/${hifi}.fastq" "$d/hifi.fastq"
    rm -rf "$d/$hifi"
  fi
  echo "DOWNLOADED $(now)  illu=$(seqlen "$d/illu_1.fastq")bp hifi=$(seqlen "$d/hifi.fastq")bp" > "$d/STATE"
  log "  [$sp] download done"
}

process_one(){
  local order="$1" sp="$2" hifi="$3" illu="$4"
  local slug; slug=$(slugof "$sp"); local d=$DATA/$slug r=$RES/$slug; mkdir -p "$r"
  echo "ASSEMBLING $(now)" > "$r/STATE"

  # ngs45 on Illumina  (wall-clock timed)
  local ngs45_sec=NA easy45_sec=NA t0
  if [ -s "$d/illu_1.fastq" ] && [ ! -s "$r/ngs45_out/nrDNA_45S.fasta" ]; then
    log "  [$sp] ngs45 (Illumina)"
    t0=$(date +%s)
    $LP timeout 4h $NGS45 run -1 "$d/illu_1.fastq" -2 "$d/illu_2.fastq" -o "$r/ngs45_out" \
        -t "$T" --subsample 200000 --bait-rounds 2 --call-variants -v > "$r/ngs45.log" 2>&1 \
        || log "  [$sp] ngs45 FAILED (see ngs45.log)"
    ngs45_sec=$(( $(date +%s) - t0 ))
    printf '%s\tIllumina/ngs45\t%s\t%s\n' "$sp" "$illu" "$ngs45_sec" >> "$BENCH/timing.tsv"
  fi
  # easy45 on HiFi  (wall-clock timed)
  if [ -s "$d/hifi.fastq" ] && [ ! -s "$r/easy45_out/consensus.fasta" ]; then
    log "  [$sp] easy45 (HiFi)"
    t0=$(date +%s)
    $LP timeout 4h $EASY45 run -i "$d/hifi.fastq" -o "$r/easy45_out" -t "$T" -v > "$r/easy45.log" 2>&1 \
        || log "  [$sp] easy45 FAILED (see easy45.log)"
    easy45_sec=$(( $(date +%s) - t0 ))
    printf '%s\tHiFi/easy45\t%s\t%s\n' "$sp" "$hifi" "$easy45_sec" >> "$BENCH/timing.tsv"
  fi
  # GenBank ITS reference (best-effort, retries)
  if [ ! -s "$r/its_ref.fasta" ]; then
    for a in 1 2 3; do
      esearch -db nuccore -query "\"$sp\"[Organism] AND internal transcribed spacer[Title] AND 400:900[SLEN]" </dev/null 2>/dev/null \
        | efetch -format fasta 2>/dev/null | head -n 800 > "$r/its_ref.fasta"
      [ -s "$r/its_ref.fasta" ] && break; sleep 15
    done
  fi

  # ---- comparisons ----
  local ngs_unit=$r/ngs45_out/nrDNA_45S.fasta  ngs_its=$r/ngs45_out/its.fasta
  local ez_cons=$r/easy45_out/consensus.fasta  ez_its=$r/easy45_out/its.fasta
  local ref=$r/its_ref.fasta
  local ngs_ulen ngs_ilen ez_ulen ez_ilen ribo
  ngs_ulen=$(seqlen "$ngs_unit"); ngs_ilen=$(seqlen "$ngs_its")
  ez_ulen=$(seqlen "$ez_cons");   ez_ilen=$(seqlen "$ez_its")
  ribo=$( [ -s "$r/ngs45_out/summary.tsv" ] && tail -1 "$r/ngs45_out/summary.tsv" | awk -F'\t' '{print $NF}' || echo NA )
  local A B C D
  A=$(best_hit "$ngs_unit" "$ez_cons")   # ngs45 unit vs easy45 consensus
  B=$(best_hit "$ngs_its"  "$ez_its")    # ngs45 ITS  vs easy45 ITS
  C=$(best_hit "$ngs_its"  "$ref")       # ngs45 ITS  vs GenBank
  D=$(best_hit "$ez_its"   "$ref")       # easy45 ITS vs GenBank
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$order" "$sp" "$illu" "$hifi" \
    "$ngs_ulen" "$ngs_ilen" "$ribo" "$ez_ulen" "$ez_ilen" \
    "$(echo $A|cut -d' ' -f1)/$(echo $A|cut -d' ' -f3)mm" \
    "$(echo $B|cut -d' ' -f1)" "$(echo $C|cut -d' ' -f1)" "$(echo $D|cut -d' ' -f1)" \
    "$ngs45_sec" "$easy45_sec" "OK" >> "$TABLE"
  echo "DONE $(now)" > "$r/STATE"
  log "  [$sp] compared: ngs45_unit=$ngs_ulen ez_unit=$ez_ulen  unit_pid=$(echo $A|cut -d' ' -f1)  ITS_vs_GB=$(echo $C|cut -d' ' -f1)"
}

# ============================ MAIN =========================================
log "=== overnight benchmark START ==="
mapfile -t ROWS < <(grep -vE '^DONE_SELECT|^$' "$MAN")
log "manifest: ${#ROWS[@]} species"

# headers
printf 'order\tspecies\tillu_srr\thifi_srr\tngs45_unit_len\tngs45_its_len\tngs45_ribo_sites\teasy45_unit_len\teasy45_its_len\tunitPID/mm\tITS_ngsVSez\tITS_ngsVSgb\tITS_ezVSgb\tngs45_sec\teasy45_sec\tstatus\n' > "$TABLE"
# per-tool wall-clock timing (NOTE: species run up to WORKERS-concurrent, so these
# are throughput timings under contention, not isolated single-job benchmarks)
printf 'species\tdata_type\taccession\twall_sec\n' > "$BENCH/timing.tsv"

# ---- Overlapped: download everything in parallel; assemble each species as
#      soon as its own inputs land (up to $WORKERS assemblies at a time) ----
log "--- launching all downloads in parallel ---"
for row in "${ROWS[@]}"; do
  IFS=$'\t' read -r order sp hifi hm illu im <<< "$row"
  download_one "$sp" "$hifi" "$illu" &
done

log "--- assembling as data lands (max $WORKERS concurrent) ---"
APIDS=()
for row in "${ROWS[@]}"; do
  IFS=$'\t' read -r order sp hifi hm illu im <<< "$row"
  while [ "$(alive_count ${APIDS[@]+"${APIDS[@]}"})" -ge "$WORKERS" ]; do sleep 10; done
  { wait_for_data "$sp" "$hifi"; process_one "$order" "$sp" "$hifi" "$illu"; } &
  APIDS+=("$!")
done
wait
log "--- all species processed ---"

# ---- final summary ----
{ echo "ngs45 overnight benchmark — finished $(now)"; echo; column -t -s$'\t' "$TABLE"; } > "$BENCH/SUMMARY.txt"
touch "$BENCH/OVERNIGHT_DONE"
log "=== overnight benchmark DONE -> $TABLE ==="
