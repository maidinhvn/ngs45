#!/bin/bash
# Wait until the concurrent benchmark round is fully finished, THEN run the
# isolated (sequential, contention-free) timing pass. Launched in the background
# now; it idles cheaply until the concurrent round completes.
set -u
BENCH=/path/to/ngs245/bench
CLOG=$BENCH/chain.log
log(){ echo "[$(date '+%F %T')] $*" | tee -a "$CLOG"; }

log "chain started; waiting for concurrent round (OVERNIGHT_DONE)..."
while [ ! -f "$BENCH/OVERNIGHT_DONE" ]; do sleep 60; done
log "OVERNIGHT_DONE seen; draining any lingering assembly processes..."
while ps -u vutrinh -o args 2>/dev/null \
      | grep -E 'ngs45\.cli|/easy45 run|s2_spades|spades\.py|bowtie2-align|minimap2' \
      | grep -v grep | grep -q .; do sleep 30; done
log "machine clear."
# 1) HiFi rescue: swap short-read HiFi -> long-read HiFi for species whose easy45 failed
log "step 1/2: HiFi rescue (long-read runs for failed easy45)"
bash "$BENCH/run_hifi_rescue.sh" >> "$BENCH/rescue.stdout" 2>&1
# 2) isolated timing + final table (uses rescued HiFi)
log "step 2/2: ISOLATED timing pass + final table"
bash "$BENCH/run_isolated_timing.sh" >> "$BENCH/isolated.stdout" 2>&1
log "ALL DONE -> $BENCH/benchmark_final.tsv , $BENCH/timing_isolated.tsv"
