#!/bin/bash
# Definitive benchmark pass:
#   1) acquire long-read HiFi subsets (ENA streaming) for species that need them
#   2) isolated, sequential ngs45 + easy45 with clean /usr/bin/time timing and the
#      final comparison table (benchmark_final.tsv).
# Illumina for all species is already downloaded under bench/data/<sp>/.
set -u
BENCH=/path/to/ngs245/bench
L=$BENCH/definitive.log
log(){ echo "[$(date '+%F %T')] $*" | tee -a "$L"; }

log "STEP 1/2: HiFi subset download (ENA)"
bash "$BENCH/download_hifi_subset.sh" >> "$BENCH/hifi_subset.stdout" 2>&1
log "STEP 2/2: isolated timing + final table"
bash "$BENCH/run_isolated_timing.sh"   >> "$BENCH/isolated.stdout" 2>&1
log "DEFINITIVE DONE -> $BENCH/benchmark_final.tsv , $BENCH/timing_isolated.tsv, $BENCH/BENCHMARK_FINAL.txt"
