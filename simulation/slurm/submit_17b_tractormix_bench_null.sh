#!/bin/bash
## Build the Tractor-Mix Type-I (null) benchmark manifest and submit the
## array job. 3 traits x N_SEEDS seeds = 30 tasks by default.
##
## Usage:
##   MODE=common  bash submit_17b_tractormix_bench_null.sh
##   MODE=lowfreq bash submit_17b_tractormix_bench_null.sh
##
## Optional env overrides:
##   N_SEEDS   (default 10)        seeds per trait (matches generate_null_pheno.R)
##   TIMEOUT   (default 20:00:00)  per-task wall limit
##   MEM       (default 32G)
##   TMIX_SIF  (default /data/wzhougroup/lhu/tools/tractormix.sif)

set -euo pipefail
MODE=${MODE:?MODE (common|lowfreq) must be exported}

N_SEEDS=${N_SEEDS:-10}
TIMEOUT=${TIMEOUT:-20:00:00}
MEM=${MEM:-32G}

BASE="${FELIX_SIM_BASE:?Set FELIX_SIM_BASE to the simulation directory before submitting}"
export FELIX_SIM_BASE
MANIFEST_DIR=${BASE}/data/${MODE}/tractormix
mkdir -p "$MANIFEST_DIR"
MANIFEST=${MANIFEST_DIR}/bench_manifest_null.tsv

TRAITS=(quant bin10 bin01)

## Write manifest: idx<TAB>trait<TAB>seed
printf "idx\ttrait\tseed\n" > "$MANIFEST"
idx=0
for tr in "${TRAITS[@]}"; do
    for s in $(seq 1 "$N_SEEDS"); do
        idx=$((idx + 1))
        printf "%d\t%s\t%d\n" "$idx" "$tr" "$s" >> "$MANIFEST"
    done
done

N_TASKS=$idx
echo "MODE=${MODE}"
echo "Wrote manifest: $MANIFEST  (${N_TASKS} tasks)"

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SCRIPT=${SCRIPT_DIR}/17b_tractormix_bench_null.sh

FORCE_RERUN=${FORCE_RERUN:-0}

JOB=$(sbatch --parsable \
    --array=1-${N_TASKS} \
    --time=${TIMEOUT} \
    --mem=${MEM} \
    --export=ALL,MODE=${MODE},MANIFEST=${MANIFEST},FORCE_RERUN=${FORCE_RERUN} \
    "${SCRIPT}")
echo "Submitted job ${JOB} (array=1-${N_TASKS}, time=${TIMEOUT}, mem=${MEM}, FORCE_RERUN=${FORCE_RERUN})"
