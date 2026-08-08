#!/bin/bash
## Wrapper for 14_saige_alt.sh that handles the cluster's array-size limit
## (default 1000) by chaining multiple sbatch submissions, each with a
## different START_BATCH_OFFSET. Chunks run back-to-back via afterany
## dependency, so the queue is well-behaved even with many phenos.
##
## Usage:
##   MODE=common  bash submit_14_saige_alt.sh
##   MODE=lowfreq bash submit_14_saige_alt.sh
##
## Optional env overrides:
##   BATCH_SIZE   (default 10)   phenos per array task
##   ARRAY_LIMIT  (default 1000) cluster max array size

set -euo pipefail
MODE=${MODE:?MODE (common|lowfreq) must be exported}

BATCH_SIZE=${BATCH_SIZE:-10}
ARRAY_LIMIT=${ARRAY_LIMIT:-1000}

BASE="${FELIX_SIM_BASE:?Set FELIX_SIM_BASE to the simulation directory before submitting}"
export FELIX_SIM_BASE
PHENO_DIR=${BASE}/data/${MODE}/pheno/alt

if [[ ! -d $PHENO_DIR ]]; then
    echo "Pheno dir not found: $PHENO_DIR" >&2
    exit 1
fi

## Count pheno files directly (not via manifest, so we don't depend on its
## format). 14_saige_alt.sh uses the same glob.
N_PHENOS=$(cd "$PHENO_DIR" && ls *.tsv 2>/dev/null \
           | grep -vE '^(manifest|causal)\.' | wc -l)
if (( N_PHENOS == 0 )); then
    echo "No pheno files in $PHENO_DIR" >&2
    exit 1
fi
N_BATCHES=$(( (N_PHENOS + BATCH_SIZE - 1) / BATCH_SIZE ))

echo "MODE=${MODE}"
echo "Manifest phenos: ${N_PHENOS}"
echo "Batches needed:  ${N_BATCHES}  (BATCH_SIZE=${BATCH_SIZE})"
echo "Array limit:     ${ARRAY_LIMIT}"

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SCRIPT=${SCRIPT_DIR}/14_saige_alt.sh

PREV_JOB=""
CHUNK=0
while (( CHUNK * ARRAY_LIMIT < N_BATCHES )); do
    OFFSET=$(( CHUNK * ARRAY_LIMIT ))
    REMAINING=$(( N_BATCHES - OFFSET ))
    NTASKS=$REMAINING
    if (( NTASKS > ARRAY_LIMIT )); then NTASKS=$ARRAY_LIMIT; fi

    DEP=""
    if [[ -n "$PREV_JOB" ]]; then DEP="--dependency=afterany:${PREV_JOB}"; fi

    JOB=$(sbatch --parsable \
        --array=1-${NTASKS} \
        $DEP \
        --export=ALL,MODE=${MODE},BATCH_SIZE=${BATCH_SIZE},START_BATCH_OFFSET=${OFFSET} \
        "${SCRIPT}")

    printf "  chunk %d: array=1-%d offset=%d job=%s%s\n" \
        "$CHUNK" "$NTASKS" "$OFFSET" "$JOB" \
        "$( [[ -n "$DEP" ]] && echo " (afterany:${PREV_JOB})" )"

    PREV_JOB=$JOB
    CHUNK=$((CHUNK + 1))
done

echo "Submitted ${CHUNK} chunk(s) covering ${N_BATCHES} batches."
echo "Final job id: ${PREV_JOB}"
