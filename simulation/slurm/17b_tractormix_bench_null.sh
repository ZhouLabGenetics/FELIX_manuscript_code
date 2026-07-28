#!/bin/bash
#SBATCH --job-name=3way_tmix_null
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --output=/data/wzhougroup/lhu/saige_tractor/simulation/3way/log/%x_%A_%a.out
#SBATCH --error=/data/wzhougroup/lhu/saige_tractor/simulation/3way/log/%x_%A_%a.err

## Run one Tractor-Mix Type-I (null) benchmark task. Wall time + peak RSS
## are captured by /usr/bin/time -v into bench.log inside the TAG output
## dir. SLURM --time wall is set by the submit wrapper (default 20h).
##
## Run via submit_17b_tractormix_bench_null.sh, which sets MANIFEST +
## array size.

set -euo pipefail
MODE=${MODE:?MODE must be exported}
MANIFEST=${MANIFEST:?MANIFEST must be exported}
module load singularity

BASE=/data/wzhougroup/lhu/saige_tractor/simulation/3way
TMIX_SIF=${TMIX_SIF:-/data/wzhougroup/lhu/tools/tractormix.sif}
SING="singularity exec --bind /data/wzhougroup/lhu:/data/wzhougroup/lhu --home /data/wzhougroup/lhu"

## Pull row $SLURM_ARRAY_TASK_ID from the manifest (skip header).
ROW=$(awk -v i="${SLURM_ARRAY_TASK_ID}" 'NR>1 && $1==i' "$MANIFEST")
if [[ -z "$ROW" ]]; then
    echo "No row ${SLURM_ARRAY_TASK_ID} in $MANIFEST" >&2
    exit 1
fi
IFS=$'\t' read -r IDX TRAIT SEED <<< "$ROW"
echo "idx=${IDX} trait=${TRAIT} seed=${SEED}"

TAG=${TRAIT}_seed$(printf "%02d" "$SEED")
OUT_DIR=${BASE}/data/${MODE}/tractormix_out/null/${TAG}
mkdir -p "$OUT_DIR"
BENCH_LOG=${OUT_DIR}/bench.log

## To re-run with the patched timing methodology (time outside container),
## set FORCE_RERUN=1 when submitting; otherwise existing results are kept.
if [[ -s ${OUT_DIR}/result.tsv && "${FORCE_RERUN:-0}" != "1" ]]; then
    echo "Already done: ${TAG} -- skipping (set FORCE_RERUN=1 to redo)."
    exit 0
fi
if [[ "${FORCE_RERUN:-0}" == "1" ]]; then
    rm -f ${OUT_DIR}/result.tsv ${OUT_DIR}/null.rda ${OUT_DIR}/bench.log
fi

## /usr/bin/time -v on the HOST wraps the singularity invocation so we
## also capture container startup -- matches the FELIX bench (17c)
## measurement methodology, which times the singularity exec from outside.
## TM_NCORE is consumed by run_tractormix_null.R's TractorMix.score call.
/usr/bin/time -v -o "$BENCH_LOG" \
    $SING $TMIX_SIF \
        Rscript ${BASE}/scripts/R/run_tractormix_null.R \
            "$MODE" "$TRAIT" "$SEED"

echo "Done: ${TAG}"
echo "Bench log: $BENCH_LOG"
