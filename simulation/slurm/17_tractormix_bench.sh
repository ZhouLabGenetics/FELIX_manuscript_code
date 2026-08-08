#!/bin/bash
#SBATCH --job-name=3way_tmix
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --output=log/%x_%A_%a.out
#SBATCH --error=log/%x_%A_%a.err

## Run one Tractor-Mix benchmark task (one row of MANIFEST). Wall time and
## peak RSS are captured by /usr/bin/time -v into bench.log inside the TAG
## output dir. The SLURM --time wall is set by the submit wrapper (default
## 20h), so bin01 runs that don't finish are reported as DNF by the aggregator.
##
## Run via submit_17_tractormix_bench.sh, which sets MANIFEST + array size.

set -euo pipefail
MODE=${MODE:?MODE must be exported}
MANIFEST=${MANIFEST:?MANIFEST must be exported}
module load singularity

BASE="${FELIX_SIM_BASE:?Set FELIX_SIM_BASE to the simulation directory before submitting}"
export FELIX_SIM_BASE
TMIX_SIF=${TMIX_SIF:-/data/wzhougroup/lhu/tools/tractormix.sif}
SING="singularity exec --bind /data/wzhougroup/lhu:/data/wzhougroup/lhu --home /data/wzhougroup/lhu"

## Pull row $SLURM_ARRAY_TASK_ID from the manifest (skip header).
ROW=$(awk -v i="${SLURM_ARRAY_TASK_ID}" 'NR>1 && $1==i' "$MANIFEST")
if [[ -z "$ROW" ]]; then
    echo "No row ${SLURM_ARRAY_TASK_ID} in $MANIFEST" >&2
    exit 1
fi
IFS=$'\t' read -r IDX TRAIT SCEN BETA SEED <<< "$ROW"
echo "idx=${IDX} trait=${TRAIT} scen=${SCEN} beta=${BETA} seed=${SEED}"

BETA_TAG=$(awk -v b="$BETA" 'BEGIN{ printf "beta%03d", int(b*100 + 0.5) }')
TAG=${TRAIT}_${SCEN}_${BETA_TAG}_seed$(printf "%02d" "$SEED")
OUT_DIR=${BASE}/data/${MODE}/tractormix_out/alt/${TAG}
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
## TM_NCORE is consumed by run_tractormix.R's TractorMix.score call.
/usr/bin/time -v -o "$BENCH_LOG" \
    $SING $TMIX_SIF \
        Rscript ${BASE}/R/run_tractormix.R \
            "$MODE" "$TRAIT" "$SCEN" "$BETA" "$SEED"

echo "Done: ${TAG}"
echo "Bench log: $BENCH_LOG"
