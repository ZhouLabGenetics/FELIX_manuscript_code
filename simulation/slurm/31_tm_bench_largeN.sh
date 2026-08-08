#!/bin/bash
#SBATCH --job-name=3way_tm_bench_largeN
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=log/%x_%A_%a.out
#SBATCH --error=log/%x_%A_%a.err

## Tractor-Mix large-N type-I benchmark on the bench_N<...>_P<...> data.
## Wraps /usr/bin/time -v OUTSIDE the singularity invocation (matches the
## SAIGE bench measurement methodology -- captures container startup too).
##
## Array layout: 15 tasks = 5 N's x 3 traits.
##   id  1..3 : N=10000,  trait in {quant, bin10, bin01}  -- fast baseline
##   id  4..6 : N=50000,  trait in {quant, bin10, bin01}
##   id  7..9 : N=100000, trait in {quant, bin10, bin01}
##   id 10..12: N=150000, trait in {quant, bin10, bin01}  -- likely OOM
##   id 13..15: N=200000, trait in {quant, bin10, bin01}  -- expected OOM
##
## TM materialises sparse K as dense at N>=150k (180+ GB) and either OOMs or
## takes impractically long. That failure point IS the benchmark headline.
##
## Submit:
##   sbatch --array=1-3   --mem=16G  --time=01:00:00  31_tm_bench_largeN.sh    # N=10k (fast)
##   sbatch --array=4-6   --mem=64G  --time=04:00:00  31_tm_bench_largeN.sh    # N=50k (safe)
##   sbatch --array=7-9   --mem=128G --time=08:00:00  31_tm_bench_largeN.sh    # N=100k (borderline)
##   sbatch --array=10-12 --mem=256G --time=24:00:00  31_tm_bench_largeN.sh    # N=150k (expect OOM)
##   sbatch --array=13-15 --mem=512G --time=24:00:00  31_tm_bench_largeN.sh    # N=200k (expect OOM)

set -euo pipefail
module load singularity

BASE="${FELIX_SIM_BASE:?Set FELIX_SIM_BASE to the simulation directory before submitting}"
export FELIX_SIM_BASE
TMIX_SIF=${TMIX_SIF:-/data/wzhougroup/lhu/tools/tractormix.sif}
SING="singularity exec --bind /data/wzhougroup/lhu:/data/wzhougroup/lhu --home /data/wzhougroup/lhu"

TASK=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID must be set (use sbatch --array=...)}
N_IDX=$(( (TASK - 1) / 3 ))      # 0,1,2,3,4
T_IDX=$(( (TASK - 1) % 3 ))      # 0,1,2

case "$N_IDX" in
  0) BENCH_N=10000  ;;
  1) BENCH_N=50000  ;;
  2) BENCH_N=100000 ;;
  3) BENCH_N=150000 ;;
  4) BENCH_N=200000 ;;
  *) echo "Bad N_IDX=$N_IDX" >&2; exit 1 ;;
esac
case "$T_IDX" in
  0) TRAIT=quant ;;
  1) TRAIT=bin10 ;;
  2) TRAIT=bin01 ;;
esac
BENCH_P=${BENCH_P:-1000}
BENCH_DIR=${BASE}/data/bench_N${BENCH_N}_P${BENCH_P}

export BENCH_DIR BENCH_N BENCH_P
TM_NCORE=${TM_NCORE:-1}; export TM_NCORE

echo "==== TM bench: N=${BENCH_N}  TRAIT=${TRAIT}  ===="
TAG=${TRAIT}_seed01
OUT_DIR=${BENCH_DIR}/tractormix_out/null/${TAG}
mkdir -p "$OUT_DIR"
BENCH_LOG=${OUT_DIR}/bench.log

if [[ -s ${OUT_DIR}/result.tsv && "${FORCE_RERUN:-0}" != "1" ]]; then
    echo "Already done: ${TAG} -- skipping (set FORCE_RERUN=1 to redo)."
    exit 0
fi
if [[ "${FORCE_RERUN:-0}" == "1" ]]; then
    rm -f ${OUT_DIR}/result.tsv ${OUT_DIR}/null.rda ${OUT_DIR}/bench.log
fi

/usr/bin/time -v -o "$BENCH_LOG" \
    $SING $TMIX_SIF \
        Rscript ${BASE}/R/run_tractormix_largeN.R "$TRAIT"

echo "==== Done.  bench log: $BENCH_LOG ===="
