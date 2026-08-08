#!/bin/bash
## Build the Tractor-Mix benchmark manifest and submit the array job.
##
## The manifest lists 9 conditions x N_SEEDS seeds = 45 tasks by default,
## with mid-of-grid betas per (trait, scenario). Edit the CONDITIONS block
## below if you want to expand the benchmark.
##
## Usage:
##   MODE=common bash submit_17_tractormix_bench.sh
##
## Optional env overrides:
##   N_SEEDS      (default 5)    seeds per condition
##   TIMEOUT      (default 20:00:00) per-task wall limit (DNF cutoff for bin01)
##   MEM          (default 32G)
##   TMIX_SIF     (default /data/wzhougroup/lhu/tools/tractormix.sif)

set -euo pipefail
MODE=${MODE:?MODE (common|lowfreq) must be exported}

N_SEEDS=${N_SEEDS:-5}
TIMEOUT=${TIMEOUT:-20:00:00}
MEM=${MEM:-32G}

BASE="${FELIX_SIM_BASE:?Set FELIX_SIM_BASE to the simulation directory before submitting}"
export FELIX_SIM_BASE
MANIFEST_DIR=${BASE}/data/${MODE}/tractormix
mkdir -p "$MANIFEST_DIR"
MANIFEST=${MANIFEST_DIR}/bench_manifest.tsv

## Mid-of-grid betas per (trait, scenario), matching the BETA_GRIDS in
## generate_alt_pheno.R. Update if you re-tune the grids.
##   trait    scen     beta
read -r -d '' CONDITIONS <<'EOF' || true
quant   shared   0.30
quant   afr      0.30
quant   hetero   0.30
bin10   shared   1.20
bin10   afr      3.00
bin10   hetero   1.80
bin01   shared   1.80
bin01   afr      3.00
bin01   hetero   3.00
EOF

## Write manifest: idx<TAB>trait<TAB>scen<TAB>beta<TAB>seed
printf "idx\ttrait\tscen\tbeta\tseed\n" > "$MANIFEST"
idx=0
while IFS=$' \t' read -r tr sc bt; do
    [[ -z "$tr" ]] && continue
    for s in $(seq 1 "$N_SEEDS"); do
        idx=$((idx + 1))
        printf "%d\t%s\t%s\t%s\t%d\n" "$idx" "$tr" "$sc" "$bt" "$s" >> "$MANIFEST"
    done
done <<< "$CONDITIONS"

N_TASKS=$idx
echo "MODE=${MODE}"
echo "Wrote manifest: $MANIFEST  (${N_TASKS} tasks)"

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SCRIPT=${SCRIPT_DIR}/17_tractormix_bench.sh

JOB=$(sbatch --parsable \
    --array=1-${N_TASKS} \
    --time=${TIMEOUT} \
    --mem=${MEM} \
    --export=ALL,MODE=${MODE},MANIFEST=${MANIFEST} \
    "${SCRIPT}")
echo "Submitted job ${JOB} (array=1-${N_TASKS}, time=${TIMEOUT}, mem=${MEM})"
