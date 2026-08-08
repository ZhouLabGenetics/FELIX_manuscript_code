#!/bin/bash
#SBATCH --job-name=3way_agg_res
#SBATCH --time=01:00:00
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=log/%x_%A.out
#SBATCH --error=log/%x_%A.err

## Aggregate FELIX output into Type-I error and Power tables.
##   MODE=common  sbatch 16_aggregate_results.sh
##   MODE=lowfreq sbatch 16_aggregate_results.sh

set -euo pipefail
MODE=${MODE:?MODE (common|lowfreq) must be exported}
module load singularity

BASE="${FELIX_SIM_BASE:?Set FELIX_SIM_BASE to the simulation directory before submitting}"
export FELIX_SIM_BASE
RTOOLS=/data/wzhougroup/lhu/tools/rtools_latest.sif
SING="singularity exec --bind /data/wzhougroup/lhu:/data/wzhougroup/lhu --home /data/wzhougroup/lhu"

$SING $RTOOLS Rscript ${BASE}/R/aggregate_type1.R "${MODE}"
$SING $RTOOLS Rscript ${BASE}/R/aggregate_power.R "${MODE}"
