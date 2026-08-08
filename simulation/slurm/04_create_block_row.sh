#!/bin/bash
#SBATCH --job-name=3way_blockrow
#SBATCH --time=08:00:00
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --output=log/%x_%A_%a.out
#SBATCH --error=log/%x_%A_%a.err
#SBATCH --array=1-20

## Reorganise 50 (Ind+Rel) blocks into 20 block_rows for downstream testing.
##   MODE=common  sbatch 04_create_block_row.sh
##   MODE=lowfreq sbatch 04_create_block_row.sh
##
## NB: for lowfreq we simulated 2M variants. Bump --array to 1-40 if you want
## the full 2M worth of block_rows; the downstream stages always iterate
## Block_row1..20 anyway.

set -euo pipefail
MODE=${MODE:?MODE (common|lowfreq) must be exported}
module load singularity

BASE="${FELIX_SIM_BASE:?Set FELIX_SIM_BASE to the simulation directory before submitting}"
export FELIX_SIM_BASE
RTOOLS=/data/wzhougroup/lhu/tools/rtools_latest.sif
SING="singularity exec --bind /data/wzhougroup/lhu:/data/wzhougroup/lhu --home /data/wzhougroup/lhu"

$SING $RTOOLS Rscript ${BASE}/R/create_Block_row.R "${SLURM_ARRAY_TASK_ID}" "${MODE}"
