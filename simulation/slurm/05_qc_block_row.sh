#!/bin/bash
#SBATCH --job-name=3way_qcbr
#SBATCH --time=04:00:00
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --output=log/%x_%A_%a.out
#SBATCH --error=log/%x_%A_%a.err
#SBATCH --array=1-20

## MAF-based QC per block_row.
##   MODE=common  sbatch 05_qc_block_row.sh
##   MODE=lowfreq sbatch 05_qc_block_row.sh

set -euo pipefail
MODE=${MODE:?MODE (common|lowfreq) must be exported}
module load singularity

BASE="${FELIX_SIM_BASE:?Set FELIX_SIM_BASE to the simulation directory before submitting}"
export FELIX_SIM_BASE
RTOOLS=/data/wzhougroup/lhu/tools/rtools_latest.sif
SING="singularity exec --bind /data/wzhougroup/lhu:/data/wzhougroup/lhu --home /data/wzhougroup/lhu"

$SING $RTOOLS Rscript ${BASE}/R/qc_Block_row.R "${SLURM_ARRAY_TASK_ID}" "${MODE}"
