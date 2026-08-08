#!/bin/bash
#SBATCH --job-name=3way_ind
#SBATCH --time=12:00:00
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --output=log/%x_%A_%a.out
#SBATCH --error=log/%x_%A_%a.err
#SBATCH --array=1-50

## 50 array tasks x 100 samples = 5000 independent admixed individuals.
##   MODE=common  sbatch 02_simu_ind.sh
##   MODE=lowfreq sbatch 02_simu_ind.sh

set -euo pipefail
MODE=${MODE:?MODE (common|lowfreq) must be exported}
module load singularity

BASE="${FELIX_SIM_BASE:?Set FELIX_SIM_BASE to the simulation directory before submitting}"
export FELIX_SIM_BASE
RTOOLS=/data/wzhougroup/lhu/tools/rtools_latest.sif
SING="singularity exec --bind /data/wzhougroup/lhu:/data/wzhougroup/lhu --home /data/wzhougroup/lhu"

$SING $RTOOLS Rscript ${BASE}/R/Pedigree3way_Ind.R "${SLURM_ARRAY_TASK_ID}" "${MODE}"
