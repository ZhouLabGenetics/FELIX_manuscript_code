#!/bin/bash
#SBATCH --job-name=3way_plink
#SBATCH --time=01:00:00
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=log/%x_%A.out
#SBATCH --error=log/%x_%A.err

## Build a pruned PLINK fileset (~3000 markers) used by SAIGE createSparseGRM.
##   MODE=common  sbatch 07_make_plink_pruned.sh
##   MODE=lowfreq sbatch 07_make_plink_pruned.sh

set -euo pipefail
MODE=${MODE:?MODE (common|lowfreq) must be exported}
module load singularity

BASE="${FELIX_SIM_BASE:?Set FELIX_SIM_BASE to the simulation directory before submitting}"
export FELIX_SIM_BASE
RTOOLS=/data/wzhougroup/lhu/tools/rtools_latest.sif
SING="singularity exec --bind /data/wzhougroup/lhu:/data/wzhougroup/lhu --home /data/wzhougroup/lhu"

$SING $RTOOLS Rscript ${BASE}/R/make_plink_pruned.R "${MODE}"
