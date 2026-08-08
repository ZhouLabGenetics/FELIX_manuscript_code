#!/bin/bash
#SBATCH --job-name=3way_sparseGRM
#SBATCH --time=02:00:00
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=log/%x_%A.out
#SBATCH --error=log/%x_%A.err

## Build the sparse GRM for FELIX using SAIGE's createSparseGRM.R.
## Output files: <mode>/plink/sparseGRM_0.125_2000.sparseGRM.mtx
##               <mode>/plink/sparseGRM_0.125_2000.sparseGRM.mtx.sampleIDs.txt
##   MODE=common  sbatch 08_create_sparse_grm.sh
##   MODE=lowfreq sbatch 08_create_sparse_grm.sh

set -euo pipefail
MODE=${MODE:?MODE (common|lowfreq) must be exported}
module load singularity 
BASE="${FELIX_SIM_BASE:?Set FELIX_SIM_BASE to the simulation directory before submitting}"
export FELIX_SIM_BASE
SAIGE=/data/wzhougroup/lhu/tools/saige_151.sif
SING="singularity exec --bind /data/wzhougroup/lhu:/data/wzhougroup/lhu --home /data/wzhougroup/lhu"

PLINK_PREFIX=${BASE}/data/${MODE}/plink/pruned
OUT_PREFIX=${BASE}/data/${MODE}/plink/sparseGRM

$SING $SAIGE createSparseGRM.R \
    --plinkFile=${PLINK_PREFIX} \
    --nThreads=4 \
    --outputPrefix=${OUT_PREFIX} \
    --numRandomMarkerforSparseKin=2000 \
    --relatednessCutoff=0.125
