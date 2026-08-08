#!/bin/bash
#SBATCH --job-name=3way_null_pheno
#SBATCH --time=01:00:00
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --output=log/%x_%A.out
#SBATCH --error=log/%x_%A.err

## Generate null phenotypes (10 seeds x 3 trait types) for Type-I error tests.
##   MODE=common  sbatch 10_gen_null_pheno.sh
##   MODE=lowfreq sbatch 10_gen_null_pheno.sh

set -euo pipefail
MODE=${MODE:?MODE (common|lowfreq) must be exported}
module load singularity

BASE="${FELIX_SIM_BASE:?Set FELIX_SIM_BASE to the simulation directory before submitting}"
export FELIX_SIM_BASE
RTOOLS=/data/wzhougroup/lhu/tools/rtools_latest.sif
SING="singularity exec --bind /data/wzhougroup/lhu:/data/wzhougroup/lhu --home /data/wzhougroup/lhu"

$SING $RTOOLS Rscript ${BASE}/R/generate_null_pheno.R "${MODE}"
