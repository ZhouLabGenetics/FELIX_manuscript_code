#!/bin/bash
#SBATCH --job-name=3way_tmix_prep
#SBATCH --time=00:30:00
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --output=log/%x_%A.out
#SBATCH --error=log/%x_%A.err

## One-time prep: convert Block_row1/G{Afr,Eur,Nat}.tsv to Tractor-Mix
## hapdose layout under data/<MODE>/tractormix/hapdose/.
##
## Uses the rtools container (just needs data.table) -- Tractor-Mix container
## not required at this step.
##
##   MODE=common  sbatch 15_make_tractormix_hapdose.sh
##   MODE=lowfreq sbatch 15_make_tractormix_hapdose.sh

set -euo pipefail
MODE=${MODE:?MODE (common|lowfreq) must be exported}
module load singularity

BASE="${FELIX_SIM_BASE:?Set FELIX_SIM_BASE to the simulation directory before submitting}"
export FELIX_SIM_BASE
RTOOLS=/data/wzhougroup/lhu/tools/rtools_latest.sif
SING="singularity exec --bind /data/wzhougroup/lhu:/data/wzhougroup/lhu --home /data/wzhougroup/lhu"

$SING $RTOOLS Rscript ${BASE}/R/make_tractormix_hapdose.R "${MODE}"
