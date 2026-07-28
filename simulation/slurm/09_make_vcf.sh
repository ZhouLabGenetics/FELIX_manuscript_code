#!/bin/bash
#SBATCH --job-name=3way_vcf
#SBATCH --time=04:00:00
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --output=/data/wzhougroup/lhu/saige_tractor/simulation/3way/log/%x_%A_%a.out
#SBATCH --error=/data/wzhougroup/lhu/saige_tractor/simulation/3way/log/%x_%A_%a.err
#SBATCH --array=1-20

## Emit one bgzipped VCF per block_row with the 3-way dosage format
## DS1:DS2:DS3:ANC1:ANC2:ANC3:DSALL that FELIX step2 consumes.
##   MODE=common  sbatch 09_make_vcf.sh
##   MODE=lowfreq sbatch 09_make_vcf.sh

set -euo pipefail
MODE=${MODE:?MODE (common|lowfreq) must be exported}
module load singularity

BASE=/data/wzhougroup/lhu/saige_tractor/simulation/3way
RTOOLS=/data/wzhougroup/lhu/tools/rtools_latest.sif
SING="singularity exec --bind /data/wzhougroup/lhu:/data/wzhougroup/lhu --home /data/wzhougroup/lhu"

$SING $RTOOLS Rscript ${BASE}/scripts/R/make_vcf.R "${SLURM_ARRAY_TASK_ID}" "${MODE}"
