#!/bin/bash
#SBATCH --job-name=3way_phased_vcf
#SBATCH --time=04:00:00
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --output=log/%x_%A_%a.out
#SBATCH --error=log/%x_%A_%a.err

## Emit one bgzipped PHASED VCF per block_row with FORMAT=GT:HA
## (phased genotype + per-haplotype ancestry). Synthesizes phase from the
## existing Block_row<No>/{GAfr,GEur,GNat,LAafr,LAeur,LAnat}.tsv tables --
## no re-simulation required. The DS1/DS2/DS3/ANC1/ANC2/ANC3 fields of the
## unphased VCF can be reconstructed losslessly from GT+HA downstream.
##
## Default: chr1 only (matches the chr1 footprint used for Tractor-Mix /
## SAIGE benchmarks). To do all 20 chromosomes, pass --array=1-20 at submit:
##
##   MODE=common  sbatch                   09b_make_phased_vcf.sh    # chr1 only
##   MODE=common  sbatch --array=1-20      09b_make_phased_vcf.sh    # all chr
##   MODE=lowfreq sbatch                   09b_make_phased_vcf.sh

set -euo pipefail
MODE=${MODE:?MODE (common|lowfreq) must be exported}
module load singularity

BASE="${FELIX_SIM_BASE:?Set FELIX_SIM_BASE to the simulation directory before submitting}"
export FELIX_SIM_BASE
RTOOLS=/data/wzhougroup/lhu/tools/rtools_latest.sif
SING="singularity exec --bind /data/wzhougroup/lhu:/data/wzhougroup/lhu --home /data/wzhougroup/lhu"

## SLURM_ARRAY_TASK_ID is undefined when --array isn't passed -> default chr1.
NO=${SLURM_ARRAY_TASK_ID:-1}

$SING $RTOOLS Rscript ${BASE}/R/make_phased_vcf.R "${NO}" "${MODE}"
