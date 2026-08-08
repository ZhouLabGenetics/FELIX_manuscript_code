#!/bin/bash
#SBATCH --job-name=3way_step2_null
#SBATCH --time=06:00:00
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=log/%x_%A_%a.out
#SBATCH --error=log/%x_%A_%a.err
#SBATCH --array=1-30

## FELIX step2 for NULL phenotypes.
## For each of the 30 null pheno files, loop through block_row VCFs (chr1..20)
## and emit per-chromosome association files used for Type-I error aggregation.
##
##   MODE=common  sbatch 13_saige_step2_null.sh
##   MODE=lowfreq sbatch 13_saige_step2_null.sh

set -euo pipefail
MODE=${MODE:?MODE (common|lowfreq) must be exported}
module load singularity

BASE="${FELIX_SIM_BASE:?Set FELIX_SIM_BASE to the simulation directory before submitting}"
export FELIX_SIM_BASE
SAIGET=/data/wzhougroup/lhu/tools/FELIX_latest.sif
SING="singularity exec --bind /data/wzhougroup/lhu:/data/wzhougroup/lhu --home /data/wzhougroup/lhu"

PHENO_DIR=${BASE}/data/${MODE}/pheno/null
mapfile -t PHENOS < <(ls ${PHENO_DIR}/*.tsv | sort)
PHENO=${PHENOS[$((SLURM_ARRAY_TASK_ID-1))]}
TAG=$(basename ${PHENO} .tsv)

OUT_DIR=${BASE}/data/${MODE}/saige_out/null/${TAG}
GMMAT=${OUT_DIR}/step1.rda
VARR=${OUT_DIR}/step1.varianceRatio.txt

VCF_DIR=${BASE}/data/${MODE}/vcf

for CHR in $(seq 1 20); do
    VCF=${VCF_DIR}/chr${CHR}.vcf.gz
    if [[ ! -f ${VCF} ]]; then
        echo "Missing ${VCF}, skipping" >&2
        continue
    fi

    $SING $SAIGET step2_SPAtests.R \
        --vcfFile=${VCF} \
        --vcfFileIndex=${VCF}.csi \
        --vcfField=DS1 \
        --LOCO=FALSE \
        --AlleleOrder=ref-first \
        --SAIGEOutputFile=${OUT_DIR}/chr${CHR}.SAIGE.txt \
        --chrom=${CHR} \
        --minMAF=0 \
        --minMAC=0.5 \
        --GMMATmodelFile=${GMMAT} \
        --varianceRatioFile=${VARR} \
        --is_admixed=TRUE \
        --number_of_ancestry=3
done
