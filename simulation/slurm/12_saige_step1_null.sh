#!/bin/bash
#SBATCH --job-name=3way_step1_null
#SBATCH --time=04:00:00
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=/data/wzhougroup/lhu/saige_tractor/simulation/3way/log/%x_%A_%a.out
#SBATCH --error=/data/wzhougroup/lhu/saige_tractor/simulation/3way/log/%x_%A_%a.err
#SBATCH --array=1-30

## FELIX step1 for NULL phenotypes.
## 3 trait types (quant, bin01, bin10) x 10 seeds = 30 pheno files.
## Each array task fits the null model for one pheno file and writes
##   data/<mode>/saige_out/null/<pheno_tag>/step1.rda
##   data/<mode>/saige_out/null/<pheno_tag>/step1.varianceRatio.txt
##
##   MODE=common  sbatch 12_saige_step1_null.sh
##   MODE=lowfreq sbatch 12_saige_step1_null.sh

set -euo pipefail
MODE=${MODE:?MODE (common|lowfreq) must be exported}
module load singularity

BASE=/data/wzhougroup/lhu/saige_tractor/simulation/3way
SAIGET=/data/wzhougroup/lhu/tools/saigetractor148.sif
SING="singularity exec --bind /data/wzhougroup/lhu:/data/wzhougroup/lhu --home /data/wzhougroup/lhu"

PHENO_DIR=${BASE}/data/${MODE}/pheno/null
mapfile -t PHENOS < <(ls ${PHENO_DIR}/*.tsv | sort)
PHENO=${PHENOS[$((SLURM_ARRAY_TASK_ID-1))]}
TAG=$(basename ${PHENO} .tsv)
TRAIT_PREFIX=${TAG%%_seed*}
case "$TRAIT_PREFIX" in
  quant) TRAIT=quantitative ;;
  bin01|bin10) TRAIT=binary ;;
  *) echo "Unknown trait prefix: $TRAIT_PREFIX" >&2; exit 1 ;;
esac

OUT_DIR=${BASE}/data/${MODE}/saige_out/null/${TAG}
mkdir -p ${OUT_DIR}

PLINK_PREFIX=${BASE}/data/${MODE}/plink/pruned
SPARSE_MTX=${BASE}/data/${MODE}/plink/sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx
SPARSE_IDS=${SPARSE_MTX}.sampleIDs.txt

$SING $SAIGET step1_fitNULLGLMM.R \
    --plinkFile=${PLINK_PREFIX} \
    --useSparseGRMtoFitNULL=TRUE \
    --sparseGRMFile=${SPARSE_MTX} \
    --sparseGRMSampleIDFile=${SPARSE_IDS} \
    --phenoFile=${PHENO} \
    --phenoCol=PHENO \
    --covarColList=AFR,NAT \
    --sampleIDColinphenoFile=IID \
    --traitType=${TRAIT} \
    --IsOverwriteVarianceRatioFile=TRUE \
    --isCateVarianceRatio=FALSE \
    --outputPrefix=${OUT_DIR}/step1
