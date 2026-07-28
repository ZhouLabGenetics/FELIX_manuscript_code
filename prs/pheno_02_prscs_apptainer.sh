#!/bin/bash
# PRS-CS posterior effect sizes (1000G EUR LD panel, default shrinkage), per
# chromosome array job. Env: PHENO_ID, ANCESTRY, N_GLOBAL, N_LOCAL.

#SBATCH --job-name=prscs
#SBATCH --partition=normal
#SBATCH --array=1-44
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=20G
#SBATCH --time=10:00:00
#SBATCH --output=/data/wzhougroup/lhu/saige_tractor/prs_pipeline/logs/prscs_%A_%a.out
#SBATCH --error=/data/wzhougroup/lhu/saige_tractor/prs_pipeline/logs/prscs_%A_%a.err

set -euo pipefail

: "${PHENO_ID:?Set PHENO_ID env variable}"
: "${ANCESTRY:?Set ANCESTRY (AFR or EUR)}"
: "${N_GLOBAL:?Set N_GLOBAL}"
: "${N_LOCAL:?Set N_LOCAL}"

BASE_DIR="/data/wzhougroup/lhu/saige_tractor/prs_pipeline"
SHARED_DIR="${BASE_DIR}/shared"
OUT_DIR="${BASE_DIR}/${PHENO_ID}/${ANCESTRY}"

PRSCS_SCRIPT="/data/wzhougroup/lhu/tools/PRScs/PRScs.py"
APPTAINER_SIF="/data/wzhougroup/lhu/tools/python3.sif"
BIM_PREFIX="${SHARED_DIR}/ukb_imp_hm3_all"

if [[ "${ANCESTRY}" == "AFR" ]]; then
    PRSCS_REF_DIR="/data/wzhougroup/lhu/tools/ref_ld/ldblk_ukbb_afr"
elif [[ "${ANCESTRY}" == "EUR" ]]; then
    PRSCS_REF_DIR="/data/wzhougroup/lhu/tools/ref_ld/ldblk_1kg_eur"
else
    echo "ERROR: ANCESTRY must be AFR or EUR (got: ${ANCESTRY})" >&2
    exit 1
fi

module load Apptainer/1.4.2-1.el9

# Array layout: 1-22 = global, 23-44 = local
TASK_ID=${SLURM_ARRAY_TASK_ID}
if [ ${TASK_ID} -le 22 ]; then
    CHR=${TASK_ID}
    METHOD="global"
    SST_FILE="${OUT_DIR}/sumstats/global_prscs_input.txt"
    N_GWAS=${N_GLOBAL}
    OUT_PREFIX="${OUT_DIR}/prscs/global/global"
else
    CHR=$((TASK_ID - 22))
    METHOD="local"
    SST_FILE="${OUT_DIR}/sumstats/local_prscs_input.txt"
    N_GWAS=${N_LOCAL}
    OUT_PREFIX="${OUT_DIR}/prscs/local/local"
fi

echo "PRS-CS: ANC=${ANCESTRY} METHOD=${METHOD} CHR=${CHR} N=${N_GWAS}"
echo "  sst:  ${SST_FILE}"
echo "  ref:  ${PRSCS_REF_DIR}"
echo "  bim:  ${BIM_PREFIX}"

apptainer exec --home /data/wzhougroup/lhu \
  --bind /data/wzhougroup/lhu:/data/wzhougroup/lhu \
  "${APPTAINER_SIF}" \
  python "${PRSCS_SCRIPT}" \
    --ref_dir="${PRSCS_REF_DIR}" \
    --bim_prefix="${BIM_PREFIX}" \
    --sst_file="${SST_FILE}" \
    --n_gwas="${N_GWAS}" \
    --chrom="${CHR}" \
    --out_dir="${OUT_PREFIX}"

echo "Done: ANC=${ANCESTRY} METHOD=${METHOD} CHR=${CHR}"
