#!/bin/bash
# Evaluate prediction R2 in UK Biobank (local vs global, per ancestry). Env:
# PHENO_ID, PHENO_COL.

#SBATCH --job-name=eval_r2
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=01:00:00
#SBATCH --output=/data/wzhougroup/lhu/saige_tractor/prs_pipeline/logs/eval_%j.out
#SBATCH --error=/data/wzhougroup/lhu/saige_tractor/prs_pipeline/logs/eval_%j.err

set -euo pipefail

: "${PHENO_ID:?Set PHENO_ID}"
: "${PHENO_COL:?Set PHENO_COL}"

BASE_DIR="/data/wzhougroup/lhu/saige_tractor/prs_pipeline"
OUT_DIR="${BASE_DIR}/${PHENO_ID}"

R_SIF="/data/wzhougroup/lhu/tools/rtools_latest.sif"
R_BIND="/data/wzhougroup/lhu:/data/wzhougroup/lhu"

UKB_PHENO="${UKB_PHENO:-${BASE_DIR}/ukb_eid_${PHENO_COL}.tsv}"

if [ ! -f "${UKB_PHENO}" ]; then
    echo "ERROR: phenotype file not found: ${UKB_PHENO}" >&2
    echo "  Set UKB_PHENO env var to override default path." >&2
    exit 1
fi

mkdir -p "${OUT_DIR}/results"
module load Apptainer/1.4.2-1.el9

echo "==========================================="
echo "Step 4: evaluate R² (4-way) — PHENO=${PHENO_ID} (${PHENO_COL})"
echo "==========================================="

# Sanity: all 4 sscore files exist
for f in \
    "${OUT_DIR}/AFR/scores/local_prs.sscore" \
    "${OUT_DIR}/AFR/scores/global_prs.sscore" \
    "${OUT_DIR}/EUR/scores/local_prs.sscore" \
    "${OUT_DIR}/EUR/scores/global_prs.sscore"
do
    if [ ! -s "${f}" ]; then
        echo "ERROR: missing or empty: ${f}" >&2
        exit 1
    fi
done

apptainer exec --bind "${R_BIND}" "${R_SIF}" \
    Rscript "${BASE_DIR}/scripts/evaluate_r2_4way.R" \
        "${OUT_DIR}/AFR/scores/local_prs.sscore" \
        "${OUT_DIR}/AFR/scores/global_prs.sscore" \
        "${OUT_DIR}/EUR/scores/local_prs.sscore" \
        "${OUT_DIR}/EUR/scores/global_prs.sscore" \
        "${UKB_PHENO}" \
        "${PHENO_COL}" \
        "${OUT_DIR}/results"

echo ""
echo "=== Results: ${OUT_DIR}/results/r2_4way.tsv ==="
cat "${OUT_DIR}/results/r2_4way.tsv"
