#!/bin/bash
# Orchestrate the PRS pipeline for one phenotype (SLURM chain):
# prep sumstats -> PRS-CS -> intersect -> score -> evaluate R2, per ancestry.

set -euo pipefail

if [ $# -ne 6 ]; then
    cat <<EOF >&2
Usage: $0 <PHENO_ID> <PHENO_COL> <N_AFR_GLOBAL> <N_AFR_LOCAL> <N_EUR_GLOBAL> <N_EUR_LOCAL>

Arguments:
  PHENO_ID       Phenotype id (e.g., 3027114)
  PHENO_COL      Phenotype column name in UKB pheno tsv (e.g., TotalCholesterol)
  N_AFR_GLOBAL   Effective N from AFR global GWAS
  N_AFR_LOCAL    Effective N from AFR local-ancestry-aware GWAS
  N_EUR_GLOBAL   Effective N from EUR global GWAS
  N_EUR_LOCAL    Effective N from EUR local-ancestry-aware GWAS
EOF
    exit 1
fi

PHENO_ID=$1
PHENO_COL=$2
N_AFR_GLOBAL=$3
N_AFR_LOCAL=$4
N_EUR_GLOBAL=$5
N_EUR_LOCAL=$6

BASE_DIR="/data/wzhougroup/lhu/saige_tractor/prs_pipeline"
SCRIPT_DIR="${BASE_DIR}"   # pheno_*.sh live at top level; R scripts live in ${BASE_DIR}/scripts/

mkdir -p "${BASE_DIR}/logs"
mkdir -p "${BASE_DIR}/${PHENO_ID}/AFR/logs" "${BASE_DIR}/${PHENO_ID}/EUR/logs"
mkdir -p "${BASE_DIR}/${PHENO_ID}/results"

echo "==========================================="
echo "Phenotype: ${PHENO_ID} (${PHENO_COL})"
echo "  N AFR: global=${N_AFR_GLOBAL}, local=${N_AFR_LOCAL}"
echo "  N EUR: global=${N_EUR_GLOBAL}, local=${N_EUR_LOCAL}"
echo "==========================================="

# ---- Step 1: prep sumstats (per ancestry) — fast, interactive ----
echo ""
echo ">>> Step 1: prep sumstats (AFR + EUR, interactive)"
for ANC in AFR EUR; do
    ANCESTRY=${ANC} PHENO_ID=${PHENO_ID} \
        bash "${SCRIPT_DIR}/pheno_01_prep_sumstats.sh"
done

# ---- Step 2: PRS-CS array ----
echo ""
echo ">>> Step 2: PRS-CS"
JID_AFR_PRSCS=$(sbatch --parsable \
    --export=ALL,PHENO_ID=${PHENO_ID},ANCESTRY=AFR,N_GLOBAL=${N_AFR_GLOBAL},N_LOCAL=${N_AFR_LOCAL} \
    "${SCRIPT_DIR}/pheno_02_prscs_apptainer.sh")
echo "  AFR PRS-CS array: ${JID_AFR_PRSCS}"

JID_EUR_PRSCS=$(sbatch --parsable \
    --export=ALL,PHENO_ID=${PHENO_ID},ANCESTRY=EUR,N_GLOBAL=${N_EUR_GLOBAL},N_LOCAL=${N_EUR_LOCAL} \
    "${SCRIPT_DIR}/pheno_02_prscs_apptainer.sh")
echo "  EUR PRS-CS array: ${JID_EUR_PRSCS}"

# ---- Step 2.5: intersect (per ancestry, depends on step 2) ----
echo ""
echo ">>> Step 2.5: intersect"
JID_AFR_INT=$(sbatch --parsable \
    --dependency=afterok:${JID_AFR_PRSCS} \
    --export=ALL,PHENO_ID=${PHENO_ID},ANCESTRY=AFR \
    "${SCRIPT_DIR}/pheno_02p5_intersect_prscs.sh")
echo "  AFR intersect: ${JID_AFR_INT}"

JID_EUR_INT=$(sbatch --parsable \
    --dependency=afterok:${JID_EUR_PRSCS} \
    --export=ALL,PHENO_ID=${PHENO_ID},ANCESTRY=EUR \
    "${SCRIPT_DIR}/pheno_02p5_intersect_prscs.sh")
echo "  EUR intersect: ${JID_EUR_INT}"

# ---- Step 3: score ----
echo ""
echo ">>> Step 3: score"
JID_AFR_SCORE=$(sbatch --parsable \
    --dependency=afterok:${JID_AFR_INT} \
    --export=ALL,PHENO_ID=${PHENO_ID},ANCESTRY=AFR \
    "${SCRIPT_DIR}/pheno_03_score.sh")
echo "  AFR score: ${JID_AFR_SCORE}"

JID_EUR_SCORE=$(sbatch --parsable \
    --dependency=afterok:${JID_EUR_INT} \
    --export=ALL,PHENO_ID=${PHENO_ID},ANCESTRY=EUR \
    "${SCRIPT_DIR}/pheno_03_score.sh")
echo "  EUR score: ${JID_EUR_SCORE}"

# ---- Step 4: evaluate R² (depends on both step-3 jobs) ----
echo ""
echo ">>> Step 4: evaluate R²"
JID_EVAL=$(sbatch --parsable \
    --dependency=afterok:${JID_AFR_SCORE}:${JID_EUR_SCORE} \
    --export=ALL,PHENO_ID=${PHENO_ID},PHENO_COL=${PHENO_COL} \
    "${SCRIPT_DIR}/pheno_04_evaluate.sh")
echo "  Eval: ${JID_EVAL}"

echo ""
echo "==========================================="
echo "All jobs submitted for ${PHENO_ID}."
echo ""
echo "Job IDs:"
echo "  Step 2 PRS-CS:     AFR=${JID_AFR_PRSCS}  EUR=${JID_EUR_PRSCS}"
echo "  Step 2.5 intersect: AFR=${JID_AFR_INT}   EUR=${JID_EUR_INT}"
echo "  Step 3 score:      AFR=${JID_AFR_SCORE}  EUR=${JID_EUR_SCORE}"
echo "  Step 4 evaluate:   ${JID_EVAL}"
echo ""
echo "Monitor:    squeue -u \$USER"
echo "Cancel all: scancel ${JID_AFR_PRSCS} ${JID_EUR_PRSCS} ${JID_AFR_INT} ${JID_EUR_INT} ${JID_AFR_SCORE} ${JID_EUR_SCORE} ${JID_EVAL}"
echo "Final result: ${BASE_DIR}/${PHENO_ID}/results/r2_4way.tsv"
echo "==========================================="
