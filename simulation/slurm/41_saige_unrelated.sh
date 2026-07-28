#!/bin/bash
#SBATCH --job-name=3way_saige_unr
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=/data/wzhougroup/lhu/saige_tractor/simulation/3way/log/%x_%A_%a.out
#SBATCH --error=/data/wzhougroup/lhu/saige_tractor/simulation/3way/log/%x_%A_%a.err

## FELIX on the unrelated benchmark cohort.
## Step1 (sparse-GRM null + variance ratio) then step2 chr1 association on
## the VCF. Hybrid format is produced by 40_simu_unrelated.sh at
## ${UNR_DIR}/hybrid/ -- a separate step2 invocation can read it once we
## confirm the hybrid-input step2 syntax.
##
## Array layout: 3 traits x N_SEEDS seeds.
##   id 1..N_SEEDS                : quant
##   id N_SEEDS+1..2*N_SEEDS      : bin10
##   id 2*N_SEEDS+1..3*N_SEEDS    : bin01
##
## Submit (N_SEEDS=10 default => array 1-30 per cohort):
##   sbatch --array=1-30 --mem=16G --time=02:00:00 41_saige_unrelated.sh                  # 100% unrelated
##   UNREL_PCT=75  sbatch --array=1-30 --mem=16G --time=02:00:00 41_saige_unrelated.sh    # 75/25
##   UNREL_PCT=50  sbatch --array=1-30 --mem=16G --time=02:00:00 41_saige_unrelated.sh    # 50/50
##   UNREL_PCT=25  sbatch --array=1-30 --mem=16G --time=02:00:00 41_saige_unrelated.sh    # 25/75

set -euo pipefail
module load apptainer 2>/dev/null || module load singularity 2>/dev/null || true

BASE=/data/wzhougroup/lhu/saige_tractor/simulation/3way
SAIGET=${SAIGET:-/data/wzhougroup/lhu/tools/saigetractor_1.4.9-tractor-hybrid.1.sif}
APPT="apptainer exec --bind /data/wzhougroup/lhu:/data/wzhougroup/lhu --home /data/wzhougroup/lhu"

UNR_N=${UNR_N:-10000}
UNR_P=${UNR_P:-50000}
UNREL_PCT=${UNREL_PCT:-100}
UNREL_TAG=$(printf "%03d" "$UNREL_PCT")
UNR_DIR=${UNR_DIR:-${BASE}/data/unr_N${UNR_N}_P${UNR_P}_unrel${UNREL_TAG}}
N_SEEDS=${N_SEEDS:-10}

TASK=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID must be set (use sbatch --array=...)}
T_IDX=$(( (TASK - 1) / N_SEEDS ))
S_IDX=$(( (TASK - 1) % N_SEEDS + 1 ))

case "$T_IDX" in
  0) TRAIT=quant ; TTYPE=quantitative ;;
  1) TRAIT=bin10 ; TTYPE=binary       ;;
  2) TRAIT=bin01 ; TTYPE=binary       ;;
  *) echo "Bad T_IDX=$T_IDX (only 0..2 valid for 3 traits)" >&2; exit 1 ;;
esac

TAG=$(printf "%s_seed%02d" "$TRAIT" "$S_IDX")
echo "==== SAIGE unr: $TAG  (UNR_N=$UNR_N UNR_P=$UNR_P) ===="

PHENO=${UNR_DIR}/pheno/null/${TAG}.tsv
[[ -f "$PHENO" ]] || { echo "Pheno missing: $PHENO" >&2; exit 1; }

OUT_DIR=${UNR_DIR}/saige_out/null/${TAG}
mkdir -p "$OUT_DIR"

PLINK_PREFIX=${UNR_DIR}/plink/pruned
SPARSE_MTX=$(ls ${UNR_DIR}/plink/sparseGRM*.mtx 2>/dev/null | head -1 || true)
[[ -f "$SPARSE_MTX" ]] || { echo "Sparse GRM missing under ${UNR_DIR}/plink/" >&2; exit 1; }
SPARSE_IDS=${SPARSE_MTX}.sampleIDs.txt
VCF=${UNR_DIR}/vcf/chr1.vcf.gz
[[ -f "$VCF"        ]] || { echo "VCF missing: $VCF" >&2; exit 1; }

echo "  PLINK=$PLINK_PREFIX"
echo "  SPARSE_MTX=$SPARSE_MTX"
echo "  PHENO=$PHENO"
echo "  VCF=$VCF"
echo "  OUT_DIR=$OUT_DIR"

## ---- step1: null + variance ratio --------------------------------------
echo "==== step1 ===="
$APPT $SAIGET step1_fitNULLGLMM.R \
    --plinkFile=${PLINK_PREFIX} \
    --useSparseGRMtoFitNULL=TRUE \
    --sparseGRMFile=${SPARSE_MTX} \
    --sparseGRMSampleIDFile=${SPARSE_IDS} \
    --phenoFile=${PHENO} \
    --phenoCol=PHENO \
    --covarColList=AFR,NAT \
    --sampleIDColinphenoFile=IID \
    --traitType=${TTYPE} \
    --IsOverwriteVarianceRatioFile=TRUE \
    --isCateVarianceRatio=FALSE \
    --outputPrefix=${OUT_DIR}/step1

GMMAT=${OUT_DIR}/step1.rda
VARR=${OUT_DIR}/step1.varianceRatio.txt
[[ -s "$GMMAT" && -s "$VARR" ]] || { echo "step1 outputs missing" >&2; exit 1; }

## ---- step2: chr1 association (VCF input) -------------------------------
echo "==== step2 (chr1, VCF) ===="
$APPT $SAIGET step2_SPAtests.R \
    --vcfFile=${VCF} \
    --vcfFileIndex=${VCF}.csi \
    --vcfField=DS1 \
    --LOCO=FALSE \
    --AlleleOrder=ref-first \
    --SAIGEOutputFile=${OUT_DIR}/chr1.SAIGE.txt \
    --chrom=1 \
    --minMAF=0 \
    --minMAC=0.5 \
    --GMMATmodelFile=${GMMAT} \
    --varianceRatioFile=${VARR} \
    --is_admixed=TRUE \
    --number_of_ancestry=3

echo "==== Done.  $TAG ===="
