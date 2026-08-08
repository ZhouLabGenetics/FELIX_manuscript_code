#!/bin/bash
#SBATCH --job-name=3way_saige_bench_largeN
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=log/%x_%A_%a.out
#SBATCH --error=log/%x_%A_%a.err

## FELIX large-N benchmark on the bench_N<...>_P<...> data.
## Runs step1 (sparse-GRM null + variance ratio) then step2 chr1 with
## /usr/bin/time -v wrapping each invocation (matches 31_tm_bench_largeN.sh
## measurement methodology).
##
## Array layout: 15 tasks = 5 N's x 3 traits, indexed identically to 31_*.
##   id  1..3 : N=10000   id  4..6 : N=50000   id  7..9 : N=100000
##   id 10..12: N=150000  id 13..15: N=200000
##
## Submit:
##   sbatch --array=1-15 --mem=64G --time=08:00:00 32_saige_bench_largeN.sh

set -euo pipefail
module load singularity

BASE="${FELIX_SIM_BASE:?Set FELIX_SIM_BASE to the simulation directory before submitting}"
export FELIX_SIM_BASE
SAIGET=${SAIGET:-/data/wzhougroup/lhu/tools/FELIX_latest.sif}
SING="singularity exec --bind /data/wzhougroup/lhu:/data/wzhougroup/lhu --home /data/wzhougroup/lhu"

TASK=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID must be set (use sbatch --array=...)}
N_IDX=$(( (TASK - 1) / 3 ))
T_IDX=$(( (TASK - 1) % 3 ))

case "$N_IDX" in
  0) BENCH_N=10000  ;;
  1) BENCH_N=50000  ;;
  2) BENCH_N=100000 ;;
  3) BENCH_N=150000 ;;
  4) BENCH_N=200000 ;;
esac
case "$T_IDX" in
  0) TRAIT=quant ; TTYPE=quantitative ;;
  1) TRAIT=bin10 ; TTYPE=binary       ;;
  2) TRAIT=bin01 ; TTYPE=binary       ;;
esac
BENCH_P=${BENCH_P:-1000}
BENCH_DIR=${BASE}/data/bench_N${BENCH_N}_P${BENCH_P}

echo "==== SAIGE bench: N=${BENCH_N}  TRAIT=${TRAIT}  ===="
TAG=${TRAIT}_seed01
PHENO=${BENCH_DIR}/pheno/null/${TAG}.tsv
[[ -f "$PHENO" ]] || { echo "Pheno missing: $PHENO" >&2; exit 1; }

OUT_DIR=${BENCH_DIR}/saige_bench/null/${TAG}
mkdir -p "$OUT_DIR"
STEP1_BENCH=${OUT_DIR}/step1.bench.log
STEP2_BENCH=${OUT_DIR}/step2_chr1.bench.log

PLINK_PREFIX=${BENCH_DIR}/plink/pruned
SPARSE_MTX=${BENCH_DIR}/plink/sparseGRM_relatednessCutoff_0.125_${BENCH_P}_randomMarkersUsed.sparseGRM.mtx
if [[ ! -f "$SPARSE_MTX" ]]; then
    ## createSparseGRM names the file with the numRandomMarker actually used;
    ## fall back to whatever sparseGRM*.mtx exists if the exact name differs.
    SPARSE_MTX=$(ls ${BENCH_DIR}/plink/sparseGRM*.mtx 2>/dev/null | head -1 || true)
fi
[[ -f "$SPARSE_MTX" ]] || { echo "Sparse GRM missing under ${BENCH_DIR}/plink/" >&2; exit 1; }
SPARSE_IDS=${SPARSE_MTX}.sampleIDs.txt
VCF=${BENCH_DIR}/vcf/chr1.vcf.gz
[[ -f "$VCF"        ]] || { echo "VCF missing: $VCF" >&2; exit 1; }

echo "  PLINK=$PLINK_PREFIX"
echo "  SPARSE_MTX=$SPARSE_MTX"
echo "  PHENO=$PHENO"
echo "  OUT_DIR=$OUT_DIR"

## ---- step1: null + variance ratio --------------------------------------
echo "==== step1 ===="
/usr/bin/time -v -o "$STEP1_BENCH" \
    $SING $SAIGET step1_fitNULLGLMM.R \
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

## ---- step2: chr1 association ------------------------------------------
echo "==== step2 (chr1) ===="
/usr/bin/time -v -o "$STEP2_BENCH" \
    $SING $SAIGET step2_SPAtests.R \
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

echo "==== Done.  step1=${STEP1_BENCH}  step2=${STEP2_BENCH} ===="
