#!/bin/bash
#SBATCH --job-name=3way_saige_bench
#SBATCH --time=06:00:00
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=/data/wzhougroup/lhu/saige_tractor/simulation/3way/log/%x_%j.out
#SBATCH --error=/data/wzhougroup/lhu/saige_tractor/simulation/3way/log/%x_%j.err

## FELIX head-to-head benchmark vs Tractor-Mix.
## Runs step1 (variance components) + step2 (chr1 only, matching the chr1
## hapdose footprint Tractor-Mix consumes) under /usr/bin/time -v. Writes
## to a separate saige_bench/ tree so the existing saige_out/ is untouched.
##
## Per-step bench logs (wall, peak RSS) are captured separately so we can
## quote both the amortised step1 cost and the per-chromosome step2 cost.
##
## Usable two ways:
##   (a) sbatch:
##       MODE=common TRAIT=quant SEED=1 sbatch 17c_saige_bench.sh
##       MODE=common TRAIT=quant SEED=1 KIND=alt SCEN=shared BETA=0.30 \
##           sbatch 17c_saige_bench.sh
##   (b) interactive (e.g. srun --pty bash, then):
##       MODE=common TRAIT=quant SEED=1 bash 17c_saige_bench.sh
##
## Env:
##   MODE   common | lowfreq          (required)
##   TRAIT  quant  | bin10 | bin01    (required)
##   SEED   1-10                      (required)
##   KIND   null (default) | alt
##   SCEN   shared | hetero           (required if KIND=alt)
##   BETA   numeric                   (required if KIND=alt)

set -euo pipefail
MODE=${MODE:?MODE (common|lowfreq) must be exported}
TRAIT=${TRAIT:?TRAIT (quant|bin10|bin01) must be exported}
SEED=${SEED:?SEED (1-10) must be exported}
KIND=${KIND:-null}
module load singularity 2>/dev/null || true

BASE=/data/wzhougroup/lhu/saige_tractor/simulation/3way
SAIGET=${SAIGET:-/data/wzhougroup/lhu/tools/saigetractor148.sif}
SING="singularity exec --bind /data/wzhougroup/lhu:/data/wzhougroup/lhu --home /data/wzhougroup/lhu"

case "$TRAIT" in
  quant)        TTYPE=quantitative ;;
  bin10|bin01)  TTYPE=binary ;;
  *) echo "Unknown TRAIT: $TRAIT" >&2; exit 1 ;;
esac

if [[ "$KIND" == "null" ]]; then
    TAG=${TRAIT}_seed$(printf "%02d" "$SEED")
    PHENO=${BASE}/data/${MODE}/pheno/null/${TAG}.tsv
elif [[ "$KIND" == "alt" ]]; then
    SCEN=${SCEN:?SCEN must be exported when KIND=alt}
    BETA=${BETA:?BETA must be exported when KIND=alt}
    BETA_TAG=$(awk -v b="$BETA" 'BEGIN{ printf "beta%03d", int(b*100 + 0.5) }')
    TAG=${TRAIT}_${SCEN}_${BETA_TAG}_seed$(printf "%02d" "$SEED")
    PHENO=${BASE}/data/${MODE}/pheno/alt/${TAG}.tsv
else
    echo "Unknown KIND: $KIND (expected null or alt)" >&2; exit 1
fi
[[ -f "$PHENO" ]] || { echo "Pheno not found: $PHENO" >&2; exit 1; }

OUT_DIR=${BASE}/data/${MODE}/saige_bench/${KIND}/${TAG}
mkdir -p "$OUT_DIR"
STEP1_BENCH=${OUT_DIR}/step1.bench.log
STEP2_BENCH=${OUT_DIR}/step2_chr1.bench.log

PLINK_PREFIX=${BASE}/data/${MODE}/plink/pruned
SPARSE_MTX=${BASE}/data/${MODE}/plink/sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx
SPARSE_IDS=${SPARSE_MTX}.sampleIDs.txt
VCF=${BASE}/data/${MODE}/vcf/chr1.vcf.gz

[[ -f "$VCF"        ]] || { echo "VCF missing: $VCF" >&2; exit 1; }
[[ -f "$SPARSE_MTX" ]] || { echo "Sparse GRM missing: $SPARSE_MTX" >&2; exit 1; }

echo "==> KIND=${KIND}  MODE=${MODE}  TRAIT=${TRAIT}  SEED=${SEED}  TAG=${TAG}"
echo "    PHENO=${PHENO}"
echo "    OUT_DIR=${OUT_DIR}"

## ---- step1 (null model + variance ratio) -------------------------------
## /usr/bin/time -v wraps the singularity invocation on the HOST, so it
## captures the full container lifetime + child R process.
echo "==> step1 (variance components)"
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

## ---- step2 chr1 (association test) -------------------------------------
echo "==> step2 (chr1 only, matching Tractor-Mix footprint)"
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

echo "==> Done: ${TAG}"
echo "    step1 bench: ${STEP1_BENCH}"
echo "    step2 bench: ${STEP2_BENCH}"
echo "    step2 out:   ${OUT_DIR}/chr1.SAIGE.txt"
