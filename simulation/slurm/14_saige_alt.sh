#!/bin/bash
#SBATCH --job-name=3way_alt
#SBATCH --time=05:00:00
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=/data/wzhougroup/lhu/saige_tractor/simulation/3way/log/%x_%A_%a.out
#SBATCH --error=/data/wzhougroup/lhu/saige_tractor/simulation/3way/log/%x_%A_%a.err

## Combined FELIX step1 + step2 for ALT phenotypes (power).
##
## Reads pheno filenames from data/<MODE>/pheno/alt/manifest.tsv (column "file")
## so only manifest-listed phenos are run. Skips any pheno whose SAIGE output
## already exists (allows incremental reruns of new conditions only).
##
## BATCH_SIZE phenos per array task. With ~10 min per pheno, BATCH_SIZE=10
## gives ~100 min per task (well inside the 5-hour wall).
##
## START_BATCH_OFFSET allows chunked submission past the cluster's array
## limit (typically 1000). See submit_14_saige_alt.sh for the wrapper.
##
##   MODE=common  sbatch --array=1-N 14_saige_alt.sh
##   MODE=common  sbatch --array=1-1000 --export=ALL,MODE=common,START_BATCH_OFFSET=1000 14_saige_alt.sh

set -euo pipefail
MODE=${MODE:?MODE (common|lowfreq) must be exported}
module load singularity

BATCH_SIZE=${BATCH_SIZE:-10}
START_BATCH_OFFSET=${START_BATCH_OFFSET:-0}

BASE=/data/wzhougroup/lhu/saige_tractor/simulation/3way
SAIGET=/data/wzhougroup/lhu/tools/FELIX_latest.sif
SING="singularity exec --bind /data/wzhougroup/lhu:/data/wzhougroup/lhu --home /data/wzhougroup/lhu"

PHENO_DIR=${BASE}/data/${MODE}/pheno/alt

## Glob pheno files directly — manifest format doesn't matter for SAIGE.
## Filenames encode (trait, scen, beta, seed) so the aggregator can read
## TAGs back. causal.txt and manifest.tsv are excluded.
mapfile -t PHENOS < <(cd "${PHENO_DIR}" && ls *.tsv 2>/dev/null \
                      | grep -vE '^(manifest|causal)\.' | sort)
N_PHENOS=${#PHENOS[@]}
if (( N_PHENOS == 0 )); then
    echo "No pheno files found in ${PHENO_DIR}" >&2
    exit 1
fi

PLINK_PREFIX=${BASE}/data/${MODE}/plink/pruned
SPARSE_MTX=${BASE}/data/${MODE}/plink/sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx
SPARSE_IDS=${SPARSE_MTX}.sampleIDs.txt
CHR=1
VCF=${BASE}/data/${MODE}/vcf/chr${CHR}.vcf.gz

GLOBAL_BATCH=$(( SLURM_ARRAY_TASK_ID - 1 + START_BATCH_OFFSET ))
START=$(( GLOBAL_BATCH * BATCH_SIZE ))

if (( START >= N_PHENOS )); then
    echo "GLOBAL_BATCH=${GLOBAL_BATCH} START=${START} >= N_PHENOS=${N_PHENOS} — nothing to do."
    exit 0
fi

echo "MODE=${MODE}  N_PHENOS=${N_PHENOS}  BATCH_SIZE=${BATCH_SIZE}  GLOBAL_BATCH=${GLOBAL_BATCH}  START=${START}"

for OFFSET in $(seq 0 $((BATCH_SIZE - 1))); do
    IDX=$((START + OFFSET))
    if (( IDX >= N_PHENOS )); then break; fi

    PHENO_FILE=${PHENOS[$IDX]}
    PHENO=${PHENO_DIR}/${PHENO_FILE}
    TAG=$(basename "${PHENO}" .tsv)
    TRAIT_PREFIX=${TAG%%_*}
    case "$TRAIT_PREFIX" in
      quant)         TRAIT=quantitative ;;
      bin01|bin10)   TRAIT=binary ;;
      *) echo "Unknown trait prefix: $TRAIT_PREFIX (TAG=$TAG)" >&2; continue ;;
    esac

    OUT_DIR=${BASE}/data/${MODE}/saige_out/alt/${TAG}

    ## Skip when a non-empty SAIGE step2 output is already present.
    if [[ -s ${OUT_DIR}/chr${CHR}.SAIGE.txt ]]; then
        echo "=== [$((OFFSET+1))/$BATCH_SIZE] skip (already done): ${TAG} ==="
        continue
    fi

    mkdir -p "${OUT_DIR}"

    echo "=== [$((OFFSET+1))/$BATCH_SIZE] step1: ${TAG} (${TRAIT}) ==="
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

    echo "=== [$((OFFSET+1))/$BATCH_SIZE] step2: ${TAG} ==="
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
        --GMMATmodelFile=${OUT_DIR}/step1.rda \
        --varianceRatioFile=${OUT_DIR}/step1.varianceRatio.txt \
        --is_admixed=TRUE \
        --number_of_ancestry=3

    echo "=== [$((OFFSET+1))/$BATCH_SIZE] done: ${TAG} ==="
done
