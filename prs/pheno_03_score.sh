#!/bin/bash
# Score PRS in UK Biobank with PLINK2 and sum across chromosomes. Env:
# PHENO_ID, ANCESTRY.

#SBATCH --job-name=prs_score
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=/data/wzhougroup/lhu/saige_tractor/prs_pipeline/logs/prs_score_%j.out
#SBATCH --error=/data/wzhougroup/lhu/saige_tractor/prs_pipeline/logs/prs_score_%j.err

set -euo pipefail

: "${PHENO_ID:?Set PHENO_ID}"
: "${ANCESTRY:?Set ANCESTRY (AFR or EUR)}"

BASE_DIR="/data/wzhougroup/lhu/saige_tractor/prs_pipeline"
SHARED_DIR="${BASE_DIR}/shared"
OUT_DIR="${BASE_DIR}/${PHENO_ID}/${ANCESTRY}"
SCORE_DIR="${OUT_DIR}/scores"

PLINK2_SIF="/data/wzhougroup/lhu/tools/plink2_alpha2.3_jan2020.sif"
PLINK2_BIND="/data/wzhougroup/lhu:/data/wzhougroup/lhu"

mkdir -p "${SCORE_DIR}"
module load Apptainer/1.4.2-1.el9

cd "${SCORE_DIR}"

echo "==========================================="
echo "Step 3: PRS scoring — PHENO=${PHENO_ID}, ANCESTRY=${ANCESTRY}"
echo "==========================================="

if [[ "${ANCESTRY}" == "AFR" ]]; then
    PGEN_DIR="${SHARED_DIR}/afr_imputed_data"

    # Sanity: per-chr score files exist
    for chr in {1..22}; do
        for MODE in local global; do
            f="${SCORE_DIR}/${MODE}_chr${chr}.txt"
            if [ ! -s "${f}" ]; then
                echo "ERROR: ${f} missing or empty — run step 2.5 first." >&2
                exit 1
            fi
        done
    done

    for MODE in local global; do
        echo ""
        echo ">>> ${MODE} (AFR, per-chr loop)"

        for chr in {1..22}; do
            OUT_PFX="${SCORE_DIR}/${MODE}_chr${chr}"
            PFILE="${PGEN_DIR}/ukb_afr_chr${chr}"
            SCORE_FILE="${SCORE_DIR}/${MODE}_chr${chr}.txt"

            # Only skip if the existing sscore is NEWER than its score-file input —
            # otherwise the cache is stale (e.g., posteriors were regenerated).
            if [ -f "${OUT_PFX}.sscore" ] && [ "${OUT_PFX}.sscore" -nt "${SCORE_FILE}" ]; then
                echo "  chr${chr} already done (sscore newer than score input), skipping"
                continue
            fi
            if [ -f "${OUT_PFX}.sscore" ]; then
                echo "  chr${chr} stale sscore, regenerating"
                rm -f "${OUT_PFX}.sscore"
            fi

            apptainer exec --bind "${PLINK2_BIND}" "${PLINK2_SIF}" \
                plink2 \
                  --pfile "${PFILE}" \
                  --rm-dup exclude-mismatch \
                  --score "${SCORE_FILE}" 2 4 6 cols=+scoresums \
                  --out "${OUT_PFX}" \
                  --memory 14000 \
                2>&1 | tee "${OUT_PFX}.plinklog"

            if [ ! -f "${OUT_PFX}.sscore" ]; then
                echo "ERROR: ${MODE} chr${chr} produced no .sscore" >&2
                exit 1
            fi
            echo "  chr${chr} done"
        done

        # Combine 22 per-chr sscore files into one
        echo "  Combining ${MODE} across chromosomes..."
        head -1 "${MODE}_chr1.sscore" > "${MODE}_prs.sscore"
        awk -F'\t' '
        FNR == 1 { next }
        {
            key = $1 SUBSEP $2
            score[key] += $NF
            if (!(key in seen)) {
                seen[key] = 1
                order[++n] = key
                nf = NF
                for (j = 1; j < nf; j++) {
                    prefix[key] = (j==1 ? $j : prefix[key] "\t" $j)
                }
            }
        }
        END {
            for (i = 1; i <= n; i++) {
                k = order[i]
                printf "%s\t%.6g\n", prefix[k], score[k]
            }
        }' "${MODE}"_chr*.sscore >> "${MODE}_prs.sscore"

        N=$(tail -n +2 "${MODE}_prs.sscore" | wc -l)
        echo "  ${MODE}_prs.sscore: ${N} samples"
    done

elif [[ "${ANCESTRY}" == "EUR" ]]; then
    PFILE="${SHARED_DIR}/training-white_gate"   # single all-chromosome pfile

    # Sanity: intersected (all-chr) files exist
    for MODE in local global; do
        f="${SCORE_DIR}/${MODE}_intersect.txt"
        if [ ! -s "${f}" ]; then
            echo "ERROR: ${f} missing or empty — run step 2.5 first." >&2
            exit 1
        fi
    done

    for MODE in local global; do
        echo ""
        echo ">>> ${MODE} (EUR, single pfile)"

        OUT_PFX="${SCORE_DIR}/${MODE}_prs"
        SCORE_FILE="${SCORE_DIR}/${MODE}_intersect.txt"

        if [ -f "${OUT_PFX}.sscore" ] && [ "${OUT_PFX}.sscore" -nt "${SCORE_FILE}" ]; then
            echo "  ${MODE} already done (sscore newer than score input), skipping"
            continue
        fi
        if [ -f "${OUT_PFX}.sscore" ]; then
            echo "  ${MODE} stale sscore, regenerating"
            rm -f "${OUT_PFX}.sscore"
        fi

        apptainer exec --bind "${PLINK2_BIND}" "${PLINK2_SIF}" \
            plink2 \
              --pfile "${PFILE}" \
              --rm-dup exclude-mismatch \
              --score "${SCORE_FILE}" 2 4 6 cols=+scoresums \
              --out "${OUT_PFX}" \
              --memory 14000 \
            2>&1 | tee "${OUT_PFX}.plinklog"

        if [ ! -f "${OUT_PFX}.sscore" ]; then
            echo "ERROR: ${MODE} (EUR) produced no .sscore" >&2
            exit 1
        fi
        N=$(tail -n +2 "${OUT_PFX}.sscore" | wc -l)
        echo "  ${OUT_PFX}.sscore: ${N} samples"
    done

else
    echo "ERROR: ANCESTRY must be AFR or EUR (got: ${ANCESTRY})" >&2
    exit 1
fi

echo "==========================================="
echo "DONE: PHENO=${PHENO_ID}, ANCESTRY=${ANCESTRY}"
echo "  Output: ${SCORE_DIR}/{local,global}_prs.sscore"
echo "==========================================="
