#!/bin/bash
# Intersect local and global PRS-CS posteriors to a common variant set, split
# per chromosome. Env: PHENO_ID, ANCESTRY.

#SBATCH --job-name=intersect_prscs
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=/data/wzhougroup/lhu/saige_tractor/prs_pipeline/logs/intersect_%j.out
#SBATCH --error=/data/wzhougroup/lhu/saige_tractor/prs_pipeline/logs/intersect_%j.err

set -euo pipefail

: "${PHENO_ID:?Set PHENO_ID}"
: "${ANCESTRY:?Set ANCESTRY (AFR or EUR)}"

BASE_DIR="/data/wzhougroup/lhu/saige_tractor/prs_pipeline"
OUT_DIR="${BASE_DIR}/${PHENO_ID}/${ANCESTRY}"
SCORE_DIR="${OUT_DIR}/scores"

mkdir -p "${SCORE_DIR}"
cd "${SCORE_DIR}"

echo "==========================================="
echo "Step 2.5: intersect PRS-CS posteriors"
echo "  PHENO_ID=${PHENO_ID}  ANCESTRY=${ANCESTRY}"
echo "==========================================="

# A: concatenate per-chr posteriors for each mode
for MODE in local global; do
    SRC_DIR="${OUT_DIR}/prscs/${MODE}"
    OUT_FILE="${SCORE_DIR}/${MODE}_all_chr.txt"

    if [ ! -d "${SRC_DIR}" ]; then
        echo "ERROR: ${SRC_DIR} not found — did step 2 finish for ${MODE}?" >&2
        exit 1
    fi

    echo ">>> Concatenating ${MODE} posteriors..."
    cat "${SRC_DIR}"/${MODE}_pst_eff_a1_b0.5_phi*_chr*.txt > "${OUT_FILE}"
    N=$(wc -l < "${OUT_FILE}")
    echo "    ${MODE}: ${N} SNPs in ${OUT_FILE}"
    if [ "${N}" -eq 0 ]; then
        echo "ERROR: ${MODE} posterior file empty." >&2
        exit 1
    fi
done

# B: intersect on SNP ID
echo ">>> Intersecting on SNP ID..."
awk '{print $2}' local_all_chr.txt  | sort -u > local_snps.txt
awk '{print $2}' global_all_chr.txt | sort -u > global_snps.txt
comm -12 local_snps.txt global_snps.txt > common_snps.txt

N_LOCAL=$(wc -l < local_snps.txt)
N_GLOBAL=$(wc -l < global_snps.txt)
N_COMMON=$(wc -l < common_snps.txt)
echo "    local-only:  ${N_LOCAL}"
echo "    global-only: ${N_GLOBAL}"
echo "    intersect:   ${N_COMMON}"

if [ "${N_COMMON}" -eq 0 ]; then
    echo "ERROR: no SNPs in common." >&2
    exit 1
fi

# C: filter both posterior files to common SNP set
echo ">>> Writing intersected posterior files..."
for MODE in local global; do
    awk 'NR==FNR{a[$1]=1; next} ($2 in a)' \
        common_snps.txt "${MODE}_all_chr.txt" \
        > "${MODE}_intersect.txt"
    N=$(wc -l < "${MODE}_intersect.txt")
    echo "    ${MODE}_intersect.txt: ${N} SNPs"
done

# Sanity: identical SNP sets
diff <(awk '{print $2}' local_intersect.txt  | sort) \
     <(awk '{print $2}' global_intersect.txt | sort) > /tmp/intersect_diff_$$ || true
if [ -s /tmp/intersect_diff_$$ ]; then
    echo "WARNING: intersected files do not have identical SNP sets:"
    head /tmp/intersect_diff_$$
fi
rm -f /tmp/intersect_diff_$$

# D: split per-chr (always produced; AFR scoring needs it, EUR scoring ignores)
echo ">>> Splitting intersected files per chromosome..."
for MODE in local global; do
    for chr in {1..22}; do
        awk -v c=$chr '$1==c' "${MODE}_intersect.txt" > "${MODE}_chr${chr}.txt"
    done
    echo "    ${MODE}: chr1..chr22 split files written"
done

rm -f local_snps.txt global_snps.txt

echo "==========================================="
echo "DONE: PHENO=${PHENO_ID}, ANCESTRY=${ANCESTRY}"
echo "==========================================="
