#!/bin/bash
# pt_predict.sh — one P+T PRS: clump a predictor's GWAS on an LD panel, then score the
# clumped loci on a validation-ancestry target. plink1.9 for both (validated dialect).
#   bash pt_predict.sh <TRAIT> <VAL> <PRED> <LDREF>
#     TRAIT : 3024929 | BMI | height | ...
#     VAL   : EUR AFR EAS CSA        (validation target)
#     PRED  : meta | felix | matched_<VAL> | felix_tract_<VAL>   (predictor_inputs/<TRAIT>/<PRED>.txt)
#     LDREF : eur10k | prop10k | AFR | EAS | CSA                 (clumping panel; a keeps/<LDREF>.keep)
# Genotypes come from two beds only: ldset (all LD panels + AFR/EAS/CSA targets, subset by
# --keep) and eurval (the EUR validation target).
source /humgen/atgu1/fin/lhu/projects/saige_tractor/prs/scripts/config.sh
TRAIT="$1"; VAL="$2"; PRED="$3"; LDREF="$4"
SUM="${BASE}/predictor_inputs/${TRAIT}/${PRED}.txt"      # SNP A1 A2 BETA P (SNP=CHR:POS)
OUT="${BASE}/results/${TRAIT}/${VAL}"; mkdir -p "${OUT}"
TAG="${PRED}__${LDREF}_${CLUMP_P1}"
[ -s "${SUM}" ] || { echo "  [skip] no predictor file ${SUM}"; exit 0; }

LDKEEP="${BASE}/keeps/${LDREF}.keep"                     # LD panel = ldset filtered to this keep
if [ "${VAL}" = "EUR" ]; then TGTPFX="${BASE}/geno_bed/eurval/chr"; TGTKEEP=""    # eurval already IS the EUR target
else TGTPFX="${BASE}/geno_bed/ldset/chr"; TGTKEEP="${BASE}/keeps/${VAL}.keep"; fi

for c in $(seq 1 22); do
  LD="${BASE}/geno_bed/ldset/chr${c}"; TGT="${TGTPFX}${c}"
  [ -s "${LD}.bed" ] && [ -s "${TGT}.bed" ] || { echo "  chr${c}: missing bed, skip"; continue; }
  # 1) clump predictor GWAS on the LD panel (ldset --keep <LDREF>) -> index SNPs at p<=CLUMP_P1
  "${PLINK19}" --bfile "${LD}" --keep "${LDKEEP}" --clump "${SUM}" \
    --clump-snp-field SNP --clump-field P \
    --clump-p1 "${CLUMP_P1}" --clump-p2 1 --clump-r2 "${R2}" --clump-kb "${KB}" \
    --out "${OUT}/${TAG}_c${c}" 2>&1 | tail -1 || true
  [ -s "${OUT}/${TAG}_c${c}.clumped" ] || { echo "  chr${c}: no clumps"; continue; }
  awk 'NR==1{for(i=1;i<=NF;i++) if($i=="SNP")k=i; next} (NF && $k!=""){print $k}' \
    "${OUT}/${TAG}_c${c}.clumped" > "${OUT}/${TAG}_c${c}.idx"
  # 2) score those loci on the validation target
  "${PLINK19}" --bfile "${TGT}" ${TGTKEEP:+--keep "${TGTKEEP}"} --extract "${OUT}/${TAG}_c${c}.idx" \
    --score "${SUM}" 1 2 4 header sum --out "${OUT}/${TAG}_c${c}" 2>&1 | tail -1 || true
done

# 3) sum SCORESUM across chromosomes -> one sscore
fs=$(ls "${OUT}/${TAG}_c"*.profile 2>/dev/null || true)
if [ -z "${fs}" ]; then echo "  [${TAG}] no scored loci (no genome-wide hits at ${CLUMP_P1})"; exit 0; fi
awk 'FNR==1{ii=sc=0; for(i=1;i<=NF;i++){if($i=="IID")ii=i; if($i=="SCORESUM")sc=i}; next}
     {s[$ii]+=$sc; seen[$ii]=1}
     END{print "IID\tSCORE_SUM"; for(k in seen) printf "%s\t%.6g\n", k, s[k]}' ${fs} \
  > "${OUT}/${TAG}.sscore"
echo "  -> ${OUT}/${TAG}.sscore"
rm -f "${OUT}/${TAG}_c"*.profile "${OUT}/${TAG}_c"*.clumped "${OUT}/${TAG}_c"*.idx "${OUT}/${TAG}_c"*.nosex "${OUT}/${TAG}_c"*.log
