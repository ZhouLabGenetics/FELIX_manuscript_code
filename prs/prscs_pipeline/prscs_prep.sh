#!/bin/bash
# prscs_prep.sh — PREP for the (later) PRS-CS stage; safe to run now while the bgen
# extraction is still going. PRS-CS matches sumstats to its LD panel by rsID, but our
# predictor score files are keyed by CHR:POS. This maps them to rsID (via the PRS-CS HM3
# snpinfo) and writes PRS-CS-format sumstats "SNP A1 A2 BETA P" under prscs_inputs/.
source "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)/config.sh"
cd "${BASE}"
mkdir -p prscs_inputs

# CHR:POS -> rsID map from the PRS-CS HM3 snpinfo (cols SNP=rsID, CHR, BP)
MAP="hm3/chrpos2rsid.txt"
[ -s "${MAP}" ] || awk 'NR==1{for(i=1;i<=NF;i++){if($i=="SNP")s=i;if($i=="CHR")c=i;if($i=="BP")b=i} next}
                        (s&&c&&b){print $c":"$b"\t"$s}' "${SNPINFO}" > "${MAP}"
echo "chrpos->rsID map: $(wc -l < "${MAP}") HM3 variants"

for TRAIT in ${TRAITS}; do
  for pf in predictor_inputs/${TRAIT}/*.txt; do
    [ -e "${pf}" ] || continue
    nm=$(basename "${pf}" .txt); out="prscs_inputs/${TRAIT}/${nm}.txt"; mkdir -p "$(dirname "${out}")"
    awk 'NR==FNR{r[$1]=$2; next}
         FNR==1{print "SNP\tA1\tA2\tBETA\tP"; next}
         ($1 in r){print r[$1]"\t"$2"\t"$3"\t"$4"\t"$5}' "${MAP}" "${pf}" > "${out}"
  done
  echo "  ${TRAIT}: $(ls prscs_inputs/${TRAIT}/ 2>/dev/null | wc -l) predictor files -> rsID-keyed"
done
# SANITY: predictor CHR:POS (hg19, lifted) must match snpinfo (hg19). A high mapping rate
# confirms the builds align; near-0 means a build mismatch (predictors not lifted / snpinfo
# not hg19) -> stop and check.
IN="predictor_inputs/BMI/meta.txt"; OUT="prscs_inputs/BMI/meta.txt"
[ -s "${OUT}" ] && echo "SANITY BMI/meta: mapped $(($(wc -l < "${OUT}")-1)) / $(($(wc -l < "${IN}")-1)) variants to rsID (expect most)"
echo "PRS-CS inputs under ${BASE}/prscs_inputs/  (rsID-keyed)."

# --- rsID bim for PRScs --bim_prefix (target SNP list, from snpinfo; hg19) ---
awk 'NR==1{for(i=1;i<=NF;i++) n[$i]=i; next}
     {print $n["CHR"], $n["SNP"], 0, $n["BP"], $n["A1"], $n["A2"]}' OFS='\t' \
     "${SNPINFO}" > "${PRSCS_BIM}.bim"
echo "PRScs bim: $(wc -l < "${PRSCS_BIM}.bim") HM3 rsIDs -> ${PRSCS_BIM}.bim"

# --- combos to run: (trait target predictor); felix+meta always, matched_<target> if present ---
: > prscs_combos.txt
for TRAIT in ${TRAITS}; do
  for TGT in ${ANCS}; do
    for PRED in felix meta "matched_${TGT}"; do
      [ -s "prscs_inputs/${TRAIT}/${PRED}.txt" ] && echo "${TRAIT} ${TGT} ${PRED}" >> prscs_combos.txt
    done
  done
done
echo "prscs_combos.txt: $(wc -l < prscs_combos.txt) (trait,target,predictor) combos"
echo "NEXT: qsub -t 1-\$(( $(wc -l < prscs_combos.txt) * 22 )) 05_prscs.qsub   (n_gwas from ${NGWAS})"
