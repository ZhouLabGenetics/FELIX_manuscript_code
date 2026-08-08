#!/bin/bash
# run_pt.sh — all P+T PRS for ONE trait: for each validation ancestry, score the 3
# predictors. Multi-ancestry predictors (meta, felix) are clumped on BOTH eur10k and
# prop10k; the ancestry-matched predictor is clumped on its matched panel.
#   bash run_pt.sh <TRAIT>            # CLUMP_P1=5e-8 default; CLUMP_P1=5e-6 to switch
source "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)/config.sh"
TRAIT="$1"
echo "=================  ${TRAIT}  (thresh ${CLUMP_P1})  ================="
for VAL in ${ANCS}; do
  echo "--- validation ${VAL} ---"
  for PRED in meta felix; do                       # multi-ancestry: two LD panels
    for LDREF in eur10k prop10k; do
      bash "${SCRIPTS_DIR}/pt_predict.sh" "${TRAIT}" "${VAL}" "${PRED}" "${LDREF}"
    done
  done
  # single-ancestry predictors clump on the ancestry-matched panel; EUR uses the 10k EUR
  # panel (the full 451k EUR is only a scoring target, not an LD ref)
  LDM="${VAL}"; [ "${VAL}" = "EUR" ] && LDM="eur10k"
  bash "${SCRIPTS_DIR}/pt_predict.sh" "${TRAIT}" "${VAL}" "matched_${VAL}" "${LDM}"
  # supplementary global-vs-local (tract) story: FELIX ancestry-tract beta on matched panel
  bash "${SCRIPTS_DIR}/pt_predict.sh" "${TRAIT}" "${VAL}" "felix_tract_${VAL}" "${LDM}"
done
echo "done ${TRAIT} -> ${BASE}/results/${TRAIT}/<VAL>/*.sscore"
