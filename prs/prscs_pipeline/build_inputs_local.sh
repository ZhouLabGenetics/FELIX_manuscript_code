#!/bin/bash
# build_inputs_local.sh — build all predictor score files LOCALLY, then rsync the
# small outputs to the new cluster (instead of the ~46 GB of raw sumstats).
set -euo pipefail
ABA=/Users/lhu/Desktop/5pop_st/aba_sumstat
SAIGE=/Users/lhu/Desktop/5pop_st/saige_sumstat
OUT=/Users/lhu/Desktop/5pop_st/prs/predictor_inputs
HERE="$(cd "$(dirname "$0")" && pwd)"
TRAITS="3006923 3007070 3009744 3013721 3022192 3024929 3027114 3028288 3035995 BMI height"
for t in ${TRAITS}; do
  python3 "${HERE}/make_predictor_inputs.py" "${t}" "${ABA}" "${SAIGE}" "${OUT}"
done
echo
echo "Built under ${OUT}. Transfer with:"
echo "  rsync -avz ${OUT}/ lhu@<newcluster>:/humgen/atgu1/fin/lhu/projects/saige_tractor/prs/predictor_inputs/"
