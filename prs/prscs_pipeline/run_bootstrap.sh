#!/bin/bash
# run_bootstrap.sh — incremental R2 (covariate-adjusted) + paired bootstrap CI per trait x
# validation ancestry. Compares all predictors against FELIX (baseline felix__eur10k):
# a NEGATIVE "difference vs baseline" with CI below 0 means FELIX predicts better.
#   bash run_bootstrap.sh                    # all 11 traits
#   bash run_bootstrap.sh 3006923            # ONE trait (used by 04_bootstrap.qsub array)
#   B=1000 CLUMP_P1=5e-6 bash run_bootstrap.sh 3006923
# EUR is ~400k people, so B*fit dominates; parallelize with 04_bootstrap.qsub (one trait/task).
source "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)/config.sh"
# load R-4.1 (base R). `use` returns nonzero when already loaded, which would trip the
# config's `set -e` and exit silently -> disable errexit/nounset just for this block.
set +eu; source /broad/software/scripts/useuse; use R-4.1; set -eu

B="${B:-2000}"
PHENO="${BASE}/pheno/subset_ukbb_biomarkers.tsv"
COVARBGZ="${BASE}/covar/final_samples.txt.bgz"
COVAR="${BASE}/covar/covars.tsv"
MAPRAW="${BASE}/sample_id_mapping.txt"
GMAP="${BASE}/geno2eid.txt"
SEL="${1:-${TRAITS}}"                         # optional single trait (array), else all
# insurance-run knobs: evaluate EUR on a subset + write to separate files, without
# touching the full-EUR run's outputs.
EUR_KEEP="${EUR_KEEP:-}"                       # if set, restrict EUR to this keep file (e.g. 10k)
ONLY_VAL="${ONLY_VAL:-}"                       # if set, only this ancestry (e.g. EUR)
OUTSUF="${OUTSUF:-}"                           # output filename suffix (e.g. _eur10k)
DUMP="${DUMP:-}"                                # if set, also dump per-replicate draws (for violins)

# covariates: decompress once (atomic, so parallel tasks can't corrupt it)
[ -s "${COVAR}" ] || { gzip -cd "${COVARBGZ}" > "${COVAR}.tmp.$$" && mv "${COVAR}.tmp.$$" "${COVAR}"; }

# eid<->genotype-id map (atomic; usually IDENTITY -> empty file -> no --map)
# any one sscore for this tag, PIPE-FREE (`ls | head` gives ls SIGPIPE=141 under set -e/pipefail)
SAMPLE=""; for f in "${BASE}"/results/*/*/*_${CLUMP_P1}.sscore; do [ -e "${f}" ] && { SAMPLE="${f}"; break; }; done
[ -n "${SAMPLE}" ] || { echo "no *_${CLUMP_P1}.sscore anywhere -- run scoring first (P+T 03 / PRS-CS 06)"; exit 1; }
if [ ! -f "${GMAP}" ]; then
  python3 "${SCRIPTS_DIR}/build_geno2eid.py" "${SAMPLE}" "${PHENO}" "${MAPRAW}" "${GMAP}.tmp.$$"
  mv "${GMAP}.tmp.$$" "${GMAP}"
fi
MAPFLAG=""; [ -s "${GMAP}" ] && MAPFLAG="--map=${GMAP}"

for TRAIT in ${SEL}; do
  COL=$(col_of "${TRAIT}")
  for VAL in ${ONLY_VAL:-${ANCS}}; do
    D="${BASE}/results/${TRAIT}/${VAL}"; [ -d "${D}" ] || continue
    args=()
    for f in "${D}"/*_${CLUMP_P1}.sscore; do
      [ -e "${f}" ] || continue
      nm=$(basename "${f}"); nm=${nm%_${CLUMP_P1}.sscore}; args+=("${nm}=${f}")
    done
    [ ${#args[@]} -ge 1 ] || { echo "${TRAIT}/${VAL}: no sscores"; continue; }
    BL="${BASELINE:-felix__eur10k}"     # P+T default felix__eur10k; PRS-CS set BASELINE=felix
    printf '%s\n' "${args[@]}" | grep -q "^${BL}=" || BL="${args[0]%%=*}"
    KEEPARG=""; [ "${VAL}" = "EUR" ] && [ -n "${EUR_KEEP}" ] && KEEPARG="--keep=${EUR_KEEP}"
    echo "=== ${TRAIT} / ${VAL} (${COL}) baseline=${BL} B=${B}${KEEPARG:+ (EUR subset)} ==="
    ${RUN_R} "${SCRIPTS_DIR}/bootstrap_r2.R" \
      --pheno="${PHENO}" --pheno_id=eid --pheno_col="${COL}" ${MAPFLAG} ${KEEPARG} \
      --covar="${COVAR}" --covar_id=s --covar_cols="${COVAR_COLS}" \
      --baseline="${BL}" --B="${B}" --seed=1 ${DUMP:+--dump_boot=1} \
      --out="${D}/incR2_${CLUMP_P1}${OUTSUF}.tsv" "${args[@]}"
  done
done
echo "DONE (${SEL}) -> results/<trait>/<anc>/incR2_${CLUMP_P1}.tsv (+ _diff.tsv)"
