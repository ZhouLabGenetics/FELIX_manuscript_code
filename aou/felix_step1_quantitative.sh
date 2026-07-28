#!/bin/bash
# FELIXassoc step 1: fit the null GLMM for a quantitative trait (All of Us workbench).
set -o pipefail
set -o errexit

PLINK_BFILE="${PLINK_FILES%.*}"
LOGFILE="${OUT_DIR}/step1_felix_runinfo_${phen_col}.log"

"${SAIGE_STEP1}" \
  --plinkFile="${PLINK_BFILE}" \
  --useSparseGRMtoFitNULL=TRUE \
  --sparseGRMFile="${sparseGRM}" \
  --sparseGRMSampleIDFile="${sparseGRM_IDlist}" \
  --phenoFile="${pheno_file}" \
  --phenoCol="${phen_col}" \
  --invNormalize TRUE \
  --covarColList=PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,age,sex,age2,age_sex,age2_sex \
  --qCovarColList=sex \
  --sampleIDColinphenoFile=IID \
  --traitType=quantitative \
  --IsOverwriteVarianceRatioFile=TRUE \
  --isCateVarianceRatio=FALSE \
  --outputPrefix="${OUT_DIR}"/step1Out_"${phen_col}" \
  --maxiter=100
