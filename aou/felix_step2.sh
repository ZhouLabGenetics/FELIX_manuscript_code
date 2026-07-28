#!/bin/bash
# FELIXassoc step 2: homogeneous + heterogeneous score tests from the FELIXla
# input and the step-1 null model, per chromosome (All of Us workbench).
set -o pipefail
set -o errexit

FELIXLA_INPUT="${FELIXLA_FILES%.*}"
GMMAT_FILE="${step1OUT_DIR}/step1Out_${pheno}.rda"
VAR_FILE="${step1OUT_DIR}/step1Out_${pheno}.varianceRatio.txt"
OUT_PREFIX="${OUT_DIR}/step2Out_felix_${pheno}_chr${chr}"
LOGFILE="${OUT_DIR}/step2_felix_runinfo_${pheno}_${chr}.log"

/bin/time -v step2_SPAtests.R \
  --tractorHybridPrefix="${FELIXLA_INPUT}" \
  --LOCO=FALSE \
  --SAIGEOutputFile="${OUT_PREFIX}" \
  --chrom="${chr}" \
  --minMAF=0 \
  --minMAC=1 \
  --number_of_ancestry=5 \
  --pvalcutoff_of_haplotype=0.000005 \
  --is_admixed=TRUE \
  --is_imputed_data=FALSE \
  --is_Firth_beta=TRUE \
  --GMMATmodelFile="${GMMAT_FILE}" \
  --varianceRatioFile="${VAR_FILE}" \
  --maxMAF_in_groupTest=0.5 \
  --MACCutoff_to_CollapseUltraRare=3
