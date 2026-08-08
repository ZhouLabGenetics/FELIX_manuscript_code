#!/bin/bash
# config.sh — shared paths/config for the UKB multi-ancestry PRS pipeline on the NEW
# cluster (UGER/qsub). Every script does `source config.sh`. Edit here, once.
set -euo pipefail

BASE=.
SCRIPTS_DIR="${BASE}/scripts"          # all pipeline scripts live here (absolute; qsub-safe)
TOOLS=/humgen/atgu1/fin/lhu/tools
PLINK2="${TOOLS}/plink2"                              # v2.00a6LM (Aug 2024) — bgen->bed
PLINK19="${TOOLS}/plink"                             # v1.9.0-b.7.11 — clump + score (validated dialect)
BGENIX="${TOOLS}/bgen/build/apps/bgenix"            # v1.1.7 — index-based HM3 pre-extraction from bgen
SIF="${TOOLS}/rtools.sif"                   # (no longer used for R)
RUN_R="Rscript"
LOGDIR="${BASE}/log"

# --- inputs already on the cluster -----------------------------------------
GENODIR=~/UKB/genotype   # for_grm .fam = per-ancestry sample lists
BGENDIR=~/ukbb/imputed_v3                          # ukb_imp_chr${c}_v3.bgen (+ .bgi)
BGEN_SAMPLE=~/ukb.sample
SAMPLE_PLINK="${BASE}/ukb_imp.sample"                 # plink2-ready .sample (2 header lines); built by 01_prep_refs.sh
PRSREF="${TOOLS}/PRS_ref"
SNPINFO="${PRSREF}/ldblk_ukbb_eur/snpinfo_ukbb_hm3"     # HM3 variant list (CHR SNP BP A1 A2 ...)

# --- the 4 UKB validation ancestries ---------------------------------------
ANCS="EUR AFR EAS CSA"
fam_of(){ case "$1" in
  EUR) echo "${GENODIR}/ukb.EUR.for_grm.1.pruned.1e7.plink.fam";;
  AFR) echo "${GENODIR}/ukb.AFR.for_grm.1.pruned.plink.fam";;
  EAS) echo "${GENODIR}/ukb.EAS.for_grm.1.pruned.plink.fam";;
  CSA) echo "${GENODIR}/ukb.CSA.for_grm.1.pruned.plink.fam";;
esac; }

# --- LD reference panels ----------------------------------------------------
EUR10K=10000
# prop10k = proportional to the pooled-discovery ancestry sizes from height's FELIX
# N_haplo (anc1=AFR 98909, anc2=EAS 13711, anc3=EUR 310264, anc4=AMR 33639, anc5=SAS 13319).
# AMR excluded (not built in UKB); renormalized over EUR/AFR/EAS/CSA(=SAS) x 10000.
PROP="EUR:7113 AFR:2268 EAS:314 CSA:305"

# --- traits (11 quantitative) + their column name in the phenotype file -----
TRAITS="3006923 3007070 3009744 3013721 3022192 3024929 3027114 3028288 3035995 BMI height"
col_of(){ case "$1" in
  3006923) echo AlanineAminotransferase;; 3007070) echo HDL;; 3009744) echo MCHC;;
  3013721) echo AspartateAminotransferase;; 3022192) echo Triglycerides;; 3024929) echo platelets;;
  3027114) echo TotalCholesterol;; 3028288) echo LDL;; 3035995) echo AlkalinePhosphatase;;
  BMI) echo BMI;; height) echo height;;
esac; }

# --- covariates -------------------------------------------------------------
COVAR_COLS="PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,age,sex,age_sex,age2,age2_sex"

CLUMP_P1="${CLUMP_P1:-5e-8}"   # single P+T threshold (5e-8 default, 5e-6 switch)
R2=0.1; KB=250                 # clumping window

# --- PRS-CS ------------------------------------------------------------------
# python3.9 (scipy+h5py) via: set +eu; source /broad/software/scripts/useuse; use Python-3.9; set -eu
PRSCS="${TOOLS}/PRScs/PRScs.py"
PRSCS_BIM="${BASE}/geno_bed/prscs_hm3"     # rsID bim for --bim_prefix (built by prscs_prep.sh)
NGWAS="${BASE}/n_gwas_prscs.tsv"           # trait <tab> predictor <tab> n_gwas
SST_DIR="${BASE}/prscs_inputs"             # rsID-keyed sumstats (SNP A1 A2 BETA P)
# validation ancestry -> PRS-CS UKBB LD panel (CSA uses the SAS panel)
panel_of(){ case "$1" in EUR) echo eur;; AFR) echo afr;; EAS) echo eas;; CSA) echo sas;; esac; }
# n_gwas lookup for (trait, predictor) from NGWAS
ngwas_of(){ awk -v t="$1" -v p="$2" '$1==t && $2==p{print $3; exit}' "${NGWAS}"; }

# --- efficiency knobs (disk/memory) ----------------------------------------
PLINK2_MEM=20000               # cap plink2 workspace (MB) UNDER h_vmem, else OOM-killed
# EUR validation cohort. "" = full 451,509 (what we want for real numbers). Set to an
# integer only if you deliberately want a subsample. bgenix HM3 pre-extraction keeps the
# intermediate small so full EUR is affordable.
EUR_VALID_N=""
