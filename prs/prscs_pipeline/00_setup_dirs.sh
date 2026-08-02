#!/bin/bash
# 00_setup_dirs.sh — create the directory tree on the new cluster. Run once from the
# scripts dir (where config.sh lives).
source /humgen/atgu1/fin/lhu/projects/saige_tractor/prs/scripts/config.sh
mkdir -p "${BASE}"/{predictor_inputs,scripts,pheno,covar,hm3,keeps,geno_bed,results,log}
# geno_bed holds: full/ (transient per-chr import, deleted after subsetting),
# ldset/ (all LD panels + AFR/EAS/CSA targets, ~35k), eurval/ (EUR validation target)
for g in full ldset eurval; do mkdir -p "${BASE}/geno_bed/${g}"; done
echo "tree ready under ${BASE}"
ls -la "${BASE}"
