#!/bin/bash
#SBATCH --job-name=3way_simu_unr
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=log/%x_%A.out
#SBATCH --error=log/%x_%A.err

## Build a 3-way admixed benchmark cohort with mixed relatedness.
## N individuals (UNREL_PCT% unrelated, rest in FAMSIZE=10 sib families),
## P common markers; all the inputs FELIX needs (PLINK + sparse GRM
## + VCF) plus the new hybrid format.
##
## Config (env vars; defaults shown):
##   UNR_N     = 10000
##   UNR_P     = 50000
##   UNREL_PCT = 100   (integer 0..100; 100 => all unrelated)
##   UNR_DIR   = ${BASE}/data/unr_N${UNR_N}_P${UNR_P}_unrel${UNREL_PCT}
##   N_SEEDS   = 10
##
## To bash-loop all four mixtures in one interactive session, see end of file.
##
## Pipeline:
##   (a) simu_unrelated.R       -> Block_row1, Admprop, SNP_AF, pheno/null/*
##   (b) make_vcf.R 1 unr       -> vcf/chr1.vcf.gz (DS1/DS2/DS3/ANC fields)
##   (c) make_plink_pruned.R unr -> plink/pruned.{bed,bim,fam}
##   (d) createSparseGRM.R       -> plink/sparseGRM*.mtx
##   (e) tractor_dosage_vcf_to_hybrid -> hybrid/chr1.*
##
## All steps single-threaded (1-CPU cluster policy). Cluster uses apptainer;
## the binary accepts the same syntax as singularity.
##
## Submit (single cohort):
##   sbatch --mem=20G --time=02:00:00 40_simu_unrelated.sh                  # 100% unrelated
##   UNREL_PCT=50 sbatch --mem=20G --time=02:00:00 40_simu_unrelated.sh     # 50/50
##
## Interactive 4-mixture loop in a 20G shell (no sbatch):
##   srun --mem=20G --cpus-per-task=1 --time=04:00:00 --pty bash
##   cd ${FELIX_SIM_BASE}/slurm
##   for PCT in 100 75 50 25; do
##     UNREL_PCT=$PCT bash ./40_simu_unrelated.sh
##   done

set -euo pipefail
module load apptainer 2>/dev/null || module load singularity 2>/dev/null || true

BASE="${FELIX_SIM_BASE:?Set FELIX_SIM_BASE to the simulation directory before submitting}"
export FELIX_SIM_BASE
RTOOLS=/data/wzhougroup/lhu/tools/rtools_latest.sif
SAIGE=/data/wzhougroup/lhu/tools/saige_151.sif
HYBRID_SIF=/data/wzhougroup/lhu/tools/FELIX_latest.sif
APPT="apptainer exec --bind /data/wzhougroup/lhu:/data/wzhougroup/lhu --home /data/wzhougroup/lhu"

UNR_N=${UNR_N:-10000}
UNR_P=${UNR_P:-50000}
UNREL_PCT=${UNREL_PCT:-100}
UNREL_TAG=$(printf "%03d" "$UNREL_PCT")
UNR_DIR=${UNR_DIR:-${BASE}/data/unr_N${UNR_N}_P${UNR_P}_unrel${UNREL_TAG}}
N_SEEDS=${N_SEEDS:-10}

export UNR_N UNR_P UNR_DIR UNREL_PCT
mkdir -p "$UNR_DIR"

echo "==== unr config ===="
echo "  UNR_N     = $UNR_N"
echo "  UNR_P     = $UNR_P"
echo "  UNREL_PCT = $UNREL_PCT"
echo "  UNR_DIR   = $UNR_DIR"
echo "  N_SEEDS   = $N_SEEDS"

## ---- (a) simulator ---------------------------------------------------
echo "==== [a] simu_unrelated.R ===="
$APPT $RTOOLS Rscript ${BASE}/R/simu_unrelated.R $N_SEEDS

## ---- (b) dosage VCF for SAIGE step2 --------------------------------
echo "==== [b] make_vcf.R 1 unr ===="
$APPT $RTOOLS Rscript ${BASE}/R/make_vcf.R 1 unr

## ---- (c) PLINK pruned for sparse GRM -------------------------------
echo "==== [c] make_plink_pruned.R unr ===="
$APPT $RTOOLS Rscript ${BASE}/R/make_plink_pruned.R unr

## ---- (d) sparse GRM for SAIGE step1 --------------------------------
## Unrelated cohort => sparseGRM essentially diagonal at relatednessCutoff=0.125.
echo "==== [d] createSparseGRM.R ===="
PLINK_PREFIX=${UNR_DIR}/plink/pruned
OUT_PREFIX=${UNR_DIR}/plink/sparseGRM
N_RAND=2000
$APPT $SAIGE createSparseGRM.R \
    --plinkFile=${PLINK_PREFIX} \
    --nThreads=1 \
    --outputPrefix=${OUT_PREFIX} \
    --numRandomMarkerforSparseKin=${N_RAND} \
    --relatednessCutoff=0.125

## ---- (e) hybrid-format conversion ----------------------------------
## tractor_dosage_vcf_to_hybrid signature:
##   <input.vcf.gz> <n_ancestries> <sample_chunks> <out_prefix>
## sample_chunks = ceil(N / 32). If this turns out to be floor instead, set
## CHUNK manually via env var CHUNK_OVERRIDE.
echo "==== [e] tractor_dosage_vcf_to_hybrid ===="
CHUNK=${CHUNK_OVERRIDE:-$(( (UNR_N + 31) / 32 ))}
echo "  CHUNK = $CHUNK  (ceil($UNR_N / 32); override via CHUNK_OVERRIDE)"
mkdir -p ${UNR_DIR}/hybrid
$APPT $HYBRID_SIF bash -c \
    "cd ${UNR_DIR} && tractor_dosage_vcf_to_hybrid vcf/chr1.vcf.gz 3 ${CHUNK} hybrid/chr1"

echo "==== Done.  UNR_DIR=${UNR_DIR} ===="
