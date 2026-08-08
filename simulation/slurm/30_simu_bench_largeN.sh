#!/bin/bash
#SBATCH --job-name=3way_simu_bench
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=log/%x_%A_%a.out
#SBATCH --error=log/%x_%A_%a.err

## Build everything needed to benchmark Tractor-Mix vs FELIX at one
## sample size in {10000, 50000, 100000, 150000, 200000}, with P=1000 markers
## and FAMSIZE=10 full sibs. SLURM array index selects N:
##
##   SLURM_ARRAY_TASK_ID=1 -> N=10000   (small baseline; TM expected to finish fast)
##   SLURM_ARRAY_TASK_ID=2 -> N=50000
##   SLURM_ARRAY_TASK_ID=3 -> N=100000
##   SLURM_ARRAY_TASK_ID=4 -> N=150000
##   SLURM_ARRAY_TASK_ID=5 -> N=200000
##
## Kinship is ONLY written as sparse RDS (block-diagonal, ~16 MB regardless
## of N). The dense NxN Kinship.tsv is intentionally NOT written: at N=100k
## R's copy semantics during fwrite peak at ~3x the matrix size (~240 GB),
## which is what OOM-killed the previous N=100k run. TM's runner densifies
## the sparse RDS in-memory at run time (cheapest form TM can possibly use,
## and the right thing to put on the clock for the benchmark).
##
## All steps are single-threaded to match cluster's 1-CPU policy.
##
## Pipeline per N:
##   (a) simu_bench_largeN.R       -> Block_row1, Admprop, SNP_AF,
##                                     kinship_sparse.rds,
##                                     pheno/null/{quant,bin10,bin01}_seed01.tsv
##   (b) make_vcf.R 1 bench        -> vcf/chr1.vcf.gz (DS1/DS2/DS3/ANC fields)
##   (c) make_tractormix_hapdose.R bench
##                                 -> tractormix/hapdose/chr1.anc{0,1,2}.dosage.txt
##   (d) make_plink_pruned.R bench -> plink/pruned.{bed,bim,fam}
##   (e) createSparseGRM.R         -> plink/sparseGRM*.mtx (+ sampleIDs)
##
## Memory peak per N is now dominated by Block_row1 fwrites of the 7 genotype
## matrices (each ~N*P*8 bytes as data.table during write). 32 GB is plenty
## at all N tested here.
##
## Submit (run one at a time):
##   sbatch --array=1 --mem=16G --time=01:00:00 30_simu_bench_largeN.sh   # N=10k
##   sbatch --array=2 --mem=32G --time=03:00:00 30_simu_bench_largeN.sh   # N=50k
##   sbatch --array=3 --mem=32G --time=04:00:00 30_simu_bench_largeN.sh   # N=100k
##   sbatch --array=4 --mem=32G --time=06:00:00 30_simu_bench_largeN.sh   # N=150k
##   sbatch --array=5 --mem=32G --time=08:00:00 30_simu_bench_largeN.sh   # N=200k

set -euo pipefail
module load singularity

BASE="${FELIX_SIM_BASE:?Set FELIX_SIM_BASE to the simulation directory before submitting}"
export FELIX_SIM_BASE
RTOOLS=/data/wzhougroup/lhu/tools/rtools_latest.sif
SAIGE=/data/wzhougroup/lhu/tools/saige_151.sif
SING="singularity exec --bind /data/wzhougroup/lhu:/data/wzhougroup/lhu --home /data/wzhougroup/lhu"

case "${SLURM_ARRAY_TASK_ID:-1}" in
  1) BENCH_N=10000  ;;
  2) BENCH_N=50000  ;;
  3) BENCH_N=100000 ;;
  4) BENCH_N=150000 ;;
  5) BENCH_N=200000 ;;
  *) echo "Unknown SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"; exit 1 ;;
esac
BENCH_P=${BENCH_P:-1000}
BENCH_DIR=${BASE}/data/bench_N${BENCH_N}_P${BENCH_P}

export BENCH_DIR BENCH_N BENCH_P
mkdir -p "$BENCH_DIR"

echo "==== bench config ===="
echo "  BENCH_N    = $BENCH_N"
echo "  BENCH_P    = $BENCH_P"
echo "  BENCH_DIR  = $BENCH_DIR"

## ---- (a) simulator ---------------------------------------------------
echo "==== [a] simu_bench_largeN.R ===="
$SING $RTOOLS Rscript ${BASE}/R/simu_bench_largeN.R

## ---- (b) dosage VCF for SAIGE step2 ---------------------------------
echo "==== [b] make_vcf.R 1 bench ===="
$SING $RTOOLS Rscript ${BASE}/R/make_vcf.R 1 bench

## ---- (c) Tractor-Mix hapdose ----------------------------------------
echo "==== [c] make_tractormix_hapdose.R bench ===="
$SING $RTOOLS Rscript ${BASE}/R/make_tractormix_hapdose.R bench

## ---- (d) PLINK pruned for sparse GRM --------------------------------
echo "==== [d] make_plink_pruned.R bench ===="
$SING $RTOOLS Rscript ${BASE}/R/make_plink_pruned.R bench

## ---- (e) sparse GRM for SAIGE step1 ---------------------------------
echo "==== [e] createSparseGRM.R ===="
PLINK_PREFIX=${BENCH_DIR}/plink/pruned
OUT_PREFIX=${BENCH_DIR}/plink/sparseGRM
## At P=1000 markers we cap numRandomMarkerforSparseKin to that available.
N_RAND=$(( BENCH_P < 2000 ? BENCH_P : 2000 ))
$SING $SAIGE createSparseGRM.R \
    --plinkFile=${PLINK_PREFIX} \
    --nThreads=1 \
    --outputPrefix=${OUT_PREFIX} \
    --numRandomMarkerforSparseKin=${N_RAND} \
    --relatednessCutoff=0.125

echo "==== Done. BENCH_DIR=${BENCH_DIR} ===="
