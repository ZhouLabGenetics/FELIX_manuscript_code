#!/bin/bash
#SBATCH --job-name=3way_simu_admix
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=/data/wzhougroup/lhu/saige_tractor/simulation/3way/log/%x_%A_%a.out
#SBATCH --error=/data/wzhougroup/lhu/saige_tractor/simulation/3way/log/%x_%A_%a.err

## Full pipeline for one all-unrelated 3-way admixture cohort.
## Runs:
##   (a) simu_admix_50k.R       -> Block_row1/{GTot,GAfr,GEur,GNat,LAafr,LAeur,LAnat}.tsv
##                                 Admprop.tsv, SNP_AF.tsv,
##                                 kinship_sparse.rds,
##                                 plink/pruned.{bed,bim,fam},
##                                 pheno/null/{quant,bin10,bin01}_seed01..NN.tsv
##   (b) make_vcf.R 1 admix     -> vcf/chr1.vcf.gz (DS1:DS2:DS3:ANC1:ANC2:ANC3:DSALL)
##                                 vcf/chr1.vcf.gz.csi
##   (c) tractor_dosage_vcf_to_hybrid
##                              -> hybrid/chr1.ancblock.{bin,idx,mks}
##                                 hybrid/chr1.common.geno.bin
##                                 hybrid/chr1.common.variant.{idx,mks}
##                                 hybrid/chr1.rare.carrier.bin
##                                 hybrid/chr1.rare.variant.{idx,mks}
##                                 hybrid/chr1.{meta,samples}
##
## Scenario selection -- either of these works:
##   SLURM_ARRAY_TASK_ID = 1|2|3  (under sbatch --array=1-3)
##   ADMIX_SCENARIO env var        (interactive: 2way_50_50 | 2way_25_75 | 3way_20_30_50)
##
## SUBMIT (cluster batch):
##   sbatch --array=1-3 --mem=48G --time=02:00:00 50_simu_admix_50k.sh
##
## INTERACTIVE, 3 PARALLEL SESSIONS (open 3 separate terminals):
##   # terminal 1
##   srun --mem=48G --cpus-per-task=1 --time=04:00:00 --pty bash
##   cd /data/wzhougroup/lhu/saige_tractor/simulation/3way/scripts/slurm
##   ADMIX_SCENARIO=2way_50_50    bash ./50_simu_admix_50k.sh
##   # terminal 2
##   srun --mem=48G --cpus-per-task=1 --time=04:00:00 --pty bash
##   cd /data/wzhougroup/lhu/saige_tractor/simulation/3way/scripts/slurm
##   ADMIX_SCENARIO=2way_25_75    bash ./50_simu_admix_50k.sh
##   # terminal 3
##   srun --mem=48G --cpus-per-task=1 --time=04:00:00 --pty bash
##   cd /data/wzhougroup/lhu/saige_tractor/simulation/3way/scripts/slurm
##   ADMIX_SCENARIO=3way_20_30_50 bash ./50_simu_admix_50k.sh
##
## INTERACTIVE, 3 SCENARIOS SEQUENTIALLY in one shell:
##   srun --mem=48G --cpus-per-task=1 --time=06:00:00 --pty bash
##   cd /data/wzhougroup/lhu/saige_tractor/simulation/3way/scripts/slurm
##   for SC in 2way_50_50 2way_25_75 3way_20_30_50; do
##     ADMIX_SCENARIO=$SC bash ./50_simu_admix_50k.sh
##   done
##
## Memory: persistent ~28 GB during sim (7 P*N int matrices); write peak ~32 GB.
##         make_vcf.R streaming peak ~30 GB. --mem=48G is comfortable.
## Wall  : sim ~3-5 min, Block_row writes ~5 min, make_vcf ~5-8 min,
##         hybrid converter ~1-3 min -> ~15-25 min per scenario.

set -euo pipefail
module load apptainer 2>/dev/null || module load singularity 2>/dev/null || true

BASE=/data/wzhougroup/lhu/saige_tractor/simulation/3way
RTOOLS=/data/wzhougroup/lhu/tools/rtools_latest.sif
HYBRID_SIF=${HYBRID_SIF:-/data/wzhougroup/lhu/tools/FELIX_latest.sif}
APPT="apptainer exec --bind /data/wzhougroup/lhu:/data/wzhougroup/lhu --home /data/wzhougroup/lhu"

ADMIX_N=${ADMIX_N:-20000}
ADMIX_P=${ADMIX_P:-50000}
N_SEEDS=${N_SEEDS:-2}

## Resolve scenario from SLURM array OR ADMIX_SCENARIO env.
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  case "$SLURM_ARRAY_TASK_ID" in
    1) ADMIX_SCENARIO=2way_50_50    ;;
    2) ADMIX_SCENARIO=2way_25_75    ;;
    3) ADMIX_SCENARIO=3way_20_30_50 ;;
    *) echo "Bad SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID (use 1, 2, or 3)" >&2; exit 1 ;;
  esac
fi
if [[ -z "${ADMIX_SCENARIO:-}" ]]; then
  echo "ADMIX_SCENARIO env var must be set when running outside sbatch --array." >&2
  echo "  e.g.  ADMIX_SCENARIO=2way_50_50 bash $0" >&2
  exit 1
fi
case "$ADMIX_SCENARIO" in
  2way_50_50|2way_25_75|3way_20_30_50) ;;
  *) echo "Bad ADMIX_SCENARIO=$ADMIX_SCENARIO" >&2; exit 1 ;;
esac

ADMIX_DIR=${ADMIX_DIR:-${BASE}/data/admix_${ADMIX_SCENARIO}_N${ADMIX_N}_P${ADMIX_P}}
export ADMIX_DIR ADMIX_N ADMIX_P ADMIX_SCENARIO
mkdir -p "$ADMIX_DIR"

echo "==== admix config ===="
echo "  ADMIX_SCENARIO = $ADMIX_SCENARIO"
echo "  ADMIX_N        = $ADMIX_N"
echo "  ADMIX_P        = $ADMIX_P"
echo "  ADMIX_DIR      = $ADMIX_DIR"
echo "  N_SEEDS        = $N_SEEDS"

## (a) simulator -- writes Block_row1/, kinship_sparse.rds, plink/pruned.{bed,bim,fam}, phenos
echo "==== [a] simu_admix_50k.R ===="
$APPT $RTOOLS Rscript ${BASE}/scripts/R/simu_admix_50k.R $N_SEEDS

## (b) make_vcf -- writes vcf/chr1.vcf.gz with DS1/DS2/DS3/ANC1/ANC2/ANC3/DSALL
echo "==== [b] make_vcf.R 1 admix ===="
$APPT $RTOOLS Rscript ${BASE}/scripts/R/make_vcf.R 1 admix

## (c) hybrid format -- writes hybrid/chr1.* (all 11 files the user listed)
echo "==== [c] tractor_dosage_vcf_to_hybrid ===="
CHUNK=${CHUNK_OVERRIDE:-$(( (ADMIX_N + 31) / 32 ))}
echo "  CHUNK = $CHUNK  (ceil($ADMIX_N / 32); override via CHUNK_OVERRIDE)"
mkdir -p ${ADMIX_DIR}/hybrid
$APPT $HYBRID_SIF bash -c \
    "cd ${ADMIX_DIR} && tractor_dosage_vcf_to_hybrid vcf/chr1.vcf.gz 3 ${CHUNK} hybrid/chr1"

echo "==== Done.  ADMIX_DIR=${ADMIX_DIR} ===="
