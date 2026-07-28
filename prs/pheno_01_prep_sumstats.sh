#!/bin/bash
# Prep global and local sumstats for PRS-CS (lift to GRCh37, HapMap3 match,
# strand-aware allele match, optional GC correction). Env: PHENO_ID, ANCESTRY.

#SBATCH --job-name=prep_sumstats
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=3:00:00
#SBATCH --output=/data/wzhougroup/lhu/saige_tractor/prs_pipeline/logs/step1_prep_%j.out
#SBATCH --error=/data/wzhougroup/lhu/saige_tractor/prs_pipeline/logs/step1_prep_%j.err

set -euo pipefail

: "${PHENO_ID:?Set PHENO_ID env variable}"
: "${ANCESTRY:?Set ANCESTRY (AFR or EUR)}"

if [[ "${ANCESTRY}" != "AFR" && "${ANCESTRY}" != "EUR" ]]; then
    echo "ERROR: ANCESTRY must be AFR or EUR (got: ${ANCESTRY})" >&2
    exit 1
fi

BASE_DIR="/data/wzhougroup/lhu/saige_tractor/prs_pipeline"
SHARED_DIR="${BASE_DIR}/shared"
OUT_DIR="${BASE_DIR}/${PHENO_ID}/${ANCESTRY}"

# Per-ancestry inputs
# Local-sumstat effect/pvalue columns are looked up by NAME (header-based)
# so the script is robust to column-order shifts in the upstream file.
if [[ "${ANCESTRY}" == "AFR" ]]; then
    GLOBAL_SUMSTAT="/data/wzhougroup/lhu/saige_tractor/aou/100pheno/allbyall/processed/${PHENO_ID}_AFR.tsv"
    BETA_COL="BETA_c_anc1"
    P_COL="p.value_c_anc1"
else
    GLOBAL_SUMSTAT="/data/wzhougroup/lhu/saige_tractor/aou/100pheno/allbyall/processed/${PHENO_ID}_EUR.tsv"
    BETA_COL="BETA_c_anc3"
    P_COL="p.value_c_anc3"
fi

LOCAL_SUMSTAT="/data/wzhougroup/lhu/saige_tractor/aou/100pheno/merged_felix_pheno_${PHENO_ID}.txt.gz"
CHAIN_FILE="/data/wzhougroup/lhu/tte_ml/gate/prs/preproc/hg38ToHg19.over.chain"
BIM_FILE="${SHARED_DIR}/ukb_imp_hm3_all.bim"

LIFTOVER_SIF="/data/wzhougroup/lhu/tools/bcftools-liftover_latest.sif"
LIFTOVER_BIND="/data/wzhougroup/lhu:/data/wzhougroup/lhu"

module load Apptainer/1.4.2-1.el9

mkdir -p "${OUT_DIR}"/{sumstats,prscs/global,prscs/local,scores,logs}
cd "${OUT_DIR}/sumstats"

echo "==========================================="
echo "Step 1: prep sumstats"
echo "  PHENO_ID=${PHENO_ID}  ANCESTRY=${ANCESTRY}"
echo "  GLOBAL_SUMSTAT=${GLOBAL_SUMSTAT}"
echo "  LOCAL_SUMSTAT=${LOCAL_SUMSTAT}"
echo "  BETA col: ${BETA_COL}"
echo "  P col:    ${P_COL}"
echo "==========================================="

# ---- A: extract columns from local sumstat (by header name) ----
# Single zcat | awk pass: find column positions from the header, then
# stream the body emitting CHR, POS, Allele1, Allele2, BETA, P in that order.
echo "=== Extracting columns from local sumstat ==="
zcat "${LOCAL_SUMSTAT}" | awk -v beta="${BETA_COL}" -v p="${P_COL}" '
BEGIN { FS = OFS = "\t" }
NR == 1 {
    for (i = 1; i <= NF; i++) {
        if ($i == "CHR")     c_chr  = i
        if ($i == "POS")     c_pos  = i
        if ($i == "Allele1") c_a1   = i
        if ($i == "Allele2") c_a2   = i
        if ($i == beta)      c_beta = i
        if ($i == p)         c_p    = i
    }
    missing = ""
    if (!c_chr)  missing = missing " CHR"
    if (!c_pos)  missing = missing " POS"
    if (!c_a1)   missing = missing " Allele1"
    if (!c_a2)   missing = missing " Allele2"
    if (!c_beta) missing = missing " " beta
    if (!c_p)    missing = missing " " p
    if (missing != "") {
        printf "ERROR: missing column(s) in local sumstat:%s\n", missing > "/dev/stderr"
        exit 1
    }
    print "CHR", "POS", "Allele1", "Allele2", beta, p
    next
}
{ print $c_chr, $c_pos, $c_a1, $c_a2, $c_beta, $c_p }
' > local_extracted.tsv
echo "  $(wc -l < local_extracted.tsv) lines (incl header)"

# ---- B: liftover hg38 -> hg19 ----
echo "=== Liftover hg38 -> hg19 ==="

awk 'NR>1 && $5!="NA" && $7!="NA" {
    print "chr"$1, $2-1, $2, $1":"$2
}' OFS='\t' "${GLOBAL_SUMSTAT}" > global_hg38.bed

awk 'NR>1 {
    print "chr"$1, $2-1, $2, $1":"$2
}' OFS='\t' local_extracted.tsv > local_hg38.bed

/data/wzhougroup/lhu/tools/liftOver global_hg38.bed "${CHAIN_FILE}" global_hg19.bed global_unmapped.bed
echo "  Global: $(wc -l < global_hg38.bed) -> $(wc -l < global_hg19.bed) lifted"

/data/wzhougroup/lhu/tools/liftOver local_hg38.bed "${CHAIN_FILE}" local_hg19.bed local_unmapped.bed
echo "  Local:  $(wc -l < local_hg38.bed) -> $(wc -l < local_hg19.bed) lifted"

# ---- C: build SNP mappings ----
echo "=== Building SNP mappings ==="

awk '
NR==FNR {
    gsub("chr","",$1)
    hg38[$1":"$3] = $4
    next
}
{
    key = $1":"$4
    if (key in hg38)
        print hg38[key], $2, toupper($5), toupper($6)
}
' OFS='\t' global_hg19.bed "${BIM_FILE}" > global_mapping.txt
echo "  Global mapping: $(wc -l < global_mapping.txt) variants"

awk '
NR==FNR {
    gsub("chr","",$1)
    hg38[$1":"$3] = $4
    next
}
{
    key = $1":"$4
    if (key in hg38)
        print hg38[key], $2, toupper($5), toupper($6)
}
' OFS='\t' local_hg19.bed "${BIM_FILE}" > local_mapping.txt
echo "  Local mapping:  $(wc -l < local_mapping.txt) variants"

# ---- D: format PRS-CS input ----
echo "=== Formatting PRS-CS input ==="

# Global sumstat cols: 1=CHR 2=POS 3=A1 4=A2 5=BETA 6=SE 7=P
awk '
NR==FNR {
    map_snp[$1] = $2; map_a1[$1] = $3; map_a2[$1] = $4
    next
}
FNR==1 { print "SNP\tA1\tA2\tBETA\tP"; next }
{
    key = $1":"$2
    if (!(key in map_snp)) next
    beta = $5; pval = $7
    if (beta == "NA" || pval == "NA" || beta == "" || pval == "") next

    sa1 = toupper($3); sa2 = toupper($4)
    ba1 = map_a1[key]; ba2 = map_a2[key]

    if ((sa1 == ba1 && sa2 == ba2) || (sa1 == ba2 && sa2 == ba1)) {
        print map_snp[key], sa2, sa1, beta, pval
        next
    }
    split("A:T,T:A,C:G,G:C", pairs, ",")
    for (p in pairs) { split(pairs[p], kv, ":"); comp[kv[1]] = kv[2] }
    ca1 = comp[sa1]; ca2 = comp[sa2]
    if ((ca1 == ba1 && ca2 == ba2) || (ca1 == ba2 && ca2 == ba1)) {
        print map_snp[key], sa2, sa1, beta, pval
    }
}
' OFS='\t' global_mapping.txt "${GLOBAL_SUMSTAT}" > global_prscs_input.txt

N_GLOBAL=$(awk 'NR>1' global_prscs_input.txt | wc -l)
echo "  Global PRS-CS input: ${N_GLOBAL} SNPs"

# Local extracted cols: 1=CHR 2=POS 3=A1 4=A2 5=BETA 6=P
awk '
NR==FNR {
    map_snp[$1] = $2; map_a1[$1] = $3; map_a2[$1] = $4
    next
}
FNR==1 { print "SNP\tA1\tA2\tBETA\tP"; next }
{
    key = $1":"$2
    if (!(key in map_snp)) next
    beta = $5; pval = $6
    if (beta == "NA" || pval == "NA" || beta == "" || pval == "") next

    sa1 = toupper($3); sa2 = toupper($4)
    ba1 = map_a1[key]; ba2 = map_a2[key]

    if ((sa1 == ba1 && sa2 == ba2) || (sa1 == ba2 && sa2 == ba1)) {
        print map_snp[key], sa2, sa1, beta, pval
        next
    }
    split("A:T,T:A,C:G,G:C", pairs, ",")
    for (p in pairs) { split(pairs[p], kv, ":"); comp[kv[1]] = kv[2] }
    ca1 = comp[sa1]; ca2 = comp[sa2]
    if ((ca1 == ba1 && ca2 == ba2) || (ca1 == ba2 && ca2 == ba1)) {
        print map_snp[key], sa2, sa1, beta, pval
    }
}
' OFS='\t' local_mapping.txt local_extracted.tsv > local_prscs_input.txt

N_LOCAL=$(awk 'NR>1' local_prscs_input.txt | wc -l)
echo "  Local  PRS-CS input: ${N_LOCAL} SNPs"

# ---- E: Genomic control correction ----
# Compute lambda_GC = median(chi-sq from p) / 0.4549 for each sumstat.
# If lambda > 1, divide all chi-sq by lambda and recompute p-values.
# Originals saved as *_raw.txt; corrected versions overwrite the canonical
# file names so pheno_02 needs no changes.
echo "=== Genomic control correction ==="
PRSCS_APPTAINER_SIF="/data/wzhougroup/lhu/tools/python3.sif"
apptainer exec --home /data/wzhougroup/lhu \
  --bind /data/wzhougroup/lhu:/data/wzhougroup/lhu \
  "${PRSCS_APPTAINER_SIF}" \
  python3 - <<'PYEOF'
import os
import scipy.stats as st

def gc_correct(path, label):
    with open(path) as f:
        header = next(f)
        rows = [line.rstrip("\n").split("\t") for line in f]
    pvals = []
    for r in rows:
        try:
            p = float(r[-1])
            if p <= 0:   p = 1e-300
            if p >  1:   p = 1.0
        except ValueError:
            p = 0.5
        pvals.append(p)
    chi2 = [st.chi2.isf(p, 1) for p in pvals]
    median_chi2 = sorted(chi2)[len(chi2)//2]
    lam = median_chi2 / 0.4549

    raw_path = path.replace(".txt", "_raw.txt")
    os.rename(path, raw_path)

    with open(path, "w") as f:
        f.write(header)
        for r, ch2 in zip(rows, chi2):
            if lam > 1.0:
                ch2_corr = ch2 / lam
                new_p = max(st.chi2.sf(ch2_corr, 1), 1e-300)
                r[-1] = f"{new_p:.6e}"
            f.write("\t".join(r) + "\n")

    action = "corrected" if lam > 1.0 else "no correction (lambda<=1)"
    print(f"  {label:6s}: lambda_GC = {lam:.4f}   [{action}]", flush=True)
    print(f"    raw  -> {raw_path}", flush=True)
    print(f"    used -> {path}", flush=True)

gc_correct("global_prscs_input.txt", "GLOBAL")
gc_correct("local_prscs_input.txt",  "LOCAL")
PYEOF

echo "=== Done: PHENO=${PHENO_ID}, ANCESTRY=${ANCESTRY} ==="
