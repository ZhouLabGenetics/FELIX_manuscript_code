#!/bin/bash
# liftover_predictors.sh — the AoU predictor score files are hg38 (All of Us is GRCh38),
# but the UKB bgen + PRS-CS snpinfo/LD panels are hg19. Lift the predictor SNP keys
# (CHR:POS) hg38 -> hg19 so they match the genotype bim IDs (and the hg19 PRS-CS refs).
# Confirmed by rs11240767: UKB bgen 728951 (hg19) vs predictors 793571 (hg38).
# Downloads liftOver + chain if missing. Backs up originals to predictor_inputs_hg38_bak/.
#
#   bash liftover_predictors.sh
source /humgen/atgu1/fin/lhu/projects/saige_tractor/prs/scripts/config.sh
cd "${BASE}"
LO="${TOOLS}/liftOver"; CHAIN="${TOOLS}/hg38ToHg19.over.chain.gz"

# tools (login-node internet). If wget is blocked, fetch these two files manually.
[ -x "${LO}" ] || { echo "downloading liftOver"; wget -q -O "${LO}" \
  http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver && chmod +x "${LO}"; }
[ -s "${CHAIN}" ] || { echo "downloading chain"; wget -q -O "${CHAIN}" \
  http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz; }
[ -x "${LO}" ] && [ -s "${CHAIN}" ] || { echo "ERROR: need ${LO} and ${CHAIN}"; exit 1; }

[ -d predictor_inputs_hg38_bak ] || cp -r predictor_inputs predictor_inputs_hg38_bak
mkdir -p lift

# 1) unique hg38 CHR:POS across ALL predictor files -> BED (0-based start, name=hg38 key)
awk 'FNR>1{print $1}' predictor_inputs/*/*.txt | sort -u > lift/hg38_snps.txt
awk -F: '{print "chr"$1"\t"($2-1)"\t"$2"\t"$1":"$2}' lift/hg38_snps.txt > lift/hg38.bed
"${LO}" lift/hg38.bed "${CHAIN}" lift/hg19.bed lift/unmapped.bed
awk '{c=$1; sub(/^chr/,"",c); print $4"\t"c":"$3}' lift/hg19.bed > lift/hg38_to_hg19.txt
echo "lifted $(wc -l < lift/hg38_to_hg19.txt) of $(wc -l < lift/hg38_snps.txt) unique SNPs"

# 2) rewrite each predictor file's SNP hg38->hg19 (drop variants that did not lift)
for f in predictor_inputs/*/*.txt; do
  awk 'NR==FNR{m[$1]=$2; next}
       FNR==1{print; next}
       ($1 in m){$1=m[$1]; print}' OFS='\t' lift/hg38_to_hg19.txt "${f}" > "${f}.tmp" && mv "${f}.tmp" "${f}"
done
echo "predictor_inputs are now hg19. Sanity vs genotypes:"
awk 'NR==FNR{g[$2];next} FNR>1{t++; if($1 in g)o++} END{printf "  chr22 overlap %d/%d = %.1f%%\n",o,t,100*o/t}' \
  geno_bed/ldset/chr22.bim predictor_inputs/BMI/meta.txt 2>/dev/null || echo "  (run after 02 finishes chr22)"
