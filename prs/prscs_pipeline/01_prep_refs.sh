#!/bin/bash
# 01_prep_refs.sh — build (a) per-chromosome HM3 CHR:POS extract lists, and (b) the
# sample keep lists for every group we will extract: the 4 ancestries (full) plus the
# two multi-ancestry LD panels eur10k and prop10k. Pure awk/shuf; runs in seconds.
source config.sh
cd "${BASE}"

# (0) rebuild a plink2-ready Oxford .sample from the header-less ...rmFirstTwoLines
# file, which is a SINGLE column of sample IDs. plink2 wants ID_1 ID_2 missing + a type
# row, and uses ID_2 as IID -> we write "id id 0" per sample. (If the source ever has
# >=3 cols we take its first three as ID_1 ID_2 missing.)
{ echo "ID_1 ID_2 missing"; echo "0 0 0";
  awk 'NF==1{print $1,$1,0; next} NF>=3{print $1,$2,$3; next} {print $1,$2,0}' "${BGEN_SAMPLE}"
} > "${SAMPLE_PLINK}"
echo "SAMPLE_PLINK: ${SAMPLE_PLINK}  ($(($(wc -l < "${SAMPLE_PLINK}") - 2)) samples)"
echo "  first rows:"; head -4 "${SAMPLE_PLINK}" | sed 's/^/    /'

# (a) HM3 rsID list (BUILD-INDEPENDENT) from the PRS-CS snpinfo. We select bgen variants
# by rsID, then let plink derive CHR:POS from the bgen's OWN hg19 positions (step 02).
# This sidesteps the snpinfo BP build (it is NOT hg19) vs the hg19 predictors/bgen.
echo "building HM3 rsID list from ${SNPINFO}"
awk 'NR==1{for(i=1;i<=NF;i++) if($i=="SNP"||$i=="rsid"||$i=="ID"||$i=="rsID") s=i; next} (s){print $s}' \
  "${SNPINFO}" > hm3/hm3_rsids.txt
echo "hm3/hm3_rsids.txt: $(wc -l < hm3/hm3_rsids.txt) HM3 rsIDs"

# (b) full per-ancestry keep lists from the for_grm .fam files. Use IID for BOTH columns
# (FID=IID=eid) so they match our .sample, whose ID_1=ID_2=eid. The for_grm .fam has FID=0,
# which would fail plink2 --keep (it matches FID AND IID). Using $2,$2 works regardless.
for a in ${ANCS}; do
  awk '{print $2, $2}' "$(fam_of "$a")" > "keeps/${a}.keep"
  echo "keeps/${a}.keep: $(wc -l < "keeps/${a}.keep") samples"
done

# seeded random sampler (reproducible subsets)
seeded(){ openssl enc -aes-256-ctr -pass pass:"$1" -nosalt </dev/zero 2>/dev/null; }

# eur10k = 10k random EUR
shuf --random-source=<(seeded 42) -n "${EUR10K}" "keeps/EUR.keep" > "keeps/eur10k.keep"
echo "keeps/eur10k.keep: $(wc -l < keeps/eur10k.keep) samples"

# prop10k = proportional mix across ancestries (from PROP in config.sh)
: > "keeps/prop10k.keep"
for kv in ${PROP}; do
  a="${kv%%:*}"; n="${kv##*:}"
  shuf --random-source=<(seeded "43${a}") -n "${n}" "keeps/${a}.keep" >> "keeps/prop10k.keep"
done
echo "keeps/prop10k.keep: $(wc -l < keeps/prop10k.keep) samples (target ~10000)"

# ldset = every LD panel + the small (AFR/EAS/CSA) targets -> one small bed (~35k)
sort -u keeps/AFR.keep keeps/EAS.keep keeps/CSA.keep keeps/eur10k.keep keeps/prop10k.keep > keeps/ldset.keep
# eurval = EUR validation target (subsampled for the sanity check; empty EUR_VALID_N = all)
if [ -n "${EUR_VALID_N}" ]; then
  shuf --random-source=<(seeded 44) -n "${EUR_VALID_N}" keeps/EUR.keep > keeps/eurval.keep
else
  cp keeps/EUR.keep keeps/eurval.keep
fi
# needed = everything we import from the bgen (once per chromosome)
sort -u keeps/ldset.keep keeps/eurval.keep > keeps/needed.keep
echo "ldset $(wc -l <keeps/ldset.keep) | eurval $(wc -l <keeps/eurval.keep) | needed $(wc -l <keeps/needed.keep) samples"
