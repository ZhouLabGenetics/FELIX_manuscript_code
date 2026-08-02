#!/usr/bin/env python3
# make_predictor_inputs.py — build P+T score files for all predictors from the AoU
# sumstats. Output "SNP A1 A2 BETA P" where SNP=CHR:POS (hg19; ABA & FELIX share coords)
# and A1=effect allele=Allele2 (SAIGE). Pure python3 stdlib. Run LOCALLY; rsync outputs.
#   python3 make_predictor_inputs.py <trait_id> <aba_dir> <saige_dir> <out_dir>
#
# Predictors written per trait:
#   meta.txt              AllxAll multi-ancestry meta  (ABA <id>_META.tsv, BETA)
#   felix.txt             FELIX full-cohort            (BETA_c_ancALL)   [multi-ancestry]
#   matched_<ANC>.txt     AllxAll ancestry-matched     (ABA <id>_<ANC>) [single-ancestry]
#   felix_tract_<ANC>.txt FELIX ancestry-TRACT beta    (BETA_c_anc<k>)  [supp. global-vs-local]
import sys, os, gzip, csv

trait, aba_dir, saige_dir, out_dir = sys.argv[1:5]

# validation ancestry -> ABA per-ancestry file suffix (UKB CSA matched to ABA SAS)
MATCH = {"EUR": "EUR", "AFR": "AFR", "EAS": "EAS", "CSA": "SAS"}
# validation ancestry -> FELIX anc index (anc1=AFR anc2=EAS anc3=EUR anc4=AMR anc5=SAS)
TRACT = {"AFR": "1", "EAS": "2", "EUR": "3", "CSA": "5"}

def felix_path(t):
    if t in ("BMI", "height"):
        return os.path.join(saige_dir, f"merged_FELIX_{t}.txt.gz")
    return os.path.join(saige_dir, f"merged_FELIX_pheno_{t}.txt.gz")

def opent(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)

def write_scores(rows, path):
    seen = set(); n = 0
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as o:
        o.write("SNP\tA1\tA2\tBETA\tP\n")
        for snp, a1, a2, beta, p in rows:
            if snp in seen or beta in ("", "NA", "nan") or p in ("", "NA", "nan"):
                continue
            seen.add(snp); n += 1
            o.write(f"{snp}\t{a1}\t{a2}\t{beta}\t{p}\n")
    return n

def aba_rows(path):
    with opent(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            yield f"{row['CHR']}:{row['POS']}", row['Allele2'], row['Allele1'], row['BETA'], row['Pvalue']

def felix_multi(path, specs):
    """Single pass over the big FELIX file, writing several beta columns at once.
    specs: list of (outpath, beta_col, p_col)."""
    w = {}; seen = {}; cnt = {}
    for outp, _, _ in specs:
        os.makedirs(os.path.dirname(outp), exist_ok=True)
        fh = open(outp, "w"); fh.write("SNP\tA1\tA2\tBETA\tP\n")
        w[outp] = fh; seen[outp] = set(); cnt[outp] = 0
    with opent(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            snp = f"{row['CHR']}:{row['POS']}"; a1 = row['Allele2']; a2 = row['Allele1']
            for outp, bc, pc in specs:
                b = row.get(bc, ""); p = row.get(pc, "")
                if snp in seen[outp] or b in ("", "NA", "nan") or p in ("", "NA", "nan"):
                    continue
                seen[outp].add(snp); cnt[outp] += 1
                w[outp].write(f"{snp}\t{a1}\t{a2}\t{b}\t{p}\n")
    for fh in w.values(): fh.close()
    return cnt

od = os.path.join(out_dir, trait)

# predictor 1: ABA META
meta = os.path.join(aba_dir, f"{trait}_META.tsv")
print(f"{trait} meta: {write_scores(aba_rows(meta), os.path.join(od,'meta.txt'))} SNPs"
      if os.path.exists(meta) else f"{trait} meta: MISSING {meta}")

# predictors 3 + tract: FELIX ancALL and the 4 ancestry tracts, in ONE pass
fp = felix_path(trait)
if os.path.exists(fp):
    specs = [(os.path.join(od, "felix.txt"), "BETA_c_ancALL", "p.value_c_ancALL")]
    for anc, k in TRACT.items():
        specs.append((os.path.join(od, f"felix_tract_{anc}.txt"), f"BETA_c_anc{k}", f"p.value_c_anc{k}"))
    for outp, c in felix_multi(fp, specs).items():
        print(f"{trait} {os.path.basename(outp)[:-4]}: {c} SNPs")
else:
    print(f"{trait} felix: MISSING {fp}")

# predictor 2: ABA ancestry-matched, where available
for val, suf in MATCH.items():
    ap = os.path.join(aba_dir, f"{trait}_{suf}.tsv")
    if os.path.exists(ap):
        print(f"{trait} matched_{val} (ABA {suf}): {write_scores(aba_rows(ap), os.path.join(od,f'matched_{val}.txt'))} SNPs")
    else:
        print(f"{trait} matched_{val}: no ABA {suf} file (skip)")
