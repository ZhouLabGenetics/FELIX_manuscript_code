#!/usr/bin/env python3
# make_prscs_supp_tables.py — clean, manuscript-ready supplementary tables from the PRS-CS
# aggregate outputs (summary_<tag>.tsv + diff_vs_felix_<tag>.tsv).
#   python3 make_prscs_supp_tables.py <results_dir> <tag>   # tag=PRSCS (or 5e-8 for P+T)
# Writes, in <results_dir>:
#   supp_table_<tag>_incR2.tsv          incremental R2 (%) + 95% CI, per trait x ancestry x method
#   supp_table_<tag>_felix_advantage.tsv  Delta-R2 (FELIX - comparator) + 95% CI + bootstrap P
import sys, os, csv
RES, TAG = sys.argv[1], sys.argv[2]

TRAIT = {"3006923":"Alanine aminotransferase","3007070":"HDL cholesterol","3009744":"MCHC",
         "3013721":"Aspartate aminotransferase","3022192":"Triglycerides","3024929":"Platelet count",
         "3027114":"Total cholesterol","3028288":"LDL cholesterol","3035995":"Alkaline phosphatase",
         "BMI":"BMI","height":"Height"}
ANC_ORD = {"EUR":0,"AFR":1,"EAS":2,"CSA":3}
def method_name(pred):
    if pred=="felix": return "FELIX ancALL"
    if pred=="meta":  return "All by All META"
    if pred.startswith("matched"): return "All by All ancestrally-matched"
    return pred
def pct(x): return f"{100*float(x):.3f}"

# ---- Table 1: incremental R2 ----
rows=[]
with open(os.path.join(RES,f"summary_{TAG}.tsv")) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        if r["predictor"]=="felix" and r["ldref"] not in ("","eur10k"): continue   # keep the single FELIX/META
        if r["predictor"]=="meta"  and r["ldref"] not in ("","eur10k"): continue
        if r["predictor"].startswith("felix_tract"): continue                       # tract is a separate supp fig
        rows.append((TRAIT.get(r["trait"],r["trait"]), r["ancestry"], method_name(r["predictor"]),
                     int(r["N"]), pct(r["incR2"]), pct(r["CI_low"]), pct(r["CI_high"])))
MORD={"FELIX ancALL":0,"All by All META":1,"All by All ancestrally-matched":2}
rows.sort(key=lambda x:(ANC_ORD.get(x[1],9), x[0], MORD.get(x[2],9)))
with open(os.path.join(RES,f"supp_table_{TAG}_incR2.tsv"),"w") as o:
    o.write("Trait\tValidation ancestry\tMethod\tN\tIncremental R2 (%)\t95% CI lower (%)\t95% CI upper (%)\n")
    for x in rows: o.write("\t".join(map(str,x))+"\n")
print(f"wrote supp_table_{TAG}_incR2.tsv ({len(rows)} rows)")

# ---- Table 2: FELIX advantage (FELIX - comparator) ----
drows=[]
with open(os.path.join(RES,f"diff_vs_felix_{TAG}.tsv")) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        comp = r["comparison"].split(" - ")[0]            # baseline is felix; comparison = "<X> - felix"
        if comp.startswith("felix"): continue             # skip felix-vs-felix / tract
        # FELIX - comparator = -(comparator - felix); flip CI bounds
        d = -float(r["d_estimate"]); lo=-float(r["CI_high"]); hi=-float(r["CI_low"])
        p = float(r["p_boot"]); p2 = 2*min(p, 1-p)         # two-sided bootstrap P for the difference
        drows.append((TRAIT.get(r["trait"],r["trait"]), r["ancestry"], method_name(comp),
                      f"{100*d:.3f}", f"{100*lo:.3f}", f"{100*hi:.3f}", f"{p2:.3g}"))
drows.sort(key=lambda x:(ANC_ORD.get(x[1],9), x[0], x[2]))
with open(os.path.join(RES,f"supp_table_{TAG}_felix_advantage.tsv"),"w") as o:
    o.write("Trait\tValidation ancestry\tComparator\tDelta R2 FELIX-comparator (%)\t95% CI lower (%)\t95% CI upper (%)\tBootstrap P\n")
    for x in drows: o.write("\t".join(map(str,x))+"\n")
print(f"wrote supp_table_{TAG}_felix_advantage.tsv ({len(drows)} rows)")
