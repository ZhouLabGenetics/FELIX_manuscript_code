#!/usr/bin/env python3
# Combine incremental R2 (Stable 9) + FELIX advantage (Stable 10) into ONE table.
# One row per (trait, ancestry, method); the Delta-R2 columns are blank for FELIX (the reference).
import sys, csv, os
RES, TAG = sys.argv[1], sys.argv[2]
TRAIT={"3006923":"Alanine aminotransferase","3007070":"HDL cholesterol","3009744":"MCHC",
 "3013721":"Aspartate aminotransferase","3022192":"Triglycerides","3024929":"Platelet count",
 "3027114":"Total cholesterol","3028288":"LDL cholesterol","3035995":"Alkaline phosphatase","BMI":"BMI","height":"Height"}
ANC={"EUR":0,"AFR":1,"EAS":2,"CSA":3}; MORD={"FELIX ancALL":0,"All by All META":1,"All by All ancestrally-matched":2}
def mname(p): return "FELIX ancALL" if p=="felix" else "All by All META" if p=="meta" else "All by All ancestrally-matched" if p.startswith("matched") else p
def pc(x): return f"{100*float(x):.3f}"
adv={}
for r in csv.DictReader(open(os.path.join(RES,f"diff_vs_felix_{TAG}.tsv")),delimiter="\t"):
    comp=r["comparison"].split(" - ")[0]
    if comp.startswith("felix"): continue
    p=float(r["p_boot"])
    adv[(r["trait"],r["ancestry"],mname(comp))]=(f"{-100*float(r['d_estimate']):.3f}",
        f"{-100*float(r['CI_high']):.3f}", f"{-100*float(r['CI_low']):.3f}", f"{2*min(p,1-p):.3g}")
rows=[]
for r in csv.DictReader(open(os.path.join(RES,f"summary_{TAG}.tsv")),delimiter="\t"):
    if r["predictor"].startswith("felix_tract"): continue
    if r["predictor"] in ("felix","meta") and r["ldref"] not in ("","eur10k"): continue
    m=mname(r["predictor"]); a=adv.get((r["trait"],r["ancestry"],m),("","","",""))
    rows.append([TRAIT.get(r["trait"],r["trait"]), r["ancestry"], m, r["N"], pc(r["incR2"]),
                 pc(r["CI_low"]), pc(r["CI_high"]), *a])
rows.sort(key=lambda x:(ANC[x[1]], x[0], MORD.get(x[2],9)))
out=os.path.join(RES,f"supp_table_{TAG}_combined.tsv")
with open(out,"w") as o:
    o.write("Trait\tValidation ancestry\tMethod\tN\tIncremental R2 (%)\t95% CI lower (%)\t95% CI upper (%)\t"
            "Delta R2 FELIX minus method (%)\tDelta 95% CI lower (%)\tDelta 95% CI upper (%)\tBootstrap P\n")
    for x in rows: o.write("\t".join(map(str,x))+"\n")
print("wrote",out,f"({len(rows)} rows)")
