#!/usr/bin/env python3
# Main Figure 4: European PRS prediction, Local (FELIXassoc) vs Global
# (All-by-All). (a) grouped R2 bars per trait; (b) relative change per trait.
import os, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PRS_DIR = os.path.dirname(os.path.abspath(__file__))
OUT     = os.path.join(PRS_DIR, "Manuscript_Figure4_PRS")
PREFER_ADJUSTED = False
RAW, ADJ = "results/r2_4way.tsv", "results/r2_4way_adjusted_EUR.tsv"

C_LOCAL  = "#C77BA6"   # FELIXassoc (softened #AA3377)
C_GLOBAL = "#6699C4"   # All by All  (softened #0077BB)

TRAITS = [("3006923","ALT"),("3007070","HDL cholesterol"),("3009744","MCHC"),
          ("3013721","AST"),("3022192","Triglycerides"),("3024929","Platelet count"),
          ("3027114","Total cholesterol"),("3028288","LDL cholesterol"),
          ("3035995","ALP"),("BMI","BMI"),("height","Height")]

def read_eur(idd):
    fa, fr = os.path.join(PRS_DIR, idd, ADJ), os.path.join(PRS_DIR, idd, RAW)
    f = fa if (PREFER_ADJUSTED and os.path.exists(fa)) else fr
    g = l = None
    for row in csv.reader(open(f), delimiter="\t"):
        if not row or row[0] == "": continue
        if row[0] == "EUR-global": g = float(row[2])
        if row[0] == "EUR-local":  l = float(row[2])
    return g, l

name = [t[1] for t in TRAITS]
glob = np.array([read_eur(t[0])[0] for t in TRAITS])
loc  = np.array([read_eur(t[0])[1] for t in TRAITS])
pct  = 100 * (loc - glob) / glob
nimp = int(np.sum(loc > glob))
adjusted = PREFER_ADJUSTED and os.path.exists(os.path.join(PRS_DIR, TRAITS[0][0], ADJ))
r2lab = ("Incremental $R^2$" if adjusted else "Prediction $R^2$") + " (European validation)"

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 13,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "axes.linewidth": 1.1, "svg.fonttype": "none"})

fig = plt.figure(figsize=(15.2, 6.2))
gs = fig.add_gridspec(1, 2, width_ratios=[1.5, 1.18], wspace=0.30)

# ---- panel a: grouped bars, sorted by Local R2 (desc) ----------------------
oa = np.argsort(-loc)
axA = fig.add_subplot(gs[0, 0])
x = np.arange(len(oa)); w = 0.4
axA.bar(x - w/2, loc[oa],  w, color=C_LOCAL,  label="Local (FELIXassoc)",  zorder=3)
axA.bar(x + w/2, glob[oa], w, color=C_GLOBAL, label="Global (All by All)", zorder=3)
axA.set_xticks(x); axA.set_xticklabels([name[i] for i in oa], rotation=32,
                                        ha="right", fontsize=11.5, fontweight="bold")
axA.set_ylabel(r2lab, fontsize=13.5)
axA.set_ylim(0, max(loc.max(), glob.max()) * 1.08)
axA.tick_params(axis="y", labelsize=11.5)
axA.margins(x=0.01)
axA.legend(frameon=False, fontsize=12.5, loc="upper right", ncol=1,
           handlelength=1.1, handleheight=1.1)
axA.text(-0.02, 1.06, "a", transform=axA.transAxes, fontsize=18, fontweight="bold")

# ---- panel b: relative change, sorted (desc) -> monotone win story ---------
ob = np.argsort(pct)             # ascending so largest gain on top after barh
axB = fig.add_subplot(gs[0, 1])
cols = [C_LOCAL if pct[i] > 0 else C_GLOBAL for i in ob]
axB.barh(np.arange(len(ob)), pct[ob], color=cols, zorder=3, height=0.68)
axB.axvline(0, color="0.35", linewidth=1.0, zorder=2)
axB.set_yticks(np.arange(len(ob))); axB.set_yticklabels([name[i] for i in ob],
                                                        fontsize=11.5, fontweight="bold")
axB.set_xlabel("Relative change in $R^2$ (%)\nLocal vs Global", fontsize=13)
for j, i in enumerate(ob):                       # value labels just past bar ends
    axB.text(pct[i] + (0.6 if pct[i] >= 0 else -0.6), j, f"{pct[i]:+.1f}",
             va="center", ha="left" if pct[i] >= 0 else "right",
             fontsize=10.5, color="0.15")
axB.set_xlim(pct.min() - 3.2, pct.max() + 3.4)   # headroom so labels never hit the axis
axB.text(-0.02, 1.06, "b", transform=axB.transAxes, fontsize=18, fontweight="bold")

fig.suptitle("Full-cohort inclusion improves European polygenic prediction",
             fontsize=18, fontweight="bold", x=0.02, ha="left", y=1.02)

fig.savefig(OUT + ".pdf", bbox_inches="tight")
fig.savefig(OUT + ".png", bbox_inches="tight", dpi=400)
print(f"wrote {OUT}.pdf/.png   (FELIXassoc improves {nimp}/{len(TRAITS)} EUR traits)")
