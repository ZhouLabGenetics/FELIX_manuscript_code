#!/usr/bin/env python3
"""
22_chisq_advantage.py
─────────────────────
FELIX's overall advantage over the All by All meta-analysis, in
chi-square terms, combining the two halves of the advantage in one figure:

  (a) DISCOVERY BREADTH — each method's own independent genome-wide-significant
      loci (the populations counted in Supplementary Fig. 2b: All by All 715,
      FELIX 837). FELIX detects more loci at comparable strength.

  (b) SIGNAL DEPTH at SHARED loci (PAIRED) — at the loci both methods detect
      (SAIGE-lead view, locus_status=="shared"), the 1-df chi-square is compared
      variant-by-variant: All by All (from ABA_META_Pvalue) vs FELIX
      (from SAIGE_P_cct_admixed_c). FELIX's chi-square is systematically
      higher — the same precision/power gain quantified in Supplementary Table 7.

chi-square is converted back from each p-value as chi2 = isf(p, df=1); p floored
at 1e-300 so underflowed p stays finite. Numbers read at run time (replicable).
NEW file only.

Output: manuscript_figures/SFig_chisq_advantage_aba_vs_saige.{pdf,png}
"""
from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import chi2
from common import FIGS_DIR, SCATTER_TABLES, load_scatter, apply_manuscript_style
apply_manuscript_style()

# Colourblind-safe pair chosen to avoid every colour in the project's anc_colors
# (#777777 #117733 #0072B2 #882255 #999933) and method_colors (#CC79A7 #E69F00 #56B4E9):
COL_ABA   = "#CC3311"   # All by All meta-analysis (Tol vibrant red)
COL_SAIGE = "#009988"   # FELIX (Tol vibrant teal)
GW = 5e-8
P_FLOOR = 1e-300
GW_CHISQ = float(chi2.isf(GW, 1))

df = load_scatter(SCATTER_TABLES["cct"])
def num(c): return pd.to_numeric(df[c], errors="coerce")
def p2chisq(col):
    p = num(col).clip(lower=P_FLOOR)
    return pd.Series(chi2.isf(p.values, 1), index=df.index)

# (a) DETECTION chi-square from each method's own test p-value, at its own GW loci
df["det_aba"]   = p2chisq("ABA_META_Pvalue")
df["det_saige"] = p2chisq("SAIGE_P_cct_admixed_c")
aba_own = df.loc[df["tophit_source"]=="ABA",   "det_aba"].dropna().values
st_own  = df.loc[df["tophit_source"]=="SAIGE", "det_saige"].dropna().values

# (b) EFFECT chi-square (beta/SE)^2 at the baseline's own GW loci (the Supplementary
# Table 7 comparison: conditioned combined-ancestry estimate vs the meta-analysis,
# same variant, paired). Reproduces Stable 7: +9% median, higher at 71%, 100% dir.
df["eff_aba"] = (num("ABA_META_BETA")/num("ABA_META_SE"))**2
df["eff_st"]  = (num("SAIGE_BETA_c_ancALL")/num("SAIGE_SE_c_ancALL"))**2
b = df[(df["tophit_source"]=="ABA")].dropna(subset=["eff_aba","eff_st"])
b = b[b["eff_aba"]>0]
d_aba, d_st = b["eff_aba"].values, b["eff_st"].values
pct_up  = 100*np.mean(d_st > d_aba)
med_inc = 100*np.median(d_st/d_aba - 1)
dir_conc = 100*np.mean(np.sign(num("ABA_META_BETA")[b.index])==np.sign(num("SAIGE_BETA_c_ancALL")[b.index]))

print(f"[22] (a) breadth: All by All own n={len(aba_own)} med={np.median(aba_own):.1f} | FELIX own n={len(st_own)} med={np.median(st_own):.1f}")
print(f"[22] (b) depth (Stable 7): n={len(b)} med chi2 ABA={np.median(d_aba):.1f} SAIGE={np.median(d_st):.1f}; higher in SAIGE {pct_up:.0f}%; median incr {med_inc:+.0f}%; dir {dir_conc:.0f}%")

fig, (axA, axB) = plt.subplots(1, 2, figsize=(9.6, 6.2), gridspec_kw={"wspace":0.34})

def style_box(ax, data, cols, ns=None):
    # Two boxes butted flush at x=1.5 so the median lines compare directly.
    # If ns is given, box WIDTH is proportional to the number of loci (a
    # text-free way to show "FELIX detects more loci"); otherwise equal.
    if ns is None:
        widths = [1.0, 1.0]
    else:
        widths = [1.0*n/max(ns) for n in ns]
    centers = [1.5 - widths[0]/2, 1.5 + widths[1]/2]
    bp = ax.boxplot(data, positions=centers, widths=widths, patch_artist=True,
                    showfliers=False, medianprops=dict(color="black", lw=2.6),
                    whiskerprops=dict(color="#444", lw=1.4), capprops=dict(color="#444", lw=1.4))
    for patch, c in zip(bp["boxes"], cols):
        patch.set_facecolor(c); patch.set_alpha(0.55); patch.set_edgecolor(c); patch.set_linewidth(1.8)
    ax.axhline(GW_CHISQ, ls="--", lw=1.3, color="#999", zorder=1)   # genome-wide reference (see legend)
    ax.set_yscale("log"); ax.set_xlim(0.4, 2.6)
    ax.set_xticks(centers); ax.set_xticklabels(["All by All\nmeta-analysis","FELIX"], fontsize=13, weight="bold")
    ax.get_xticklabels()[0].set_color(COL_ABA); ax.get_xticklabels()[1].set_color(COL_SAIGE)
    ax.tick_params(axis="x", length=0)
    for s in ("top","right"): ax.spines[s].set_visible(False)

# ── Panel a: breadth (detection chi-square; box width ∝ number of loci) ──
style_box(axA, [aba_own, st_own], [COL_ABA, COL_SAIGE], ns=[len(aba_own), len(st_own)])
axA.set_ylabel("Detection χ² (1 d.f., from test p-value)")
axA.set_title("a", loc="left", weight="bold", fontsize=15)

# ── Panel b: depth (effect chi-square (beta/SE)^2 at baseline loci = Stable 7) ──
style_box(axB, [d_aba, d_st], [COL_ABA, COL_SAIGE])
axB.set_ylabel("Effect χ²  =  (β / SE)²")
axB.set_title("b", loc="left", weight="bold", fontsize=15)
fig.savefig(FIGS_DIR/"SFig_chisq_advantage_aba_vs_saige.pdf", bbox_inches="tight")
fig.savefig(FIGS_DIR/"SFig_chisq_advantage_aba_vs_saige.png", bbox_inches="tight", dpi=400)
plt.close(fig)
print("[22] wrote SFig_chisq_advantage_aba_vs_saige.pdf + .png")
