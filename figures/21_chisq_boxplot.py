#!/usr/bin/env python3
"""
21_chisq_boxplot.py
───────────────────
Overall signal-strength comparison, All by All vs SAIGE-Tractor, as the boxplot
analog of Supplementary Fig. 2b. Where 2b counts each method's independent
genome-wide-significant loci, this plots the DISTRIBUTION of their signal: the
1-degree-of-freedom chi-square converted back from each locus's p-value
(chi2 = isf(p, df=1)), pooled across the 24 analysed phenotypes.

  left  box : All by All meta-analysis  — chi-square at its GW-significant loci
              (from ABA_META_Pvalue)
  right box : SAIGE-Tractor combined (CCT) test — chi-square at its GW-significant
              loci (from SAIGE_P_cct_admixed_c)

Two side-by-side boxplots; log y-axis (top-locus chi-square spans orders of
magnitude); genome-wide chi-square threshold marked. Numbers read at run time
(replicable). NEW file only.

Output: manuscript_figures/SFig_chisq_boxplot_aba_vs_saige.{pdf,png}
"""
from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import chi2
from common import FIGS_DIR, SCATTER_TABLES, load_scatter, safe_float, apply_manuscript_style
apply_manuscript_style()

COL_ABA   = "#0072B2"   # All by All (deep blue)
COL_SAIGE = "#CC79A7"   # SAIGE-Tractor (brand reddish purple)
GW = 5e-8
P_FLOOR = 1e-300        # floor to keep chi-square finite for underflowed p (p==0)
GW_CHISQ = float(chi2.isf(GW, 1))   # ~29.72

df = load_scatter(SCATTER_TABLES["cct"])

def to_chisq(pseries):
    p = pd.to_numeric(pseries, errors="coerce").astype(float)
    p = p.clip(lower=P_FLOOR)
    return pd.Series(chi2.isf(p.values, 1), index=p.index)

df["chisq_aba"]   = to_chisq(df["ABA_META_Pvalue"])
df["chisq_saige"] = to_chisq(df["SAIGE_P_cct_admixed_c"])

# Each method's own independent genome-wide-significant loci = its tophit list
# (tophit_source), the same populations counted in Supplementary Fig. 2b.
# Shared loci are double-listed (one row per source method), so each method's box
# uses its own source rows: SAIGE -> 837, All by All -> 715.
aba   = df.loc[df["tophit_source"] == "ABA",   "chisq_aba"].dropna().values
saige = df.loc[df["tophit_source"] == "SAIGE", "chisq_saige"].dropna().values

p_aba_src   = pd.to_numeric(df.loc[df["tophit_source"]=="ABA","ABA_META_Pvalue"], errors="coerce")
p_saige_src = pd.to_numeric(df.loc[df["tophit_source"]=="SAIGE","SAIGE_P_cct_admixed_c"], errors="coerce")
print(f"[21] All by All tophits n={len(aba)} ({(p_aba_src<GW).sum()} with p<5e-8)  median chi2={np.median(aba):.1f}")
print(f"[21] SAIGE-Tractor tophits n={len(saige)} ({(p_saige_src<GW).sum()} with p<5e-8)  median chi2={np.median(saige):.1f}")
print(f"[21] phenotypes with loci: {df.loc[df['tophit_source'].isin(['ABA','SAIGE']),'phenotype'].nunique()} of 24")

fig, ax = plt.subplots(figsize=(6.8, 6.4))
data = [aba, saige]
cols = [COL_ABA, COL_SAIGE]
labels = [f"All by All\nmeta-analysis", f"SAIGE-Tractor\ncombined test"]
positions = [1, 2]

bp = ax.boxplot(data, positions=positions, widths=0.55, patch_artist=True,
                showfliers=False, medianprops=dict(color="black", lw=2.2),
                whiskerprops=dict(color="#555", lw=1.4),
                capprops=dict(color="#555", lw=1.4))
for patch, c in zip(bp["boxes"], cols):
    patch.set_facecolor(c); patch.set_alpha(0.35); patch.set_edgecolor(c); patch.set_linewidth(1.8)

# jittered points
rng = np.random.default_rng(0)
for pos, d, c in zip(positions, data, cols):
    x = pos + rng.uniform(-0.16, 0.16, size=len(d))
    ax.scatter(x, d, s=10, color=c, alpha=0.28, edgecolor="none", zorder=3)

ax.axhline(GW_CHISQ, ls="--", lw=1.6, color="#333", zorder=2)
ax.text(2.48, GW_CHISQ*1.05, "genome-wide\n(χ² = 29.7)", ha="right", va="bottom",
        fontsize=10.5, color="#333")

ax.set_yscale("log")
ax.set_xticks(positions); ax.set_xticklabels(labels, fontsize=14, weight="bold")
bp["boxes"][0].axes.get_xticklabels()[0].set_color(COL_ABA)
ax.get_xticklabels()[1].set_color(COL_SAIGE)
ax.set_ylabel("Association strength,  χ² (1 d.f.) from p-value")
ax.set_title("Overall signal strength across 24 phenotypes",
             loc="left", weight="bold", fontsize=14)
# annotate N and median above each box
for pos, d, c in zip(positions, data, cols):
    ax.text(pos, d.max()*1.35, f"n = {len(d)}\nmedian χ² = {np.median(d):.0f}",
            ha="center", va="bottom", fontsize=11, color=c, weight="bold")
ax.set_ylim(top=max(aba.max(), saige.max())*3.0)
for s in ("top", "right"): ax.spines[s].set_visible(False)

fig.savefig(FIGS_DIR / "SFig_chisq_boxplot_aba_vs_saige.pdf", bbox_inches="tight")
fig.savefig(FIGS_DIR / "SFig_chisq_boxplot_aba_vs_saige.png", bbox_inches="tight", dpi=400)
plt.close(fig)
print(f"[21] wrote SFig_chisq_boxplot_aba_vs_saige.pdf + .png")
