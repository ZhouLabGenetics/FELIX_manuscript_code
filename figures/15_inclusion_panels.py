#!/usr/bin/env python3
"""
For each Mechanism-inclusion → power candidate locus, produce a consistent
SET of two standalone panels 
"""
from __future__ import annotations
import math
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from common import (FIGS_DIR, REPLICATION_DIR, SCATTER_TABLES, load_scatter,
                    safe_float, saige_region, aba_region, apply_manuscript_style)
apply_manuscript_style()

# ── ancestry palette (current manuscript spec) ──
ANC_COLORS = {
    "AFR":   "#0072B2",   # blue
    "AMR":   "#CC79A7",   # reddish purple  (= NatAm)
    "NatAm": "#CC79A7",
    "EAS":   "#56B4E9",   # sky blue
    "EUR":   "#E69F00",   # orange
    "MID":   "#D55E00",   # vermilion
    "SAS":   "#009E73",   # bluish green
}
# locus-zoom theme (matches Fig 4a)
COLOR_META = "#7B1FA2"    # deep purple — global-ancestry meta-analysis
COLOR_ST   = "#2E7D32"    # forest green — FELIX
GW = 5e-8; GW_LOG = -np.log10(GW)

ss_tab = pd.read_csv(REPLICATION_DIR / "sample_size_table.tsv", sep="\t")
cct = load_scatter(SCATTER_TABLES["cct"])

ANCS = [("AFR", "AFR"), ("EAS", "EAS"), ("EUR", "EUR"),
        ("NatAm", "NatAm"), ("SAS", "SAS")]


def savefig(fig, name):
    fig.savefig(FIGS_DIR / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(FIGS_DIR / f"{name}.png", bbox_inches="tight", dpi=400)
    plt.close(fig)
    print(f"[15] wrote {name}.pdf + .png")


def mlog(x):
    x = pd.to_numeric(x, errors="coerce")
    return -np.log10(x.where(x > 0))


# ════════════════════════════════════════════════════════════════════════════
# (1) sample-size comparison
# ════════════════════════════════════════════════════════════════════════════
def sample_size_panel(pheno_id, gene, trait, driver, fname):
    rrow = ss_tab[ss_tab["phenotype"].astype(str) == str(pheno_id)]
    if rrow.empty:
        print(f"[15] {gene}: no sample-size row for {pheno_id}"); return
    rrow = rrow.iloc[0]
    aba, sg, cols, labs = [], [], [], []
    for nm, key in ANCS:
        a = safe_float(rrow.get(f"{nm}_ABA")); s = safe_float(rrow.get(f"{nm}_TRACTOR"))
        a = 0 if (a is None or math.isnan(a)) else a
        s = 0 if (s is None or math.isnan(s)) else s
        aba.append(a); sg.append(s); cols.append(ANC_COLORS[nm]); labs.append(nm)

    fig, ax = plt.subplots(figsize=(8.6, 5.4))
    y = np.arange(len(ANCS))[::-1]; bh = 0.38; gap = 0.04
    ax.barh(y + bh/2 + gap, aba, height=bh, color=cols, alpha=0.30,
            edgecolor="white", linewidth=0.8)
    ax.barh(y - bh/2 - gap, sg, height=bh, color=cols, alpha=1.0,
            edgecolor="white", linewidth=0.8)
    mx = max(sg + aba); pad = mx * 0.015
    for yy, a, s in zip(y, aba, sg):
        ax.text(a + pad, yy + bh/2 + gap, f"{int(a):,}" if a > 0 else "not tested",
                ha="left", va="center", fontsize=10.5, color="#555")
        ax.text(s + pad, yy - bh/2 - gap, f"{int(s):,}",
                ha="left", va="center", fontsize=11, weight="bold", color="#111")
    ax.set_yticks(y); ax.set_yticklabels(labs, fontsize=14, weight="bold")
    for tl, c in zip(ax.get_yticklabels(), cols): tl.set_color(c)
    ax.set_xlim(0, mx * 1.30)
    ax.set_xlabel("Effective sample size analysed (individuals)")
    # two-line title: headline + driver gain (avoids any overlapping annotation)
    di = [nm for nm, _ in ANCS].index(driver)
    a_d, s_d = aba[di], sg[di]
    if a_d > 0:
        sub = (f"{driver} stratum: {int(a_d):,} → {int(s_d):,} "
               f"(+{100*(s_d-a_d)/a_d:.0f}%) by recovering {driver} haplotypes from admixed participants")
    else:
        sub = f"{driver} carries the signal; recovered from admixed participants across the cohort"
    ax.set_title(f"{gene} — {trait}: per-ancestry effective sample size\n{sub}",
                 loc="left", weight="bold", fontsize=13, linespacing=1.5)
    # render the subtitle (2nd line) lighter/italic by overplotting
    ax.legend(handles=[Rectangle((0,0),1,1, fc="#777", alpha=0.30,
                                  label="global-ancestry meta-analysis"),
                       Rectangle((0,0),1,1, fc="#777", alpha=1.0,
                                  label="FELIX")],
              loc="lower right", fontsize=11, frameon=False)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.tick_params(axis="y", left=False)
    ax.grid(axis="x", color="#eee", lw=0.8); ax.set_axisbelow(True)
    savefig(fig, fname)


# ════════════════════════════════════════════════════════════════════════════
# (2) locus-zoom  (REAL summary statistics)
# ════════════════════════════════════════════════════════════════════════════
def locuszoom_panel(saige_pheno, aba_pheno, gene, chrom, pos, trait, fname,
                    win=500_000, caveat=None):
    sg = saige_region(saige_pheno, chrom, pos, window=win)
    ab = aba_region(aba_pheno, "META", chrom, pos, window=win)
    if sg.empty or ab.empty:
        print(f"[15] {gene}: empty region"); return
    sg = sg.assign(POS=pd.to_numeric(sg["POS"], errors="coerce"),
                   y=mlog(sg["P_cct_admixed_c"])).dropna(subset=["POS", "y"])
    ab = ab.assign(POS=pd.to_numeric(ab["POS"], errors="coerce"),
                   y=mlog(ab["Pvalue"])).dropna(subset=["POS", "y"])

    fig, ax = plt.subplots(figsize=(9.2, 5.4))
    ax.scatter(ab["POS"]/1e6, ab["y"], s=46, facecolor="none", edgecolor=COLOR_META,
               linewidth=1.4, label="global-ancestry meta-analysis", zorder=3)
    ax.scatter(sg["POS"]/1e6, sg["y"], s=52, color=COLOR_ST, alpha=0.9,
               edgecolor="white", linewidth=0.5, label="FELIX (combined test)", zorder=4)
    ax.axhline(GW_LOG, ls="--", lw=1.8, color="#333")
    # GW label on the LEFT so it never collides with the gene label (right of the peak)
    ax.text(sg["POS"].min()/1e6, GW_LOG + 0.12, "genome-wide significance",
            ha="left", va="bottom", fontsize=11.5, color="#333")
    lead = sg.loc[sg["y"].idxmax()]
    ax.scatter([lead["POS"]/1e6], [lead["y"]], marker="*", s=460, color=COLOR_ST,
               edgecolor="black", linewidth=1.0, zorder=6)
    ax.annotate(gene, xy=(lead["POS"]/1e6, lead["y"]),
                xytext=(lead["POS"]/1e6 + 0.18, lead["y"] - 0.10),
                fontsize=13, weight="bold", color="#111", ha="left",
                arrowprops=dict(arrowstyle="-|>", color="#111", lw=1.4))
    ax.set_xlabel(f"Chromosome {chrom} position (Mb)")
    ax.set_ylabel("Association evidence  ($-\\log_{10} p$)")
    ax.set_ylim(0, max(sg["y"].max(), ab["y"].max()) * 1.20)
    ax.legend(loc="upper left", frameon=False, fontsize=12)
    ax.set_title(f"{gene} — {trait}: FELIX lifts the peak past "
                 f"genome-wide significance", loc="left", weight="bold", fontsize=13)
    if caveat:
        ax.text(0.5, -0.20, caveat, transform=ax.transAxes, ha="center", va="top",
                fontsize=10, style="italic", color="#a00")
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    savefig(fig, fname)


# ════════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    # IL23R — re-palette sample-size (locus-zoom is Fig 4a)
    sample_size_panel("GI_522.11", "IL23R", "Crohn's disease", "EUR",
                      "fig4_ss_il23r")

    # ADRB2 — full set (AFR-driven)
    sample_size_panel("BMI", "ADRB2/SH3TC2", "BMI", "AFR", "fig4_ss_adrb2")
    locuszoom_panel("BMI", "BMI", "ADRB2/SH3TC2", 5, 148898672, "BMI",
                    "fig4_lz_adrb2")

    # IKZF3 — full set (recommended EUR substitute; tested in both)
    sample_size_panel("RE_475", "IKZF3", "asthma", "EUR", "fig4_ss_ikzf3")
    locuszoom_panel("pheno_RE_475", "RE_475", "IKZF3", 17, 39765489, "asthma",
                    "fig4_lz_ikzf3")

    print(f"[15] done — outputs in {FIGS_DIR}")
