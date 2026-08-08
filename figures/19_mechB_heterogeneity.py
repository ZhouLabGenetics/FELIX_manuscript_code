#!/usr/bin/env python3
"""
19_mechB_heterogeneity.py
─────────────────────────
Mechanism B (cross-ancestry marginal heterogeneity): a variant whose per-ancestry
effects point in OPPOSITE directions is cancelled by any shared-effect combination
(the All-by-All fixed-effect meta-analysis and FELIX's homogeneous test),
and is recovered only by FELIX's heterogeneous (K-df) test and the combined
(CCT) test. Two worked examples:

  B1  APOC1 (APOE/APOC1) — alkaline phosphatase   (chr19:44,919,689)
  B2  HPR / TXNL4B       — LDL cholesterol         (chr16:72,080,103)

Visualization (the three-beat story, per example):
  LEFT  forest plot   — per-ancestry effect (beta, 95% CI) for FELIX
                        (filled) and the All-by-All strata (open), coloured by
                        ancestry. African positive, European negative — the sign
                        flip is the WHY. The grey diamond at ~0 is the All-by-All
                        fixed-effect meta estimate: the cancellation made visible.
  RIGHT test ladder   — -log10 p for the four combination choices (All-by-All meta,
                        FELIX homogeneous / heterogeneous / CCT) against the
                        genome-wide line. Only het and CCT cross it: the WHAT.

All numbers are read at run time from the scatter CCT table (replicable).

Output (NEW files only — never overwrites fig4C1/fig4C2):
  manuscript_figures/fig4_mechB_heterogeneity.{pdf,png}     (combined 2x2)
  manuscript_figures/fig4_mechB1_apoc1.{pdf,png}            (B1 standalone)
  manuscript_figures/fig4_mechB2_hpr.{pdf,png}              (B2 standalone)
"""
from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from common import (FIGS_DIR, SCATTER_TABLES, load_scatter, safe_float,
                    apply_manuscript_style)
apply_manuscript_style()

# ancestry palette (Okabe–Ito, matches the rest of Figure 4)
CB = {"AFR": "#D55E00", "EAS": "#56B4E9", "EUR": "#0072B2",
      "NatAm": "#E69F00", "SAS": "#CC79A7"}
ANC = [("AFR", "anc1", "ABA_AFR"), ("EAS", "anc2", "ABA_EAS"),
       ("EUR", "anc3", "ABA_EUR"), ("NatAm", "anc4", "ABA_AMR"),
       ("SAS", "anc5", "ABA_SAS")]
COL_META   = "#7a7a7a"     # neutral grey — shared-effect combinations that cancel
COL_RECOV  = "#CC79A7"     # FELIX brand — tests that recover the signal
GW = 5e-8; GW_LOG = -np.log10(GW)

df = load_scatter(SCATTER_TABLES["cct"])


def row(gene, ph, pos):
    sub = df[(df["phenotype"].astype(str) == str(ph)) &
             (df["Gene"].astype(str).str.contains(gene, na=False, regex=False))]
    if pos is not None and not sub.empty:
        sub = sub.iloc[(sub["Pos"] - pos).abs().argsort()[:1]]
    return sub.iloc[0] if not sub.empty else None


def mlog(p, cap=24.0):
    v = safe_float(p)
    if v is None or v <= 0:
        return cap          # p == 0 (underflow) -> display cap
    return min(-np.log10(v), cap)


def pfmt(p):
    v = safe_float(p)
    if v is None:
        return "NA"
    if v == 0:
        return "<1e-320"
    return f"{v:.0e}" if v < 1e-3 else f"{v:.2g}"


# ── forest panel: per-ancestry opposite-direction effects + meta cancellation ──
def forest(ax, r, trait, title):
    recs = []
    for nm, suf, ab in ANC:
        sb = safe_float(r.get(f"SAIGE_BETA_c_{suf}")); sse = safe_float(r.get(f"SAIGE_SE_c_{suf}"))
        sp = safe_float(r.get(f"SAIGE_p.value_c_{suf}"))
        ab_b = safe_float(r.get(f"{ab}_BETA")); ab_se = safe_float(r.get(f"{ab}_SE"))
        if sb is None or sse is None:
            continue
        recs.append((nm, sb, sse, sp, ab_b, ab_se, CB[nm]))
    recs.sort(key=lambda t: t[1])                      # negative (EUR) at bottom -> positive (AFR) at top
    ys = np.arange(len(recs))
    for y, (nm, sb, sse, sp, ab_b, ab_se, col) in zip(ys, recs):
        # All-by-All stratum (open square), slightly below
        if ab_b is not None and ab_se is not None:
            ax.plot([ab_b - 1.96*ab_se, ab_b + 1.96*ab_se], [y - 0.16]*2,
                    color=col, lw=2.2, alpha=0.8, solid_capstyle="round", zorder=3)
            ax.scatter([ab_b], [y - 0.16], s=110, facecolor="white", edgecolor=col,
                       linewidth=2.0, marker="s", zorder=4)
        # FELIX (filled circle), slightly above
        ax.plot([sb - 1.96*sse, sb + 1.96*sse], [y + 0.16]*2,
                color=col, lw=4.5, solid_capstyle="round", zorder=5)
        ax.scatter([sb], [y + 0.16], s=240, color=col, edgecolor="white",
                   linewidth=1.8, zorder=6)
        star = " ★" if (sp is not None and sp < GW) else ""
        ax.text(sb + (0.006 if sb >= 0 else -0.006), y + 0.46, f"{nm}{star}",
                ha=("left" if sb >= 0 else "right"), va="bottom",
                fontsize=12, weight="bold", color=col)
    ax.axvline(0, color="#333", lw=1.6, zorder=2)
    # All-by-All fixed-effect meta estimate — the cancellation, as a grey diamond at ~0
    mb = safe_float(r.get("ABA_META_BETA")); mse = safe_float(r.get("ABA_META_SE"))
    mp = r.get("ABA_META_Pvalue")
    if mb is not None:
        yj = -0.9
        if mse is not None:
            ax.plot([mb - 1.96*mse, mb + 1.96*mse], [yj]*2, color=COL_META, lw=3.0,
                    solid_capstyle="round", zorder=5)
        ax.scatter([mb], [yj], s=300, marker="D", color=COL_META, edgecolor="white",
                   linewidth=1.4, zorder=6)
        ax.text(mb, yj - 0.42, f"All by All meta (fixed effect): p = {pfmt(mp)}  →  cancels",
                ha="center", va="top", fontsize=10.5, color=COL_META, weight="bold")
    # direction cues
    xl = ax.get_xlim()
    ax.text(xl[1]*0.75, len(recs)+0.15, "African  +  →", ha="center", va="center",
            fontsize=11, color=CB["AFR"], weight="bold", style="italic")
    ax.text(xl[0]*0.75, len(recs)+0.15, "←  European  −", ha="center", va="center",
            fontsize=11, color=CB["EUR"], weight="bold", style="italic")
    ax.set_yticks([]); ax.set_ylim(-1.7, len(recs)+0.7)
    ax.set_xlabel(f"Per-ancestry effect on {trait}  (β, 95% CI)")
    ax.set_title(title, loc="left", weight="bold", fontsize=12.5)
    ax.legend(handles=[
        Line2D([], [], marker="o", color="#555", ls="none", markersize=11, label="FELIX"),
        Line2D([], [], marker="s", color="none", markerfacecolor="white",
               markeredgecolor="#555", markeredgewidth=2.0, ls="none", markersize=10,
               label="All by All (stratum)")],
        loc="lower right", frameon=False, fontsize=10)
    for s in ("top", "right", "left"): ax.spines[s].set_visible(False)


# ── ladder panel: which combination recovers the signal ──────────────────────
def ladder(ax, r):
    tests = [("All by All\nmeta-analysis", r.get("ABA_META_Pvalue")),
             ("FELIX\nhomogeneous", r.get("SAIGE_P_hom_admixed_c")),
             ("FELIX\nheterogeneous", r.get("SAIGE_P_het_admixed_c")),
             ("FELIX\nCCT (combined)", r.get("SAIGE_P_cct_admixed_c"))]
    ys = np.arange(len(tests))[::-1]                   # meta on top
    vals = [mlog(p) for _, p in tests]
    cols = [COL_RECOV if v >= GW_LOG else COL_META for v in vals]
    ax.barh(ys, vals, height=0.62, color=cols, edgecolor="white", linewidth=0.8, zorder=3)
    for y, (lab, p), v, c in zip(ys, tests, vals, cols):
        xt = v + max(vals)*0.02
        capped = (safe_float(p) == 0)
        ax.text(xt, y, ("≥ " if capped else "") + f"p = {pfmt(p)}", ha="left", va="center",
                fontsize=10.5, color=c, weight="bold")
    ax.axvline(GW_LOG, ls="--", lw=1.8, color="#333", zorder=2)
    ax.text(GW_LOG + max(vals)*0.012, 0.4, "genome-wide", ha="left", va="center",
            fontsize=9.5, color="#333", rotation=90)
    ax.set_yticks(ys); ax.set_yticklabels([t[0] for t in tests], fontsize=11)
    ax.set_xlim(0, max(vals)*1.30)
    ax.set_xlabel("Association evidence  ($-\\log_{10} p$)")
    ax.set_title("only the heterogeneous / combined test recovers it",
                 loc="left", weight="bold", fontsize=12.5, pad=12)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.tick_params(axis="y", left=False)


def savefig(fig, name):
    fig.savefig(FIGS_DIR / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(FIGS_DIR / f"{name}.png", bbox_inches="tight", dpi=400)
    plt.close(fig)
    print(f"[19] wrote {name}.pdf + .png")


EXAMPLES = [
    ("B1", "APOC1", "3035995", 44919689, "alkaline phosphatase",
     "Mechanism B1   APOC1 (APOE/APOC1) — alkaline phosphatase"),
    ("B2", "HPR", "3028288", 72080103, "LDL cholesterol",
     "Mechanism B2   HPR / TXNL4B — LDL cholesterol"),
]


def standalone(tag, gene, ph, pos, trait, title, fname):
    r = row(gene, ph, pos)
    if r is None:
        print(f"[19] {gene}/{ph}: not found"); return
    fig, axes = plt.subplots(1, 2, figsize=(15.0, 5.6),
                             gridspec_kw={"width_ratios": [1.35, 1.0], "wspace": 0.30})
    forest(axes[0], r, trait, title)
    ladder(axes[1], r)
    savefig(fig, fname)


def combined():
    fig, axes = plt.subplots(2, 2, figsize=(15.2, 11.2),
                             gridspec_kw={"width_ratios": [1.35, 1.0],
                                          "wspace": 0.30, "hspace": 0.42})
    for i, (tag, gene, ph, pos, trait, title) in enumerate(EXAMPLES):
        r = row(gene, ph, pos)
        if r is None:
            print(f"[19] {gene}/{ph}: not found"); continue
        forest(axes[i, 0], r, trait, title)
        ladder(axes[i, 1], r)
        # echo the numbers for the manuscript / sanity check
        print(f"[19] {tag} {gene}: meta p={pfmt(r.get('ABA_META_Pvalue'))}, "
              f"hom p={pfmt(r.get('SAIGE_P_hom_admixed_c'))}, "
              f"het p={pfmt(r.get('SAIGE_P_het_admixed_c'))}, "
              f"cct p={pfmt(r.get('SAIGE_P_cct_admixed_c'))}")
    savefig(fig, "fig4_mechB_heterogeneity")


if __name__ == "__main__":
    combined()
    standalone(*EXAMPLES[0], "fig4_mechB1_apoc1")
    standalone(*EXAMPLES[1], "fig4_mechB2_hpr")
    print(f"[19] done — outputs in {FIGS_DIR}")
