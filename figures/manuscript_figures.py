#!/usr/bin/env python3
"""
04_manuscript_figures.py
────────────────────────
Publication-quality PDF + PNG figures emphasizing SAIGE-Tractor's
advantages over the AllxAll meta pipeline.

Each figure is written to `manuscript_figures/` as both `.pdf` (vector,
for the manuscript) and `.png` (raster, for slides/preview).

Figures
  Fig 1  sample_size_panel        — local vs global N per pheno × ancestry
  Fig 2  discovery_counts         — shared / SAIGE-only / ABA-only per pheno
  Fig 3  shared_pvalue_scatter    — CCT, HET, HOM panels of -log10(p) SAIGE vs ABA
  Fig 4  per_ancestry_scatter     — 5-panel grid (anc1..5 vs AFR/EAS/EUR/AMR/SAS)
  Fig 5  la_conditioning_boost    — top 20 LA-conditioning unmasking events
  Fig 6  het_vs_hom_differential  — Δ -log10(p) between HET and HOM in SAIGE-only
  Fig 7  replication_barcode      — SAIGE-only loci replicating in ABA
  Fig 8  saige_only_overview      — categorical reasons for SAIGE-only loci
"""
from __future__ import annotations
import json, math
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Patch

from common import (
    FIGS_DIR, REPLICATION_DIR, ANCS, ANC_NAMES, PHENO_LABELS,
    SCATTER_TABLES, load_scatter, safe_float,
    apply_manuscript_style, NAME_SAIGE, NAME_ABA, NAME_ABA_META, NAME_SAIGE_MEGA,
    COLOR_SAIGE, COLOR_ABA, COLOR_SHARED, ANC_COLORS,
)
apply_manuscript_style()

def savefig(fig, name):
    pdf = FIGS_DIR / f"{name}.pdf"
    png = FIGS_DIR / f"{name}.png"
    fig.savefig(pdf); fig.savefig(png)
    plt.close(fig)
    print(f"[04] wrote {pdf.name} + {png.name}")

# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 1 — Sample-size panel
# ═══════════════════════════════════════════════════════════════════════════
def fig1_sample_size():
    ss = pd.read_csv(REPLICATION_DIR / "sample_size_table_long.tsv", sep="\t")
    ss["ABA_N"]     = pd.to_numeric(ss["ABA_N"],     errors="coerce").fillna(0)
    ss["TRACTOR_N"] = pd.to_numeric(ss["TRACTOR_N"], errors="coerce").fillna(0)
    phenos = (ss.groupby("phenotype")["ABA_N"].sum()
                .sort_values(ascending=True).index.tolist())
    if not phenos: return
    labels = [PHENO_LABELS.get(p, p)[:34] for p in phenos]

    fig, axes = plt.subplots(1, 2, figsize=(15, max(5.5, 0.36*len(phenos))),
                              gridspec_kw={"wspace":0.45})

    # LEFT: All by All total N stacked per ancestry (global). Includes MID.
    ax = axes[0]
    ancs_plot = ["AFR","EAS","EUR","NatAm","SAS","MID"]
    bottom = np.zeros(len(phenos))
    for anc in ancs_plot:
        vals = []
        for p in phenos:
            row = ss[(ss["phenotype"]==p) & (ss["ancestry"]==anc)]
            v = row["ABA_N"].iloc[0] if not row.empty else 0
            vals.append(0 if pd.isna(v) else v)
        ax.barh(labels, vals, left=bottom, color=ANC_COLORS[anc], label=anc,
                edgecolor="white", linewidth=0.4, height=0.75)
        bottom = bottom + np.array(vals)
    ax.set_xlabel(f"{NAME_ABA} sample size — N$_{{case}}$ + N$_{{ctrl}}$ per global-ancestry stratum")
    ax.set_title(f"a   Global-ancestry sample size ({NAME_ABA})", weight="bold")
    ax.set_xscale("symlog", linthresh=1000)
    ax.grid(axis="x", color="#eee", linewidth=0.8)
    ax.set_axisbelow(True)
    ax.legend(loc="lower right", title="Ancestry", framealpha=0.95)

    # RIGHT: SAIGE-Tractor individual-equivalent N (haplotypes ÷ 2) per local ancestry
    ax = axes[1]
    bottom = np.zeros(len(phenos))
    for anc in ["AFR","EAS","EUR","NatAm","SAS"]:
        vals = []
        for p in phenos:
            row = ss[(ss["phenotype"]==p) & (ss["ancestry"]==anc)]
            v = row["TRACTOR_N"].iloc[0] if not row.empty else 0
            vals.append(0 if pd.isna(v) else v)
        ax.barh(labels, vals, left=bottom, color=ANC_COLORS[anc], label=anc,
                edgecolor="white", linewidth=0.4, height=0.75)
        bottom = bottom + np.array(vals)
    ax.set_xlabel(f"{NAME_SAIGE} effective sample size — individual-equivalent N per local-ancestry test")
    ax.set_title(f"b   Local-ancestry effective N ({NAME_SAIGE})", weight="bold")
    ax.set_xscale("symlog", linthresh=1000)
    ax.grid(axis="x", color="#eee", linewidth=0.8); ax.set_axisbelow(True)

    fig.suptitle(f"Effective sample sizes per ancestry — global ({NAME_ABA}) vs local ({NAME_SAIGE})",
                 x=0.05, ha="left", y=1.02, weight="bold")
    savefig(fig, "fig1_sample_size_panel")

# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 2 — Discovery counts per phenotype
# ═══════════════════════════════════════════════════════════════════════════
def fig2_discovery_counts():
    df = load_scatter(SCATTER_TABLES["cct"])
    rows = []
    for ph, g in df.groupby("phenotype"):
        n_shared     = int((g[g["tophit_source"]=="SAIGE"]["locus_status"]=="shared").sum())
        n_saige_only = int((g[g["tophit_source"]=="SAIGE"]["locus_status"]=="SAIGE_only").sum())
        n_aba_only   = int((g[g["tophit_source"]=="ABA"  ]["locus_status"]=="ABA_only").sum())
        rows.append((ph, n_shared, n_saige_only, n_aba_only))
    d = pd.DataFrame(rows, columns=["pheno","shared","saige_only","aba_only"])
    d["total"] = d["shared"] + d["saige_only"] + d["aba_only"]
    d = d.sort_values("total", ascending=True)
    labels = [PHENO_LABELS.get(str(p), str(p))[:28] for p in d["pheno"]]

    fig, ax = plt.subplots(figsize=(10, max(5, 0.35*len(d))))
    y = np.arange(len(d))
    ax.barh(y, d["shared"], color=COLOR_SHARED, label="Shared", height=0.78)
    ax.barh(y, d["saige_only"], left=d["shared"], color=COLOR_SAIGE,
            label=f"{NAME_SAIGE} only", height=0.78)
    ax.barh(y, d["aba_only"], left=d["shared"]+d["saige_only"], color=COLOR_ABA,
            label=f"{NAME_ABA} only", height=0.78)
    ax.set_yticks(y); ax.set_yticklabels(labels)
    ax.set_xlabel("Number of independent loci (CCT-conditioned vs meta-analysis)")
    ax.set_title(f"Locus discovery per phenotype — {NAME_SAIGE} CCT vs {NAME_ABA_META}",
                 loc="left", weight="bold")
    ax.legend(loc="lower right", framealpha=0.95)
    ax.grid(axis="x", color="#eee", linewidth=0.8); ax.set_axisbelow(True)
    savefig(fig, "fig2_discovery_counts")

# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 3 — Shared p-value scatter, CCT/HET/HOM panels
# ═══════════════════════════════════════════════════════════════════════════
def _scatter_panel(ax, df, label_y, color_pos, color_neg):
    sub = df[df["locus_status"]=="shared"]
    s = sub[sub["tophit_source"]=="SAIGE"][["SAIGE_pvalue","ABA_pvalue","Gene","phenotype"]].copy()
    s["x"] = -np.log10(pd.to_numeric(s["SAIGE_pvalue"], errors="coerce").clip(1e-320))
    s["y"] = -np.log10(pd.to_numeric(s["ABA_pvalue"],   errors="coerce").clip(1e-320))
    s = s.dropna(subset=["x","y"])
    if s.empty:
        ax.text(0.5,0.5,"no shared loci",transform=ax.transAxes,ha="center")
        return
    maxv = max(s["x"].max(), s["y"].max(), 9.0) * 1.05   # ≥9 so 5e-8 line is in range
    diff = s["x"] - s["y"]
    colors = np.where(diff>1, color_pos, np.where(diff<-1, color_neg, "#666"))
    ax.scatter(s["x"], s["y"], c=colors, s=42, alpha=0.75, edgecolor="white", linewidth=0.4)
    ax.plot([0,maxv],[0,maxv], ls="--", c="#888", lw=1.2)
    ax.axvline(7.3, c="#c44", lw=1.0, alpha=0.5); ax.axhline(7.3, c="#c44", lw=1.0, alpha=0.5)
    ax.set_xlim(0, maxv); ax.set_ylim(0, maxv)
    ax.set_xlabel(f"{NAME_SAIGE}  $-\\log_{{10}}(p)$")
    ax.set_ylabel(label_y)
    ax.set_aspect("equal")
    n_s = int((diff > 1).sum()); n_a = int((diff < -1).sum())
    ax.text(0.04, 0.97,
            f"{NAME_SAIGE} stronger: {n_s}\n{NAME_ABA} stronger:   {n_a}",
            transform=ax.transAxes, va="top", fontsize=12, family="monospace",
            bbox=dict(facecolor="white", edgecolor="#cfd8dc", alpha=0.95, pad=6))

def fig3_shared_pvalue_scatter():
    cct = load_scatter(SCATTER_TABLES["cct"])
    het = load_scatter(SCATTER_TABLES["het"])
    hom = load_scatter(SCATTER_TABLES["hom"])
    fig, axes = plt.subplots(1, 3, figsize=(20, 7))
    aba_lab = f"{NAME_ABA_META}  $-\\log_{{10}}(p)$"
    _scatter_panel(axes[0], cct, aba_lab, COLOR_SAIGE, COLOR_ABA)
    axes[0].set_title(f"a   CCT (combined) vs {NAME_ABA_META}",
                       loc="left", weight="bold")
    _scatter_panel(axes[1], hom, aba_lab, COLOR_SAIGE, COLOR_ABA)
    axes[1].set_title(f"b   HOM (homogeneous) vs {NAME_ABA_META}",
                       loc="left", weight="bold")
    _scatter_panel(axes[2], het, aba_lab, COLOR_SAIGE, COLOR_ABA)
    axes[2].set_title(f"c   HET (heterogeneous) vs {NAME_ABA_META}",
                       loc="left", weight="bold")
    fig.suptitle(f"Shared-locus signal strength — {NAME_SAIGE} combiners vs {NAME_ABA_META}",
                 x=0.06, ha="left", y=1.02, weight="bold")
    savefig(fig, "fig3_shared_pvalue_scatter")

# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 4 — Per-ancestry scatter (5-panel grid)
# ═══════════════════════════════════════════════════════════════════════════
def fig4_per_ancestry_scatter():
    keys  = ["anc_afr","anc_eas","anc_eur","anc_amr","anc_sas"]
    names = ["AFR","EAS","EUR","NatAm (AMR)","SAS"]
    # Give SAS a narrower column because it has very few hits
    fig = plt.figure(figsize=(24, 6))
    gs  = fig.add_gridspec(1, 5, width_ratios=[1, 1, 1, 1, 0.55], wspace=0.45)
    axes = [fig.add_subplot(gs[0, i]) for i in range(5)]
    for ax, key, name in zip(axes, keys, names):
        try:
            df = load_scatter(SCATTER_TABLES[key])
        except FileNotFoundError:
            ax.set_visible(False); continue
        _scatter_panel(ax, df,
                       f"{NAME_ABA}  $-\\log_{{10}}(p)$",
                       ANC_COLORS.get(name.split()[0], "#c0392b"), "#888")
        ax.set_title(name, loc="center", weight="bold", fontsize=20)
    fig.suptitle(f"Per-ancestry signal strength — {NAME_SAIGE} ancestry-specific test vs {NAME_ABA} per-ancestry GWAS",
                 x=0.06, ha="left", y=1.04, weight="bold")
    savefig(fig, "fig4_per_ancestry_scatter")

# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 5 — LA-conditioning boost waterfall
# ═══════════════════════════════════════════════════════════════════════════
def fig5_la_boost():
    df = load_scatter(SCATTER_TABLES["cct"])
    so = df[(df["tophit_source"]=="SAIGE") & (df["locus_status"]=="SAIGE_only")].copy()
    p_marg = pd.to_numeric(so["SAIGE_P_cct_admixed"],   errors="coerce").clip(1e-320)
    p_cond = pd.to_numeric(so["SAIGE_P_cct_admixed_c"], errors="coerce").clip(1e-320)
    so["boost"] = -np.log10(p_cond) + np.log10(p_marg)
    so = so.dropna(subset=["boost"]).sort_values("boost", ascending=False).head(25)
    if so.empty: return
    fig, ax = plt.subplots(figsize=(11, max(4, 0.34*len(so))))
    labels = [f"{str(r['Gene']).split(';')[0][:20]}  ({PHENO_LABELS.get(str(r['phenotype']), str(r['phenotype']))[:18]})"
              for _, r in so.iterrows()]
    y = np.arange(len(so))
    ax.barh(y, so["boost"], color="#1f77b4", edgecolor="white")
    ax.set_yticks(y); ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel("Local-ancestry conditioning boost  (Δ orders of magnitude on $-\\log_{10}p$)")
    ax.set_title(f"Local-ancestry conditioning unmasks suppressed signals — top 25 {NAME_SAIGE}-only loci",
                 loc="left", weight="bold")
    ax.grid(axis="x", color="#eee", linewidth=0.8); ax.set_axisbelow(True)
    for spine in ("top","right"): ax.spines[spine].set_visible(False)
    savefig(fig, "fig5_la_conditioning_boost")

# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 6 — HET vs HOM differential among SAIGE-only loci
# ═══════════════════════════════════════════════════════════════════════════
def fig6_het_hom():
    df = load_scatter(SCATTER_TABLES["cct"])
    so = df[(df["tophit_source"]=="SAIGE") & (df["locus_status"]=="SAIGE_only")].copy()
    p_het = pd.to_numeric(so["SAIGE_P_het_admixed_c"], errors="coerce").clip(1e-320)
    p_hom = pd.to_numeric(so["SAIGE_P_hom_admixed_c"], errors="coerce").clip(1e-320)
    so["delta"] = -np.log10(p_het) + np.log10(p_hom)   # >0 ⇒ HET more sig.
    so = so.dropna(subset=["delta"])
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.hist(so["delta"], bins=40, color="#888", edgecolor="white")
    ax.axvline(0, c="#222", lw=1)
    ax.axvline(2, c="#e65100", ls="--", lw=1, label="Δ = ±2")
    ax.axvline(-2, c="#01579b", ls="--", lw=1)
    n_het = int((so["delta"] >  2).sum()); n_hom = int((so["delta"] < -2).sum())
    ax.text(0.98, 0.95,
            f"HET ≫ HOM (Δ>2): {n_het}\nHOM ≫ HET (Δ<−2): {n_hom}",
            transform=ax.transAxes, ha="right", va="top",
            bbox=dict(facecolor="white", edgecolor="#ddd", pad=4))
    ax.set_xlabel("Δ $-\\log_{10}(p)$  (HET − HOM)")
    ax.set_ylabel(f"{NAME_SAIGE}-only loci")
    ax.set_title(f"Heterogeneous (HET) vs homogeneous (HOM) combiner in {NAME_SAIGE}-only loci",
                 loc="left", weight="bold")
    ax.legend(loc="upper right")
    for spine in ("top","right"): ax.spines[spine].set_visible(False)
    savefig(fig, "fig6_het_vs_hom")

# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 7 — Replication of SAIGE-Tractor unique loci in ABA
# ═══════════════════════════════════════════════════════════════════════════
def fig7_replication_barcode():
    f = REPLICATION_DIR / "replication_saige_only_in_aba.tsv"
    if not f.exists(): return
    df = pd.read_csv(f, sep="\t")
    if df.empty: return
    df["neglog10_aba"] = -np.log10(pd.to_numeric(df["ABA_META_p_window"], errors="coerce").clip(1e-320))
    df["neglog10_saige"] = -np.log10(pd.to_numeric(df["SAIGE_p"], errors="coerce").clip(1e-320))
    df = df.sort_values("neglog10_saige", ascending=False).head(60)
    fig, ax = plt.subplots(figsize=(11, max(5, 0.27*len(df))))
    y = np.arange(len(df))
    # SAIGE bar (positive direction)
    ax.barh(y, df["neglog10_saige"], color=COLOR_SAIGE,
            label=f"{NAME_SAIGE} CCT  $-\\log_{{10}}(p)$",
            height=0.72, alpha=0.92)
    # ABA bar (mirrored to the left)
    ax.barh(y, -df["neglog10_aba"], color=COLOR_ABA,
            label=f"{NAME_ABA_META} best in ±500 kb  $-\\log_{{10}}(p)$",
            height=0.72, alpha=0.92)
    ax.axvline(7.3, c="#c44", lw=1.0, alpha=0.5)
    ax.axvline(-7.3, c="#c44", lw=1.0, alpha=0.5)
    ax.axvline(5.3, c="#888", ls=":", lw=1.0, alpha=0.5)
    ax.axvline(-5.3, c="#888", ls=":", lw=1.0, alpha=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r.Gene[:14]} — {PHENO_LABELS.get(str(r.phenotype), str(r.phenotype))[:16]}"
                        for r in df.itertuples()], fontsize=11)
    ax.invert_yaxis()
    xt = ax.get_xticks()
    ax.set_xticklabels([f"{abs(int(v))}" for v in xt])
    ax.set_xlabel(f"$-\\log_{{10}}(p)$    ← {NAME_ABA_META} (window)        {NAME_SAIGE} CCT →")
    ax.set_title(f"Locus-level replication of {NAME_SAIGE} unique top hits in {NAME_ABA_META} (±500 kb window)",
                 loc="left", weight="bold")
    ax.legend(loc="lower right", framealpha=0.95)
    for spine in ("top","right"): ax.spines[spine].set_visible(False)
    savefig(fig, "fig7_replication_barcode")

# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 8 — SAIGE-only overview: reasons / categories
# ═══════════════════════════════════════════════════════════════════════════
def _classify_saige_only(r):
    betas = []
    for s,_,_,_ in ANCS:
        b = safe_float(r.get(f"SAIGE_BETA_c_{s}"))
        se= safe_float(r.get(f"SAIGE_SE_c_{s}"))
        if b is not None and se and se>0: betas.append(b)
    n_pos = sum(b>0 for b in betas); n_neg = sum(b<0 for b in betas)
    has_opp = n_pos>=1 and n_neg>=1
    p_marg = safe_float(r.get("SAIGE_P_cct_admixed"))
    p_cond = safe_float(r.get("SAIGE_P_cct_admixed_c"))
    la_boost = 0
    if p_marg and p_cond and p_marg>0 and p_cond>0:
        la_boost = -math.log10(p_cond) + math.log10(p_marg)
    p_het = safe_float(r.get("SAIGE_P_het_admixed_c"))
    p_hom = safe_float(r.get("SAIGE_P_hom_admixed_c"))
    het_adv = False
    if p_het and p_hom and p_hom>0:
        het_adv = -math.log10(p_het) + math.log10(p_hom) > 5
    sig_ancs = [s for s,_,_,_ in ANCS
                if (p:=safe_float(r.get(f"SAIGE_p.value_c_{s}"))) and p<1e-4]
    rng = max(betas)-min(betas) if betas else 0
    if has_opp and rng>0.03:               return "Opposing directions"
    if la_boost >= 5:                      return "LA conditioning boost"
    if het_adv:                            return "Effect heterogeneity"
    if len(sig_ancs) == 1:                 return "Single-ancestry signal"
    if len(betas) <= 1:                    return "Single-ancestry signal"
    return "LA refined resolution"

def fig8_reasons_pie():
    df = load_scatter(SCATTER_TABLES["cct"])
    so = df[(df["tophit_source"]=="SAIGE") & (df["locus_status"]=="SAIGE_only")].copy()
    so["reason"] = so.apply(_classify_saige_only, axis=1)
    cnt = so["reason"].value_counts()
    fig, ax = plt.subplots(figsize=(7,6))
    colors = ["#d62728","#1f77b4","#ff7f0e","#9467bd","#2ca02c","#8c564b"]
    wedges, _, _ = ax.pie(cnt.values, labels=None, autopct="%d", startangle=80,
                          colors=colors[:len(cnt)],
                          wedgeprops=dict(edgecolor="white", linewidth=2))
    ax.legend(wedges, [f"{k}  ({v})" for k,v in cnt.items()],
              loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False)
    ax.set_title(f"Why {NAME_SAIGE} finds these loci ({cnt.sum()} {NAME_SAIGE}-only loci)",
                 loc="left", weight="bold")
    savefig(fig, "fig8_saige_only_reasons")

# ───────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    fig1_sample_size()
    fig2_discovery_counts()
    fig3_shared_pvalue_scatter()
    fig4_per_ancestry_scatter()
    fig5_la_boost()
    fig6_het_hom()
    fig7_replication_barcode()
    fig8_reasons_pie()
    print("[04] all figures written to manuscript_figures/")
