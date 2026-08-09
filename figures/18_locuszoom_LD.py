#!/usr/bin/env python3
"""
18_locuszoom_LD.py
"""
from __future__ import annotations
import subprocess, shutil
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from common import (ROOT, FIGS_DIR, saige_region, aba_region, apply_manuscript_style)
apply_manuscript_style()

# ── 1000 Genomes GRCh38 high-coverage phased panel ──
G1K_BASE = ("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
            "1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV")
def g1k_vcf(chrom):
    return (f"{G1K_BASE}/1kGP_high_coverage_Illumina.chr{chrom}"
            f".filtered.SNV_INDEL_SV_phased_panel.vcf.gz")

REFDIR = ROOT / "ref_1000g"
REGDIR = REFDIR / "regions"; REGDIR.mkdir(parents=True, exist_ok=True)
PLINK  = shutil.which("plink") or "/Users/lhu/miniconda3/bin/plink"
TABIX  = shutil.which("tabix") or "/Users/lhu/miniconda3/bin/tabix"
BGZIP  = shutil.which("bgzip") or "/Users/lhu/miniconda3/bin/bgzip"

# driving-ancestry -> 1000G superpopulation sample list
DRIVER_SP = {"EUR": "EUR", "AFR": "AFR", "NatAm": "AMR", "AMR": "AMR",
             "EAS": "EAS", "SAS": "SAS"}

GW = 5e-8; GW_LOG = -np.log10(GW)
COLOR_META = "#7B1FA2"          # deep purple — global-ancestry meta-analysis track
LEAD_COLOR = "#7B1FA2"          # lead variant diamond
# classic LocusZoom r^2 rainbow
LD_BINS  = [0.8, 0.6, 0.4, 0.2, 0.0]
LD_COLS  = ["#D7191C", "#FDAE61", "#4DAF4A", "#74ADD1", "#313695"]
LD_LABS  = ["0.8 – 1.0", "0.6 – 0.8", "0.4 – 0.6", "0.2 – 0.4", "0.0 – 0.2"]

def r2_color(r2):
    if r2 is None or (isinstance(r2, float) and np.isnan(r2)):
        return "#BDBDBD"
    for thr, col in zip(LD_BINS, LD_COLS):
        if r2 >= thr:
            return col
    return LD_COLS[-1]

def mlog(x):
    x = pd.to_numeric(x, errors="coerce")
    return -np.log10(x.where(x > 0))

def savefig(fig, name):
    fig.savefig(FIGS_DIR / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(FIGS_DIR / f"{name}.png", bbox_inches="tight", dpi=400)
    plt.close(fig)
    print(f"[18] wrote {name}.pdf + .png")


def region_r2(gene, chrom, lead_pos, driver, win=500_000, force=False):
    """Stream the 1000G region for the driving ancestry, build a plink bfile,
    and compute r^2 to the SNP nearest the lead. Returns (pos->r2 dict,
    lead_bp, n_samples). Cached under ref_1000g/regions/<gene>."""
    tag = REGDIR / gene.lower().replace("/", "_")
    sp  = DRIVER_SP[driver]
    keep = REFDIR / f"samples_{sp}.txt"
    bed  = Path(f"{tag}.bed")
    if force or not bed.exists():
        lo, hi = max(1, lead_pos - win), lead_pos + win
        vcf = g1k_vcf(chrom)
        print(f"[18] {gene}: streaming 1000G chr{chrom}:{lo}-{hi} ({sp} founders)…")
        # tabix region -> bgzip
        with open(f"{tag}.vcf.gz", "wb") as out:
            p1 = subprocess.Popen([TABIX, "-h", vcf, f"chr{chrom}:{lo}-{hi}"],
                                  stdout=subprocess.PIPE)
            p2 = subprocess.Popen([BGZIP], stdin=p1.stdout, stdout=out)
            p1.stdout.close(); p2.communicate()
        # keep file FID=IID
        kp = f"{tag}.keep"
        pd.read_csv(keep, header=None)[0].to_frame().assign(b=lambda d: d[0]) \
            .to_csv(kp, sep=" ", header=False, index=False)
        subprocess.run([PLINK, "--vcf", f"{tag}.vcf.gz", "--double-id",
                        "--keep", kp, "--snps-only", "just-acgt",
                        "--biallelic-only", "strict", "--make-bed",
                        "--out", str(tag), "--allow-extra-chr"],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    bim = pd.read_csv(f"{tag}.bim", sep="\t", header=None,
                      names=["chr", "id", "cm", "bp", "a1", "a2"])
    # nearest 1000G SNP to the lead position = LD anchor
    anchor = bim.iloc[(bim["bp"] - lead_pos).abs().argmin()]
    subprocess.run([PLINK, "--bfile", str(tag), "--r2",
                    "--ld-snp", str(anchor["id"]), "--ld-window", "99999",
                    "--ld-window-kb", str(int(2 * win / 1000)),
                    "--ld-window-r2", "0", "--out", str(tag),
                    "--allow-extra-chr"],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    ld = pd.read_csv(f"{tag}.ld", sep=r"\s+")
    r2 = dict(zip(ld["BP_B"].astype(int), ld["R2"].astype(float)))
    r2[int(anchor["bp"])] = 1.0
    nfam = sum(1 for _ in open(f"{tag}.fam"))
    return r2, int(anchor["bp"]), nfam


def lz_LD_panel(saige_pheno, aba_pheno, gene, chrom, pos, trait, driver,
                fname, win=500_000, caveat=None):
    sg = saige_region(saige_pheno, chrom, pos, window=win)
    ab = aba_region(aba_pheno, "META", chrom, pos, window=win)
    if sg.empty:
        print(f"[18] {gene}: empty SAIGE region"); return
    sg = sg.assign(POS=pd.to_numeric(sg["POS"], errors="coerce"),
                   y=mlog(sg["P_cct_admixed_c"])).dropna(subset=["POS", "y"])
    if not ab.empty:
        ab = ab.assign(POS=pd.to_numeric(ab["POS"], errors="coerce"),
                       y=mlog(ab["Pvalue"])).dropna(subset=["POS", "y"])
    lead = sg.loc[sg["y"].idxmax()]
    lead_pos = int(lead["POS"])

    r2map, anchor_bp, nref = region_r2(gene.split("/")[0], chrom, lead_pos, driver, win=win)
    sg["col"] = sg["POS"].astype(int).map(r2map).map(r2_color)
    if not ab.empty:
        ab["col"] = ab["POS"].astype(int).map(r2map).map(r2_color)

    fig, ax = plt.subplots(figsize=(9.4, 5.6))
    # Both summary statistics overlaid, BOTH coloured by r^2 to the lead.
    # All by All meta-analysis = hollow circles (edge = r^2); drawn first/under.
    if not ab.empty:
        ax.scatter(ab["POS"]/1e6, ab["y"], s=52, facecolor="none",
                   edgecolors=ab["col"], linewidth=1.7, zorder=3)
    # FELIX combined test = filled circles (fill = r^2).
    ax.scatter(sg["POS"]/1e6, sg["y"], s=60, c=sg["col"], edgecolor="white",
               linewidth=0.5, zorder=4)
    ax.axhline(GW_LOG, ls="--", lw=1.8, color="#333")
    ax.text(sg["POS"].min()/1e6, GW_LOG + 0.12, "genome-wide significance",
            ha="left", va="bottom", fontsize=11.5, color="#333")
    # lead variant
    ax.scatter([lead["POS"]/1e6], [lead["y"]], marker="D", s=210,
               color=LEAD_COLOR, edgecolor="black", linewidth=1.1, zorder=7)
    ax.annotate(gene.split("/")[0], xy=(lead["POS"]/1e6, lead["y"]),
                xytext=(lead["POS"]/1e6 + 0.16, lead["y"] - 0.10),
                fontsize=13, weight="bold", color="#111", ha="left",
                arrowprops=dict(arrowstyle="-|>", color="#111", lw=1.4))
    ax.set_xlabel(f"Chromosome {chrom} position (Mb)")
    ax.set_ylabel("Association evidence  ($-\\log_{10} p$)")
    ymax = max(sg["y"].max(), ab["y"].max() if not ab.empty else 0)
    ax.set_ylim(0, ymax * 1.20)
    ax.set_title(f"{gene} — {trait}: FELIX vs All by All, "
                 f"variants coloured by LD ($r^2$) to the lead",
                 loc="left", weight="bold", fontsize=13)

    # r^2 legend (LocusZoom rainbow) + method legend
    ld_handles = [Line2D([], [], marker="o", color="none", markerfacecolor=c,
                         markeredgecolor="white", markersize=11, label=l)
                  for c, l in zip(LD_COLS, LD_LABS)]
    leg1 = ax.legend(handles=ld_handles, title="$r^2$ to lead",
                     loc="upper right", frameon=False, fontsize=11,
                     title_fontsize=12, handletextpad=0.3, labelspacing=0.25)
    ax.add_artist(leg1)
    meth = [Line2D([], [], marker="D", color="none", markerfacecolor=LEAD_COLOR,
                   markeredgecolor="black", markersize=11, label="lead variant"),
            Line2D([], [], marker="o", color="none", markerfacecolor="#777",
                   markeredgecolor="white", markersize=10,
                   label="FELIX (combined test)"),
            Line2D([], [], marker="o", color="none", markerfacecolor="none",
                   markeredgecolor="#777", markersize=10, markeredgewidth=1.7,
                   label="All by All (meta-analysis)")]
    ax.legend(handles=meth, loc="upper left", frameon=False, fontsize=11)
    ax.text(0.99, -0.16, f"LD: 1000G GRCh38 high-coverage, {DRIVER_SP[driver]} "
            f"founders (n={nref})", transform=ax.transAxes, ha="right",
            va="top", fontsize=9.5, color="#777")
    if caveat:
        ax.text(0.5, -0.24, caveat, transform=ax.transAxes, ha="center",
                va="top", fontsize=10, style="italic", color="#a00")
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    savefig(fig, fname)


if __name__ == "__main__":
    # IL23R — Crohn's (EUR-driven)
    lz_LD_panel("pheno_GI_522.11", "GI_522.11", "IL23R", 1, 67240275,
                "Crohn's disease", "EUR", "fig4_lz_il23r_LD")
    # IKZF3 — asthma (EUR-driven; 17q21)
    lz_LD_panel("pheno_RE_475", "RE_475", "IKZF3", 17, 39765489,
                "asthma", "EUR", "fig4_lz_ikzf3_LD")
    # ADRB2/SH3TC2 — BMI (AFR-driven)
    lz_LD_panel("BMI", "BMI", "ADRB2/SH3TC2", 5, 148898672,
                "BMI", "AFR", "fig4_lz_adrb2_LD")

    print(f"[18] done — outputs in {FIGS_DIR}")
