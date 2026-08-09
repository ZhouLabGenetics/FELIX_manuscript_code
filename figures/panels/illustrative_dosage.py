#!/usr/bin/env python3
""" illustrative haplotype -> ancestry-specific dosage layout"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrow, Circle

COL = {"AFR": "#9B4E8E", "EUR": "#D8CB8B"}
DARK = "#222222"

# --- worked example ---
snps = [0.12, 0.30, 0.50, 0.68, 0.86]
# local-ancestry segments (start, end, ancestry) for each haplotype
segM = [(0.00, 0.40, "EUR"), (0.40, 0.72, "AFR"), (0.72, 1.00, "EUR")]
segP = [(0.00, 0.55, "AFR"), (0.55, 1.00, "EUR")]
ancM = ["EUR", "EUR", "AFR", "AFR", "EUR"]   # ancestry at each SNP on maternal hap
ancP = ["AFR", "AFR", "AFR", "EUR", "EUR"]
alM  = [1, 0, 1, 1, 0]                        # allele on maternal (1=alt, 0=ref)
alP  = [0, 1, 1, 0, 1]
DS = {"AFR": [], "EUR": []}
ANC = {"AFR": [], "EUR": []}
for i in range(5):
    for k in ("AFR", "EUR"):
        DS[k].append(alM[i] * (ancM[i] == k) + alP[i] * (ancP[i] == k))
        ANC[k].append((ancM[i] == k) + (ancP[i] == k))

fig, ax = plt.subplots(figsize=(15.5, 4.6))
BAR_H = 0.60
yM, yP = 4.55, 3.75          # maternal / paternal haplotype bars
xL = 0.0

def draw_hap(y, segs, ancs, alleles, label):
    ax.add_patch(Rectangle((0, y), 1.0, BAR_H, facecolor="none", edgecolor="grey", lw=0.8, zorder=4))
    for x0, x1, a in segs:
        ax.add_patch(Rectangle((x0, y), x1 - x0, BAR_H, facecolor=COL[a], edgecolor="white",
                               lw=0.6, zorder=3))
    ax.text(-0.02, y + BAR_H / 2, label, fontsize=12, fontweight="bold", ha="right", va="center",
            color="#333")
    for i, xs in enumerate(snps):
        filled = alleles[i] == 1
        ax.scatter([xs], [y + BAR_H / 2], s=150, marker="o", zorder=6,
                   facecolor=(DARK if filled else "white"), edgecolor=DARK, linewidths=1.4)

# SNP guide lines + labels
for i, xs in enumerate(snps):
    ax.plot([xs, xs], [1.05, yM + BAR_H + 0.06], color="#bbb", lw=0.8, ls=(0, (2, 2)), zorder=1)
    ax.text(xs, yM + BAR_H + 0.14, f"s{i+1}", fontsize=10.5, ha="center", va="bottom", color="#555")

draw_hap(yM, segM, ancM, alM, "maternal")
draw_hap(yP, segP, ancP, alP, "paternal")

# deconvolution arrow (no text)
ax.annotate("", xy=(0.5, 3.15), xytext=(0.5, 3.55),
            arrowprops=dict(arrowstyle="-|>", color="#444", lw=2.2))

# ancestry-specific dosage tracks
def draw_track(ybase, k, label):
    ax.text(-0.02, ybase + 0.30, label, fontsize=12, fontweight="bold", ha="right", va="center",
            color=COL[k] if k != "EUR" else "#9a7d1f")
    ax.plot([0, 1.0], [ybase, ybase], color="#ccc", lw=0.8, zorder=1)
    for i, xs in enumerate(snps):
        v = DS[k][i]
        if v > 0:
            ax.add_patch(Rectangle((xs - 0.016, ybase), 0.032, v * 0.42, facecolor=COL[k],
                                   edgecolor="white", lw=0.5, zorder=3))
        ax.text(xs, ybase - 0.14, str(v), fontsize=10.5, ha="center", va="top",
                fontweight="bold", color=(COL[k] if k != "EUR" else "#9a7d1f"))

draw_track(2.35, "AFR", "AFR dosage")
draw_track(1.15, "EUR", "EUR dosage")

# highlight the worked SNP (s3: both haplotypes AFR, both carry alt -> DS_AFR = 2)
ax.add_patch(Rectangle((snps[2] - 0.03, 1.0), 0.06, yM + BAR_H - 0.9, facecolor="none",
                       edgecolor="#AA3377", lw=1.4, ls=(0, (3, 2)), zorder=7))

# (title, equation and DS definition text removed - visual-only per request)

# manual legend row (inside axes, no clipping): ancestry swatches + allele glyphs
yl = 0.45
ax.add_patch(Rectangle((0.02, yl - 0.05), 0.03, 0.16, facecolor=COL["AFR"], edgecolor="grey", lw=0.5))
ax.text(0.06, yl + 0.03, "AFR", fontsize=11, va="center")
ax.add_patch(Rectangle((0.12, yl - 0.05), 0.03, 0.16, facecolor=COL["EUR"], edgecolor="grey", lw=0.5))
ax.text(0.16, yl + 0.03, "EUR", fontsize=11, va="center")
ax.text(0.005, yl + 0.03, "local ancestry:", fontsize=10.5, style="italic", color="#666",
        va="center", ha="right")
ax.scatter([0.44], [yl + 0.03], s=150, marker="o", facecolor=DARK, edgecolor=DARK, linewidths=1.4)
ax.text(0.46, yl + 0.03, "alt allele (1)", fontsize=11, va="center")
ax.scatter([0.60], [yl + 0.03], s=150, marker="o", facecolor="white", edgecolor=DARK, linewidths=1.4)
ax.text(0.62, yl + 0.03, "ref allele (0)", fontsize=11, va="center")

ax.set_xlim(-0.12, 1.02)
ax.set_ylim(0.15, 5.45)
ax.axis("off")
plt.tight_layout()
fig.savefig("figures_R/png/fig1_panelB_dosage_deconvolution.png", dpi=400, bbox_inches="tight",
            facecolor="white")
fig.savefig("figures_R/pdf/fig1_panelB_dosage_deconvolution.pdf", bbox_inches="tight",
            facecolor="white")
print("wrote fig1_panelB_dosage_deconvolution (png+pdf)")
