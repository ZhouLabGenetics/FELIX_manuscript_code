#!/usr/bin/env python3
"""Figure 1 Panel A - admixed genomes for 3 individuals, ORIGINAL diploid layout
(2 haplotype bars per individual, one chromosome shown 0-100%), but with a REALISTIC
segment-length distribution (broken-stick recombination breakpoints -> natural mix of
long and short tracts) and plausible ancestry mixtures. Drop-in replacement.
Colours sampled from the current figure (locked ancestry palette)."""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch
import matplotlib.patheffects as pe

rng = np.random.default_rng(11)
COL = {"EUR": "#D8CB8B", "AMR": "#DB8F75", "AFR": "#9B4E8E", "EAS": "#6699CC", "SAS": "#8ACCC9"}
LEG = ["EUR", "AMR", "AFR", "EAS", "SAS"]

def hap_segments(props, mean_breaks):
    """broken-stick breakpoints -> exponential-ish segment lengths (mix of long & short)."""
    nb = 0 if mean_breaks == 0 else rng.poisson(mean_breaks)
    bps = np.sort(rng.uniform(0, 1, nb))
    edges = np.concatenate([[0.0], bps, [1.0]])
    ancs = rng.choice(list(props), size=len(edges) - 1, p=list(props.values()))
    return list(zip(edges[:-1], edges[1:], ancs))

INDIV = [
    ("Indiv 1   ·   Global ancestry = EUR",                 {"EUR": 1.0}, 0),
    ("Indiv 2   ·   Global cluster = AMR (NatAm)",          {"EUR": .50, "AMR": .42, "AFR": .08}, 12),
    ("Indiv 3   ·   Between global clusters — excluded", {"AFR": .40, "EUR": .35, "AMR": .20, "EAS": .05}, 16),
]

fig, ax = plt.subplots(figsize=(12.5, 5.2))
BAR_H = 0.62
HAP_GAP = 0.16
IND_GAP = 1.15
top = 0.0
blocks = []
for title, props, mb in INDIV:
    yM = top - BAR_H
    yP = yM - HAP_GAP - BAR_H
    for y in (yM, yP):
        # frame
        ax.add_patch(Rectangle((0, y), 1.0, BAR_H, facecolor="none",
                               edgecolor="grey", linewidth=0.6, zorder=3))
        for x0, x1, a in hap_segments(props, mb):
            ax.add_patch(Rectangle((x0, y), x1 - x0, BAR_H, facecolor=COL[a],
                                   edgecolor="white", linewidth=0.5, zorder=2))
    # individual title
    ax.text(0.0, top + 0.16, title, fontsize=13.5, fontweight="bold", va="bottom", ha="left")
    # diploid brace
    ymid = (yM + BAR_H + yP) / 2
    ax.annotate("", xy=(-0.045, yP - 0.02), xytext=(-0.045, yM + BAR_H + 0.02),
                arrowprops=dict(arrowstyle="-", color="#555", lw=1.4))
    ax.text(-0.055, ymid, "diploid", fontsize=10.5, style="italic", color="#555",
            va="center", ha="right")
    blocks.append(yP)
    top = yP - IND_GAP

bottom = blocks[-1]
ax.set_xlim(-0.135, 1.015)
ax.set_ylim(bottom - 0.35, 0.55)
# x axis
ax.set_xticks([0, .25, .5, .75, 1.0])
ax.set_xticklabels(["0", "25", "50", "75", "100"], fontsize=11)
ax.set_xlabel("chromosome position (%)", fontsize=13)
ax.set_yticks([])
for s in ("left", "right", "top"):
    ax.spines[s].set_visible(False)
ax.spines["bottom"].set_bounds(0, 1.0)
ax.tick_params(axis="x", length=4)

# legend
handles = [Rectangle((0, 0), 1, 1, facecolor=COL[a], edgecolor="grey", linewidth=0.5) for a in LEG]
ax.legend(handles, LEG, title="Local ancestry", ncol=5, loc="upper center",
          bbox_to_anchor=(0.5, -0.13), frameon=False, fontsize=12, title_fontsize=12,
          handlelength=1.3, columnspacing=1.6)

plt.tight_layout()
out = "figures_R"
fig.savefig(f"{out}/png/fig1_panelA_admixed_genomes.png", dpi=400, bbox_inches="tight",
            facecolor="white")
fig.savefig(f"{out}/pdf/fig1_panelA_admixed_genomes.pdf", bbox_inches="tight",
            facecolor="white")
print("wrote fig1_panelA_admixed_genomes (png+pdf)")
