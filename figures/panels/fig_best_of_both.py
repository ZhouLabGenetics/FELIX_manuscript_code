#!/usr/bin/env python3
"""Best-of-both-worlds graphical-abstract panel"""
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle

plt.rcParams["font.family"] = "DejaVu Sans"

# ---- palette (locked) ------------------------------------------------------
AFR = "#9B4E8E"; EUR = "#D8CB8B"
FELIX_C = "#332288"        # indigo — distinct from AFR purple and from all ancestries
FOUND = "#0E7C7B"          # detected (deep teal; ✓ shape also encodes)
MISS  = "#C0504D"          # missed (muted red; ✗ shape also encodes) — clearly not teal
INK = "#1E1E1E"; MUT = "#5C5C5C"; CHR = "#E4E4EA"; CHR_E = "#AEAEB8"

fig, ax = plt.subplots(figsize=(20.0, 7.8))
ax.set_xlim(0, 20.0); ax.set_ylim(0, 7.8); ax.axis("off")   # match ratio -> circles stay round

Y_TOP, Y_BOT = 3.95, 1.85                # the two locus rows
COL = {"hom": 9.4, "het": 13.2, "felix": 17.0}   # spread across full width (banner)
HEAD_Y = 6.0

# ---- helpers ---------------------------------------------------------------
def chromosome(yC):
    """a simple but realistic ideogram: two rounded arms, a centromere, faint bands."""
    x0, x1, h = 0.95, 4.55, 0.36
    xcen = 2.15                          # asymmetric centromere (short p arm, long q arm)
    for a, b in [(x0, xcen - 0.13), (xcen + 0.13, x1)]:
        ax.add_patch(FancyBboxPatch((a, yC - h/2), b - a, h,
                     boxstyle="round,pad=0,rounding_size=0.18",
                     fc=CHR, ec=CHR_E, lw=1.6, zorder=3))
        # faint chromosome bands
        n = 3 if (b - a) < 1.0 else 5
        for i in range(1, n):
            xb = a + (b - a) * i / n
            ax.plot([xb, xb], [yC - h/2 + 0.05, yC + h/2 - 0.05],
                    color="#C9C9D2", lw=1.4, zorder=4)
    ax.add_patch(Circle((xcen, yC), 0.15, fc="#9A9AA6", ec="white", lw=1.2, zorder=5))
    return yC + h/2

def lollipop(x, ytop, ych, n):
    ax.plot([x, x], [ytop, ych - 0.2], color="#7C7C88", lw=2.4, zorder=4)
    ax.add_patch(Circle((x, ych), 0.23, fc=INK, ec="white", lw=1.6, zorder=6))
    ax.text(x, ych, n, ha="center", va="center", fontsize=13.5, fontweight="bold",
            color="white", zorder=7)

def badge(x, y, n):
    ax.add_patch(Circle((x, y), 0.22, fc=INK, ec="none", zorder=6))
    ax.text(x, y, n, ha="center", va="center", fontsize=13.5, fontweight="bold",
            color="white", zorder=7)

def eff_arrows(xbase, yc, ancestry_specific):
    ax.plot([xbase, xbase], [yc - 0.5, yc + 0.4], color="#CFCFCF", lw=1.6, zorder=4)
    lenAFR, lenEUR = (0.88, 0.4) if ancestry_specific else (0.72, 0.72)
    for anc, col, yoff, ln in [("AFR", AFR, 0.24, lenAFR), ("EUR", EUR, -0.34, lenEUR)]:
        y = yc + yoff
        ax.text(xbase - 0.95, y, anc, ha="right", va="center", fontsize=13,
                fontweight="bold", color=col, zorder=5)
        ax.add_patch(FancyArrowPatch((xbase, y), (xbase + ln, y), arrowstyle="-|>",
                     mutation_scale=22, lw=7, color=col, shrinkA=0, shrinkB=0, zorder=5))

def mark(x, y, ok, note=None):
    r = 0.46
    ax.add_patch(Circle((x, y), r, fc=(FOUND if ok else "white"),
                 ec=(FOUND if ok else MISS), lw=(0 if ok else 4.2), zorder=6))
    ax.text(x, y, "✓" if ok else "✗", ha="center", va="center",
            fontsize=40, fontweight="bold", zorder=7,
            color=("white" if ok else MISS))
    if note:
        ax.text(x, y - 0.72, note, ha="center", va="top", fontsize=14,
                color=INK, fontweight="bold", zorder=6)

def header(x, lines, fc, fg, hero=False):
    w, h = 2.3, 1.08
    ax.add_patch(FancyBboxPatch((x - w/2, HEAD_Y - h/2), w, h,
                 boxstyle="round,pad=0.02,rounding_size=0.16",
                 fc=fc, ec=(FELIX_C if hero else "#C4C4C4"),
                 lw=(2.8 if hero else 1.6), zorder=3))
    ax.text(x, HEAD_Y, lines, ha="center", va="center",
            fontsize=(19 if hero else 16), fontweight="bold", color=fg,
            zorder=4, linespacing=1.0)

# ---- title -----------------------------------------------------------------
ax.text(0.18, 7.45, "Best of both worlds", ha="left", va="center",
        fontsize=24, fontweight="bold", color=INK)

# ---- ONE genome, two kinds of locus (left) ---------------------------------
GX0, GX1, GY0, GY1 = 0.4, 5.1, 0.55, 6.95
ax.add_patch(FancyBboxPatch((GX0, GY0), GX1 - GX0, GY1 - GY0,
             boxstyle="round,pad=0.02,rounding_size=0.2",
             fc="#FAFAFA", ec="#9A9A9A", lw=2.2, zorder=2))
gc = (GX0 + GX1) / 2
ax.text(gc, GY1 - 0.34, "One individual's genome", ha="center", va="center",
        fontsize=15.5, fontweight="bold", color=INK, zorder=5)
top_arm = chromosome(5.5)
lollipop(1.55, top_arm, 6.05, "1")
lollipop(3.5, top_arm, 6.05, "2")
for yc, n, title, aspec in [(Y_TOP, "1", "Shared-effect locus", False),
                            (Y_BOT, "2", "Ancestry-specific locus", True)]:
    badge(1.1, yc + 0.72, n)
    ax.text(2.95, yc + 0.72, title, ha="center", va="center", fontsize=13.5,
            fontweight="bold", color=INK, zorder=5)
    eff_arrows(2.95, yc, aspec)

# ---- flow rails ------------------------------------------------------------
for y in (Y_TOP, Y_BOT):
    ax.plot([GX1 + 0.05, COL["felix"] + 0.6], [y, y], color="#DEDEDE", lw=2.6, zorder=1)
    ax.add_patch(FancyArrowPatch((GX1 - 0.02, y), (GX1 + 0.55, y), arrowstyle="-|>",
                 mutation_scale=20, lw=2.6, color="#B9B9BF", zorder=2))
for x in COL.values():
    ax.plot([x, x], [0.75, 5.35], color="#F0F0F0", lw=1.8, zorder=0)

# ---- test headers ----------------------------------------------------------
header(COL["hom"], "Homogeneous\ntest", "#ECECEC", INK)
header(COL["het"], "Heterogeneous\ntest", "#ECECEC", INK)
header(COL["felix"], "FELIXassoc", FELIX_C, "white", hero=True)

# ---- marks -----------------------------------------------------------------
mark(COL["hom"], Y_TOP, True)
mark(COL["hom"], Y_BOT, False, "averaged away")
mark(COL["het"], Y_TOP, False, "underpowered")
mark(COL["het"], Y_BOT, True)
mark(COL["felix"], Y_TOP, True)
mark(COL["felix"], Y_BOT, True)

# highlight the FELIXassoc column = catches both
ax.add_patch(FancyBboxPatch((COL["felix"] - 0.78, Y_BOT - 0.82),
             1.56, (Y_TOP - Y_BOT) + 1.64,
             boxstyle="round,pad=0.02,rounding_size=0.2",
             fc="none", ec=FELIX_C, lw=2.6, ls=(0, (5, 3)), alpha=0.9, zorder=5))

plt.tight_layout()
out = "figures_R"
fig.savefig(f"{out}/png/fig_best_of_both.png", dpi=600, bbox_inches="tight", facecolor="white")
fig.savefig(f"{out}/pdf/fig_best_of_both.pdf", bbox_inches="tight", facecolor="white")
print("wrote fig_best_of_both (png+pdf)")
