#!/usr/bin/env python3
"""Figure 1, Panel B — sub-panels 2, 3, 4 regenerated as a consistent, clean set
(style matching the user's reference). Writes three separate PNGs/PDFs:
  fig1_panel2_storage, fig1_panel3_tests, fig1_panel4_outputs
Panel 2: FELIXla storage (fixes old 'HAPLAI-Packed' name).
Panel 3: the two complementary score tests (H0's) -> Cauchy combination.
Panel 4: what FELIXassoc returns per variant, as a symbolic (no-number) table with an
         aggregate ALL row (= homogeneous / standard-GWAS) and the combined-test box.
Grounded in src/Main.cpp outputs (P_hom, P_het, P_cct; per-ancestry AF/N; conditioned)."""
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle, Circle, Ellipse

plt.rcParams["font.family"] = "DejaVu Sans"

# consistent 5-ancestry palette + accents
ANC = [("AFR", "#D55E00"), ("EAS", "#56B4E9"), ("EUR", "#0077BB"),
       ("AMR", "#E69F00"), ("SAS", "#CC79A7")]
INK="#1E1E1E"; MUT="#6A6A6A"; GRID="#D8D8D8"
BROWN="#8A6D1F"; PURPLE="#6A2C6E"; GREEN="#2E6B3E"
NAVY="#26315C"; ORANGE="#C4571F"; INDIGO="#332288"; BLUE="#2C6FB0"; GOLD="#C99A2E"
AFRp="#9B4E8E"; EURk="#D8CB8B"                      # local-ancestry tract colours (panel 2)

def card(ax, W, H, title, tcolor):
    ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off")
    ax.text(0.12, H-0.32, title, ha="left", va="center", fontsize=20,
            fontweight="bold", color=tcolor)
    ax.add_patch(FancyBboxPatch((0.14, 0.18), W-0.28, H-0.94,
                 boxstyle="round,pad=0.02,rounding_size=0.14",
                 fc="white", ec=tcolor, lw=2.6, zorder=1))

def save(fig, name):
    fig.savefig(f"figures_R/png/{name}.png", dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(f"figures_R/pdf/{name}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig); print("wrote", name)

# ============================ PANEL 2: storage ==============================
def panel2():
    W,H=5.6,7.2; fig,ax=plt.subplots(figsize=(W,H)); card(ax,W,H,"2. Scalable & compressed storage",BROWN)
    # --- dosage mini-illustration (top) ---
    x0,x1=0.55,5.05; ytop,ybot=6.05,5.62
    segs=[(0,1.4,EURk),(1.4,2.2,AFRp),(3.6,1.0,EURk)]   # maternal tracts
    for (a,w,c) in segs: ax.add_patch(Rectangle((x0+a*0.72,ytop),w*0.72,0.34,fc=c,ec="white",lw=1.2,zorder=3))
    segs2=[(0,0.9,AFRp),(0.9,2.7,EURk),(3.6,1.0,AFRp)]
    for (a,w,c) in segs2: ax.add_patch(Rectangle((x0+a*0.72,ybot),w*0.72,0.34,fc=c,ec="white",lw=1.2,zorder=3))
    ax.text(x0-0.05,ytop+0.17,"mat",ha="right",va="center",fontsize=9,color=MUT)
    ax.text(x0-0.05,ybot+0.17,"pat",ha="right",va="center",fontsize=9,color=MUT)
    for i,xx in enumerate([0.7,1.5,2.3,3.1,3.9]):
        for yy,fill in [(ytop+0.17,i%2==0),(ybot+0.17,i%2==1)]:
            ax.add_patch(Circle((x0+xx*0.72,yy),0.055,fc=(INK if fill else "white"),ec=INK,lw=1.1,zorder=4))
    ax.text((x0+x1)/2,6.55,"phased haplotypes  →  per-ancestry dosage",ha="center",fontsize=11,style="italic",color=MUT)
    # --- compress cylinders ---
    def cyl(cx,cy,w,h,fc):
        ax.add_patch(Rectangle((cx-w/2,cy-h/2),w,h,fc=fc,ec="#7a7a7a",lw=1.4,zorder=3))
        ax.add_patch(Ellipse((cx,cy+h/2),w,0.28,fc=fc,ec="#7a7a7a",lw=1.4,zorder=4))
        ax.add_patch(Ellipse((cx,cy-h/2),w,0.28,fc=fc,ec="#7a7a7a",lw=1.4,zorder=2))
    cyl(1.55,3.6,1.5,1.7,"#ECECEC"); ax.text(1.55,2.35,"Raw genotype\nVCF",ha="center",va="center",fontsize=11,fontweight="bold",color=INK)
    cyl(4.15,3.7,0.95,1.0,"#CFE0F2"); ax.text(4.15,2.9,"FELIXla",ha="center",va="center",fontsize=12.5,fontweight="bold",color=BLUE)
    ax.add_patch(FancyArrowPatch((2.5,3.7),(3.5,3.75),arrowstyle="-|>",mutation_scale=22,lw=3,color=BROWN,zorder=5))
    ax.text(3.0,4.15,"compress",ha="center",fontsize=11,fontweight="bold",color=BROWN)
    # --- reruns box ---
    ax.add_patch(FancyBboxPatch((0.9,0.75),3.8,0.9,boxstyle="round,pad=0.02,rounding_size=0.12",
                 fc="white",ec=BROWN,lw=2.2,zorder=3))
    ax.text(2.8,1.2,"single-pass streaming I/O",ha="center",va="center",fontsize=13,fontweight="bold",color=INK)
    save(fig,"fig1_panel2_storage")

# ============================ PANEL 3: score tests ==========================
def panel3():
    W,H=6.2,7.2; fig,ax=plt.subplots(figsize=(W,H)); card(ax,W,H,"3. Two complementary score tests",PURPLE)
    chipy=6.05; xs=[0.95,1.95,2.95,3.95,4.95]
    for (nm,c),x in zip(ANC,xs):
        ax.add_patch(FancyBboxPatch((x-0.4,chipy-0.22),0.8,0.44,boxstyle="round,pad=0.02,rounding_size=0.1",
                     fc=c,ec="none",zorder=3))
        ax.text(x,chipy,nm,ha="center",va="center",fontsize=11,fontweight="bold",color="white",zorder=4)
    # fan lines
    for x in xs:
        ax.plot([x,1.65],[chipy-0.24,4.55],color="#BFD3E8",lw=1.1,zorder=1)
        ax.plot([x,4.45],[chipy-0.24,4.55],color="#EAD9A8",lw=1.1,zorder=1)
    # test boxes (wider gap so the two titles don't collide)
    ax.add_patch(FancyBboxPatch((0.45,3.15),2.4,1.35,boxstyle="round,pad=0.02,rounding_size=0.12",
                 fc="#EAF2FA",ec=BLUE,lw=2.2,zorder=3))
    ax.text(1.65,4.16,"Homogeneous",ha="center",fontsize=12.5,fontweight="bold",color=BLUE)
    ax.text(1.65,3.86,"test (1 d.o.f.)",ha="center",fontsize=11,color=BLUE)
    ax.text(1.65,3.48,r"$H_0:\ \beta_{\mathrm{summed}}=0$",ha="center",fontsize=11.5,color=INK)
    ax.add_patch(FancyBboxPatch((3.25,3.15),2.4,1.35,boxstyle="round,pad=0.02,rounding_size=0.12",
                 fc="#FBF4E0",ec=GOLD,lw=2.2,zorder=3))
    ax.text(4.45,4.16,"Heterogeneous",ha="center",fontsize=12.5,fontweight="bold",color="#9A7A1E")
    ax.text(4.45,3.86,"test (K d.o.f.)",ha="center",fontsize=11,color="#9A7A1E")
    ax.text(4.45,3.48,r"$H_0:\ \beta_{\mathrm{AFR}}=\cdots=\beta_{\mathrm{SAS}}=0$",ha="center",fontsize=10.5,color=INK)
    ax.add_patch(FancyArrowPatch((1.65,3.13),(2.75,2.15),arrowstyle="-|>",mutation_scale=18,lw=2.2,color=BLUE,zorder=4))
    ax.add_patch(FancyArrowPatch((4.45,3.13),(3.35,2.15),arrowstyle="-|>",mutation_scale=18,lw=2.2,color=GOLD,zorder=4))
    ax.add_patch(FancyBboxPatch((1.2,1.15),3.8,0.95,boxstyle="round,pad=0.02,rounding_size=0.12",
                 fc="white",ec=ORANGE,lw=2.4,zorder=3))
    ax.text(3.1,1.62,r"$p_{\mathrm{CCT}}=\mathrm{CCT}(p_{\mathrm{hom}},\,p_{\mathrm{het}})$",
            ha="center",va="center",fontsize=15,fontweight="bold",color=ORANGE,zorder=4)
    save(fig,"fig1_panel3_tests")

# ============================ PANEL 4: outputs ==============================
def panel4():
    W,H=5.9,7.2; fig,ax=plt.subplots(figsize=(W,H)); card(ax,W,H,"4. Method outputs",GREEN)
    ax.text(0.42,6.18,"per-ancestry estimates at every variant",fontsize=12.5,style="italic",color=INK)
    L,R=0.42,5.5; colc={"b":3.05,"se":4.05,"p":4.95}
    yhdr=5.72; rh=0.52
    # header band (navy)
    ax.add_patch(Rectangle((L,yhdr-rh/2),R-L,rh,fc=NAVY,ec="none",zorder=2))
    ax.text(0.56,yhdr,"ancestry",ha="left",va="center",fontsize=12,fontweight="bold",color="white",zorder=3)
    ax.text(colc["b"],yhdr,r"$\hat\beta$",ha="center",va="center",fontsize=15,fontweight="bold",color="white",zorder=3)
    ax.text(colc["se"],yhdr,"SE",ha="center",va="center",fontsize=13,fontweight="bold",color="white",zorder=3)
    ax.text(colc["p"],yhdr,r"$p$",ha="center",va="center",fontsize=15,fontweight="bold",color="white",zorder=3)
    y=yhdr-rh
    for i,(nm,c) in enumerate(ANC):
        if i%2==0: ax.add_patch(Rectangle((L,y-rh/2),R-L,rh,fc="#F2F5FA",ec="none",zorder=1))
        ax.add_patch(Rectangle((0.42,y-0.10),0.20,0.20,fc=c,ec="none",zorder=3))
        ax.text(0.74,y,nm,ha="left",va="center",fontsize=12,fontweight="bold",color=c,zorder=3)
        ax.text(colc["b"],y,r"$\hat\beta_{\mathrm{%s}}$"%nm,ha="center",va="center",fontsize=13,color=INK,zorder=3)
        ax.text(colc["se"],y,"SE",ha="center",va="center",fontsize=11.5,color=MUT,zorder=3)
        ax.text(colc["p"],y,r"$p$",ha="center",va="center",fontsize=12.5,color=MUT,zorder=3)
        y-=rh
    # aggregate ALL row (homogeneous / standard GWAS)
    ax.add_patch(Rectangle((L,y-rh/2),R-L,rh,fc="#E7EDE8",ec="none",zorder=1))
    ax.plot([L,R],[y+rh/2,y+rh/2],color="#8FAF97",lw=1.4,zorder=2)
    ax.text(0.74,y,"ALL",ha="left",va="center",fontsize=12,fontweight="bold",color=GREEN,zorder=3)
    ax.text(colc["b"],y,r"$\hat\beta_{\mathrm{ALL}}$",ha="center",va="center",fontsize=13,fontweight="bold",color=INK,zorder=3)
    ax.text(colc["se"],y,"SE",ha="center",va="center",fontsize=11.5,fontweight="bold",color=INK,zorder=3)
    ax.text(colc["p"],y,r"$p$",ha="center",va="center",fontsize=12.5,fontweight="bold",color=INK,zorder=3)
    # combined box
    yb=y-rh/2-0.98
    ax.add_patch(FancyBboxPatch((0.7,yb),4.6,0.72,boxstyle="round,pad=0.02,rounding_size=0.12",
                 fc="white",ec=ORANGE,lw=2.4,zorder=3))
    ax.text(3.0,yb+0.36,r"+ combined:  $p_{\mathrm{hom}}$   $p_{\mathrm{het}}$   $p_{\mathrm{CCT}}$",
            ha="center",va="center",fontsize=14,fontweight="bold",color=ORANGE,zorder=4)
    # footer: aggregate meaning + auxiliary outputs
    ax.text(0.42,yb-0.30,"ALL row = homogeneous test  (standard GWAS)",fontsize=10.5,style="italic",color=MUT)
    ax.text(0.42,yb-0.58,"also per ancestry:  allele frequency · haplotype count",fontsize=10,color=MUT)
    ax.text(0.42,yb-0.84,"conditioned estimates at local-ancestry loci",fontsize=10,color=MUT)
    save(fig,"fig1_panel4_outputs")

panel2(); panel3(); panel4()
print("done: panels 2,3,4")
