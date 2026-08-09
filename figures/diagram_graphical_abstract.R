#!/usr/bin/env Rscript
# Figure 1 components rebuilt in R style
source("scripts_R/00_theme.R"); source("scripts_R/00_people.R"); source("scripts_R/01_data.R")
suppressMessages(library(patchwork))
A <- ANC_COLORS
txt <- function(x,y,l,s=5,f="plain",c="black",h=0.5) annotate("text",x=x,y=y,label=l,size=s,fontface=f,colour=c,hjust=h)

## A — concept: one global label vs local ancestry per haplotype
seg <- function(y, cols) { n<-length(cols); w<-6/n
  lapply(seq_len(n), function(k) annotate("rect", xmin=2+(k-1)*w, xmax=2+k*w, ymin=y-0.35, ymax=y+0.35,
    fill=cols[k], colour="white", linewidth=0.5)) }
pA <- ggplot() + coord_fixed(xlim=c(0,10), ylim=c(0,10), expand=FALSE) + theme_void(base_size=18) +
  person_layer(1, 6.5, 3.2, c(A[["AFR"]],A[["EUR"]],A[["NatAm"]])) +
  txt(2.2, 8.9, "global ancestry: one label", 4.4, "bold", "#555", 0) +
  txt(2.2, 3.2, "local ancestry: per haplotype", 4.4, "bold", "#555", 0)
for (L in seg(2.2, c(A[["EUR"]],A[["EUR"]],A[["AFR"]],A[["NatAm"]]))) pA <- pA + L
for (L in seg(1.2, c(A[["AFR"]],A[["EUR"]],A[["EUR"]],A[["EUR"]]))) pA <- pA + L
pA <- pA + ggtitle("A  Haplotype-resolved ancestry")

## B — method: HOM + HET -> CCT
tb <- function(x0,x1,y0,y1,col,lab,sub) list(
  annotate("rect",xmin=x0,xmax=x1,ymin=y0,ymax=y1,fill=scales::alpha(col,0.3),colour=col,linewidth=1.3),
  txt((x0+x1)/2,(y0+y1)/2+0.5,lab,5,"bold",col), txt((x0+x1)/2,(y0+y1)/2-0.5,sub,3.4,"plain","#333"))
pB <- ggplot() + coord_fixed(xlim=c(0,10), ylim=c(0,10), expand=FALSE) + theme_void(base_size=18)
for (L in tb(0.5,4.4,6,8.4, TEST_COLORS[["HOM"]], "HOM","1 df, summed")) pB <- pB + L
for (L in tb(5.6,9.5,6,8.4, "#9a7b16", "HET","K df, per-ancestry")) pB <- pB + L
for (L in tb(2.8,7.2,2.4,4.4, TEST_COLORS[["CCT"]], "Cauchy (CCT)","adaptive power")) pB <- pB + L
pB <- pB +
  annotate("segment", x=2.4,xend=4.0,y=6,yend=4.4, colour=TEST_COLORS[["HOM"]], linewidth=1.1, arrow=arrow(length=unit(0.2,"cm"),type="closed")) +
  annotate("segment", x=7.5,xend=6.0,y=6,yend=4.4, colour="#9a7b16", linewidth=1.1, arrow=arrow(length=unit(0.2,"cm"),type="closed")) +
  ggtitle("B  Adaptive combined test")

## C — resolution: tighter SE at shared loci
dt <- load_scatter(); d <- dt[tophit_source=="ABA"]
d[, `:=`(x=as.numeric(ABA_META_SE), y=as.numeric(SAIGE_SE_c_ancALL))]
d <- d[is.finite(x)&is.finite(y)&x>0&y>0]
pC <- ggplot(d, aes(x,y)) + geom_abline(slope=1, linetype="dashed", color="grey45") +
  geom_point(size=2, alpha=0.6, color="#AA3377") + scale_x_log10() + scale_y_log10() + coord_equal() +
  labs(x="All by All SE", y="FELIXassoc SE") + ggtitle("C  Tighter effect estimates")

## D — worked example: HPR opposing per-ancestry effects
r <- get_locus(dt, "HPR", "3028288", 72080103)
dd <- data.table(anc=ANC_MAP$name,
  b=sapply(paste0("SAIGE_BETA_c_",ANC_MAP$suf), function(c) as.numeric(r[[c]])),
  se=sapply(paste0("SAIGE_SE_c_",ANC_MAP$suf), function(c) as.numeric(r[[c]])))
dd <- dd[is.finite(b)]; dd[, anc:=factor(anc, levels=rev(ANC_MAP$name))]
pD <- ggplot(dd, aes(b, anc, color=anc)) + geom_vline(xintercept=0, color="grey55") +
  geom_errorbarh(aes(xmin=b-1.96*se, xmax=b+1.96*se), height=0, linewidth=1.2) +
  geom_point(size=PTBIG) + scale_color_anc(guide="none") +
  labs(x=expression("per-ancestry effect  "*beta), y=NULL) +
  ggtitle("D  Heterogeneity recovered (HPR)")

ga <- (pA | pB) / (pC | pD) + plot_annotation(theme = theme(plot.title = element_text(size=22, face="bold")))
save_fig(ga, "fig1_graphical_abstract", width = 15, height = 12)
