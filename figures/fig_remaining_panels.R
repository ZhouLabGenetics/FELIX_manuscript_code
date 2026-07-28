#!/usr/bin/env Rscript
# Remaining 1:1 panels: per-locus mechanism forests (fig14a-21), overview_strong,
# non-LD locus-zooms (14/15), sample-size bars (15), het forests (14 C1/C2),
# ADRB2 inclusion (14 A2), and graphical-abstract panel_d (12). Shared helpers; new palette.
source("scripts_R/00_theme.R"); source("scripts_R/01_data.R")
dt <- load_scatter()
ssL <- fread(file.path(REPLIC, "sample_size_table_long.tsv"))

## ---- shared helpers ----
forest_locus <- function(gene, ph, pos, title, file, w = 9, h = 5.6) {
  r <- get_locus(dt, gene, ph, pos); if (is.null(r)) { message("miss forest ", gene); return(invisible()) }
  d <- rbindlist(lapply(1:5, function(i) { suf<-ANC_MAP$suf[i]; nm<-ANC_MAP$name[i]; aba<-ANC_MAP$aba[i]
    data.table(anc=nm, st=as.numeric(r[[paste0("SAIGE_BETA_c_",suf)]]), sse=as.numeric(r[[paste0("SAIGE_SE_c_",suf)]]),
               ab=as.numeric(r[[paste0(aba,"_BETA")]]), abse=as.numeric(r[[paste0(aba,"_SE")]])) }))
  d <- d[!is.na(st)]; if (!nrow(d)) { message("empty forest ", gene); return(invisible()) }
  d[, anc := factor(anc, levels = d[order(st)]$anc)]
  L <- rbind(d[, .(anc, method="FELIXassoc", b=st, se=sse)], d[, .(anc, method="All by All", b=ab, se=abse)])
  p <- ggplot(L, aes(b, anc, color=anc, shape=method)) + geom_vline(xintercept=0, color="grey60") +
    geom_errorbarh(aes(xmin=b-1.96*se, xmax=b+1.96*se), height=0, linewidth=1.1, position=position_dodge(0.5)) +
    geom_point(size=PTBIG, fill="white", stroke=1.2, position=position_dodge(0.5)) +
    scale_shape_manual(values=c("FELIXassoc"=16, "All by All"=22), name=NULL) + scale_color_anc(guide="none") +
    labs(x=expression("per-ancestry effect ("*beta*", 95% CI)"), y=NULL, title=title) + theme(legend.position="top")
  save_fig(p, file, width=w, height=h)
}
lz <- function(sp, ap, chr, pos, title, file) {
  sg <- saige_region(sp, chr, pos); sg[, y := mlog(P_cct_admixed_c)]; sg <- sg[is.finite(y)]
  if (!nrow(sg)) { message("miss lz ", sp); return(invisible()) }
  ab <- aba_region(ap, "META", chr, pos); if (nrow(ab)) ab[, y := mlog(Pvalue)]
  lead <- sg[which.max(y)]
  both <- rbind(sg[, .(POS, y, method="FELIXassoc")],
                if (nrow(ab)) ab[, .(POS, y, method="All by All")] else NULL)
  p <- ggplot(both, aes(POS/1e6, y, color=method)) +
    geom_point(size=2.2, alpha=0.8) +
    geom_point(data=lead, aes(POS/1e6, y), shape=18, size=8, color="#EE3377", inherit.aes=FALSE) +
    geom_hline(yintercept=GW_LOG, linetype="dashed", color="grey40") +
    scale_color_manual(values=METHOD_COLORS, name=NULL) +
    labs(x=paste0("Chromosome ", chr, " position (Mb)"), y=expression("Association  "*-log[10](p)), title=title) +
    theme(legend.position="top")
  save_fig(p, file, width=9.5, height=5.6)
}
ss_bars <- function(pheno, title, file) {
  d <- ssL[phenotype == pheno & ancestry %in% c("AFR","EAS","EUR","NatAm","SAS")]
  if (!nrow(d)) { message("miss ss ", pheno); return(invisible()) }
  m <- melt(d[, .(ancestry, `All by All`=as.numeric(ABA_N), `FELIXassoc`=as.numeric(TRACTOR_N))],
            id.vars="ancestry", variable.name="method", value.name="N")
  m[, ancestry := factor(ancestry, levels=c("AFR","EAS","EUR","NatAm","SAS"))]
  p <- ggplot(m[!is.na(N) & N > 0], aes(ancestry, N, fill=method)) +
    geom_col(position=position_dodge(0.75), width=0.7) +
    scale_fill_manual(values=METHOD_COLORS, name=NULL) +
    scale_y_continuous(labels=scales::label_number(scale_cut=scales::cut_short_scale())) +
    labs(x=NULL, y="Effective sample size", title=title) + theme(legend.position="top")
  save_fig(p, file, width=8, height=5.6)
}

## ---- per-locus mechanism forests (script 11: fig14a-21) ----
ML <- list(
  c("VEGFA","3022192",6,43791136,"fig14a_mechanism_VEGFA_triglycerides","VEGFA - triglycerides"),
  c("IL23R","GI_522.11",1,67240275,"fig14b_mechanism_IL23R_crohns","IL23R - Crohn's disease"),
  c("RCOR1","3024929",14,102663104,"fig14c_mechanism_RCOR1_platelet","RCOR1 - platelet count"),
  c("HPR","3028288",16,72080103,"fig15_mechanism_HPR_LDL","HPR/TXNL4B - LDL cholesterol"),
  c("HLA-DRA","NS_326.1",6,32445768,"fig16_mechanism_HLA_DRA_MS","HLA-DRA - multiple sclerosis"),
  c("OR51B5","282.5",11,5498240,"fig17_mechanism_OR51B5_SCA","OR51B5/HBB - sickle-cell anemia"),
  c("HFE","3009744",6,26092913,"fig18_mechanism_HFE_MCHC","HFE - MCHC"),
  c("PNPLA3","3013721",22,43928850,"fig19_mechanism_PNPLA3_AST","PNPLA3 I148M - AST"),
  c("APOC1","3035995",19,44919689,"fig20_mechanism_APOC1_ALP","APOC1 - alkaline phosphatase"),
  c("CD36","3007070",7,80671133,"fig21_mechanism_CD36_HDL","CD36 - HDL cholesterol"))
for (m in ML) forest_locus(m[1], m[2], as.numeric(m[4]), m[6], m[5])

## ---- overview_strong (6-locus facet forest) ----
strong <- data.table(gene=c("HFE","IL23R","PNPLA3","APOC1","HPR","CD36"),
  ph=c("3009744","GI_522.11","3013721","3035995","3028288","3007070"),
  pos=c(26092913,67240275,43928850,44919689,72080103,80671133),
  label=c("HFE (MCHC)","IL23R (Crohn's)","PNPLA3 (AST)","APOC1 (ALP)","HPR (LDL)","CD36 (HDL)"))
ds <- rbindlist(lapply(seq_len(nrow(strong)), function(i){ r<-get_locus(dt,strong$gene[i],strong$ph[i],strong$pos[i])
  if(is.null(r)) return(NULL)
  data.table(label=strong$label[i], anc=ANC_MAP$name,
             b=sapply(paste0("SAIGE_BETA_c_",ANC_MAP$suf), function(c) as.numeric(r[[c]])),
             se=sapply(paste0("SAIGE_SE_c_",ANC_MAP$suf), function(c) as.numeric(r[[c]]))) }))
ds <- ds[is.finite(b)]; ds[, `:=`(anc=factor(anc, levels=rev(ANC_MAP$name)), label=factor(label, levels=strong$label))]
pov <- ggplot(ds, aes(b, anc, color=anc)) + geom_vline(xintercept=0, color="grey55") +
  geom_errorbarh(aes(xmin=b-1.96*se, xmax=b+1.96*se), height=0, linewidth=1.1) + geom_point(size=PT) +
  facet_wrap(~label, ncol=3, scales="free_x") + scale_color_anc(guide="none") +
  labs(x=expression("per-ancestry effect  "*beta*" (95% CI)"), y=NULL, title="Mechanism overview (strong examples)") +
  theme(strip.background=element_blank(), strip.text=element_text(size=17, face="bold"))
save_fig(pov, "fig_mechanism_overview_strong", width=15, height=9)

## ---- non-LD locus-zooms (14 A1, 15) ----
lz("pheno_GI_522.11","GI_522.11",1,67240275,"IL23R - Crohn's disease","fig4A1_il23r_locuszoom")
lz("BMI","BMI",5,148898672,"ADRB2/SH3TC2 - BMI","fig4_lz_adrb2")
lz("pheno_RE_475","RE_475",17,39765489,"IKZF3 - asthma","fig4_lz_ikzf3")
lz("height","height",6,44831920,"SUPT3H - height","fig4_lz_supt3h")

## ---- sample-size bars (15) ----
ss_bars("GI_522.11","IL23R - Crohn's disease","fig4_ss_il23r")
ss_bars("BMI","ADRB2 - BMI","fig4_ss_adrb2")
ss_bars("RE_475","IKZF3 - asthma","fig4_ss_ikzf3")
ss_bars("height","SUPT3H - height","fig4_ss_supt3h")

## ---- het forests (14 C1/C2), ADRB2 inclusion (14 A2), panel_d (12) ----
forest_locus("APOC1","3035995",44919689,"APOC1 - alkaline phosphatase","fig4C1_apoc1_heterogeneity")
forest_locus("HPR","3028288",72080103,"HPR/TXNL4B - LDL cholesterol","fig4C2_hpr_heterogeneity")
forest_locus("ADRB2","BMI",148898672,"ADRB2/SH3TC2 - BMI (inclusion)","fig4A2_adrb2_inclusion_AFR")
forest_locus("HPR","3028288",72080103,"HPR/TXNL4B - LDL cholesterol","panel_d_hpr")
lz("pheno_GI_522.11","GI_522.11",1,67240275,"IL23R - Crohn's disease","panel_d_il23r")
