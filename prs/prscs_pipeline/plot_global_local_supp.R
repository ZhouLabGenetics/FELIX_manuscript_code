#!/usr/bin/env Rscript
# plot_global_local_supp.R — SUPPLEMENTARY figure: FELIXassoc ancestry-resolved (local) vs
# All by All ancestry-matched (global) incremental R2, EUR and AFR targets, per trait.
# Reads summary_<tag>.tsv (uses felix_tract_<ANC> = local, matched_<ANC> = global).
#   Rscript plot_global_local_supp.R <results_dir> <tag> [outstem]
suppressMessages({library(ggplot2)})
a <- commandArgs(trailingOnly=TRUE); RES<-a[1]; TAG<-a[2]
OUT <- ifelse(length(a)>=3, a[3], paste0("fig_supp_global_local_",TAG))
NAME <- c("3006923"="ALT","3007070"="HDL","3009744"="MCHC","3013721"="AST",
          "3022192"="Triglycerides","3024929"="Platelets","3027114"="Total chol.",
          "3028288"="LDL","3035995"="Alk. phos.","BMI"="BMI","height"="Height")
d <- read.table(file.path(RES,paste0("summary_",TAG,".tsv")), header=TRUE, sep="\t",
                check.names=FALSE, stringsAsFactors=FALSE)
rows <- list()
for(anc in c("EUR","AFR")){
  g <- d[d$ancestry==anc & d$predictor==paste0("matched_",anc), ]
  l <- d[d$ancestry==anc & d$predictor==paste0("felix_tract_",anc), ]
  for(x in list(list(g,"Global (All by All)"), list(l,"Local (FELIXassoc)"))){
    dd<-x[[1]]; if(!nrow(dd)) next
    rows[[length(rows)+1]] <- data.frame(anc=anc, trait=NAME[dd$trait], grp=x[[2]],
                                          incR2=dd$incR2*100, lo=dd$CI_low*100, hi=dd$CI_high*100)
  }
}
df <- do.call(rbind, rows)
df$anc   <- factor(df$anc, levels=c("EUR","AFR"))
df$grp   <- factor(df$grp, levels=c("Global (All by All)","Local (FELIXassoc)"))
df$trait <- factor(df$trait, levels=unname(NAME))
pal <- c("Global (All by All)"="#6699C4","Local (FELIXassoc)"="#C77BA6")
p <- ggplot(df, aes(trait, incR2, fill=grp)) +
  geom_col(position=position_dodge(0.8), width=0.72) +
  geom_errorbar(aes(ymin=lo, ymax=hi), position=position_dodge(0.8), width=0.25, linewidth=0.4) +
  scale_fill_manual(values=pal, name=NULL) +
  facet_wrap(~ anc, ncol=1, scales="free_y") +
  labs(x=NULL, y="Incremental R² (%)") +
  theme_minimal(base_size=19) +
  theme(legend.position="top", legend.text=element_text(size=17),
        axis.title.y=element_text(size=21), axis.text.y=element_text(size=15),
        axis.text.x=element_text(angle=40, hjust=1, size=15),
        strip.text=element_text(face="bold", size=21), panel.grid.major.x=element_blank())
ggsave(paste0(OUT,".pdf"), p, width=13, height=11)
ggsave(paste0(OUT,".png"), p, width=13, height=11, dpi=300)
cat("wrote", paste0(OUT,".pdf/.png"), "\n")
