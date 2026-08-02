#!/usr/bin/env Rscript
# plot_violin_bytrait.R — LAYOUT B: x = trait (human-readable), 3 dodged violins of the
# 2000 bootstrap incR2 draws per method (FELIX ancALL / All by All META / All by All
# ancestrally-matched), 4 ancestry panels, SHARED y. Needs the *_draws.tsv.gz from a
# bootstrap run with DUMP=1. Works for P+T (tag=5e-8) and PRS-CS (tag=PRSCS).
#   Rscript plot_violin_bytrait.R <results_dir> <tag> [outstem]
suppressMessages({library(ggplot2)})
a <- commandArgs(trailingOnly=TRUE); RES <- a[1]; TAG <- a[2]
OUT <- ifelse(length(a)>=3, a[3], paste0("fig_violin_bytrait_", TAG))
NAME <- c("3006923"="ALT","3007070"="HDL","3009744"="MCHC","3013721"="AST",
          "3022192"="Triglycerides","3024929"="Platelets","3027114"="Total chol.",
          "3028288"="LDL","3035995"="Alk. phos.","BMI"="BMI","height"="Height")
ANCS <- c("EUR","AFR","EAS","CSA"); TR <- names(NAME)
MCOL <- c("FELIX ancALL"="#AA3377","All by All META"="#0077BB",
          "All by All ancestrally-matched"="#EEE685")
pick <- function(cols, pats){ for(p in pats){ h<-grep(p,cols,value=TRUE); if(length(h)) return(h[1]) }; NA }

rows <- list()
for(t in TR) for(anc in ANCS){
  f <- file.path(RES,t,anc,paste0("incR2_",TAG,"_draws.tsv.gz"))
  if(!file.exists(f)) f <- file.path(RES,t,anc,paste0("incR2_",TAG,"_eur10k_draws.tsv.gz"))
  if(!file.exists(f)) next
  d <- read.table(gzfile(f), header=TRUE, sep="\t", check.names=FALSE)
  cf <- pick(colnames(d), c("^felix__eur10k$","^felix$"))
  cm <- pick(colnames(d), c("^meta__eur10k$","^meta$"))
  cx <- pick(colnames(d), c("^matched_"))
  for(m in list(c(cf,"FELIX ancALL"), c(cm,"All by All META"), c(cx,"All by All ancestrally-matched"))){
    if(is.na(m[1])) next
    rows[[length(rows)+1]] <- data.frame(trait=NAME[[t]], ancestry=anc, method=m[2], val=d[[m[1]]]*100)
  }
}
if(!length(rows)) stop("no *_draws.tsv.gz found under ", RES, " for tag ", TAG,
                       " -- re-run bootstrap with DUMP=1")
df <- do.call(rbind, rows)
df$ancestry <- factor(df$ancestry, levels=ANCS)
df$method   <- factor(df$method,   levels=names(MCOL))
df$trait    <- factor(df$trait,    levels=unname(NAME))

p <- ggplot(df, aes(trait, val, fill=method)) +
  geom_violin(position=position_dodge(0.8), width=0.82, linewidth=0.3, scale="width") +
  scale_fill_manual(values=MCOL, name=NULL) +
  facet_wrap(~ ancestry, ncol=2, scales="fixed") +          # SHARED y across all 4 panels
  labs(x=NULL, y="Incremental R² (%)") +
  theme_minimal(base_size=20) +
  theme(legend.position="top",
        legend.text=element_text(size=18),
        panel.grid.major.x=element_blank(),
        axis.title.y=element_text(size=22),
        axis.text.y=element_text(size=17),
        axis.text.x=element_text(angle=40, hjust=1, size=16),
        strip.text=element_text(face="bold", size=22))
ggsave(paste0(OUT,".pdf"), p, width=17, height=11)
ggsave(paste0(OUT,".png"), p, width=17, height=11, dpi=300)
cat("wrote", paste0(OUT,".pdf/.png"), "\n")
