#!/usr/bin/env Rscript
# plot_violin_bymethod.R — LAYOUT A: x = 3 methods (FELIX ancALL / All by All META / All by
# All ancestrally-matched); each violin = distribution of incR2 across the 11 traits; dots =
# individual traits; 4 ancestry panels, SHARED y. Reads summary_<tag>.tsv (NO re-run).
# Works for P+T (tag=5e-8) and PRS-CS (tag=PRSCS).
#   Rscript plot_violin_bymethod.R <results_dir> <tag> [outstem]
suppressMessages({library(ggplot2)})
a <- commandArgs(trailingOnly=TRUE); RES <- a[1]; TAG <- a[2]
OUT <- ifelse(length(a)>=3, a[3], paste0("fig_violin_bymethod_", TAG))
ANCS <- c("EUR","AFR","EAS","CSA")
MCOL <- c("FELIX ancALL"="#AA3377","All by All META"="#0077BB",
          "All by All ancestrally-matched"="#EEE685")

d <- read.table(file.path(RES, paste0("summary_",TAG,".tsv")), header=TRUE, sep="\t",
                check.names=FALSE, stringsAsFactors=FALSE)
d$ldref <- as.character(d$ldref); d$ldref[is.na(d$ldref)] <- ""   # PRS-CS has empty ldref -> read as NA
d$method <- NA
d$method[d$predictor=="felix" & d$ldref %in% c("eur10k","")] <- "FELIX ancALL"
d$method[d$predictor=="meta"  & d$ldref %in% c("eur10k","")] <- "All by All META"
d$method[grepl("^matched", d$predictor)]                     <- "All by All ancestrally-matched"
d <- d[!is.na(d$method), ]
d$incR2 <- d$incR2*100
d$ancestry <- factor(d$ancestry, levels=ANCS)
d$method   <- factor(d$method,   levels=names(MCOL))

p <- ggplot(d, aes(method, incR2, fill=method)) +
  geom_violin(width=0.85, scale="width", linewidth=0.3, alpha=0.9) +
  geom_jitter(width=0.12, height=0, size=1.9, alpha=0.8, color="grey20") +
  scale_fill_manual(values=MCOL, name=NULL) +
  facet_wrap(~ ancestry, ncol=2, scales="fixed") +          # SHARED y across all 4 panels
  labs(x=NULL, y="Incremental R² (%)") +
  theme_minimal(base_size=20) +
  theme(legend.position="top",
        legend.text=element_text(size=18),
        axis.title.y=element_text(size=22),
        axis.text.y=element_text(size=17),
        axis.text.x=element_blank(),
        panel.grid.major.x=element_blank(),
        strip.text=element_text(face="bold", size=22))
ggsave(paste0(OUT,".pdf"), p, width=12, height=10)
ggsave(paste0(OUT,".png"), p, width=12, height=10, dpi=300)
cat("wrote", paste0(OUT,".pdf/.png"), "\n")
