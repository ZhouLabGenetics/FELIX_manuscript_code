#!/usr/bin/env Rscript
# plot_results.R — three figures from summary_<thresh>.tsv (base R, no ggplot dep):
#  1. main_predictors : FELIX vs META vs matched incR2, faceted by ancestry (eur10k LD)
#  2. ld_sensitivity  : felix/meta incR2 with eur10k vs prop10k LD (scatter, y=x line)
#  3. tract_global_local : matched (global AllxAll) vs FELIX tract (local), per ancestry
#   Rscript plot_results.R <results_dir> <thresh>
a <- commandArgs(trailingOnly=TRUE); res <- a[1]; thr <- a[2]
d <- read.table(file.path(res, sprintf("summary_%s.tsv", thr)), header=TRUE, sep="\t",
                stringsAsFactors=FALSE, check.names=FALSE)
ancs <- c("EUR","AFR","EAS","CSA"); ancs <- ancs[ancs %in% d$ancestry]
eb <- function(x, lo, hi, col) arrows(x, lo, x, hi, angle=90, code=3, length=0.02, col=col)

## 1. main predictors (eur10k LD for felix/meta; matched uses its own panel)
pdf(file.path(res, sprintf("fig_main_predictors_%s.pdf", thr)), width=11, height=8)
par(mfrow=c(2,2), mar=c(7,4,3,1))
cols <- c(felix="#C0392B", meta="#2C6FB0", matched="#7A7A7A")
for(an in ancs){
  s <- d[d$ancestry==an,]
  fe <- s[s$predictor=="felix"  & s$ldref=="eur10k",]
  me <- s[s$predictor=="meta"   & s$ldref=="eur10k",]
  ma <- s[grepl("^matched", s$predictor),]
  tr <- sort(unique(c(fe$trait, me$trait, ma$trait)))
  if(!length(tr)){ plot.new(); title(an); next }
  yr <- range(0, s$CI_high[s$predictor %in% c("felix","meta") & s$ldref=="eur10k" | grepl("^matched",s$predictor)], na.rm=TRUE)
  plot(NA, xlim=c(0.5,length(tr)+0.5), ylim=yr, xaxt="n", xlab="", ylab="incremental R2", main=paste("validation:",an))
  axis(1, at=seq_along(tr), labels=tr, las=2, cex.axis=0.7)
  put <- function(df, off, col){ if(!nrow(df)) return(); x<-match(df$trait,tr)+off
    points(x, df$incR2, pch=19, col=col); eb(x, df$CI_low, df$CI_high, col) }
  put(fe,-0.15,cols["felix"]); put(me,0,cols["meta"]); put(ma,0.15,cols["matched"])
  abline(h=0, col="grey80")
  legend("topright", c("FELIX ancALL","AllxAll META","AllxAll matched"), pch=19, col=cols, bty="n", cex=0.8)
}
dev.off()

## 2. LD sensitivity: eur10k (x) vs prop10k (y) for the multi-ancestry predictors
pdf(file.path(res, sprintf("fig_ld_sensitivity_%s.pdf", thr)), width=6, height=6)
w <- reshape(d[d$predictor %in% c("felix","meta") & d$ldref %in% c("eur10k","prop10k"),
               c("trait","ancestry","predictor","ldref","incR2")],
             timevar="ldref", idvar=c("trait","ancestry","predictor"), direction="wide")
lim <- range(0, w$incR2.eur10k, w$incR2.prop10k, na.rm=TRUE)
plot(w$incR2.eur10k, w$incR2.prop10k, xlim=lim, ylim=lim,
     pch=ifelse(w$predictor=="felix",19,1), col=ifelse(w$predictor=="felix","#C0392B","#2C6FB0"),
     xlab="incR2 with 10k-EUR LD", ylab="incR2 with proportional LD",
     main=sprintf("LD-reference sensitivity (%s)", thr))
abline(0,1,col="grey60",lty=2)
legend("topleft", c("FELIX ancALL","AllxAll META"), pch=c(19,1), col=c("#C0392B","#2C6FB0"), bty="n")
dev.off()

## 3. tract global-vs-local: matched (global) vs felix_tract (local), per ancestry
pdf(file.path(res, sprintf("fig_tract_global_local_%s.pdf", thr)), width=11, height=8)
par(mfrow=c(2,2), mar=c(7,4,3,1))
for(an in ancs){
  s <- d[d$ancestry==an,]
  g <- s[grepl("^matched", s$predictor),]
  l <- s[s$predictor==paste0("felix_tract_",an),]
  tr <- sort(unique(c(g$trait, l$trait)))
  if(!length(tr)){ plot.new(); title(an); next }
  yr <- range(0, g$CI_high, l$CI_high, na.rm=TRUE)
  plot(NA, xlim=c(0.5,length(tr)+0.5), ylim=yr, xaxt="n", xlab="", ylab="incremental R2",
       main=paste("validation:",an,"  global(AllxAll) vs local(FELIX tract)"))
  axis(1, at=seq_along(tr), labels=tr, las=2, cex.axis=0.7)
  put <- function(df, off, col){ if(!nrow(df)) return(); x<-match(df$trait,tr)+off
    points(x, df$incR2, pch=19, col=col); eb(x, df$CI_low, df$CI_high, col) }
  put(g,-0.1,"#7A7A7A"); put(l,0.1,"#C77BA6")
  abline(h=0, col="grey80")
  legend("topright", c("global AllxAll","local FELIX tract"), pch=19, col=c("#7A7A7A","#C77BA6"), bty="n", cex=0.8)
}
dev.off()
cat("wrote fig_main_predictors / fig_ld_sensitivity / fig_tract_global_local for", thr, "\n")
