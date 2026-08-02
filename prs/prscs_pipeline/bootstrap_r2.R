#!/usr/bin/env Rscript
# bootstrap_r2.R
# R2 (or covariate-adjusted incremental R2) with bootstrap 95% CI for one or more
# PRS score files vs a phenotype, plus the paired R2 difference vs a baseline.
#
# Plain R2:            R2 = cor(score, y)^2   (== summary(lm(y~score))$r.squared)
# Incremental R2:      R2(y ~ covars + score) - R2(y ~ covars)   [use --covar ...]
# Bootstrap: resample individuals with replacement B times (paired across all
# score files) -> percentile 95% CI, and CI + p for each (score - baseline).
#
# Usage (raw R2, AFR, restricted to the Pan-UKB AFR keep set):
#   Rscript bootstrap_r2.R \
#     --pheno=subset_ukbb_biomarkers.tsv --pheno_id=eid --pheno_col=TotalCholesterol \
#     --map=sample_id_mapping.txt --keep=shared/afr_keep.txt \
#     --baseline=global --B=2000 --seed=1 --out=3027114/results/r2_afr_boot.tsv \
#     global=3027114/AFR/scores/global_prs.sscore \
#     local=3027114/AFR/scores/local_prs.sscore \
#     partial=3027114/AFR/scores/partial_prs.mapped.sscore
#
# Incremental R2: add
#     --covar=covars_by_eid.tsv --covar_id=eid \
#     --covar_cols=PC1,PC2,...,PC20,age,sex,age_sex,age2,age2_sex
#   The covariate file MUST be joinable to the phenotype id (eid) BEFORE the
#   eid->sscore map. See the runbook: the final_samples.txt.bgz covariates are in
#   a separate Pan-UKB id space and must first be bridged to eid.

args <- commandArgs(trailingOnly = TRUE)
getf <- function(k, d=NULL){ h<-grep(paste0("^--",k,"="),args,value=TRUE)
  if(length(h)) sub(paste0("^--",k,"="),"",h[1]) else d }
pheno_f   <- getf("pheno");     stopifnot(!is.null(pheno_f))
pheno_id  <- getf("pheno_id","eid")
pheno_col <- getf("pheno_col"); stopifnot(!is.null(pheno_col))
out_f     <- getf("out","r2_boot.tsv")
map_f     <- getf("map"); keep_f <- getf("keep")
covar_f   <- getf("covar"); covar_id <- getf("covar_id","eid")
covar_cols<- getf("covar_cols")
baseline  <- getf("baseline","global")
B         <- as.integer(getf("B","2000")); seed <- as.integer(getf("seed","1"))
dump_boot <- getf("dump_boot")   # if set, also write the B x score bootstrap draws (for violins)
score_args<- args[grepl("=",args) & !grepl("^--",args)]
if(!length(score_args)) stop("Provide >=1 score file as name=path")
score_files <- setNames(sub("^[^=]+=","",score_args), sub("=.*$","",score_args))
set.seed(seed)

read_sscore <- function(f){ d<-read.table(f,header=TRUE,comment.char="",check.names=FALSE,sep="\t")
  sc<-grep("SCORE.*SUM$",colnames(d),value=TRUE)[1]; if(is.na(sc)) sc<-tail(colnames(d),1)
  data.frame(IID=as.character(d$IID), score=d[[sc]]) }
inv_norm <- function(x) qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))

# ---- phenotype (keyed by eid) ----
ph <- read.table(pheno_f,header=TRUE,sep="\t",check.names=FALSE)
stopifnot(pheno_id %in% colnames(ph), pheno_col %in% colnames(ph))
dat <- data.frame(IID=as.character(ph[[pheno_id]]), y=inv_norm(ph[[pheno_col]]))

# ---- optional covariates (join by eid, BEFORE the sscore map) ----
COV <- NULL
if(!is.null(covar_f)){
  stopifnot(!is.null(covar_cols))
  cc <- strsplit(covar_cols,",")[[1]]
  cv <- read.table(covar_f,header=TRUE,sep="\t",check.names=FALSE)
  miss <- setdiff(c(covar_id,cc), colnames(cv)); if(length(miss)) stop("covar cols missing: ",paste(miss,collapse=","))
  cv <- cv[,c(covar_id,cc)]; colnames(cv)[1] <- "IID"; cv$IID <- as.character(cv$IID)
  ov <- length(intersect(dat$IID, cv$IID))
  cat(sprintf("covariate join: %d of %d phenotype ids matched the covariate file\n", ov, nrow(dat)))
  if(ov==0) stop("covariate ids do not overlap the phenotype ids -- bridge them first (see runbook).")
  dat <- merge(dat, cv, by="IID"); COV <- cc
}

# ---- eid -> sscore id + keep filter ----
if(!is.null(map_f)){ m<-read.table(map_f,header=FALSE); m<-data.frame(ukb=as.character(m$V1),sid=as.character(m$V2))
  dat<-merge(dat,m,by.x="IID",by.y="ukb"); dat$IID<-dat$sid; dat$sid<-NULL
  if(!is.null(keep_f)){ k<-read.table(keep_f,header=FALSE); keep_sid<-m$sid[m$ukb %in% as.character(k$V2)] }
} else if(!is.null(keep_f)){ k<-read.table(keep_f,header=FALSE); keep_sid<-as.character(k$V2) }

# ---- merge score files onto the SAME individuals ----
for(nm in names(score_files)){ s<-read_sscore(score_files[[nm]]); colnames(s)[2]<-nm; dat<-merge(dat,s,by="IID") }
if(exists("keep_sid")) dat<-dat[dat$IID %in% keep_sid,]
dat<-dat[complete.cases(dat),]; n<-nrow(dat)
cat(sprintf("N = %d individuals; scores: %s%s\n", n, paste(names(score_files),collapse=", "),
            if(is.null(COV)) "" else sprintf("; adjusting for %d covariates (incremental R2)", length(COV))))
if(n<20) stop("too few overlapping individuals")
nms <- names(score_files)

# ---- R2 estimator (plain or incremental); base model fit ONCE per resample ----
if(!is.null(COV)) Xc <- as.matrix(cbind(1, dat[,COV,drop=FALSE]))
metric_all <- function(rows){
  y <- dat$y[rows]
  if(is.null(COV)){
    vapply(nms, function(nm) cor(dat[[nm]][rows], y)^2, numeric(1))
  } else {
    tss <- sum((y-mean(y))^2); Xb <- Xc[rows,,drop=FALSE]
    r2b <- 1 - sum(.lm.fit(Xb, y)$residuals^2)/tss            # covariate-only, once
    vapply(nms, function(nm)
      (1 - sum(.lm.fit(cbind(Xb, dat[[nm]][rows]), y)$residuals^2)/tss) - r2b, numeric(1))
  }
}
all_rows <- seq_len(n)
pt <- metric_all(all_rows); names(pt) <- nms

# ---- bootstrap ----
boot <- matrix(NA_real_, B, length(nms), dimnames=list(NULL,nms))
for(b in seq_len(B)) boot[b,] <- metric_all(sample.int(n, n, replace=TRUE))
ci <- function(v) quantile(v,c(0.025,0.975),names=FALSE)
metric <- if(is.null(COV)) "R2" else "incR2"

r2_tab <- data.frame(score=nms, N=n, estimate=pt[nms],
  CI_low=apply(boot,2,function(v)ci(v)[1])[nms], CI_high=apply(boot,2,function(v)ci(v)[2])[nms],
  metric=metric, row.names=NULL)

diff_tab <- NULL
if(baseline %in% nms){
  diff_tab <- do.call(rbind, lapply(setdiff(nms,baseline), function(nm){
    d<-boot[,nm]-boot[,baseline]
    data.frame(comparison=sprintf("%s - %s",nm,baseline), d_estimate=pt[nm]-pt[baseline],
               CI_low=ci(d)[1], CI_high=ci(d)[2], p_boot=mean(d<=0), row.names=NULL) })) }

cat(sprintf("\n== %s (bootstrap 95%% CI) ==\n", metric)); print(r2_tab,row.names=FALSE)
if(!is.null(diff_tab)){ cat("\n== difference vs baseline ==\n"); print(diff_tab,row.names=FALSE) }
dir.create(dirname(out_f),recursive=TRUE,showWarnings=FALSE)
write.table(r2_tab,out_f,sep="\t",quote=FALSE,row.names=FALSE)
if(!is.null(diff_tab)) write.table(diff_tab,sub("\\.tsv$","_diff.tsv",out_f),sep="\t",quote=FALSE,row.names=FALSE)
if(!is.null(dump_boot)){ df_f <- sub("\\.tsv$","_draws.tsv.gz",out_f)
  write.table(as.data.frame(boot), gzfile(df_f), sep="\t", quote=FALSE, row.names=FALSE)
  cat("Wrote draws:",df_f,"\n") }
cat("\nWrote:",out_f,"\n")
