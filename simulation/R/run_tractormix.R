## Run Tractor-Mix on one (mode, trait, scen, beta, seed) condition.
## Fits glmmkin null model on the existing ALT pheno, then runs
## TractorMix.score against the 3-ancestry hapdose files.
##
## Output:
##   data/<mode>/tractormix_out/alt/<TAG>/result.tsv      -- per-variant pvals
##   data/<mode>/tractormix_out/alt/<TAG>/null.rda        -- fitted glmmkin obj
##
## bench.log (wall time, peak RSS) is captured by the SLURM wrapper using
## /usr/bin/time -v, not by this script.
##
## Usage:  Rscript run_tractormix.R <mode> <trait> <scen> <beta> <seed>
##         e.g.  Rscript run_tractormix.R common quant shared 0.30 1

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5)
  stop("Usage: run_tractormix.R <mode> <trait> <scen> <beta> <seed>")
mode    <- args[1]
trait   <- args[2]
scen    <- args[3]
beta    <- as.numeric(args[4])
seed    <- as.integer(args[5])

source(file.path(Sys.getenv("FELIX_SIM_BASE", unset = stop("Set FELIX_SIM_BASE to the simulation directory before running this script")), "R", "config.R"))
set_mode(mode)
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(GMMAT)
})

TM_DIR <- Sys.getenv("TRACTORMIX_DIR", "/opt/Tractor-Mix")
source(file.path(TM_DIR, "TractorMix.score.R"))

## ---- Resolve paths -----------------------------------------------------
beta_tag <- sprintf("beta%03d", round(beta * 100))
TAG  <- sprintf("%s_%s_%s_seed%02d", trait, scen, beta_tag, seed)
pheno_file <- file.path(DATA_DIR, "pheno", "alt",
                        sprintf("%s_%s_%s_seed%02d.tsv",
                                trait, scen, beta_tag, seed))
if (!file.exists(pheno_file)) stop("Pheno not found: ", pheno_file)

hapdose <- file.path(DATA_DIR, "tractormix", "hapdose",
                     paste0("chr1.anc", 0:2, ".dosage.txt"))
stopifnot(all(file.exists(hapdose)))

out_dir <- file.path(DATA_DIR, "tractormix_out", "alt", TAG)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
result_file <- file.path(out_dir, "result.tsv")
null_rda    <- file.path(out_dir, "null.rda")

if (file.exists(result_file) && file.size(result_file) > 0) {
  cat("Already done: ", result_file, "  -- skipping.\n", sep = "")
  quit(save = "no", status = 0)
}

cat("MODE=", mode, "  TRAIT=", trait, "  SCEN=", scen,
    "  BETA=", beta, "  SEED=", seed, "  TAG=", TAG, "\n", sep = "")

## ---- Load pheno + kinship ---------------------------------------------
pheno <- fread(pheno_file, sep = "\t", header = TRUE)
stopifnot(all(c("IID", "PHENO", "AFR", "NAT") %in% names(pheno)))

K_dt <- fread(file.path(DATA_DIR, "Kinship.tsv"), sep = "\t", header = TRUE)
K    <- as.matrix(K_dt[, -1, with = FALSE])
rownames(K) <- colnames(K) <- K_dt$ID
## Align pheno row order with K row order.
pheno <- pheno[match(rownames(K), as.character(pheno$IID))]
stopifnot(identical(as.character(pheno$IID), rownames(K)))

## ---- Fit null model ----------------------------------------------------
## TractorMix.score does NOT inspect obj$family. It inspects
## as.character(obj$call) and greps for the literal strings "gaussian" or
## "binomial". That means we must pass family = gaussian() / binomial()
## LITERALLY in the call -- a variable like `family = fam` captures the
## variable name "fam" in $call and Tractor-Mix rejects it.
t0 <- Sys.time()
cat("[", format(t0), "] glmmkin null fit...\n")
if (trait == "quant") {
  nullmod <- glmmkin(fixed  = PHENO ~ AFR + NAT,
                     data   = as.data.frame(pheno),
                     id     = "IID",
                     kins   = K,
                     family = gaussian())
} else {
  nullmod <- glmmkin(fixed  = PHENO ~ AFR + NAT,
                     data   = as.data.frame(pheno),
                     id     = "IID",
                     kins   = K,
                     family = binomial())
}
t1 <- Sys.time()
cat("[", format(t1), "] null fit done (",
    round(as.numeric(t1 - t0, units = "secs"), 1), "s)\n", sep = "")

## Diagnostic: see what GMMAT actually stored in $call. Some versions of
## GMMAT normalise the call (reassign $call$family to the evaluated family
## object) which destroys the literal text TractorMix.score's grepl needs.
expected_token <- if (trait == "quant") "gaussian" else "binomial"
raw_call <- paste(as.character(nullmod$call), collapse = " | ")
cat("Raw nullmod$call: ", raw_call, "\n", sep = "")

## Bulletproof fix: overwrite $call with a hand-constructed call object that
## *definitely* contains the family token. TractorMix.score only does
## grepl(token, as.character(obj$call)), so the substring is what matters --
## the call object does not need to actually be re-runnable, just searchable.
nullmod$call <- if (trait == "quant") {
  quote(glmmkin(fixed = PHENO ~ AFR + NAT, family = gaussian()))
} else {
  quote(glmmkin(fixed = PHENO ~ AFR + NAT, family = binomial()))
}
forced_call <- paste(as.character(nullmod$call), collapse = " | ")
cat("Forced nullmod$call: ", forced_call, "\n", sep = "")
if (!grepl(expected_token, forced_call))
  stop("BUG: forced call still missing '", expected_token, "'")
cat("nullmod$call OK (contains '", expected_token, "')\n", sep = "")

save(nullmod, file = null_rda)

## ---- Score test --------------------------------------------------------
cat("[", format(Sys.time()), "] TractorMix.score...\n")
TractorMix.score(obj          = nullmod,
                 infiles      = hapdose,
                 outfiles     = result_file,
                 AC_threshold = 1,
                 n_core       = as.integer(Sys.getenv("TM_NCORE", "2")),
                 chunk_size   = 2048)
t2 <- Sys.time()
cat("[", format(t2), "] score test done (",
    round(as.numeric(t2 - t1, units = "secs"), 1), "s)\n", sep = "")
cat("Total: ", round(as.numeric(t2 - t0, units = "secs"), 1), "s\n", sep = "")
cat("Wrote: ", result_file, "\n", sep = "")
