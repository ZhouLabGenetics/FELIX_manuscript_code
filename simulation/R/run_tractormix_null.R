## Run Tractor-Mix on one (mode, trait, seed) NULL condition.
## Fits glmmkin null model on the pre-generated null pheno
## (pheno/null/<trait>_seed<NN>.tsv), then runs TractorMix.score against the
## chr1 3-ancestry hapdose files. Used to demonstrate Type-I error inflation.
##
## Output:
##   data/<mode>/tractormix_out/null/<trait>_seed<NN>/result.tsv  -- per-variant pvals
##   data/<mode>/tractormix_out/null/<trait>_seed<NN>/null.rda    -- fitted glmmkin obj
##
## bench.log (wall, peak RSS) is captured by the SLURM wrapper via
## /usr/bin/time -v, not by this script.
##
## Usage:  Rscript run_tractormix_null.R <mode> <trait> <seed>
##         e.g.  Rscript run_tractormix_null.R common quant 1

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3)
  stop("Usage: run_tractormix_null.R <mode> <trait> <seed>")
mode  <- args[1]
trait <- args[2]
seed  <- as.integer(args[3])

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
TAG <- sprintf("%s_seed%02d", trait, seed)
pheno_file <- file.path(DATA_DIR, "pheno", "null",
                        sprintf("%s_seed%02d.tsv", trait, seed))
if (!file.exists(pheno_file)) stop("Pheno not found: ", pheno_file)

hapdose <- file.path(DATA_DIR, "tractormix", "hapdose",
                     paste0("chr1.anc", 0:2, ".dosage.txt"))
stopifnot(all(file.exists(hapdose)))

out_dir <- file.path(DATA_DIR, "tractormix_out", "null", TAG)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
result_file <- file.path(out_dir, "result.tsv")
null_rda    <- file.path(out_dir, "null.rda")

if (file.exists(result_file) && file.size(result_file) > 0) {
  cat("Already done: ", result_file, "  -- skipping.\n", sep = "")
  quit(save = "no", status = 0)
}

cat("MODE=", mode, "  TRAIT=", trait, "  SEED=", seed, "  TAG=", TAG, "\n",
    sep = "")

## ---- Load pheno + kinship ---------------------------------------------
pheno <- fread(pheno_file, sep = "\t", header = TRUE)
stopifnot(all(c("IID", "PHENO", "AFR", "NAT") %in% names(pheno)))

K_dt <- fread(file.path(DATA_DIR, "Kinship.tsv"), sep = "\t", header = TRUE)
K    <- as.matrix(K_dt[, -1, with = FALSE])
rownames(K) <- colnames(K) <- K_dt$ID
pheno <- pheno[match(rownames(K), as.character(pheno$IID))]
stopifnot(identical(as.character(pheno$IID), rownames(K)))

## ---- Fit null model ----------------------------------------------------
## TractorMix.score does NOT inspect obj$family. It inspects
## as.character(obj$call) and greps for "gaussian" or "binomial". Pass
## family LITERALLY so the call object retains the searchable token.
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

## Force a $call that contains the family token (see run_tractormix.R for
## the rationale -- some GMMAT versions normalise the call and erase the
## literal "gaussian"/"binomial" string the score test greps for).
expected_token <- if (trait == "quant") "gaussian" else "binomial"
nullmod$call <- if (trait == "quant") {
  quote(glmmkin(fixed = PHENO ~ AFR + NAT, family = gaussian()))
} else {
  quote(glmmkin(fixed = PHENO ~ AFR + NAT, family = binomial()))
}
forced_call <- paste(as.character(nullmod$call), collapse = " | ")
if (!grepl(expected_token, forced_call))
  stop("BUG: forced call missing '", expected_token, "'")

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
