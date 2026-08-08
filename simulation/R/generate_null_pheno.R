## Generate null phenotypes for Type-I error evaluation.
## Produces, per mode, a pheno file per (trait_type, seed) combination.
##
## Trait types:
##   - quant  : continuous,  y = g + e
##              with g ~ N(0, sigma_g^2 * 2K), e ~ N(0, sigma_e^2 I)
##   - bin01  : binary, prevalence = 0.01
##   - bin10  : binary, prevalence = 0.10
##
## TRUE NULL: no covariate effects, no genetic effects. For quantitative
## traits the kinship random effect is kept (SAIGE models it via the sparse
## GRM). For binary traits the random effect is DROPPED so that case status
## is iid Bernoulli(prevalence) — this avoids pedigree-driven case clustering
## that creates confounding between local ancestry and case status when the
## per-ancestry case count is very small (e.g. N_case_anc1 ~ 7-9 at 1%
## prevalence with ~10% AFR admixture).
##
## Covariates AFR and NAT are written to the file so SAIGE can condition on
## them, but they have zero true effect on the phenotype.
##
## Usage:  Rscript generate_null_pheno.R <common|lowfreq>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: generate_null_pheno.R <common|lowfreq>")
mode <- args[1]

source(file.path(Sys.getenv("FELIX_SIM_BASE", unset = stop("Set FELIX_SIM_BASE to the simulation directory before running this script")), "R", "config.R"))
set_mode(mode)
library(data.table)
suppressPackageStartupMessages(library(MASS))

## --- knobs --------------------------------------------------------------
N_SEEDS      <- 10
## Variance components for the quantitative trait.
SIGMA_G2     <- 0.3
SIGMA_E2     <- 0.7
PREV_BIN01   <- 0.01
PREV_BIN10   <- 0.10
## ------------------------------------------------------------------------

pheno_dir <- file.path(DATA_DIR, "pheno", "null")
dir.create(pheno_dir, recursive = TRUE, showWarnings = FALSE)

## Load covariates (written to pheno file but with zero true effect).
adm <- fread(file.path(DATA_DIR, "Admprop.tsv"), sep = "\t", header = TRUE)
stopifnot(nrow(adm) == N_TOTAL)
iid <- adm$ID
X_afr <- adm$AFR
X_nat <- adm$NAT

## Load true pedigree kinship for quantitative random effect.
K_dt <- fread(file.path(DATA_DIR, "Kinship.tsv"), sep = "\t", header = TRUE)
K    <- as.matrix(K_dt[, -1, with = FALSE])
rownames(K) <- K_dt$ID
stopifnot(identical(rownames(K), iid))

## Random-effect covariance: 2 * K * sigma_g^2.
Sigma   <- SIGMA_G2 * 2 * K
eig     <- eigen(Sigma, symmetric = TRUE)
eig$values[eig$values < 0] <- 0
L_fac   <- eig$vectors %*% diag(sqrt(eig$values))

draw_re <- function() as.numeric(L_fac %*% rnorm(N_TOTAL))

write_pheno <- function(trait, seed, y) {
  out <- data.frame(IID = iid, PHENO = y, AFR = X_afr, NAT = X_nat)
  fn  <- file.path(pheno_dir, sprintf("%s_seed%02d.tsv", trait, seed))
  fwrite(out, fn, sep = "\t", quote = FALSE)
}

for (s in seq_len(N_SEEDS)) {
  set.seed(2026L * 1000L + s)

  g <- draw_re()
  e <- rnorm(N_TOTAL, 0, sqrt(SIGMA_E2))

  ## Quantitative: random effect + noise, no covariate effect.
  y_quant <- g + e
  write_pheno("quant", s, y_quant)

  ## Binary: iid Bernoulli(prevalence). No random effect, no covariates.
  ## This is the cleanest null for binary traits — case status is completely
  ## independent of pedigree and local ancestry.
  y_b01 <- rbinom(N_TOTAL, 1, PREV_BIN01)
  write_pheno("bin01", s, y_b01)

  y_b10 <- rbinom(N_TOTAL, 1, PREV_BIN10)
  write_pheno("bin10", s, y_b10)

  cat(sprintf("seed %02d : quant done | bin01 prev=%.3f | bin10 prev=%.3f\n",
              s, mean(y_b01), mean(y_b10)))
}

cat("All null phenotypes written under", pheno_dir, "\n")
