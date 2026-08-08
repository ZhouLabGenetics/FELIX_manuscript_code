## Generate ALT phenotypes for power evaluation.
##
## Per-(trait, scenario) beta grids let each condition be calibrated
## independently — quantitative + common saturates around 0.5, but binary low
## prevalence + AFR-only needs much larger betas to escape the floor.
##
## Existing pheno files are reused if their filename matches a (trait, scen,
## beta, seed) listed in BETA_GRIDS. New combinations are generated fresh.
## Set OVERWRITE=TRUE in the environment to force regeneration of every file.
##
## Outputs per mode:
##   data/<mode>/pheno/alt/<trait>_<scen>_beta<round(B*100)>_seed<s>.tsv
##   data/<mode>/pheno/alt/manifest.tsv
##
## Filename beta is round(B*100) zero-padded to 3 digits, so 0.10 -> beta010,
## 0.13 -> beta013, 7.50 -> beta750. The %03d encoding lets fractional betas
## (e.g. 0.13, 0.16, 0.18) coexist with the multiples-of-0.1 grid points.
##
## Manifest columns (long format, one row per pheno file):
##   trait | scenario | beta | seed | causal_snp | causal_chr | causal_pos
##         | file | beta0
##
## Usage:  Rscript generate_alt_pheno.R <common|lowfreq>
##         OVERWRITE=TRUE Rscript generate_alt_pheno.R common   # force regen

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: generate_alt_pheno.R <common|lowfreq>")
mode <- args[1]
OVERWRITE <- isTRUE(as.logical(Sys.getenv("OVERWRITE", "FALSE")))

source(file.path(Sys.getenv("FELIX_SIM_BASE", unset = stop("Set FELIX_SIM_BASE to the simulation directory before running this script")), "R", "config.R"))
set_mode(mode)
library(data.table)

## --- knobs --------------------------------------------------------------
N_SEEDS <- 100
SCENARIOS <- list(
  shared = c(afr = 1.0, eur = 1.0, nat = 1.0),
  afr    = c(afr = 1.0, eur = 0.0, nat = 0.0),
  hetero = c(afr = 0.0, eur = 0.5, nat = 1.0)
)

## Per-(mode, trait, scenario) beta grid. Tune these to keep power curves on
## the rising slope (i.e. neither flat at 0 nor saturated at 1).
## Values must be multiples of 0.01 (filenames use round(B*100)).
BETA_GRIDS <- list(
  common = list(
    quant = list(shared = c(0.10, 0.13, 0.16, 0.18, 0.20, 0.30, 0.40, 0.50),
                 afr    = c(0.10, 0.20, 0.30, 0.40, 0.50),
                 hetero = c(0.10, 0.20, 0.30, 0.40, 0.50)),
    bin10 = list(shared = c(0.40, 0.80, 1.20, 1.60, 2.00),
                 afr    = c(0.40, 0.80, 1.20, 1.60, 2.00),
                 hetero = c(0.60, 1.20, 1.80, 2.40, 3.00)),
    bin01 = list(shared = c(0.60, 1.20, 1.80, 2.40, 3.00),
                 afr    = c(0.80, 1.60, 2.40, 3.20, 4.00),
                 hetero = c(1.00, 2.00, 3.00, 4.00, 5.00))
  ),
  lowfreq = list(
    quant = list(shared = c(0.30, 0.60, 0.90, 1.20, 1.50),
                 afr    = c(0.80, 1.60, 2.40, 3.20, 4.00),
                 hetero = c(0.50, 1.00, 1.50, 2.00, 2.50)),
    bin10 = list(shared = c(0.40, 0.80, 1.20, 1.60, 2.00),
                 afr    = c(1.00, 2.00, 3.00, 4.00, 5.00),
                 hetero = c(0.60, 1.20, 1.80, 2.40, 3.00)),
    bin01 = list(shared = c(0.60, 1.20, 1.80, 2.40, 3.00),
                 afr    = c(1.50, 3.00, 4.50, 6.00, 7.50),
                 hetero = c(1.00, 2.00, 3.00, 4.00, 5.00))
  )
)

TRAITS    <- c("quant", "bin10", "bin01")
BETA_AFR  <- 0.2
BETA_NAT  <- 0.1
SIGMA_G2  <- 0.3
SIGMA_E2  <- 0.7
PREV_BIN01 <- 0.01
PREV_BIN10 <- 0.10
## ------------------------------------------------------------------------

pheno_dir <- file.path(DATA_DIR, "pheno", "alt")
dir.create(pheno_dir, recursive = TRUE, showWarnings = FALSE)

adm <- fread(file.path(DATA_DIR, "Admprop.tsv"), sep = "\t", header = TRUE)
stopifnot(nrow(adm) == N_TOTAL)
iid   <- adm$ID
X_afr <- adm$AFR
X_nat <- adm$NAT

K_dt <- fread(file.path(DATA_DIR, "Kinship.tsv"), sep = "\t", header = TRUE)
K    <- as.matrix(K_dt[, -1, with = FALSE])
rownames(K) <- K_dt$ID
stopifnot(identical(rownames(K), iid))

Sigma <- SIGMA_G2 * 2 * K
eig   <- eigen(Sigma, symmetric = TRUE)
eig$values[eig$values < 0] <- 0
L_fac <- eig$vectors %*% diag(sqrt(eig$values))
draw_re <- function() as.numeric(L_fac %*% rnorm(N_TOTAL))

## --- causal variant pick (unchanged) -----------------------------------
brdir <- file.path(DATA_DIR, "Block_row1")
read_br <- function(base) fread(file.path(brdir, paste0(base, ".tsv")), sep = "\t")
GTot <- read_br("GTot"); GAfr <- read_br("GAfr")
GEur <- read_br("GEur"); GNat <- read_br("GNat")
meta_cols <- 1:5
to_mat <- function(dt) as.matrix(dt[, -meta_cols, with = FALSE])
m_tot <- to_mat(GTot); m_afr <- to_mat(GAfr)
m_eur <- to_mat(GEur); m_nat <- to_mat(GNat)
af  <- rowMeans(m_tot) / 2
maf <- pmin(af, 1 - af)
target_maf <- if (mode == "common") 0.2 else 0.03
rv_afr <- apply(m_afr, 1, function(v) var(v) > 0)
rv_eur <- apply(m_eur, 1, function(v) var(v) > 0)
rv_nat <- apply(m_nat, 1, function(v) var(v) > 0)
ok   <- rv_afr & rv_eur & rv_nat
cand <- which(ok)
if (!length(cand)) stop("No candidate causal variant in Block_row1")
pick <- cand[which.min(abs(maf[cand] - target_maf))]
causal_meta   <- GTot[pick, 1:5]
causal_ds_afr <- m_afr[pick, ]
causal_ds_eur <- m_eur[pick, ]
causal_ds_nat <- m_nat[pick, ]
cat("Picked causal variant:", causal_meta$SNP,
    sprintf("(CHR=%s POS=%s MAF=%.3f)\n",
            causal_meta$CHR, causal_meta$POS, maf[pick]))

## Stamp the causal SNP id to a tiny side-file so aggregate_power.R can read
## it without depending on the manifest format. Causal pick is deterministic
## per mode, so this is idempotent across reruns.
writeLines(as.character(causal_meta$SNP),
           file.path(pheno_dir, "causal.txt"))

calibrate_beta0 <- function(eta, target_prev) {
  f <- function(b0) mean(plogis(eta + b0)) - target_prev
  uniroot(f, interval = c(-20, 20), tol = 1e-6)$root
}

## --- generate ----------------------------------------------------------
manifest <- list()
n_kept <- 0L
n_new  <- 0L

for (scen in names(SCENARIOS)) {
  w        <- SCENARIOS[[scen]]
  scen_idx <- which(names(SCENARIOS) == scen)

  for (tr in TRAITS) {
    grid <- BETA_GRIDS[[mode]][[tr]][[scen]]
    if (is.null(grid)) {
      cat("No grid for", tr, scen, "— skipping\n"); next
    }
    prev <- if (tr == "bin01") PREV_BIN01 else if (tr == "bin10") PREV_BIN10 else NA

    for (B in grid) {
      ba <- B * w[["afr"]]; be <- B * w[["eur"]]; bn <- B * w[["nat"]]
      genetic <- ba * causal_ds_afr + be * causal_ds_eur + bn * causal_ds_nat

      for (s in seq_len(N_SEEDS)) {
        tag   <- sprintf("%s_beta%03d_seed%02d", scen, round(B * 100), s)
        fname <- sprintf("%s_%s.tsv", tr, tag)
        fpath <- file.path(pheno_dir, fname)

        if (!OVERWRITE && file.exists(fpath)) {
          n_kept <- n_kept + 1L
          manifest[[length(manifest) + 1L]] <- data.frame(
            trait = tr, scenario = scen, beta = B, seed = s,
            causal_snp = causal_meta$SNP,
            causal_chr = causal_meta$CHR,
            causal_pos = causal_meta$POS,
            file  = fname,
            beta0 = NA_real_,
            stringsAsFactors = FALSE
          )
          next
        }

        ## Same RNG seed across traits for fixed (scen, B, s) so g is
        ## consistent — useful for sanity-checking and for re-running a
        ## subset deterministically.
        set.seed(2026L * 1e6 + scen_idx * 1e5 + round(B * 100) * 100 + s)

        g       <- draw_re()
        e       <- rnorm(N_TOTAL, 0, sqrt(SIGMA_E2))
        eta_cov <- BETA_AFR * X_afr + BETA_NAT * X_nat

        if (tr == "quant") {
          y     <- eta_cov + g + e + genetic
          beta0 <- NA_real_
        } else {
          eta_bin <- eta_cov + g + genetic
          beta0   <- calibrate_beta0(eta_bin, prev)
          y       <- rbinom(N_TOTAL, 1, plogis(eta_bin + beta0))
        }

        out <- data.frame(IID = iid, PHENO = y, AFR = X_afr, NAT = X_nat)
        fwrite(out, fpath, sep = "\t", quote = FALSE)
        n_new <- n_new + 1L

        manifest[[length(manifest) + 1L]] <- data.frame(
          trait = tr, scenario = scen, beta = B, seed = s,
          causal_snp = causal_meta$SNP,
          causal_chr = causal_meta$CHR,
          causal_pos = causal_meta$POS,
          file  = fname,
          beta0 = beta0,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

man <- rbindlist(manifest)
fwrite(man, file.path(pheno_dir, "manifest.tsv"), sep = "\t", quote = FALSE)
cat(sprintf("Manifest rows: %d  (existing reused: %d, newly written: %d)\n",
            nrow(man), n_kept, n_new))
cat("Wrote", file.path(pheno_dir, "manifest.tsv"), "\n")
