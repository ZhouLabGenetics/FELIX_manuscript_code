## Mixed-relatedness 3-way admixed simulation for FELIX benchmarking.
##
## Cohort composition controlled by UNREL_PCT (env var, integer 0..100):
##   UNREL_PCT=100 -> all unrelated (original case)
##   UNREL_PCT=75  -> 75% unrelated + 25% in FAMSIZE=10 sib families
##   UNREL_PCT=50  -> 50/50
##   UNREL_PCT=25  -> 25% unrelated + 75% related
##
## Related arm: each family draws ONE Dirichlet ancestry vector, generates
## 2 admixed parents, then FAMSIZE=10 full sibs via GenerateChild.
## Unrelated arm: each individual draws own Dirichlet, then GenerateAdm3.
##
## Outputs (under UNR_DIR):
##   Block_row1/{GTot,GAfr,GEur,GNat,LAafr,LAeur,LAnat}.tsv
##   Admprop.tsv
##   SNP_AF.tsv
##   kinship_sparse.rds        block-diag (families) + 0.5*I (unrelated)
##   pheno/null/{quant,bin10,bin01}_seed01..NN.tsv
##
## Phenotype RNG: PHENO column uses `set.seed(20260601L + seed_i)`, so the
## same noise vector is produced across cohorts. AFR/NAT covariate columns
## are cohort-specific because Admprop changes with structure.
##
## Memory: peak ~14 GB during sim (7 int matrices at N=10k, P=50k). Block_row
## writes use per-matrix rm()+gc() to keep write-phase peak under ~18 GB,
## so the whole run fits in 20 GB.
##
## Usage:
##   UNR_DIR=/.../unr_N10000_P50000_unrel100  UNR_N=10000  UNR_P=50000  \
##   UNREL_PCT=100  Rscript simu_unrelated.R [N_SEEDS]

args <- commandArgs(trailingOnly = TRUE)
n_seeds <- if (length(args) >= 1) as.integer(args[1]) else 10L
stopifnot(is.finite(n_seeds), n_seeds >= 1L, n_seeds <= 99L)

unrel_pct <- as.integer(Sys.getenv("UNREL_PCT", "100"))
stopifnot(is.finite(unrel_pct), unrel_pct >= 0L, unrel_pct <= 100L)

source(file.path(Sys.getenv("FELIX_SIM_BASE", unset = stop("Set FELIX_SIM_BASE to the simulation directory before running this script")), "R", "config.R"))
set_mode("unr")
source_utils()

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

stopifnot(N_VARIANTS >= 100L)

rel_count <- as.integer(round(N_TOTAL * (100L - unrel_pct) / 100))
if (rel_count %% FAMSIZE != 0L)
  stop("rel_count=", rel_count, " must be divisible by FAMSIZE=", FAMSIZE,
       "; adjust UNR_N or UNREL_PCT.")
n_fam     <- rel_count %/% FAMSIZE
unr_count <- N_TOTAL - rel_count

cat(sprintf("simu_unrelated  N=%d (rel=%d in %d fams, unrel=%d)  P=%d  N_SEEDS=%d\n",
            N_TOTAL, rel_count, n_fam, unr_count, N_VARIANTS, n_seeds))

set.seed(20260601L)

## ---- 1. MAFs under Balding-Nichols ------------------------------------
cat("[1/5] generating MAFs...\n")
mafs <- GenerateMafs3(N_VARIANTS,
                      fst_afr = 0.15, fst_eur = 0.10, fst_nat = 0.15,
                      runif1  = ANC_P_LOW, runif2 = ANC_P_HIGH)
snp_af <- data.table(
  SNP    = paste0("SNP", seq_len(N_VARIANTS)),
  CHR    = 1L,
  POS    = seq_len(N_VARIANTS),
  REF    = "A",
  ALT    = "G",
  MAF_ANC = mafs$maf_anc,
  MAF_AFR = mafs$maf_afr,
  MAF_EUR = mafs$maf_eur,
  MAF_NAT = mafs$maf_nat
)
fwrite(snp_af, file.path(DATA_DIR, "SNP_AF.tsv"), sep = "\t", quote = FALSE)

## ---- 2. Simulate cohort: families first, then unrelated ---------------
cat(sprintf("[2/5] simulating cohort (rel=%d / unrel=%d) ...\n",
            rel_count, unr_count))

sample_ids <- paste0("Sample", seq_len(N_TOTAL))
GTot  <- matrix(0L, N_TOTAL, N_VARIANTS)
GAfr  <- GTot; GEur <- GTot; GNat <- GTot
LAafr <- GTot; LAeur <- GTot; LAnat <- GTot
Admprop <- matrix(0,  N_TOTAL, 3, dimnames = list(NULL, c("AFR","EUR","NAT")))

t_start <- Sys.time()

## --- 2a. Related arm: N_FAM families x FAMSIZE full sibs ---
if (n_fam > 0L) {
  for (f in seq_len(n_fam)) {
    props <- as.numeric(rDirichlet(1, DIRICHLET_ALPHA))
    p1 <- GenerateAdm3(mafs$maf_afr, mafs$maf_eur, mafs$maf_nat, props)
    p2 <- GenerateAdm3(mafs$maf_afr, mafs$maf_eur, mafs$maf_nat, props)
    for (s in seq_len(FAMSIZE)) {
      i   <- (f - 1L) * FAMSIZE + s
      kid <- GenerateChild(p1, p2)
      GTot[i, ]  <- kid$hap1 + kid$hap2
      GAfr[i, ]  <- ((kid$la1 == 0L) & (kid$hap1 == 1L)) + ((kid$la2 == 0L) & (kid$hap2 == 1L))
      GEur[i, ]  <- ((kid$la1 == 1L) & (kid$hap1 == 1L)) + ((kid$la2 == 1L) & (kid$hap2 == 1L))
      GNat[i, ]  <- ((kid$la1 == 2L) & (kid$hap1 == 1L)) + ((kid$la2 == 2L) & (kid$hap2 == 1L))
      LAafr[i, ] <- (kid$la1 == 0L) + (kid$la2 == 0L)
      LAeur[i, ] <- (kid$la1 == 1L) + (kid$la2 == 1L)
      LAnat[i, ] <- (kid$la1 == 2L) + (kid$la2 == 2L)
      Admprop[i, ] <- props
    }
    if (f %% max(1L, n_fam %/% 10L) == 0L) {
      dt <- as.numeric(Sys.time() - t_start, units = "secs")
      cat(sprintf("  rel: ...family %d/%d  (%.1fs elapsed)\n", f, n_fam, dt))
    }
  }
}

## --- 2b. Unrelated arm: each indep Dirichlet + GenerateAdm3 ---
if (unr_count > 0L) {
  i_off <- rel_count
  for (k in seq_len(unr_count)) {
    i     <- i_off + k
    props <- as.numeric(rDirichlet(1, DIRICHLET_ALPHA))
    ind   <- GenerateAdm3(mafs$maf_afr, mafs$maf_eur, mafs$maf_nat, props)
    GTot[i, ]  <- ind$hap1 + ind$hap2
    GAfr[i, ]  <- ((ind$la1 == 0L) & (ind$hap1 == 1L)) + ((ind$la2 == 0L) & (ind$hap2 == 1L))
    GEur[i, ]  <- ((ind$la1 == 1L) & (ind$hap1 == 1L)) + ((ind$la2 == 1L) & (ind$hap2 == 1L))
    GNat[i, ]  <- ((ind$la1 == 2L) & (ind$hap1 == 1L)) + ((ind$la2 == 2L) & (ind$hap2 == 1L))
    LAafr[i, ] <- (ind$la1 == 0L) + (ind$la2 == 0L)
    LAeur[i, ] <- (ind$la1 == 1L) + (ind$la2 == 1L)
    LAnat[i, ] <- (ind$la1 == 2L) + (ind$la2 == 2L)
    Admprop[i, ] <- props
    if (k %% max(1L, unr_count %/% 10L) == 0L) {
      dt <- as.numeric(Sys.time() - t_start, units = "secs")
      cat(sprintf("  unrel: ...ind %d/%d  (%.1fs elapsed)\n", k, unr_count, dt))
    }
  }
}
cat(sprintf("  done sim, total %.1fs\n",
            as.numeric(Sys.time() - t_start, units = "secs")))

## ---- 3. Write Block_row1 + Admprop ------------------------------------
cat("[3/5] writing Block_row1/* and Admprop ...\n")
brdir <- file.path(DATA_DIR, "Block_row1")
dir.create(brdir, recursive = TRUE, showWarnings = FALSE)

meta_dt <- data.table(CHR = 1L,
                      POS = seq_len(N_VARIANTS),
                      SNP = paste0("SNP", seq_len(N_VARIANTS)),
                      REF = "A",
                      ALT = "G")

write_var_x_sam <- function(M, base) {
  ## M is samples x variants; transpose to variants x samples.
  out <- cbind(meta_dt, as.data.table(t(M)))
  setnames(out, c(names(meta_dt), sample_ids))
  fwrite(out, file.path(brdir, paste0(base, ".tsv")), sep = "\t", quote = FALSE)
}
## Per-matrix rm() between writes keeps peak memory below the 20 GB headroom.
write_var_x_sam(GTot,  "GTot");  rm(GTot);  invisible(gc())
write_var_x_sam(GAfr,  "GAfr");  rm(GAfr);  invisible(gc())
write_var_x_sam(GEur,  "GEur");  rm(GEur);  invisible(gc())
write_var_x_sam(GNat,  "GNat");  rm(GNat);  invisible(gc())
write_var_x_sam(LAafr, "LAafr"); rm(LAafr); invisible(gc())
write_var_x_sam(LAeur, "LAeur"); rm(LAeur); invisible(gc())
write_var_x_sam(LAnat, "LAnat"); rm(LAnat); invisible(gc())

Admprop_dt <- data.table(ID  = sample_ids,
                         AFR = Admprop[, "AFR"],
                         EUR = Admprop[, "EUR"],
                         NAT = Admprop[, "NAT"])
fwrite(Admprop_dt, file.path(DATA_DIR, "Admprop.tsv"), sep = "\t", quote = FALSE)

## ---- 4. Sparse kinship: families + identity ---------------------------
cat("[4/5] writing sparse kinship (kinship_sparse.rds) ...\n")
K_block <- matrix(0.25, FAMSIZE, FAMSIZE)
diag(K_block) <- 0.5
if (n_fam > 0L && unr_count > 0L) {
  fam_blocks <- replicate(n_fam, K_block, simplify = FALSE)
  K_sparse   <- bdiag(c(fam_blocks, list(Diagonal(x = rep(0.5, unr_count)))))
} else if (n_fam > 0L) {
  K_sparse <- bdiag(replicate(n_fam, K_block, simplify = FALSE))
} else {
  K_sparse <- Diagonal(x = rep(0.5, N_TOTAL))
}
K_sparse <- as(K_sparse, "generalMatrix")
rownames(K_sparse) <- colnames(K_sparse) <- sample_ids
saveRDS(K_sparse, file.path(DATA_DIR, "kinship_sparse.rds"))
cat("  nnz=", length(K_sparse@x), "\n", sep = "")

## ---- 5. Null phenotypes ----------------------------------------------
cat(sprintf("[5/5] generating %d null pheno seeds for {quant, bin10, bin01} ...\n", n_seeds))
pheno_dir <- file.path(DATA_DIR, "pheno", "null")
dir.create(pheno_dir, recursive = TRUE, showWarnings = FALSE)

PREV_BIN10 <- 0.10
PREV_BIN01 <- 0.01

write_pheno <- function(trait, y, seed_i) {
  out <- data.table(IID   = sample_ids,
                    PHENO = y,
                    AFR   = Admprop_dt$AFR,
                    NAT   = Admprop_dt$NAT)
  fn  <- file.path(pheno_dir, sprintf("%s_seed%02d.tsv", trait, seed_i))
  fwrite(out, fn, sep = "\t", quote = FALSE)
}

for (s in seq_len(n_seeds)) {
  set.seed(20260601L + s)
  y_quant <- rnorm(N_TOTAL)                  # null: pure noise, iid
  y_b10   <- rbinom(N_TOTAL, 1L, PREV_BIN10)
  y_b01   <- rbinom(N_TOTAL, 1L, PREV_BIN01)
  write_pheno("quant", y_quant, s)
  write_pheno("bin10", y_b10,   s)
  write_pheno("bin01", y_b01,   s)
}

cat("All done.  N=", N_TOTAL, " P=", N_VARIANTS, "  dir=", DATA_DIR, "\n", sep = "")
