## Minimal large-N benchmark simulator for FELIX vs Tractor-Mix.

source(file.path(Sys.getenv("FELIX_SIM_BASE", unset = stop("Set FELIX_SIM_BASE to the simulation directory before running this script")), "R", "config.R"))
set_mode("bench")
source_utils()

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

stopifnot(N_VARIANTS >= 100L)
stopifnot(N_TOTAL %% FAMSIZE == 0L)
N_FAM <- N_TOTAL %/% FAMSIZE

cat(sprintf("simu_bench_largeN  N=%d  P=%d  FAMSIZE=%d  N_FAM=%d\n",
            N_TOTAL, N_VARIANTS, FAMSIZE, N_FAM))

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

## ---- 2. Simulate N_FAM families of FAMSIZE full sibs ------------------
## Per family: draw Dirichlet ancestry props once (shared parents),
## generate 2 admixed parents, then generate FAMSIZE sibs via gene drop.
cat("[2/5] simulating", N_FAM, "families x", FAMSIZE, "sibs...\n")

## Output matrices: samples x variants (transpose to variants x samples when writing)
sample_ids <- paste0("Sample", seq_len(N_TOTAL))
GTot  <- matrix(0L, N_TOTAL, N_VARIANTS)
GAfr  <- GTot; GEur <- GTot; GNat <- GTot
LAafr <- GTot; LAeur <- GTot; LAnat <- GTot
Admprop <- matrix(0,  N_TOTAL, 3, dimnames = list(NULL, c("AFR","EUR","NAT")))

t_start <- Sys.time()
for (f in seq_len(N_FAM)) {
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
  if (f %% max(1L, N_FAM %/% 20L) == 0L) {
    dt <- as.numeric(Sys.time() - t_start, units = "secs")
    cat(sprintf("  ...family %d/%d  (%.1fs elapsed, ETA %.0fs)\n",
                f, N_FAM, dt, dt / f * (N_FAM - f)))
  }
}
cat(sprintf("  done sim, total %.1fs\n",
            as.numeric(Sys.time() - t_start, units = "secs")))

## ---- 3. Write Block_row1 + Admprop ------------------------------------
cat("[3/5] writing Block_row1/* and Admprop ...\n")
brdir <- file.path(DATA_DIR, "Block_row1")
dir.create(brdir, recursive = TRUE, showWarnings = FALSE)

## All Block_row1/*.tsv layouts in the main pipeline are variants x samples
## with 5 metadata columns then 1 col per sample.
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
write_var_x_sam(GTot,  "GTot")
write_var_x_sam(GAfr,  "GAfr")
write_var_x_sam(GEur,  "GEur")
write_var_x_sam(GNat,  "GNat")
write_var_x_sam(LAafr, "LAafr")
write_var_x_sam(LAeur, "LAeur")
write_var_x_sam(LAnat, "LAnat")

Admprop_dt <- data.table(ID = sample_ids,
                         AFR = Admprop[, "AFR"],
                         EUR = Admprop[, "EUR"],
                         NAT = Admprop[, "NAT"])
fwrite(Admprop_dt, file.path(DATA_DIR, "Admprop.tsv"), sep = "\t", quote = FALSE)

## Free the genotype matrices we no longer need before kinship/pheno.
rm(GTot, GAfr, GEur, GNat, LAafr, LAeur, LAnat); invisible(gc())

## ---- 4. Pedigree kinship: block-diagonal, sparse ----------------------
cat("[4/5] building pedigree kinship (sparse block-diagonal) ...\n")
## Full-sib family block: diag 0.5, off-diag 0.25.
K_block <- matrix(0.25, FAMSIZE, FAMSIZE)
diag(K_block) <- 0.5
K_sparse <- bdiag(replicate(N_FAM, K_block, simplify = FALSE))
K_sparse <- as(K_sparse, "dgCMatrix")
rownames(K_sparse) <- colnames(K_sparse) <- sample_ids
saveRDS(K_sparse, file.path(DATA_DIR, "kinship_sparse.rds"))
cat("  wrote kinship_sparse.rds  ( nnz=", length(K_sparse@x), " )\n", sep = "")

## ---- 5. Null phenotypes (1 seed per trait) ---------------------------
cat("[5/5] generating null phenotypes (1 seed each: quant, bin10, bin01) ...\n")
pheno_dir <- file.path(DATA_DIR, "pheno", "null")
dir.create(pheno_dir, recursive = TRUE, showWarnings = FALSE)

SIGMA_G2   <- 0.3
SIGMA_E2   <- 0.7
PREV_BIN10 <- 0.10
PREV_BIN01 <- 0.01

set.seed(20260601L + 1L)
## Quant random effect: per-family draw using the 10x10 block cholesky --
## avoids ever forming the dense N x N covariance.
Sigma_block <- 2 * SIGMA_G2 * K_block
L_block     <- t(chol(Sigma_block))      # lower-tri so L %*% z gives MVN draws
g <- numeric(N_TOTAL)
for (f in seq_len(N_FAM)) {
  idx     <- ((f - 1L) * FAMSIZE + 1L):(f * FAMSIZE)
  g[idx]  <- as.numeric(L_block %*% rnorm(FAMSIZE))
}
e <- rnorm(N_TOTAL, 0, sqrt(SIGMA_E2))
y_quant <- g + e
y_b10   <- rbinom(N_TOTAL, 1L, PREV_BIN10)
y_b01   <- rbinom(N_TOTAL, 1L, PREV_BIN01)

write_pheno <- function(trait, y) {
  out <- data.table(IID = sample_ids,
                    PHENO = y,
                    AFR = Admprop_dt$AFR,
                    NAT = Admprop_dt$NAT)
  fn  <- file.path(pheno_dir, sprintf("%s_seed01.tsv", trait))
  fwrite(out, fn, sep = "\t", quote = FALSE)
}
write_pheno("quant", y_quant)
write_pheno("bin10", y_b10)
write_pheno("bin01", y_b01)

cat(sprintf("  bin10 prev=%.3f  bin01 prev=%.4f\n", mean(y_b10), mean(y_b01)))
cat("All done.  N=", N_TOTAL, " P=", N_VARIANTS, "  dir=", DATA_DIR, "\n", sep = "")
