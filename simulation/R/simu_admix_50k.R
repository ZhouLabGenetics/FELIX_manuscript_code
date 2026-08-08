## All-unrelated 3-way admixture simulator with full local-ancestry tracking,
## sized for N=20k (or whatever ADMIX_N you set). Output schema matches
## simu_unrelated.R / simu_3way.R so make_vcf.R works unchanged.
##
## Three scenarios via ADMIX_SCENARIO env var:
##   2way_50_50    -- Dirichlet(alpha = c(5.0, 5.0)) over (AFR, EUR); NAT=0
##   2way_25_75    -- Dirichlet(alpha = c(2.5, 7.5)) over (AFR, EUR); NAT=0
##   3way_20_30_50 -- Dirichlet(alpha = c(2.0, 3.0, 5.0)) over (AFR, EUR, NAT)
##
## Outputs (under ADMIX_DIR):
##   Block_row1/{GTot,GAfr,GEur,GNat,LAafr,LAeur,LAnat}.tsv
##   Admprop.tsv, SNP_AF.tsv
##   kinship_sparse.rds            (0.5 * I, all unrelated)
##   plink/pruned.{bed,bim,fam}    (inline, ~3000 markers)
##   pheno/null/{quant,bin10,bin01}_seed01..NN.tsv
##     * 2-way pheno cols: IID, PHENO, AFR
##     * 3-way pheno cols: IID, PHENO, AFR, NAT
##
## Algorithm: per-individual slow-path so LA is tracked. Vectorised in
## batches of 500 individuals via:
##     u1, u2 ~ Uniform        ->  la1, la2 by threshold against cumprops
##     maf_at = M_anc[la+1, s] ->  hap1, hap2 ~ Bernoulli
##
## Memory: 7 P*N int matrices + per-batch transients. At ADMIX_N=20k, P=50k:
##   persistent  = 7 * 4 GB = 28 GB
##   per-batch   = ~0.6 GB
##   peak (first GTot.tsv write transient) ~ 32 GB.
##   Fits comfortably in --mem=48G.
##
## Wall: sim ~3-5 min, 7 TSV writes ~5 min, PLINK ~30 s -> ~10 min/scenario.
##
## Usage:
##   ADMIX_DIR=... ADMIX_N=20000 ADMIX_P=50000 ADMIX_SCENARIO=2way_50_50  \
##   Rscript simu_admix_50k.R [N_SEEDS]

args <- commandArgs(trailingOnly = TRUE)
n_seeds <- if (length(args) >= 1) as.integer(args[1]) else 2L
stopifnot(is.finite(n_seeds), n_seeds >= 1L, n_seeds <= 99L)

scenario <- Sys.getenv("ADMIX_SCENARIO", "")
stopifnot(scenario %in% c("2way_50_50", "2way_25_75", "3way_20_30_50"))

source(file.path(Sys.getenv("FELIX_SIM_BASE", unset = stop("Set FELIX_SIM_BASE to the simulation directory before running this script")), "R", "config.R"))
set_mode("admix")
source_utils()

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(genio)
})

stopifnot(N_VARIANTS >= 100L)

alpha_full <- switch(scenario,
  "2way_50_50"    = c(5.0, 5.0, 0.0),
  "2way_25_75"    = c(2.5, 7.5, 0.0),
  "3way_20_30_50" = c(2.0, 3.0, 5.0)
)
is_3way <- scenario == "3way_20_30_50"

cat(sprintf("simu_admix  scenario=%s  alpha=(%s)  N=%d  P=%d  N_SEEDS=%d\n",
            scenario, paste(alpha_full, collapse = ", "),
            N_TOTAL, N_VARIANTS, n_seeds))

set.seed(20260601L)

## ---- 1. MAFs ----------------------------------------------------------
cat("[1/6] generating MAFs ...\n")
mafs <- GenerateMafs3(N_VARIANTS,
                      fst_afr = 0.15, fst_eur = 0.10, fst_nat = 0.15,
                      runif1  = ANC_P_LOW, runif2 = ANC_P_HIGH)
M_anc <- rbind(AFR = mafs$maf_afr,
               EUR = mafs$maf_eur,
               NAT = mafs$maf_nat)

snp_af <- data.table(
  SNP = paste0("SNP", seq_len(N_VARIANTS)),
  CHR = 1L, POS = seq_len(N_VARIANTS), REF = "A", ALT = "G",
  MAF_ANC = mafs$maf_anc,
  MAF_AFR = mafs$maf_afr, MAF_EUR = mafs$maf_eur, MAF_NAT = mafs$maf_nat
)
fwrite(snp_af, file.path(DATA_DIR, "SNP_AF.tsv"), sep = "\t", quote = FALSE)

## ---- 2. Simulate cohort (batched, LA tracked) ------------------------
cat(sprintf("[2/6] simulating %d unrelated admixed individuals (BATCH=500) ...\n",
            N_TOTAL))

sample_ids <- paste0("Sample", seq_len(N_TOTAL))

## 7 persistent P*N int matrices: variants by samples.
GTot  <- matrix(0L, N_VARIANTS, N_TOTAL)
GAfr  <- matrix(0L, N_VARIANTS, N_TOTAL)
GEur  <- matrix(0L, N_VARIANTS, N_TOTAL)
GNat  <- matrix(0L, N_VARIANTS, N_TOTAL)
LAafr <- matrix(0L, N_VARIANTS, N_TOTAL)
LAeur <- matrix(0L, N_VARIANTS, N_TOTAL)
LAnat <- matrix(0L, N_VARIANTS, N_TOTAL)
Admprop <- matrix(0, N_TOTAL, 3, dimnames = list(NULL, c("AFR", "EUR", "NAT")))

BATCH     <- 500L
n_batches <- (N_TOTAL + BATCH - 1L) %/% BATCH
t_start   <- Sys.time()

for (b in seq_len(n_batches)) {
  i_lo   <- (b - 1L) * BATCH + 1L
  i_hi   <- min(b * BATCH, N_TOTAL)
  this_n <- i_hi - i_lo + 1L

  if (is_3way) {
    props_block <- rDirichlet(this_n, alpha_full)
  } else {
    p2 <- rDirichlet(this_n, alpha_full[1:2])
    props_block <- cbind(p2, NAT = 0)
  }
  Admprop[i_lo:i_hi, ] <- props_block
  cumprops <- t(apply(props_block, 1, cumsum))

  u1 <- matrix(runif(this_n * N_VARIANTS), this_n, N_VARIANTS)
  la1_block <- (u1 > cumprops[, 1]) + (u1 > cumprops[, 2])
  rm(u1)
  u2 <- matrix(runif(this_n * N_VARIANTS), this_n, N_VARIANTS)
  la2_block <- (u2 > cumprops[, 1]) + (u2 > cumprops[, 2])
  rm(u2)

  s_idx <- rep.int(seq_len(N_VARIANTS), rep.int(this_n, N_VARIANTS))
  maf_at1 <- M_anc[cbind(c(la1_block) + 1L, s_idx)]; dim(maf_at1) <- c(this_n, N_VARIANTS)
  hap1_block <- matrix(rbinom(this_n * N_VARIANTS, 1L, as.vector(maf_at1)),
                       this_n, N_VARIANTS)
  rm(maf_at1)
  maf_at2 <- M_anc[cbind(c(la2_block) + 1L, s_idx)]; dim(maf_at2) <- c(this_n, N_VARIANTS)
  hap2_block <- matrix(rbinom(this_n * N_VARIANTS, 1L, as.vector(maf_at2)),
                       this_n, N_VARIANTS)
  rm(maf_at2, s_idx)

  ## Derive + transpose-assign into the 7 persistent P*N matrices.
  GTot[,  i_lo:i_hi] <- t(hap1_block + hap2_block)
  GAfr[,  i_lo:i_hi] <- t((la1_block == 0L) * hap1_block + (la2_block == 0L) * hap2_block)
  GEur[,  i_lo:i_hi] <- t((la1_block == 1L) * hap1_block + (la2_block == 1L) * hap2_block)
  GNat[,  i_lo:i_hi] <- t((la1_block == 2L) * hap1_block + (la2_block == 2L) * hap2_block)
  LAafr[, i_lo:i_hi] <- t((la1_block == 0L) + (la2_block == 0L))
  LAeur[, i_lo:i_hi] <- t((la1_block == 1L) + (la2_block == 1L))
  LAnat[, i_lo:i_hi] <- t((la1_block == 2L) + (la2_block == 2L))
  rm(la1_block, la2_block, hap1_block, hap2_block)

  if (b %% max(1L, n_batches %/% 10L) == 0L) {
    dt <- as.numeric(Sys.time() - t_start, units = "secs")
    cat(sprintf("  ...batch %d/%d  (n=%d, %.1fs elapsed)\n",
                b, n_batches, i_hi, dt))
  }
}
cat(sprintf("  done sim, total %.1fs\n",
            as.numeric(Sys.time() - t_start, units = "secs")))

## ---- 3. Inline PLINK pruned (~3000 markers) ---------------------------
cat("[3/6] writing PLINK pruned.{bed,bim,fam} (inline) ...\n")
meta_dt <- data.table(CHR = 1L,
                      POS = seq_len(N_VARIANTS),
                      SNP = paste0("SNP", seq_len(N_VARIANTS)),
                      REF = "A", ALT = "G")
set.seed(42L)
N_MARKERS <- 3000L
n_keep <- min(N_MARKERS, N_VARIANTS)
keep_idx <- sort(sample.int(N_VARIANTS, n_keep))

geno_plink <- GTot[keep_idx, , drop = FALSE]
rownames(geno_plink) <- meta_dt$SNP[keep_idx]
colnames(geno_plink) <- sample_ids
bim <- data.frame(chr = 1L, id = meta_dt$SNP[keep_idx], posg = 0,
                  pos = meta_dt$POS[keep_idx], ref = "A", alt = "G",
                  stringsAsFactors = FALSE)
fam <- data.frame(fam = sample_ids, id = sample_ids,
                  pat = 0, mat = 0, sex = 0, pheno = -9,
                  stringsAsFactors = FALSE)
plink_dir <- file.path(DATA_DIR, "plink")
dir.create(plink_dir, recursive = TRUE, showWarnings = FALSE)
write_plink(file.path(plink_dir, "pruned"), X = geno_plink, bim = bim, fam = fam)
rm(geno_plink); invisible(gc())

## ---- 4. Write 7 Block_row1 TSVs (rm + gc between writes) -------------
cat("[4/6] writing Block_row1/*.tsv (per-matrix free) ...\n")
brdir <- file.path(DATA_DIR, "Block_row1")
dir.create(brdir, recursive = TRUE, showWarnings = FALSE)

write_var_x_sam <- function(M, base) {
  sample_cols <- asplit(M, 2L)
  names(sample_cols) <- sample_ids
  fwrite(c(as.list(meta_dt), sample_cols),
         file.path(brdir, paste0(base, ".tsv")), sep = "\t", quote = FALSE)
  rm(sample_cols); invisible(gc())
}

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

## ---- 5. Sparse kinship (0.5 * I) -------------------------------------
cat("[5/6] writing kinship_sparse.rds ...\n")
K_sparse <- Diagonal(x = rep(0.5, N_TOTAL))
K_sparse <- as(K_sparse, "generalMatrix")
rownames(K_sparse) <- colnames(K_sparse) <- sample_ids
saveRDS(K_sparse, file.path(DATA_DIR, "kinship_sparse.rds"))

## ---- 6. Null phenos --------------------------------------------------
cat(sprintf("[6/6] generating %d null pheno seeds for {quant, bin10, bin01} ...\n",
            n_seeds))
pheno_dir <- file.path(DATA_DIR, "pheno", "null")
dir.create(pheno_dir, recursive = TRUE, showWarnings = FALSE)

PREV_BIN10 <- 0.10
PREV_BIN01 <- 0.01

write_pheno <- function(trait, y, seed_i) {
  if (is_3way) {
    out <- data.table(IID = sample_ids, PHENO = y,
                      AFR = Admprop_dt$AFR, NAT = Admprop_dt$NAT)
  } else {
    out <- data.table(IID = sample_ids, PHENO = y, AFR = Admprop_dt$AFR)
  }
  fn <- file.path(pheno_dir, sprintf("%s_seed%02d.tsv", trait, seed_i))
  fwrite(out, fn, sep = "\t", quote = FALSE)
}

for (s in seq_len(n_seeds)) {
  set.seed(20260601L + s)
  y_quant <- rnorm(N_TOTAL)
  y_b10   <- rbinom(N_TOTAL, 1L, PREV_BIN10)
  y_b01   <- rbinom(N_TOTAL, 1L, PREV_BIN01)
  write_pheno("quant", y_quant, s)
  write_pheno("bin10", y_b10,   s)
  write_pheno("bin01", y_b01,   s)
}

cat("All done.  scenario=", scenario, "  N=", N_TOTAL, " P=", N_VARIANTS,
    "  dir=", DATA_DIR, "\n", sep = "")
