## Run Tractor-Mix on one (BENCH_N, trait) condition for the large-N benchmark.


args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: run_tractormix_largeN.R <quant|bin10|bin01>")
trait <- args[1]
seed  <- 1L

source(file.path(Sys.getenv("FELIX_SIM_BASE", unset = stop("Set FELIX_SIM_BASE to the simulation directory before running this script")), "R", "config.R"))
set_mode("bench")
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(GMMAT)
})

TM_DIR <- Sys.getenv("TRACTORMIX_DIR", "/opt/Tractor-Mix")
source(file.path(TM_DIR, "TractorMix.score.R"))

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

cat("BENCH_N=", N_TOTAL, "  TRAIT=", trait, "  TAG=", TAG, "\n", sep = "")

## ---- Load pheno -------------------------------------------------------
pheno <- fread(pheno_file, sep = "\t", header = TRUE)
stopifnot(all(c("IID", "PHENO", "AFR", "NAT") %in% names(pheno)))

## ---- Load kinship: dense preferred; sparse fallback materialised dense
dense_path  <- file.path(DATA_DIR, "Kinship.tsv")
sparse_path <- file.path(DATA_DIR, "kinship_sparse.rds")
if (file.exists(dense_path) && file.size(dense_path) > 0) {
  cat("[", format(Sys.time()), "] loading dense Kinship.tsv ...\n")
  K_dt <- fread(dense_path, sep = "\t", header = TRUE)
  K    <- as.matrix(K_dt[, -1, with = FALSE])
  rownames(K) <- colnames(K) <- K_dt$ID
  rm(K_dt); invisible(gc())
} else if (file.exists(sparse_path)) {
  cat("[", format(Sys.time()), "] dense Kinship.tsv missing; ",
      "loading sparse kinship and materialising dense (",
      sprintf("%.1f GB", (N_TOTAL^2 * 8) / 1024^3), ") ...\n", sep = "")
  K_sp <- readRDS(sparse_path)
  K    <- as.matrix(K_sp)
  rm(K_sp); invisible(gc())
} else {
  stop("No kinship file found (Kinship.tsv or kinship_sparse.rds) in ", DATA_DIR)
}

## Align pheno row order with K row order.
pheno <- pheno[match(rownames(K), as.character(pheno$IID))]
stopifnot(identical(as.character(pheno$IID), rownames(K)))

## ---- Fit null model ----------------------------------------------------
t0 <- Sys.time()
cat("[", format(t0), "] glmmkin null fit (this is the TM bottleneck) ...\n")
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

## Force $call to retain the family token TractorMix.score greps for.
nullmod$call <- if (trait == "quant") {
  quote(glmmkin(fixed = PHENO ~ AFR + NAT, family = gaussian()))
} else {
  quote(glmmkin(fixed = PHENO ~ AFR + NAT, family = binomial()))
}

save(nullmod, file = null_rda)

## ---- Score test --------------------------------------------------------
cat("[", format(Sys.time()), "] TractorMix.score on chr1 hapdose ...\n")
TractorMix.score(obj          = nullmod,
                 infiles      = hapdose,
                 outfiles     = result_file,
                 AC_threshold = 1,
                 n_core       = as.integer(Sys.getenv("TM_NCORE", "1")),
                 chunk_size   = 2048)
t2 <- Sys.time()
cat("[", format(t2), "] score test done (",
    round(as.numeric(t2 - t1, units = "secs"), 1), "s)\n", sep = "")
cat("Total: ", round(as.numeric(t2 - t0, units = "secs"), 1), "s\n", sep = "")
cat("Wrote: ", result_file, "\n", sep = "")
