## Build a small pruned PLINK fileset (~3000 markers) to feed SAIGE's
## createSparseGRM.R for sparse-GRM construction.
##
## Inputs : Block_row1/GTot.tsv  (post-QC, ~50k variants x 5000 samples)
## Outputs: data/<mode>/plink/pruned.{bed,bim,fam}
##
## Usage:  Rscript make_plink_pruned.R <common|lowfreq>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: make_plink_pruned.R <common|lowfreq>")
mode <- args[1]

source(file.path(Sys.getenv("FELIX_SIM_BASE", unset = stop("Set FELIX_SIM_BASE to the simulation directory before running this script")), "R", "config.R"))
set_mode(mode)
library(data.table)
suppressPackageStartupMessages(library(genio))

set.seed(42L)
N_MARKERS <- 3000L

brdir <- file.path(DATA_DIR, "Block_row1")
GTot  <- fread(file.path(brdir, "GTot.tsv"), sep = "\t")

meta_cols <- 1:5
geno <- as.matrix(GTot[, -meta_cols, with = FALSE])   # variants x samples
n_var <- nrow(geno)
n_sam <- ncol(geno)
if (n_var < N_MARKERS) {
  warning("Only ", n_var, " variants available in Block_row1; using all of them.")
  keep_idx <- seq_len(n_var)
} else {
  keep_idx <- sort(sample.int(n_var, N_MARKERS))
}

geno <- geno[keep_idx, , drop = FALSE]
meta <- GTot[keep_idx, 1:5, with = FALSE]

## genio expects a matrix where rows are variants and columns are samples.
## rownames(X) = variant id, colnames(X) = sample id.
rownames(geno) <- meta$SNP
colnames(geno) <- paste0("Sample", seq_len(n_sam))

## bim: chr, id, posg, pos, ref, alt
bim <- data.frame(
  chr  = meta$CHR,
  id   = meta$SNP,
  posg = 0,
  pos  = meta$POS,
  ref  = meta$REF,
  alt  = meta$ALT,
  stringsAsFactors = FALSE
)

## fam: one row per sample; SAIGE expects FID == IID so they match phenotype.
fam <- data.frame(
  fam = paste0("Sample", seq_len(n_sam)),
  id  = paste0("Sample", seq_len(n_sam)),
  pat = 0, mat = 0, sex = 0, pheno = -9,
  stringsAsFactors = FALSE
)

plink_dir <- file.path(DATA_DIR, "plink")
dir.create(plink_dir, recursive = TRUE, showWarnings = FALSE)
prefix <- file.path(plink_dir, "pruned")

write_plink(prefix, X = geno, bim = bim, fam = fam)

cat("Wrote", nrow(bim), "markers x", nrow(fam), "samples to", prefix, ".{bed,bim,fam}\n")
