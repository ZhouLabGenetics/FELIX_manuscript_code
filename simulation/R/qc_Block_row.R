## Filter variants in a block_row by realised total minor allele frequency.
##   common  : keep MAF in [0.05, 0.50]
##   lowfreq : keep MAF in [0.01, 0.05]
## Also drop sites with zero variance in any ancestry-specific matrix so the
## joint score test doesn't choke on degenerate columns.
##
## Usage:  Rscript qc_Block_row.R <No> <common|lowfreq>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: qc_Block_row.R <No> <common|lowfreq>")
No   <- as.integer(args[1])
mode <- args[2]

source(file.path(Sys.getenv("FELIX_SIM_BASE", unset = stop("Set FELIX_SIM_BASE to the simulation directory before running this script")), "R", "config.R"))
set_mode(mode)
library(data.table)

brdir <- file.path(DATA_DIR, paste0("Block_row", No))

read_br <- function(base) fread(file.path(brdir, paste0(base, ".tsv")), sep = "\t")
write_br <- function(dt, base) fwrite(dt, file.path(brdir, paste0(base, ".tsv")),
                                      sep = "\t", quote = FALSE, row.names = FALSE)

GTot  <- read_br("GTot")
GAfr  <- read_br("GAfr")
GEur  <- read_br("GEur")
GNat  <- read_br("GNat")
LAafr <- read_br("LAafr")
LAeur <- read_br("LAeur")
LAnat <- read_br("LAnat")

## First 5 cols are metadata (CHR, POS, SNP, REF, ALT).
meta_cols <- 1:5
geno_mat  <- as.matrix(GTot[, -meta_cols, with = FALSE])
af        <- rowMeans(geno_mat) / 2
maf       <- pmin(af, 1 - af)

var_nonzero <- function(dt) {
  m <- as.matrix(dt[, -meta_cols, with = FALSE])
  matrixStats::rowVars(m) > 0
}

if (!requireNamespace("matrixStats", quietly = TRUE)) {
  ## fallback if matrixStats is missing
  var_nonzero <- function(dt) {
    m <- as.matrix(dt[, -meta_cols, with = FALSE])
    apply(m, 1, var) > 0
  }
}

maf_ok <- maf >= MAF_LOW & maf <= MAF_HIGH
## require non-zero variance in all three ancestry-specific matrices so the
## joint 3-df score test is well-defined
keep <- maf_ok & var_nonzero(GAfr) & var_nonzero(GEur) & var_nonzero(GNat)

cat("Block_row", No, "mode", mode,
    ": input", nrow(GTot), "variants, kept", sum(keep),
    "(", round(100 * mean(keep), 1), "% )\n")

write_br(GTot[keep,],  "GTot")
write_br(GAfr[keep,],  "GAfr")
write_br(GEur[keep,],  "GEur")
write_br(GNat[keep,],  "GNat")
write_br(LAafr[keep,], "LAafr")
write_br(LAeur[keep,], "LAeur")
write_br(LAnat[keep,], "LAnat")
