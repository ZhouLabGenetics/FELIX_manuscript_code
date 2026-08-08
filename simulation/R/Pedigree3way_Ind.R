## Simulate one 100-sample block of INDEPENDENT admixed individuals.
## Array task 1..50 (50 blocks x 100 = 5000 independent samples).
##
## Usage:  Rscript Pedigree3way_Ind.R <No> <common|lowfreq>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Pedigree3way_Ind.R <No> <common|lowfreq>")
No   <- as.integer(args[1])
mode <- args[2]

source(file.path(Sys.getenv("FELIX_SIM_BASE", unset = stop("Set FELIX_SIM_BASE to the simulation directory before running this script")), "R", "config.R"))
set_mode(mode)
source_utils()
library(data.table)

## Independent samples are numbered 1..2500 across all Ind blocks.
Blockdir  <- file.path(IND_DIR, paste0("Block", No))
dir.create(Blockdir, recursive = TRUE, showWarnings = FALSE)

set.seed(10000L * ifelse(mode == "common", 1L, 2L) + No)

MAFtbl    <- fread(file.path(DATA_DIR, "SNP_AF.tsv"), sep = "\t", header = TRUE)
MAFs      <- list(maf_afr = MAFtbl$MAF_AFR,
                  maf_eur = MAFtbl$MAF_EUR,
                  maf_nat = MAFtbl$MAF_NAT)
p         <- nrow(MAFtbl)

n        <- BLOCK_N
SampleID <- paste0("Sample", ((No - 1) * n + 1):((No - 1) * n + n))

GTot    <- matrix(NA_integer_, n, p); GAfr <- GTot; GEur <- GTot; GNat <- GTot
LAafr   <- GTot; LAeur <- GTot; LAnat <- GTot
Admprop <- matrix(NA_real_, n, 3, dimnames = list(NULL, c("AFR","EUR","NAT")))

for (i in seq_len(n)) {
  props = as.numeric(rDirichlet(1, DIRICHLET_ALPHA))
  ind   = GenerateAdm3(MAFs$maf_afr, MAFs$maf_eur, MAFs$maf_nat, props)
  GTot[i, ]  = ind$hap1 + ind$hap2
  GAfr[i, ]  = GetAncestrySpecCount(ind, 0)
  GEur[i, ]  = GetAncestrySpecCount(ind, 1)
  GNat[i, ]  = GetAncestrySpecCount(ind, 2)
  LAafr[i, ] = (ind$la1 == 0) + (ind$la2 == 0)
  LAeur[i, ] = (ind$la1 == 1) + (ind$la2 == 1)
  LAnat[i, ] = (ind$la1 == 2) + (ind$la2 == 2)
  Admprop[i, ] = GetGlobAncestry3(ind)
}

rownames(GTot) = rownames(GAfr) = rownames(GEur) = rownames(GNat) = SampleID
rownames(LAafr) = rownames(LAeur) = rownames(LAnat) = SampleID
rownames(Admprop) = SampleID

## Written as variants x samples matrices (matches Tan's layout).
write_var_x_sam <- function(M, file) {
  fwrite(as.data.table(t(M)), file, sep = "\t", quote = FALSE,
         col.names = TRUE, row.names = FALSE)
}

write_var_x_sam(GTot,  file.path(Blockdir, "GTot.tsv"))
write_var_x_sam(GAfr,  file.path(Blockdir, "GAfr.tsv"))
write_var_x_sam(GEur,  file.path(Blockdir, "GEur.tsv"))
write_var_x_sam(GNat,  file.path(Blockdir, "GNat.tsv"))
write_var_x_sam(LAafr, file.path(Blockdir, "LAafr.tsv"))
write_var_x_sam(LAeur, file.path(Blockdir, "LAeur.tsv"))
write_var_x_sam(LAnat, file.path(Blockdir, "LAnat.tsv"))

write.table(Admprop, file.path(Blockdir, "Admprop.tsv"),
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)

cat("Finished Ind block", No, "for mode", mode, "\n")
