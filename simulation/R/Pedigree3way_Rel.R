## Simulate one 100-sample block of RELATED admixed individuals (10 families of 10).
## Array task 1..50 (50 blocks x 100 = 5000 related samples).
##
## Usage:  Rscript Pedigree3way_Rel.R <No> <common|lowfreq>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Pedigree3way_Rel.R <No> <common|lowfreq>")
No   <- as.integer(args[1])
mode <- args[2]

source(file.path(Sys.getenv("FELIX_SIM_BASE", unset = stop("Set FELIX_SIM_BASE to the simulation directory before running this script")), "R", "config.R"))
set_mode(mode)
source_utils()
library(data.table)
library(kinship2)

Blockdir <- file.path(REL_DIR, paste0("Block", No))
dir.create(Blockdir, recursive = TRUE, showWarnings = FALSE)

set.seed(20000L * ifelse(mode == "common", 1L, 2L) + No)

MAFtbl    <- fread(file.path(DATA_DIR, "SNP_AF.tsv"), sep = "\t", header = TRUE)
MAFs      <- list(maf_afr = MAFtbl$MAF_AFR,
                  maf_eur = MAFtbl$MAF_EUR,
                  maf_nat = MAFtbl$MAF_NAT)
p         <- nrow(MAFtbl)

n        <- BLOCK_N
## Related samples are numbered (N_IND+1)..(N_IND+N_REL), offset by N_IND.
SampleID <- paste0("Sample", N_IND + ((No - 1) * n + 1):((No - 1) * n + n))

GTot    <- matrix(NA_integer_, n, p); GAfr <- GTot; GEur <- GTot; GNat <- GTot
LAafr   <- GTot; LAeur <- GTot; LAnat <- GTot
Admprop <- matrix(NA_real_, n, 3, dimnames = list(NULL, c("AFR","EUR","NAT")))

famsize <- FAMSIZE
for (i in seq_len(n / famsize)) {
  fam = MakePedigree3(MAFs, DIRICHLET_ALPHA)
  rng = (famsize * (i - 1) + 1):(famsize * i)
  GTot[rng, ]  = fam$GenoMatTot
  GAfr[rng, ]  = fam$GenoMatAfr
  GEur[rng, ]  = fam$GenoMatEur
  GNat[rng, ]  = fam$GenoMatNat
  LAafr[rng, ] = fam$LAMat_Afr
  LAeur[rng, ] = fam$LAMat_Eur
  LAnat[rng, ] = fam$LAMat_Nat
  Admprop[rng, ] = fam$Admprop
}

rownames(GTot) = rownames(GAfr) = rownames(GEur) = rownames(GNat) = SampleID
rownames(LAafr) = rownames(LAeur) = rownames(LAnat) = SampleID
rownames(Admprop) = SampleID

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

cat("Finished Rel block", No, "for mode", mode, "\n")
