## Generate per-population allele frequencies under the Balding-Nichols model
## for the 3-way AFR/EUR/NAT simulation.
##
## Usage:  Rscript Pedigree3way_MAF.R <common|lowfreq>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Pedigree3way_MAF.R <common|lowfreq>")
mode <- args[1]

source(file.path(Sys.getenv("FELIX_SIM_BASE", unset = stop("Set FELIX_SIM_BASE to the simulation directory before running this script")), "R", "config.R"))
set_mode(mode)
source_utils()

set.seed(12345 + ifelse(mode == "common", 0, 1))

MAFs <- GenerateMafs3(
  p       = N_VARIANTS,
  fst_afr = 0.15, fst_eur = 0.10, fst_nat = 0.15,
  runif1  = ANC_P_LOW, runif2 = ANC_P_HIGH
)

out <- data.frame(
  SNP     = paste0("SNP", seq_len(N_VARIANTS)),
  MAF_AFR = MAFs$maf_afr,
  MAF_EUR = MAFs$maf_eur,
  MAF_NAT = MAFs$maf_nat,
  MAF_ANC = MAFs$maf_anc
)

write.table(out, file.path(DATA_DIR, "SNP_AF.tsv"),
            quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
cat("Mode:", mode, "   wrote", nrow(out), "variants to",
    file.path(DATA_DIR, "SNP_AF.tsv"), "\n")
