## Concatenate per-block Admprop files into a single Admprop.tsv (5000 x 3)
## with row order Sample1..Sample5000. Also writes the true pedigree kinship
## matrix (used later as a random-effect source for phenotype simulation).
##
## Usage:  Rscript aggregate_Admprop.R <common|lowfreq>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: aggregate_Admprop.R <common|lowfreq>")
mode <- args[1]

source(file.path(Sys.getenv("FELIX_SIM_BASE", unset = stop("Set FELIX_SIM_BASE to the simulation directory before running this script")), "R", "config.R"))
set_mode(mode)
source_utils()
library(data.table)

read_block <- function(dir, b) {
  f <- file.path(dir, paste0("Block", b), "Admprop.tsv")
  tbl <- read.table(f, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
  tbl
}

ind_list <- lapply(seq_len(N_BLOCK), read_block, dir = IND_DIR)
rel_list <- lapply(seq_len(N_BLOCK), read_block, dir = REL_DIR)

all_ind <- do.call(rbind, ind_list)
all_rel <- do.call(rbind, rel_list)
all_adm <- rbind(all_ind, all_rel)

stopifnot(nrow(all_adm) == N_TOTAL)
rownames(all_adm) <- paste0("Sample", 1:N_TOTAL)

out <- data.frame(ID = rownames(all_adm),
                  AFR = all_adm$AFR,
                  EUR = all_adm$EUR,
                  NAT = all_adm$NAT)
write.table(out, file.path(DATA_DIR, "Admprop.tsv"),
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
cat("Wrote", nrow(out), "rows to", file.path(DATA_DIR, "Admprop.tsv"), "\n")

## True pedigree kinship for phenotype random-effect simulation.
K <- MakeGRM3(N_IND, N_REL)   # 5000 x 5000, rownames/colnames Sample1..Sample5000
fwrite(as.data.table(K, keep.rownames = "ID"),
       file.path(DATA_DIR, "Kinship.tsv"),
       sep = "\t", quote = FALSE, row.names = FALSE)
cat("Wrote true kinship matrix to", file.path(DATA_DIR, "Kinship.tsv"), "\n")
