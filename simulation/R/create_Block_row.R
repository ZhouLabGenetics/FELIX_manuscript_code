## Reorganise per-sample blocks into "block rows" (all 5000 samples, ~50k variants each).
## Tan et al's original layout: 50 blocks x (100 samples, 1M variants)   ->
##                              20 block_rows x (5000 samples, 50k variants).
##
## Usage:  Rscript create_Block_row.R <No> <common|lowfreq>
##   No = 1..20 (block_row index)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: create_Block_row.R <No> <common|lowfreq>")
No   <- as.integer(args[1])
mode <- args[2]

source(file.path(Sys.getenv("FELIX_SIM_BASE", unset = stop("Set FELIX_SIM_BASE to the simulation directory before running this script")), "R", "config.R"))
set_mode(mode)
library(data.table)

brdir <- file.path(DATA_DIR, paste0("Block_row", No))
dir.create(brdir, recursive = TRUE, showWarnings = FALSE)

## Variant index range for this block row.
row_start <- (No - 1) * BLOCKROW_P + 1L
row_end   <- No * BLOCKROW_P

## Read a given file (GTot / GAfr / ... ) from all 25 Ind blocks + 25 Rel blocks,
## subset the rows [row_start:row_end], and cbind by sample columns in the order
## Ind Block1..25 | Rel Block1..25 (which makes Sample1..5000 line up).
concat_block_row <- function(basename) {
  ind_parts <- vector("list", N_BLOCK)
  rel_parts <- vector("list", N_BLOCK)

  for (b in seq_len(N_BLOCK)) {
    ind_file <- file.path(IND_DIR, paste0("Block", b), paste0(basename, ".tsv"))
    rel_file <- file.path(REL_DIR, paste0("Block", b), paste0(basename, ".tsv"))

    ## Tan's blocks are variants x samples matrices WITH header (sample IDs).
    ## We only need the requested row range; use fread with skip/nrows for speed.
    ## +1 for the header line.
    ind_parts[[b]] <- fread(ind_file, sep = "\t",
                            skip = row_start, nrows = BLOCKROW_P,
                            header = FALSE)
    rel_parts[[b]] <- fread(rel_file, sep = "\t",
                            skip = row_start, nrows = BLOCKROW_P,
                            header = FALSE)
  }

  ind_mat <- do.call(cbind, ind_parts)
  rel_mat <- do.call(cbind, rel_parts)
  full    <- cbind(ind_mat, rel_mat)
  setnames(full, paste0("Sample", 1:ncol(full)))
  full
}

## Build SNP metadata once (CHR, POS, SNP, REF, ALT).
## CHR is set to the block_row index so SAIGE step2 can iterate via --chrom=<No>.
SNP_meta <- data.table(
  CHR = No,
  POS = row_start:row_end,
  SNP = paste0("SNP", row_start:row_end),
  REF = "A",
  ALT = "G"
)

for (base in c("GTot", "GAfr", "GEur", "GNat", "LAafr", "LAeur", "LAnat")) {
  mat <- concat_block_row(base)
  out <- cbind(SNP_meta, mat)
  fwrite(out, file.path(brdir, paste0(base, ".tsv")),
         sep = "\t", quote = FALSE, row.names = FALSE)
  cat("Wrote", base, "for Block_row", No, "(", nrow(out), "variants )\n")
}
