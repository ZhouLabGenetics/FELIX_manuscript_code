## Convert Block_row1's per-ancestry dosage tables into Tractor-Mix hapdose
## files. Tractor-Mix expects one file per ancestry, layout:
##
##   CHR  POS  ID  REF  ALT  <sample1>  <sample2>  ...
##
## Your simulation's GAfr.tsv / GEur.tsv / GNat.tsv already have this layout
## except the 3rd metadata column is named "SNP" instead of "ID". This script
## just renames the column and writes to data/<mode>/tractormix/hapdose/.
##
## The hapdose files are phenotype-independent -- generate once per mode.
##
## Output:
##   data/<mode>/tractormix/hapdose/chr1.anc0.dosage.txt   (AFR)
##   data/<mode>/tractormix/hapdose/chr1.anc1.dosage.txt   (EUR)
##   data/<mode>/tractormix/hapdose/chr1.anc2.dosage.txt   (NAT)
##
## Convention: anc0/1/2 = AFR/EUR/NAT (same order as your VCF DS1/DS2/DS3).
##
## Usage:  Rscript make_tractormix_hapdose.R <common|lowfreq>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: make_tractormix_hapdose.R <common|lowfreq>")
mode <- args[1]

source(file.path(Sys.getenv("FELIX_SIM_BASE", unset = stop("Set FELIX_SIM_BASE to the simulation directory before running this script")), "R", "config.R"))
set_mode(mode)
library(data.table)

brdir   <- file.path(DATA_DIR, "Block_row1")
out_dir <- file.path(DATA_DIR, "tractormix", "hapdose")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_one <- function(src_base, anc_idx) {
  dt <- fread(file.path(brdir, paste0(src_base, ".tsv")), sep = "\t")
  ## Rename SNP -> ID. If the file has a different 3rd column name, force it.
  setnames(dt, names(dt)[3], "ID")
  ## Reorder to canonical 5 metadata cols + samples.
  meta <- c("CHR", "POS", "ID", "REF", "ALT")
  if (!all(meta %in% names(dt)))
    stop("Missing metadata column(s) in ", src_base, ": ",
         paste(setdiff(meta, names(dt)), collapse = ", "))
  setcolorder(dt, c(meta, setdiff(names(dt), meta)))
  out <- file.path(out_dir, sprintf("chr1.anc%d.dosage.txt", anc_idx))
  fwrite(dt, out, sep = "\t", quote = FALSE)
  cat(sprintf("Wrote %s  (%d variants, %d samples)\n",
              out, nrow(dt), ncol(dt) - 5))
}

write_one("GAfr", 0L)
write_one("GEur", 1L)
write_one("GNat", 2L)

cat("Done. Tractor-Mix hapdose files in ", out_dir, "\n", sep = "")
