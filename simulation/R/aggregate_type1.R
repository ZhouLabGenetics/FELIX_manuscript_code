## Summarise FELIX step2 output on NULL phenotypes.
##
## For each trait type (quant, bin01, bin10), concatenate p-values from ALL
## seeds into one file so the user can make QQ plots downstream.
## Also produces the threshold-based summary table.
##
## Output per mode:
##   FP/<mode>/pvals_<trait>.tsv        — all p-values, all seeds pooled
##   FP/<mode>/type1_summary.tsv        — empirical Type-I error rates
##
## Usage:  Rscript aggregate_type1.R <common|lowfreq>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: aggregate_type1.R <common|lowfreq>")
mode <- args[1]

source(file.path(Sys.getenv("FELIX_SIM_BASE", unset = stop("Set FELIX_SIM_BASE to the simulation directory before running this script")), "R", "config.R"))
set_mode(mode)
library(data.table)

in_root  <- file.path(BASE, "data", mode, "saige_out", "null")
out_root <- file.path(BASE, "FP",   mode)
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

## The 6 conditioned p-value columns from FELIX output.
P_COLS <- c("p.value_c_anc1", "p.value_c_anc2", "p.value_c_anc3",
            "P_het_admixed_c", "P_hom_admixed_c", "P_cct_admixed_c")

thresholds <- c(5e-8, 1e-5, 1e-4, 1e-3, 1e-2, 5e-2)

pheno_dirs <- list.dirs(in_root, recursive = FALSE, full.names = TRUE)
if (!length(pheno_dirs)) stop("No SAIGE outputs found under ", in_root)

## Group pheno directories by trait.
tags   <- basename(pheno_dirs)
traits <- sub("_seed.*$", "", tags)

for (tr in unique(traits)) {
  idx <- which(traits == tr)
  cat("Trait:", tr, " (", length(idx), "seeds )\n")

  chunk_list <- list()
  for (j in idx) {
    tag  <- tags[j]
    seed <- as.integer(sub("^.*_seed", "", tag))
    files <- list.files(pheno_dirs[j],
                        pattern = "^chr.*\\.SAIGE\\.txt$",
                        full.names = TRUE)
    if (!length(files)) next
    dat <- rbindlist(lapply(files, fread), fill = TRUE)

    ## Keep only the columns we care about.
    keep_cols <- intersect(c("CHR", "POS", "MarkerID", P_COLS), colnames(dat))
    if (!length(intersect(P_COLS, colnames(dat)))) {
      warning("No expected p-value columns in ", pheno_dirs[j])
      next
    }
    sub_dat <- dat[, ..keep_cols]
    sub_dat[, seed := seed]
    chunk_list[[length(chunk_list) + 1]] <- sub_dat
  }

  if (!length(chunk_list)) {
    warning("No data collected for trait ", tr)
    next
  }

  pooled <- rbindlist(chunk_list, fill = TRUE)

  ## Write the pooled p-value file (one per trait, all seeds).
  pval_file <- file.path(out_root, paste0("pvals_", tr, ".tsv"))
  fwrite(pooled, pval_file, sep = "\t", quote = FALSE)
  cat("  Wrote", nrow(pooled), "rows to", pval_file, "\n")
}

## --- Summary table: empirical Type-I error per (trait, pcol, threshold) ---
summary_rows <- list()
for (tr in unique(traits)) {
  pval_file <- file.path(out_root, paste0("pvals_", tr, ".tsv"))
  if (!file.exists(pval_file)) next
  pooled <- fread(pval_file)

  present_pcols <- intersect(P_COLS, colnames(pooled))
  seeds <- sort(unique(pooled$seed))
  for (pc in present_pcols) {
    ## Per-seed empirical rates.
    seed_rates <- pooled[, .(dummy = 1), by = seed]
    for (th in thresholds) {
      per_seed <- pooled[!is.na(get(pc)),
                         .(n_tests = .N,
                           n_sig   = sum(get(pc) < th),
                           emp_rate = mean(get(pc) < th)),
                         by = seed]
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        trait     = tr,
        pcol      = pc,
        threshold = th,
        n_seeds   = nrow(per_seed),
        mean_rate = mean(per_seed$emp_rate),
        se_rate   = sd(per_seed$emp_rate) / sqrt(nrow(per_seed)),
        stringsAsFactors = FALSE
      )
    }
  }
}

if (length(summary_rows)) {
  summ <- rbindlist(summary_rows)
  summ_file <- file.path(out_root, "type1_summary.tsv")
  fwrite(summ, summ_file, sep = "\t", quote = FALSE)
  cat("Wrote", summ_file, "\n")
}
