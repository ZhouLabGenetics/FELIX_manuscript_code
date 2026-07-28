## Summarise FELIX step2 output on ALT phenotypes into power tables.
##
## Manifest-agnostic: walks data/<mode>/saige_out/alt/<TAG>/ directly. The
## TAG name encodes (trait, scenario, beta, seed) and is parsed via regex.
## The causal SNP is read from data/<mode>/pheno/alt/causal.txt; if that
## file is missing, falls back to causal_snp column of any manifest.tsv
## (works for both the old wide and the new long format).
##
## This means you can edit BETA_GRIDS for some conditions, rerun step 11
## for only those, and aggregation still works correctly across the full
## set of SAIGE outputs that exist on disk.
##
## Output per mode:
##   Power/<mode>/causal_<trait>_<scenario>_beta<B>.tsv  -- per-seed detail
##   Power/<mode>/power_all.tsv                          -- master file
##   Power/<mode>/power_summary.tsv                      -- power + mean chisq
##
## Usage:  Rscript aggregate_power.R <common|lowfreq>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: aggregate_power.R <common|lowfreq>")
mode <- args[1]

source("/data/wzhougroup/lhu/saige_tractor/simulation/3way/scripts/R/config.R")
set_mode(mode)
library(data.table)

in_root  <- file.path(BASE, "data", mode, "saige_out", "alt")
out_root <- file.path(BASE, "Power",    mode)
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

## ---- Resolve causal SNP ID ----
causal_file   <- file.path(DATA_DIR, "pheno", "alt", "causal.txt")
manifest_file <- file.path(DATA_DIR, "pheno", "alt", "manifest.tsv")

if (file.exists(causal_file)) {
  causal_id <- trimws(readLines(causal_file)[1])
  cat("Causal SNP (from causal.txt):", causal_id, "\n")
} else if (file.exists(manifest_file)) {
  man <- fread(manifest_file)
  if (!"causal_snp" %in% colnames(man))
    stop("manifest.tsv has no causal_snp column; cannot resolve causal SNP. ",
         "Rerun generate_alt_pheno.R to write causal.txt.")
  causal_id <- as.character(man$causal_snp[1])
  cat("Causal SNP (from manifest fallback):", causal_id, "\n")
} else {
  stop("Neither causal.txt nor manifest.tsv found in ",
       file.path(DATA_DIR, "pheno", "alt"))
}

P_COLS <- c("p.value_c_anc1", "p.value_c_anc2", "p.value_c_anc3",
            "P_het_admixed_c", "P_hom_admixed_c", "P_cct_admixed_c")
P_DF <- c(p.value_c_anc1 = 1, p.value_c_anc2 = 1, p.value_c_anc3 = 1,
          P_het_admixed_c = 2, P_hom_admixed_c = 1, P_cct_admixed_c = NA)
BETA_COLS <- c("BETA_c_anc1", "BETA_c_anc2", "BETA_c_anc3", "BETA_c_ancALL")
thresholds <- c(5e-8, 1e-5, 1e-4, 1e-3, 1e-2)

## ---- Walk TAG dirs and collect causal-variant rows ----
TAG_RE <- "^(quant|bin01|bin10)_(shared|afr|hetero)_beta([0-9]+)_seed([0-9]+)$"

tag_dirs <- list.dirs(in_root, recursive = FALSE)
cat("Found", length(tag_dirs), "TAG directories under", in_root, "\n")

all_rows <- list()
n_no_match <- 0L
n_no_saige <- 0L
n_no_causal <- 0L

for (pd in tag_dirs) {
  tag <- basename(pd)
  m <- regmatches(tag, regexec(TAG_RE, tag))[[1]]
  if (!length(m)) { n_no_match <- n_no_match + 1L; next }

  trait_name <- m[2]
  scen       <- m[3]
  beta_val   <- as.numeric(m[4]) / 100
  seed_val   <- as.integer(m[5])

  f <- list.files(pd, pattern = "\\.SAIGE\\.txt$", full.names = TRUE)
  if (!length(f)) { n_no_saige <- n_no_saige + 1L; next }

  dat <- rbindlist(lapply(f, fread), fill = TRUE)
  if (!"MarkerID" %in% colnames(dat)) { n_no_saige <- n_no_saige + 1L; next }

  hit <- dat[MarkerID == causal_id]
  if (!nrow(hit)) { n_no_causal <- n_no_causal + 1L; next }

  row <- data.table(
    mode       = mode,
    trait      = trait_name,
    scenario   = scen,
    beta       = beta_val,
    seed       = seed_val,
    causal_snp = causal_id,
    CHR        = hit$CHR[1],
    POS        = hit$POS[1]
  )
  ## Coerce numeric here so a truncated SAIGE file (last line half-written
  ## -> fread infers character) doesn't poison the rbindlist downstream.
  for (pc in intersect(P_COLS, colnames(hit))) {
    v <- hit[[pc]][1]
    row[[pc]] <- if (is.numeric(v)) v else suppressWarnings(as.numeric(v))
  }
  for (bc in intersect(BETA_COLS, colnames(hit))) {
    v <- hit[[bc]][1]
    row[[bc]] <- if (is.numeric(v)) v else suppressWarnings(as.numeric(v))
  }

  all_rows[[length(all_rows) + 1L]] <- row
}

cat(sprintf("Skipped: %d unparseable TAGs, %d missing SAIGE output, %d missing causal hit\n",
            n_no_match, n_no_saige, n_no_causal))
if (!length(all_rows))
  stop("No causal variant data collected. Check ", in_root)

per_seed <- rbindlist(all_rows, fill = TRUE)

## ---- Add chi-square / -log10p columns ----
## Use a local numeric vector v rather than relying on per_seed[[pc]] being
## updated in place — earlier we saw qchisq still get a non-numeric input
## even after [[<- coercion. set() is data.table-canonical for column writes.
present_pcols <- intersect(P_COLS, colnames(per_seed))
for (pc in present_pcols) {
  v <- per_seed[[pc]]
  cat(sprintf("Column %s: typeof=%s, is.numeric=%s, n=%d\n",
              pc, typeof(v), is.numeric(v), length(v)))
  if (!is.numeric(v)) {
    n_pre <- sum(!is.na(v))
    v <- suppressWarnings(as.numeric(as.character(v)))
    n_post <- sum(!is.na(v))
    cat(sprintf("  -> coerced %s to numeric; %d values became NA.\n",
                pc, n_pre - n_post))
    set(per_seed, j = pc, value = v)
  }
  df <- P_DF[pc]
  chisq_col <- sub("^p\\.value_c_", "chisq_", sub("^P_", "chisq_", pc))
  if (!is.na(df)) {
    set(per_seed, j = chisq_col,
        value = qchisq(v, df = df, lower.tail = FALSE))
  } else {
    set(per_seed, j = paste0("neglog10p_", pc),
        value = -log10(v))
  }
}

## ---- Master file ----
fwrite(per_seed, file.path(out_root, "power_all.tsv"), sep = "\t", quote = FALSE)
cat("Wrote", nrow(per_seed), "rows to power_all.tsv\n")

## ---- Per-condition detail files ----
conditions <- unique(per_seed[, .(trait, scenario, beta)])
for (r in seq_len(nrow(conditions))) {
  tr <- conditions$trait[r]; sc <- conditions$scenario[r]; bt <- conditions$beta[r]
  sub <- per_seed[trait == tr & scenario == sc & beta == bt]
  fname <- sprintf("causal_%s_%s_beta%03d.tsv", tr, sc, round(bt * 100))
  fwrite(sub, file.path(out_root, fname), sep = "\t", quote = FALSE)
}
cat("Wrote", nrow(conditions), "per-condition detail files\n")

## ---- Power summary ----
summary_rows <- list()
for (r in seq_len(nrow(conditions))) {
  tr <- conditions$trait[r]; sc <- conditions$scenario[r]; bt <- conditions$beta[r]
  sub <- per_seed[trait == tr & scenario == sc & beta == bt]

  for (pc in present_pcols) {
    pvs <- sub[[pc]]; pvs_ok <- pvs[!is.na(pvs)]
    if (!length(pvs_ok)) next

    df <- P_DF[pc]
    if (!is.na(df)) {
      chisq_col  <- sub("^p\\.value_c_", "chisq_", sub("^P_", "chisq_", pc))
      mean_chisq <- mean(sub[[chisq_col]], na.rm = TRUE)
    } else {
      mean_chisq <- mean(-log10(pvs_ok))
    }

    for (th in thresholds) {
      summary_rows[[length(summary_rows) + 1L]] <- data.frame(
        trait      = tr, scenario = sc, beta = bt,
        pcol       = pc, threshold = th,
        n_seeds    = length(pvs_ok),
        n_hit      = sum(pvs_ok < th),
        power      = mean(pvs_ok < th),
        mean_chisq = mean_chisq,
        stringsAsFactors = FALSE
      )
    }
  }
}

summ <- rbindlist(summary_rows)
fwrite(summ, file.path(out_root, "power_summary.tsv"), sep = "\t", quote = FALSE)
cat("Wrote power_summary.tsv\n")
