## Aggregate Tractor-Mix benchmark outputs into manuscript-ready tables.
##
## Walks data/<mode>/tractormix_out/alt/<TAG>/ and joins with
## - bench.log              (GNU /usr/bin/time -v: wall sec, max RSS kB)
## - result.tsv             (per-variant p-values + ancestry effect sizes)
## - SAIGE step2 output for the same TAG (data/<mode>/saige_out/alt/<TAG>/)
## - causal.txt             (causal SNP id)
##
## Outputs under Bench/<mode>/:
##   tractormix_runtime.tsv      one row per (trait, scen, beta, seed)
##   tractormix_causal.tsv       Tractor-Mix vs SAIGE causal-SNP pvalues
##   tractormix_qq_pool.tsv      non-causal variant pvalues pooled across seeds
##
## A status field flags DNF (no result.tsv but bench.log exists) and ERR
## (neither result nor bench, indicating crash before /usr/bin/time wrote).
##
## Usage:  Rscript aggregate_tractormix.R <common|lowfreq>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: aggregate_tractormix.R <common|lowfreq>")
mode <- args[1]

source(file.path(Sys.getenv("FELIX_SIM_BASE", unset = stop("Set FELIX_SIM_BASE to the simulation directory before running this script")), "R", "config.R"))
set_mode(mode)
library(data.table)

tm_root    <- file.path(DATA_DIR, "tractormix_out", "alt")
saige_root <- file.path(DATA_DIR, "saige_out", "alt")
bench_dir  <- file.path(BASE, "Bench", mode)
dir.create(bench_dir, recursive = TRUE, showWarnings = FALSE)

causal_id <- trimws(readLines(file.path(DATA_DIR, "pheno", "alt", "causal.txt"))[1])
cat("Causal SNP:", causal_id, "\n")

TAG_RE <- "^(quant|bin01|bin10)_(shared|afr|hetero)_beta([0-9]+)_seed([0-9]+)$"
tag_dirs <- list.dirs(tm_root, recursive = FALSE)
cat("Found", length(tag_dirs), "TAG dirs under", tm_root, "\n")

## ---- Helpers ----------------------------------------------------------
parse_bench_log <- function(path) {
  ## GNU time -v format: lines like "  Elapsed (wall clock) time (h:mm:ss or m:ss): 12:34.56"
  ## and "  Maximum resident set size (kbytes): 1234567"
  if (!file.exists(path) || file.size(path) == 0)
    return(list(wall_s = NA_real_, rss_kb = NA_real_))
  L <- readLines(path, warn = FALSE)
  wall_line <- grep("Elapsed \\(wall clock\\) time", L, value = TRUE)
  rss_line  <- grep("Maximum resident set size",     L, value = TRUE)
  parse_hms <- function(s) {
    s <- sub(".*: ", "", s)
    parts <- strsplit(s, ":")[[1]]
    if (length(parts) == 3)
      as.numeric(parts[1]) * 3600 + as.numeric(parts[2]) * 60 + as.numeric(parts[3])
    else if (length(parts) == 2)
      as.numeric(parts[1]) * 60 + as.numeric(parts[2])
    else as.numeric(parts[1])
  }
  list(
    wall_s = if (length(wall_line)) parse_hms(wall_line[1]) else NA_real_,
    rss_kb = if (length(rss_line))
      as.numeric(sub(".*: ", "", rss_line[1])) else NA_real_
  )
}

## Find the p-value column. Tractor-Mix typically writes "P" (joint) and
## "Pval_anc0/1/2". We prefer the joint P for the QQ + causal scatter, and
## also retain per-ancestry pvals.
PVAL_COLS <- c("P", "Pval_anc0", "Pval_anc1", "Pval_anc2")
EFF_COLS  <- c("Eff_anc0", "Eff_anc1", "Eff_anc2",
               "SE_anc0",  "SE_anc1",  "SE_anc2")

## ---- Walk ------------------------------------------------------------
runtime_rows <- list()
causal_rows  <- list()
qq_rows      <- list()

for (pd in tag_dirs) {
  tag <- basename(pd)
  m <- regmatches(tag, regexec(TAG_RE, tag))[[1]]
  if (!length(m)) next
  trait <- m[2]; scen <- m[3]
  beta  <- as.numeric(m[4]) / 100
  seed  <- as.integer(m[5])

  bench <- parse_bench_log(file.path(pd, "bench.log"))
  result_f <- file.path(pd, "result.tsv")
  has_result <- file.exists(result_f) && file.size(result_f) > 0

  status <- if (has_result) "OK"
            else if (!is.na(bench$wall_s)) "DNF"
            else "ERR"

  runtime_rows[[length(runtime_rows) + 1L]] <- data.table(
    mode = mode, trait = trait, scenario = scen, beta = beta, seed = seed,
    status = status, wall_s = bench$wall_s, rss_kb = bench$rss_kb,
    tag = tag
  )

  if (!has_result) next
  res <- fread(result_f)
  ## Tractor-Mix uses ID for marker column; standardize to "MarkerID".
  id_col <- intersect(c("ID","MarkerID","SNP","rsID"), names(res))[1]
  if (is.na(id_col)) {
    warning("No ID column in ", result_f, "; skipping causal/qq for ", tag)
    next
  }
  setnames(res, id_col, "MarkerID")
  present_pcols <- intersect(PVAL_COLS, names(res))
  if (!length(present_pcols)) {
    warning("No pvalue cols in ", result_f, "; skipping ", tag); next
  }

  ## --- Causal hit row + paired SAIGE pvals --------------------
  hit <- res[MarkerID == causal_id]
  if (nrow(hit)) {
    saige_f <- file.path(saige_root, tag, "chr1.SAIGE.txt")
    saige_p <- NA_real_; saige_p_het <- NA_real_; saige_p_hom <- NA_real_; saige_p_cct <- NA_real_
    if (file.exists(saige_f) && file.size(saige_f) > 0) {
      sd <- tryCatch(fread(saige_f), error = function(e) NULL)
      if (!is.null(sd) && "MarkerID" %in% names(sd)) {
        sh <- sd[MarkerID == causal_id]
        if (nrow(sh)) {
          saige_p     <- suppressWarnings(as.numeric(sh$p.value_c_anc1[1]))
          saige_p_het <- suppressWarnings(as.numeric(sh$P_het_admixed_c[1]))
          saige_p_hom <- suppressWarnings(as.numeric(sh$P_hom_admixed_c[1]))
          saige_p_cct <- suppressWarnings(as.numeric(sh$P_cct_admixed_c[1]))
        }
      }
    }
    causal_rows[[length(causal_rows) + 1L]] <- data.table(
      mode = mode, trait = trait, scenario = scen, beta = beta, seed = seed,
      tractormix_P    = suppressWarnings(as.numeric(hit$P[1])),
      tractormix_Panc0 = suppressWarnings(as.numeric(hit$Pval_anc0[1])),
      tractormix_Panc1 = suppressWarnings(as.numeric(hit$Pval_anc1[1])),
      tractormix_Panc2 = suppressWarnings(as.numeric(hit$Pval_anc2[1])),
      saige_p_anc1 = saige_p, saige_p_het = saige_p_het,
      saige_p_hom = saige_p_hom, saige_p_cct = saige_p_cct
    )
  }

  ## --- Non-causal pool for QQ ----------------------------------
  nc <- res[MarkerID != causal_id, c("MarkerID", present_pcols), with = FALSE]
  nc[, `:=`(mode = mode, trait = trait, scenario = scen,
            beta = beta, seed = seed)]
  qq_rows[[length(qq_rows) + 1L]] <- nc
}

runtime <- rbindlist(runtime_rows, fill = TRUE)
fwrite(runtime, file.path(bench_dir, "tractormix_runtime.tsv"),
       sep = "\t", quote = FALSE)
cat("Wrote tractormix_runtime.tsv  (", nrow(runtime), " rows)\n", sep = "")

if (length(causal_rows)) {
  causal <- rbindlist(causal_rows, fill = TRUE)
  fwrite(causal, file.path(bench_dir, "tractormix_causal.tsv"),
         sep = "\t", quote = FALSE)
  cat("Wrote tractormix_causal.tsv  (", nrow(causal), " rows)\n", sep = "")
}
if (length(qq_rows)) {
  qq <- rbindlist(qq_rows, fill = TRUE)
  fwrite(qq, file.path(bench_dir, "tractormix_qq_pool.tsv"),
         sep = "\t", quote = FALSE)
  cat("Wrote tractormix_qq_pool.tsv (", nrow(qq), " rows)\n", sep = "")
}

## ---- Quick summary -----------------------------------------------
cat("\nStatus tally:\n")
print(runtime[, .N, by = .(trait, scenario, status)])
cat("\nMedian wall_s and max rss_kb by (trait, scen) for OK runs:\n")
print(runtime[status == "OK",
              .(median_wall_min = round(median(wall_s, na.rm = TRUE) / 60, 1),
                max_rss_gb      = round(max(rss_kb,  na.rm = TRUE) / 1e6, 2),
                n = .N),
              by = .(trait, scenario)])
