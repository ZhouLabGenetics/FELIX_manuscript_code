## QQ plots + genomic control (lambda_GC) for Tractor-Mix null runs.
## Walks data/<mode>/tractormix_out/null/<trait>_seed*/result.tsv, pools
## p-values per trait across seeds, and produces:
##   Bench/<mode>/tractormix_null_lambda.tsv              -- lambda_GC table
##   Bench/<mode>/tractormix_null_qq_<trait>.pdf          -- one QQ per trait
##
## Each QQ panel includes the joint P (if present) and the three per-ancestry
## p-values (Pval_anc0/1/2). Inflation appears as deflection above the y=x line.
##
## Optional SAIGE chr1 overlay: pass --with-saige on the command line, and
## the script will also pool data/<mode>/saige_out/null/<trait>_seed*/chr1.SAIGE.txt
## p.value_c_anc{1,2,3} + P_het/hom/cct columns for comparison.
##
## Usage:
##   Rscript plot_tractormix_qq_null.R <common|lowfreq>
##   Rscript plot_tractormix_qq_null.R <common|lowfreq> --with-saige

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1)
  stop("Usage: plot_tractormix_qq_null.R <common|lowfreq> [--with-saige]")
mode       <- args[1]
with_saige <- "--with-saige" %in% args

source(file.path(Sys.getenv("FELIX_SIM_BASE", unset = stop("Set FELIX_SIM_BASE to the simulation directory before running this script")), "R", "config.R"))
set_mode(mode)
suppressPackageStartupMessages({
  library(data.table)
})

tm_root    <- file.path(DATA_DIR, "tractormix_out", "null")
saige_root <- file.path(DATA_DIR, "saige_out",     "null")
bench_dir  <- file.path(BASE, "Bench", mode)
dir.create(bench_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(dir.exists(tm_root))

TRAITS    <- c("quant", "bin10", "bin01")
TM_PCOLS  <- c(P_joint   = "P",
               P_anc1    = "Pval_anc0",
               P_anc2    = "Pval_anc1",
               P_anc3    = "Pval_anc2")
SAIGE_PCOLS <- c(saige_p_anc1 = "p.value_c_anc1",
                 saige_p_anc2 = "p.value_c_anc2",
                 saige_p_anc3 = "p.value_c_anc3",
                 saige_p_het  = "P_het_admixed_c",
                 saige_p_hom  = "P_hom_admixed_c",
                 saige_p_cct  = "P_cct_admixed_c")

## ---- helpers -----------------------------------------------------------
lambda_gc <- function(p) {
  p <- p[is.finite(p) & p > 0 & p <= 1]
  if (!length(p)) return(NA_real_)
  qchisq(median(p), df = 1, lower.tail = FALSE) / qchisq(0.5, df = 1)
}

## Thin a vector of p-values for QQ plotting:
##   keep all P < tail_cut (tail), subsample bulk to n_bulk.
thin_for_qq <- function(p, tail_cut = 0.01, n_bulk = 2000, seed = 1) {
  p <- p[is.finite(p) & p > 0 & p <= 1]
  if (!length(p)) return(numeric(0))
  tail <- p[p < tail_cut]
  bulk <- p[p >= tail_cut]
  if (length(bulk) > n_bulk) {
    set.seed(seed)
    bulk <- sample(bulk, n_bulk)
  }
  sort(c(tail, bulk))
}

qq_points <- function(p_sorted, n_total) {
  ## Use n_total (true count, not thinned) for expected quantiles so the line
  ## reflects the actual hypothesis count.
  if (!length(p_sorted)) return(NULL)
  ranks_in_total <- (n_total - rev(seq_along(p_sorted)) + 1)
  exp_p <- (rev(seq_along(p_sorted)) - 0.5) / n_total
  data.table(obs_neglog10 = -log10(p_sorted),
             exp_neglog10 = -log10(exp_p))
}

## ---- gather TM p-values per trait ------------------------------------
gather_tm <- function(trait) {
  dirs <- list.files(tm_root, pattern = paste0("^", trait, "_seed[0-9]+$"),
                     full.names = TRUE)
  out_cols <- names(TM_PCOLS)
  pools <- setNames(replicate(length(out_cols), list(), simplify = FALSE),
                    out_cols)
  n_seeds <- 0L
  for (d in dirs) {
    f <- file.path(d, "result.tsv")
    if (!file.exists(f) || file.size(f) == 0) next
    res <- tryCatch(fread(f), error = function(e) NULL)
    if (is.null(res)) next
    n_seeds <- n_seeds + 1L
    for (i in seq_along(TM_PCOLS)) {
      col <- TM_PCOLS[[i]]
      if (col %in% names(res))
        pools[[i]] <- c(pools[[i]], list(suppressWarnings(as.numeric(res[[col]]))))
    }
  }
  pools <- lapply(pools, function(x) unlist(x, use.names = FALSE))
  list(pools = pools, n_seeds = n_seeds)
}

gather_saige <- function(trait) {
  dirs <- list.files(saige_root, pattern = paste0("^", trait, "_seed[0-9]+$"),
                     full.names = TRUE)
  out_cols <- names(SAIGE_PCOLS)
  pools <- setNames(replicate(length(out_cols), list(), simplify = FALSE),
                    out_cols)
  n_seeds <- 0L
  for (d in dirs) {
    f <- file.path(d, "chr1.SAIGE.txt")     ## chr1 only -> apples-to-apples vs TM
    if (!file.exists(f) || file.size(f) == 0) next
    sd <- tryCatch(fread(f), error = function(e) NULL)
    if (is.null(sd)) next
    n_seeds <- n_seeds + 1L
    for (i in seq_along(SAIGE_PCOLS)) {
      col <- SAIGE_PCOLS[[i]]
      if (col %in% names(sd))
        pools[[i]] <- c(pools[[i]], list(suppressWarnings(as.numeric(sd[[col]]))))
    }
  }
  pools <- lapply(pools, function(x) unlist(x, use.names = FALSE))
  list(pools = pools, n_seeds = n_seeds)
}

## ---- main loop --------------------------------------------------------
lambda_rows <- list()
TM_COLORS  <- c(P_joint = "black", P_anc1 = "#E41A1C",
                P_anc2  = "#377EB8", P_anc3 = "#4DAF4A")
S_COLORS   <- c(saige_p_anc1 = "#E41A1C", saige_p_anc2 = "#377EB8",
                saige_p_anc3 = "#4DAF4A", saige_p_het  = "purple",
                saige_p_hom  = "orange", saige_p_cct  = "black")

for (trait in TRAITS) {
  cat("==>", trait, "\n")
  tm <- gather_tm(trait)
  if (tm$n_seeds == 0L) {
    cat("   no TM result.tsv found, skipping\n"); next
  }
  cat("   TM seeds pooled:", tm$n_seeds, "\n")

  ## lambda_GC per test
  for (nm in names(tm$pools)) {
    lambda_rows[[length(lambda_rows) + 1L]] <- data.table(
      mode = mode, source = "Tractor-Mix", trait = trait,
      test = nm,
      n_pvals = sum(is.finite(tm$pools[[nm]]) & tm$pools[[nm]] > 0),
      lambda_gc = lambda_gc(tm$pools[[nm]])
    )
  }

  saige <- NULL
  if (with_saige) {
    saige <- gather_saige(trait)
    cat("   SAIGE chr1 seeds pooled:", saige$n_seeds, "\n")
    if (saige$n_seeds > 0L) {
      for (nm in names(saige$pools)) {
        lambda_rows[[length(lambda_rows) + 1L]] <- data.table(
          mode = mode, source = "FELIX", trait = trait,
          test = nm,
          n_pvals = sum(is.finite(saige$pools[[nm]]) & saige$pools[[nm]] > 0),
          lambda_gc = lambda_gc(saige$pools[[nm]])
        )
      }
    }
  }

  ## ---- QQ plot --------------------------------------------------------
  pdf_f <- file.path(bench_dir, paste0("tractormix_null_qq_", trait, ".pdf"))
  pdf(pdf_f, width = 7, height = 7, useDingbats = FALSE)
  on.exit(try(dev.off(), silent = TRUE), add = TRUE)

  ## Determine plot range from the most extreme observed point.
  max_obs <- 0
  for (nm in names(tm$pools)) {
    p <- tm$pools[[nm]]; p <- p[is.finite(p) & p > 0]
    if (length(p)) max_obs <- max(max_obs, -log10(min(p)))
  }
  if (with_saige && !is.null(saige) && saige$n_seeds > 0L) {
    for (nm in names(saige$pools)) {
      p <- saige$pools[[nm]]; p <- p[is.finite(p) & p > 0]
      if (length(p)) max_obs <- max(max_obs, -log10(min(p)))
    }
  }
  axis_max <- max(6, ceiling(max_obs) + 0.5)

  plot(NA, xlim = c(0, axis_max), ylim = c(0, axis_max),
       xlab = expression(Expected~-log[10](P)),
       ylab = expression(Observed~-log[10](P)),
       main = sprintf("%s -- %s -- null (Type-I)\n%d TM seeds pooled",
                      mode, trait, tm$n_seeds),
       las = 1)
  abline(0, 1, lty = 1, col = "grey60")

  ## TM curves: solid
  for (nm in names(tm$pools)) {
    p <- tm$pools[[nm]]
    p_finite <- p[is.finite(p) & p > 0 & p <= 1]
    if (!length(p_finite)) next
    qq <- qq_points(thin_for_qq(p_finite), length(p_finite))
    if (is.null(qq)) next
    lines(qq$exp_neglog10, qq$obs_neglog10, type = "p", pch = 16, cex = 0.4,
          col = TM_COLORS[nm])
  }

  ## SAIGE curves: triangles
  if (with_saige && !is.null(saige) && saige$n_seeds > 0L) {
    for (nm in names(saige$pools)) {
      p <- saige$pools[[nm]]
      p_finite <- p[is.finite(p) & p > 0 & p <= 1]
      if (!length(p_finite)) next
      qq <- qq_points(thin_for_qq(p_finite, seed = 2), length(p_finite))
      if (is.null(qq)) next
      lines(qq$exp_neglog10, qq$obs_neglog10, type = "p", pch = 2, cex = 0.35,
            col = S_COLORS[nm])
    }
  }

  ## Legend with lambda_GC values
  leg_txt <- character(0); leg_col <- character(0); leg_pch <- integer(0)
  for (nm in names(tm$pools)) {
    lam <- lambda_gc(tm$pools[[nm]])
    leg_txt <- c(leg_txt, sprintf("TM %s  lambda=%.3f", nm, lam))
    leg_col <- c(leg_col, TM_COLORS[nm])
    leg_pch <- c(leg_pch, 16)
  }
  if (with_saige && !is.null(saige) && saige$n_seeds > 0L) {
    for (nm in names(saige$pools)) {
      lam <- lambda_gc(saige$pools[[nm]])
      leg_txt <- c(leg_txt, sprintf("SAIGE %s  lambda=%.3f",
                                    sub("^saige_p_", "", nm), lam))
      leg_col <- c(leg_col, S_COLORS[nm])
      leg_pch <- c(leg_pch, 2)
    }
  }
  legend("topleft", legend = leg_txt, col = leg_col, pch = leg_pch,
         bty = "n", cex = 0.7)
  dev.off()
  cat("   wrote:", pdf_f, "\n")
}

## ---- write lambda table ----------------------------------------------
lambda_tbl <- rbindlist(lambda_rows, fill = TRUE)
lambda_tbl[, lambda_gc := round(lambda_gc, 4)]
lambda_f <- file.path(bench_dir, "tractormix_null_lambda.tsv")
fwrite(lambda_tbl, lambda_f, sep = "\t", quote = FALSE)
cat("\nWrote lambda table:", lambda_f, "\n")
print(lambda_tbl)
