## Shared paths / constants for the 3-way admixture simulation.
## Every R script sources this file with a fully-qualified path, e.g.
##   source("/data/wzhougroup/lhu/saige_tractor/simulation/3way/scripts/R/config.R")
## and then calls set_mode("common") or set_mode("lowfreq") to set paths
## for a given MAF regime.

BASE        <- "/data/wzhougroup/lhu/saige_tractor/simulation/3way"
SCRIPT_DIR  <- file.path(BASE, "scripts", "R")
LOG_DIR     <- file.path(BASE, "log")

## Simulation-wide constants
N_IND    <- 5000    # independent individuals
N_REL    <- 5000    # related individuals (500 families of 10)
BLOCK_N  <- 100     # samples per block
N_BLOCK  <- 50      # blocks per arm (ind and rel)
FAMSIZE  <- 10
N_TOTAL  <- N_IND + N_REL                 # 10000
N_BLOCKROWS <- 20   # see create_Block_row.R
BLOCKROW_P  <- 50000
DIRICHLET_ALPHA <- c(1.0, 5.5, 3.5)       # Latino-like AFR/EUR/NAT

set_mode <- function(mode) {
  stopifnot(mode %in% c("common", "lowfreq", "bench", "unr", "admix"))
  assign("MODE",        mode,                         envir = .GlobalEnv)

  if (mode == "admix") {
    ## Explicit-mixture admixed cohort mode. Driven by env vars:
    ##   ADMIX_DIR       -- absolute path to the data dir
    ##   ADMIX_N         -- sample size (positive int)
    ##   ADMIX_P         -- variant count (default 50000)
    ##   ADMIX_SCENARIO  -- one of: 2way_50_50 | 2way_25_75 | 3way_20_30_50
    admix_dir <- Sys.getenv("ADMIX_DIR", "")
    if (!nzchar(admix_dir))
      stop("set_mode('admix'): ADMIX_DIR env var must be set")
    admix_n <- as.integer(Sys.getenv("ADMIX_N", "0"))
    if (is.na(admix_n) || admix_n <= 0L)
      stop("set_mode('admix'): ADMIX_N env var (positive int) must be set")
    admix_p <- as.integer(Sys.getenv("ADMIX_P", "50000"))

    assign("DATA_DIR",     admix_dir,                     envir = .GlobalEnv)
    assign("IND_DIR",      file.path(admix_dir, "Ind"),   envir = .GlobalEnv)
    assign("REL_DIR",      file.path(admix_dir, "Rel"),   envir = .GlobalEnv)
    assign("BLOCKROW_DIR", admix_dir,                     envir = .GlobalEnv)
    assign("FP_COMMON_DIR", file.path(BASE, "FP",    "admix"), envir = .GlobalEnv)
    assign("POWER_DIR",     file.path(BASE, "Power", "admix"), envir = .GlobalEnv)
    assign("N_TOTAL",    admix_n,  envir = .GlobalEnv)
    assign("N_VARIANTS", admix_p,  envir = .GlobalEnv)
    assign("MAF_LOW",    0.05, envir = .GlobalEnv)
    assign("MAF_HIGH",   0.50, envir = .GlobalEnv)
    assign("ANC_P_LOW",  0.10, envir = .GlobalEnv)
    assign("ANC_P_HIGH", 0.90, envir = .GlobalEnv)
    dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)
    return(invisible(NULL))
  }

  if (mode == "unr") {
    ## Unrelated-cohort mode for FELIX benchmarking on
    ## independent admixed individuals. Driven by env vars:
    ##   UNR_DIR  -- absolute path to the data dir (e.g.
    ##               .../data/unr_N10000_P50000); created if missing.
    ##   UNR_N    -- sample size (positive int)
    ##   UNR_P    -- variant count (default 50000)
    unr_dir <- Sys.getenv("UNR_DIR", "")
    if (!nzchar(unr_dir))
      stop("set_mode('unr'): UNR_DIR env var must be set")
    unr_n <- as.integer(Sys.getenv("UNR_N", "0"))
    if (is.na(unr_n) || unr_n <= 0L)
      stop("set_mode('unr'): UNR_N env var (positive int) must be set")
    unr_p <- as.integer(Sys.getenv("UNR_P", "50000"))

    assign("DATA_DIR",     unr_dir,                       envir = .GlobalEnv)
    assign("IND_DIR",      file.path(unr_dir, "Ind"),     envir = .GlobalEnv)
    assign("REL_DIR",      file.path(unr_dir, "Rel"),     envir = .GlobalEnv)
    assign("BLOCKROW_DIR", unr_dir,                       envir = .GlobalEnv)
    assign("FP_COMMON_DIR", file.path(BASE, "FP",    "unr"), envir = .GlobalEnv)
    assign("POWER_DIR",     file.path(BASE, "Power", "unr"), envir = .GlobalEnv)
    assign("N_TOTAL",    unr_n,  envir = .GlobalEnv)
    assign("N_VARIANTS", unr_p,  envir = .GlobalEnv)
    assign("MAF_LOW",    0.05, envir = .GlobalEnv)
    assign("MAF_HIGH",   0.50, envir = .GlobalEnv)
    assign("ANC_P_LOW",  0.10, envir = .GlobalEnv)
    assign("ANC_P_HIGH", 0.90, envir = .GlobalEnv)
    dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)
    return(invisible(NULL))
  }

  if (mode == "bench") {
    ## Parametric large-N benchmark mode. Driven by env vars:
    ##   BENCH_DIR  -- absolute path to the per-N data dir (e.g.
    ##                 .../data/bench_N50000_P1000); created if missing.
    ##   BENCH_N    -- sample size (positive int)
    ##   BENCH_P    -- variant count (default 1000)
    ## Used by simu_bench_largeN.R and reused by make_vcf.R /
    ## make_tractormix_hapdose.R / make_plink_pruned.R unchanged.
    bench_dir <- Sys.getenv("BENCH_DIR", "")
    if (!nzchar(bench_dir))
      stop("set_mode('bench'): BENCH_DIR env var must be set")
    bench_n <- as.integer(Sys.getenv("BENCH_N", "0"))
    if (is.na(bench_n) || bench_n <= 0L)
      stop("set_mode('bench'): BENCH_N env var (positive int) must be set")
    bench_p <- as.integer(Sys.getenv("BENCH_P", "1000"))

    assign("DATA_DIR",     bench_dir,                    envir = .GlobalEnv)
    assign("IND_DIR",      file.path(bench_dir, "Ind"),  envir = .GlobalEnv)
    assign("REL_DIR",      file.path(bench_dir, "Rel"),  envir = .GlobalEnv)
    assign("BLOCKROW_DIR", bench_dir,                    envir = .GlobalEnv)
    assign("FP_COMMON_DIR", file.path(BASE, "FP",    "bench"), envir = .GlobalEnv)
    assign("POWER_DIR",     file.path(BASE, "Power", "bench"), envir = .GlobalEnv)
    assign("N_TOTAL",    bench_n,  envir = .GlobalEnv)
    assign("N_VARIANTS", bench_p,  envir = .GlobalEnv)
    assign("MAF_LOW",    0.05, envir = .GlobalEnv)
    assign("MAF_HIGH",   0.50, envir = .GlobalEnv)
    assign("ANC_P_LOW",  0.10, envir = .GlobalEnv)
    assign("ANC_P_HIGH", 0.90, envir = .GlobalEnv)
    dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)
    return(invisible(NULL))
  }

  assign("DATA_DIR",    file.path(BASE, "data", mode),        envir = .GlobalEnv)
  assign("IND_DIR",     file.path(BASE, "data", mode, "Ind"), envir = .GlobalEnv)
  assign("REL_DIR",     file.path(BASE, "data", mode, "Rel"), envir = .GlobalEnv)
  assign("BLOCKROW_DIR", file.path(BASE, "data", mode),       envir = .GlobalEnv)
  assign("FP_COMMON_DIR",  file.path(BASE, "FP",    mode),    envir = .GlobalEnv)
  assign("POWER_DIR",      file.path(BASE, "Power", mode),    envir = .GlobalEnv)
  if (mode == "common") {
    assign("MAF_LOW",  0.05, envir = .GlobalEnv)
    assign("MAF_HIGH", 0.50, envir = .GlobalEnv)
    assign("ANC_P_LOW",  0.10, envir = .GlobalEnv)
    assign("ANC_P_HIGH", 0.90, envir = .GlobalEnv)
    assign("N_VARIANTS", 1000000L, envir = .GlobalEnv)
  } else {
    assign("MAF_LOW",  0.01, envir = .GlobalEnv)
    assign("MAF_HIGH", 0.05, envir = .GlobalEnv)
    assign("ANC_P_LOW",  0.01, envir = .GlobalEnv)
    assign("ANC_P_HIGH", 0.10, envir = .GlobalEnv)
    ## over-simulate since the tight 1-5% filter drops many variants
    assign("N_VARIANTS", 2000000L, envir = .GlobalEnv)
  }
  dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(IND_DIR,  recursive = TRUE, showWarnings = FALSE)
  dir.create(REL_DIR,  recursive = TRUE, showWarnings = FALSE)
  invisible(NULL)
}

## Convenience loader for the core helpers.
source_utils <- function() {
  source(file.path(SCRIPT_DIR, "Utils3way.R"))
  source(file.path(SCRIPT_DIR, "Pedigree3way.R"))
}
