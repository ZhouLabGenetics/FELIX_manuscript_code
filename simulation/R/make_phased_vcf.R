## Synthesize a phased VCF from the existing Block_row<No>/ tables.
##
## The simulator stores diploid-summed counts per ancestry (GAfr/GEur/GNat
## in {0,1,2}) and per-ancestry haplotype counts (LAafr/LAeur/LAnat in
## {0,1,2}, summing to 2). True hap1/hap2 ordering from Pedigree3way_Ind.R
## is NOT preserved on disk, so we choose a CANONICAL ordering that is
## deterministic and reproducible:
##
##   Sort the two haplotypes by ancestry index ascending (AFR=0 < EUR=1 < NAT=2);
##   break ties by allele (0 before 1).
##
## Encoded as a VCF v4.2 with two phased FORMAT fields:
##
##   GT  : phased genotype  e.g. "0|1"
##   HA  : phased haplotype ancestry  e.g. "0|2"  (0=AFR, 1=EUR, 2=NAT)
##
## Per-sample field is "<g1>|<g2>:<a1>|<a2>" (= 7 chars) vs the unphased
## dosage VCF's "<DS1>:<DS2>:<DS3>:<ANC1>:<ANC2>:<ANC3>:<DSALL>" (~13 chars),
## a ~50% size reduction. The DS/ANC fields can be reconstructed losslessly
## from GT+HA downstream:
##   DS_k  = sum over haps h of  (GT_h == 1) * (HA_h == k)
##   ANC_k = sum over haps h of  (HA_h == k)
##   DSALL = sum over haps h of  GT_h
##
## Output : data/<mode>/vcf/chr<No>.phased.vcf.gz and .csi
##
## Usage:  Rscript make_phased_vcf.R <No> <common|lowfreq>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: make_phased_vcf.R <No> <common|lowfreq>")
No   <- as.integer(args[1])
mode <- args[2]

source(file.path(Sys.getenv("FELIX_SIM_BASE", unset = stop("Set FELIX_SIM_BASE to the simulation directory before running this script")), "R", "config.R"))
set_mode(mode)
library(data.table)

brdir   <- file.path(DATA_DIR, paste0("Block_row", No))
vcf_dir <- file.path(DATA_DIR, "vcf")
dir.create(vcf_dir, recursive = TRUE, showWarnings = FALSE)

read_br <- function(base) fread(file.path(brdir, paste0(base, ".tsv")), sep = "\t")
GAfr  <- read_br("GAfr")
GEur  <- read_br("GEur")
GNat  <- read_br("GNat")
LAafr <- read_br("LAafr")
LAeur <- read_br("LAeur")
LAnat <- read_br("LAnat")

meta_cols <- 1:5
meta <- GAfr[, meta_cols, with = FALSE]
sample_ids <- setdiff(colnames(GAfr), colnames(meta))
stopifnot(length(sample_ids) == N_TOTAL)

to_mat <- function(dt) {
  m <- as.matrix(dt[, -meta_cols, with = FALSE])
  storage.mode(m) <- "integer"
  m
}
g_a  <- to_mat(GAfr);  g_e  <- to_mat(GEur);  g_n  <- to_mat(GNat)
la_a <- to_mat(LAafr); la_e <- to_mat(LAeur); la_n <- to_mat(LAnat)
rm(GAfr, GEur, GNat, LAafr, LAeur, LAnat); invisible(gc())

## Free up; we no longer need the data.tables.

## ---- Decompose into (hap1, hap2) under canonical ordering --------------
## LA vector (LAa, LAe, LAn) sums to 2. Six possible ancestry pairs:
##   sa_a : (2,0,0)  AFR/AFR
##   sa_e : (0,2,0)  EUR/EUR
##   sa_n : (0,0,2)  NAT/NAT
##   mix_ae : (1,1,0)  AFR | EUR  (hap1=AFR, hap2=EUR)
##   mix_an : (1,0,1)  AFR | NAT
##   mix_en : (0,1,1)  EUR | NAT
sa_a   <- la_a == 2L
sa_e   <- la_e == 2L
sa_n   <- la_n == 2L
mix_ae <- la_a == 1L & la_e == 1L
mix_an <- la_a == 1L & la_n == 1L
mix_en <- la_e == 1L & la_n == 1L
## Coverage sanity check: every (variant, sample) must fall in exactly one bin.
covered <- sa_a + sa_e + sa_n + mix_ae + mix_an + mix_en
if (any(covered != 1L))
  stop("Bad LA decomposition: ", sum(covered != 1L),
       " cells failed sum-to-one check")
rm(covered)

dims <- dim(g_a)
h1_anc <- matrix(0L, dims[1], dims[2])
h2_anc <- matrix(0L, dims[1], dims[2])
h1_anc[sa_a]   <- 0L; h2_anc[sa_a]   <- 0L
h1_anc[sa_e]   <- 1L; h2_anc[sa_e]   <- 1L
h1_anc[sa_n]   <- 2L; h2_anc[sa_n]   <- 2L
h1_anc[mix_ae] <- 0L; h2_anc[mix_ae] <- 1L
h1_anc[mix_an] <- 0L; h2_anc[mix_an] <- 2L
h1_anc[mix_en] <- 1L; h2_anc[mix_en] <- 2L

h1_al <- matrix(0L, dims[1], dims[2])
h2_al <- matrix(0L, dims[1], dims[2])
## Same-ancestry cells: g in {0,1,2}; canonical (0-allele first) -> (0,0),(0,1),(1,1)
h1_al[sa_a]  <- as.integer(g_a[sa_a] == 2L)
h2_al[sa_a]  <- as.integer(g_a[sa_a] >= 1L)
h1_al[sa_e]  <- as.integer(g_e[sa_e] == 2L)
h2_al[sa_e]  <- as.integer(g_e[sa_e] >= 1L)
h1_al[sa_n]  <- as.integer(g_n[sa_n] == 2L)
h2_al[sa_n]  <- as.integer(g_n[sa_n] >= 1L)
## Mixed-ancestry cells: hap1 carries the lower-index ancestry's allele;
## g_k is in {0,1} there (only 1 hap per ancestry).
h1_al[mix_ae] <- g_a[mix_ae]; h2_al[mix_ae] <- g_e[mix_ae]
h1_al[mix_an] <- g_a[mix_an]; h2_al[mix_an] <- g_n[mix_an]
h1_al[mix_en] <- g_e[mix_en]; h2_al[mix_en] <- g_n[mix_en]

## ---- Round-trip sanity check (a few random cells) ----------------------
chk_n <- 1000L
set.seed(42L)
ii <- sample.int(dims[1], chk_n, replace = TRUE)
jj <- sample.int(dims[2], chk_n, replace = TRUE)
recover_g <- function(anc_target) {
  (h1_al[cbind(ii,jj)] * (h1_anc[cbind(ii,jj)] == anc_target)) +
  (h2_al[cbind(ii,jj)] * (h2_anc[cbind(ii,jj)] == anc_target))
}
recover_la <- function(anc_target) {
  (h1_anc[cbind(ii,jj)] == anc_target) + (h2_anc[cbind(ii,jj)] == anc_target)
}
stopifnot(all(recover_g(0L)  == g_a[cbind(ii,jj)]))
stopifnot(all(recover_g(1L)  == g_e[cbind(ii,jj)]))
stopifnot(all(recover_g(2L)  == g_n[cbind(ii,jj)]))
stopifnot(all(recover_la(0L) == la_a[cbind(ii,jj)]))
stopifnot(all(recover_la(1L) == la_e[cbind(ii,jj)]))
stopifnot(all(recover_la(2L) == la_n[cbind(ii,jj)]))
cat("Round-trip OK on ", chk_n, " random cells.\n", sep = "")
rm(g_a, g_e, g_n, la_a, la_e, la_n); invisible(gc())

## ---- Build per-sample strings ------------------------------------------
field <- paste0(h1_al, "|", h2_al, ":", h1_anc, "|", h2_anc)
dim(field) <- dims
rm(h1_al, h2_al, h1_anc, h2_anc); invisible(gc())
body_rows <- apply(field, 1, paste, collapse = "\t")
rm(field); invisible(gc())

## ---- VCF body + header -------------------------------------------------
chrom <- as.character(meta$CHR)
pos   <- meta$POS
id    <- meta$SNP
ref   <- meta$REF
alt   <- meta$ALT
fmt_line <- rep("GT:HA", nrow(meta))
body <- paste(chrom, pos, id, ref, alt, ".", "PASS", ".", fmt_line, body_rows,
              sep = "\t")
rm(body_rows); invisible(gc())

header_lines <- c(
  "##fileformat=VCFv4.2",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Phased genotype: hap1|hap2 alt allele indicator (0/1)\">",
  "##FORMAT=<ID=HA,Number=1,Type=String,Description=\"Phased haplotype ancestry: hap1|hap2 (0=AFR, 1=EUR, 2=NAT)\">",
  paste0("##contig=<ID=", No, ">"),
  paste(c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT", sample_ids),
        collapse = "\t")
)

vcf_path <- file.path(vcf_dir, paste0("chr", No, ".phased.vcf"))
writeLines(c(header_lines, body), vcf_path)
rm(body); invisible(gc())

vcf_gz <- paste0(vcf_path, ".gz")
if (file.exists(vcf_gz)) file.remove(vcf_gz)
rc <- system2("bgzip", args = c("-f", vcf_path))
if (rc != 0) stop("bgzip failed on ", vcf_path)
rc <- system2("tabix", args = c("-p", "vcf", "--csi", vcf_gz))
if (rc != 0) stop("tabix failed on ", vcf_gz)

cat("Wrote", vcf_gz, "(", nrow(meta), "variants,", length(sample_ids), "samples )\n")
