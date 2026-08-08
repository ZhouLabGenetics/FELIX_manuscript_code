## Convert one Block_row (post-QC) into a bgzipped VCF that FELIX
## consumes in step2. Per-sample FORMAT layout matches the 2-way example
## extended to 3-way:
##
##   FORMAT = DS1:DS2:DS3:ANC1:ANC2:ANC3:DSALL
##
##   DS1  = ancestry-specific alt dosage for AFR (= GAfr[i,j])
##   DS2  = ancestry-specific alt dosage for EUR (= GEur[i,j])
##   DS3  = ancestry-specific alt dosage for NAT (= GNat[i,j])
##   ANC1 = local-ancestry allele count for AFR  (= LAafr[i,j])
##   ANC2 = local-ancestry allele count for EUR  (= LAeur[i,j])
##   ANC3 = local-ancestry allele count for NAT  (= LAnat[i,j])
##   DSALL = DS1 + DS2 + DS3                          (= GTot[i,j])
##
## No GT field is emitted (FELIX reads dosage fields directly).
##
## Output : data/<mode>/vcf/chr<No>.vcf.gz and chr<No>.vcf.gz.csi
##          (mode = common|lowfreq plain dirs; for unr/admix/bench the mode
##           sets DATA_DIR via env vars and writes vcf/ under that DATA_DIR.)
##
## Usage:  Rscript make_vcf.R <No> <common|lowfreq|unr|admix|bench>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: make_vcf.R <No> <common|lowfreq|unr|admix|bench>")
No   <- as.integer(args[1])
mode <- args[2]

source(file.path(Sys.getenv("FELIX_SIM_BASE", unset = stop("Set FELIX_SIM_BASE to the simulation directory before running this script")), "R", "config.R"))
set_mode(mode)
library(data.table)

brdir <- file.path(DATA_DIR, paste0("Block_row", No))
vcf_dir <- file.path(DATA_DIR, "vcf")
dir.create(vcf_dir, recursive = TRUE, showWarnings = FALSE)

meta_cols <- 1:5

## Read each Block_row TSV and immediately materialise the matrix, freeing
## the data.table. Without these rm()+gc() calls peak memory ~doubles for
## large N (each fread holds the data, then as.matrix duplicates it).
read_to_mat <- function(base) {
  dt  <- fread(file.path(brdir, paste0(base, ".tsv")), sep = "\t")
  mat <- as.matrix(dt[, -meta_cols, with = FALSE])
  rm(dt); invisible(gc())
  mat
}

## Read GTot first to capture metadata + sample IDs.
GTot_dt <- fread(file.path(brdir, "GTot.tsv"), sep = "\t")
meta <- GTot_dt[, meta_cols, with = FALSE]
sample_ids <- setdiff(colnames(GTot_dt), colnames(meta))
stopifnot(length(sample_ids) == N_TOTAL)
m_all <- as.matrix(GTot_dt[, -meta_cols, with = FALSE])
rm(GTot_dt); invisible(gc())

m_ds1 <- read_to_mat("GAfr")
m_ds2 <- read_to_mat("GEur")
m_ds3 <- read_to_mat("GNat")
m_a1  <- read_to_mat("LAafr")
m_a2  <- read_to_mat("LAeur")
m_a3  <- read_to_mat("LAnat")

P_total <- nrow(meta)
chrom <- as.character(meta$CHR)
pos   <- meta$POS
id    <- meta$SNP
ref   <- meta$REF
alt   <- meta$ALT
fmt_str <- "DS1:DS2:DS3:ANC1:ANC2:ANC3:DSALL"

## Header
header_lines <- c(
  "##fileformat=VCFv4.2",
  "##FORMAT=<ID=DS1,Number=1,Type=Float,Description=\"Dosage of AFR ancestry\">",
  "##FORMAT=<ID=DS2,Number=1,Type=Float,Description=\"Dosage of EUR ancestry\">",
  "##FORMAT=<ID=DS3,Number=1,Type=Float,Description=\"Dosage of NAT ancestry\">",
  "##FORMAT=<ID=ANC1,Number=1,Type=Float,Description=\"Allele count of AFR ancestry\">",
  "##FORMAT=<ID=ANC2,Number=1,Type=Float,Description=\"Allele count of EUR ancestry\">",
  "##FORMAT=<ID=ANC3,Number=1,Type=Float,Description=\"Allele count of NAT ancestry\">",
  "##FORMAT=<ID=DSALL,Number=1,Type=Float,Description=\"Sum of dosages across ancestries\">",
  paste0("##contig=<ID=", No, ">"),
  paste(c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT", sample_ids),
        collapse = "\t")
)

## Stream the body in variant chunks so we never materialise the full
## P*N character matrix (~8 GB at N=20k, P=50k; ~20 GB at N=50k, P=50k).
vcf_path <- file.path(vcf_dir, paste0("chr", No, ".vcf"))
con <- file(vcf_path, open = "w")
writeLines(header_lines, con)

CHUNK <- 500L  # variants per chunk
for (lo in seq.int(1L, P_total, by = CHUNK)) {
  hi <- min(lo + CHUNK - 1L, P_total)
  idx <- lo:hi
  field <- paste(m_ds1[idx, ], m_ds2[idx, ], m_ds3[idx, ],
                 m_a1[idx, ],  m_a2[idx, ],  m_a3[idx, ], m_all[idx, ],
                 sep = ":")
  dim(field) <- c(length(idx), length(sample_ids))
  body_rows <- apply(field, 1, paste, collapse = "\t")
  body_chunk <- paste(chrom[idx], pos[idx], id[idx], ref[idx], alt[idx],
                      ".", "PASS", ".", fmt_str, body_rows, sep = "\t")
  writeLines(body_chunk, con)
}
close(con)

## bgzip + tabix / csi index. bgzip must be available on PATH inside the
## singularity container (bioconductor_docker ships it via htslib-tools).
vcf_gz <- paste0(vcf_path, ".gz")
if (file.exists(vcf_gz)) file.remove(vcf_gz)
rc <- system2("bgzip", args = c("-f", vcf_path))
if (rc != 0) stop("bgzip failed on ", vcf_path)
rc <- system2("tabix", args = c("-p", "vcf", "--csi", vcf_gz))
if (rc != 0) stop("tabix failed on ", vcf_gz)

cat("Wrote", vcf_gz, "(", nrow(meta), "variants,", length(sample_ids), "samples )\n")
