#!/usr/bin/env Rscript
# evaluate_r2_4way.R
# Compute R² for 4 PRS combinations (AFR/EUR x local/global) vs a phenotype.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
    stop(
        "Usage: Rscript evaluate_r2_4way.R ",
        "<AFR_local> <AFR_global> <EUR_local> <EUR_global> ",
        "<ukb_pheno.tsv> <pheno_col> <output_dir>"
    )
}

afr_local_f  <- args[1]
afr_global_f <- args[2]
eur_local_f  <- args[3]
eur_global_f <- args[4]
ukb_pheno_f  <- args[5]
pheno_col    <- args[6]
output_dir   <- args[7]

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_sscore <- function(f) {
    d <- read.table(f, header = TRUE, comment.char = "", check.names = FALSE,
                    sep = "\t")
    score_col <- grep("SCORE.*SUM$", colnames(d), value = TRUE)[1]
    if (is.na(score_col)) {
        score_col <- tail(colnames(d), 1)
        message(sprintf("  [%s] no SCORE_SUM col; using last col '%s'", f, score_col))
    }
    data.frame(IID = d$IID, score = d[[score_col]], check.names = FALSE)
}

cat("Reading sscore files...\n")
afr_local  <- read_sscore(afr_local_f)
afr_global <- read_sscore(afr_global_f)
eur_local  <- read_sscore(eur_local_f)
eur_global <- read_sscore(eur_global_f)

cat("Reading phenotype:", ukb_pheno_f, "col:", pheno_col, "\n")
pheno <- read.table(ukb_pheno_f, header = TRUE, sep = "\t", check.names = FALSE)
if (!pheno_col %in% colnames(pheno)) {
    stop(sprintf("Phenotype column '%s' not found in %s. Available: %s",
                 pheno_col, ukb_pheno_f, paste(colnames(pheno), collapse = ", ")))
}
id_col <- intersect(c("eid", "IID", "id", "ID"), colnames(pheno))[1]
if (is.na(id_col)) stop("No eid/IID/id column in phenotype file")
inv_norm <- function(x) {
    qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
}
pheno <- data.frame(IID = pheno[[id_col]], y = pheno[[pheno_col]], check.names = FALSE)
pheno$y <- inv_norm(pheno$y)

# =====================================================================
# TRANSLATE PHENOTYPE IDs
cat("Translating UKB EIDs to SSCORE IDs...\n")
map_file <- "/data/wzhougroup/lhu/saige_tractor/prs_pipeline/sample_id_mapping.txt"
id_map <- read.table(map_file, header = FALSE, stringsAsFactors = FALSE)
colnames(id_map) <- c("ukb_eid", "sscore_id")
id_map$ukb_eid <- as.character(id_map$ukb_eid)
id_map$sscore_id <- as.character(id_map$sscore_id)

pheno$IID <- as.character(pheno$IID)
pheno <- merge(pheno, id_map, by.x = "IID", by.y = "ukb_eid")
pheno$IID <- pheno$sscore_id
pheno$sscore_id <- NULL 

# =====================================================================
# MAP AND FILTER THE AFR KEEP LIST
cat("Loading and mapping AFR keep list...\n")
afr_keep_file <- "/data/wzhougroup/lhu/saige_tractor/prs_pipeline/shared/afr_keep.txt"
afr_keep_raw <- read.table(afr_keep_file, header = FALSE)
afr_keep_ukb <- as.character(afr_keep_raw$V2)

# Translate the AFR keep list UKB EIDs into SSCORE IDs using the same map_file
afr_keep_mapped <- id_map$sscore_id[id_map$ukb_eid %in% afr_keep_ukb]
# =====================================================================

compute_r2 <- function(prs, ph) {
    m <- merge(prs, ph, by = "IID")
    m <- m[!is.na(m$score) & !is.na(m$y), ]
    if (nrow(m) < 10) return(c(N = nrow(m), R2 = NA_real_))
    fit <- lm(y ~ score, data = m)
    c(N = nrow(m), R2 = summary(fit)$r.squared)
}

cat("Computing R²...\n")

# Apply the mapped AFR keep list to the AFR dataframes
afr_local_filtered  <- afr_local[afr_local$IID %in% afr_keep_mapped, ]
afr_global_filtered <- afr_global[afr_global$IID %in% afr_keep_mapped, ]

r <- rbind(
    `AFR-local`  = compute_r2(afr_local_filtered,  pheno),
    `AFR-global` = compute_r2(afr_global_filtered, pheno),
    `EUR-local`  = compute_r2(eur_local,  pheno),
    `EUR-global` = compute_r2(eur_global, pheno)
)

out_f <- file.path(output_dir, "r2_4way.tsv")
write.table(r, out_f, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

cat("\nResults:\n")
print(r)
cat("\nSaved to:", out_f, "\n")
