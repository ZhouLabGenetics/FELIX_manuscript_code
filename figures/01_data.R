# 01_data.R — shared data loaders (mirror scripts/common.py). source after 00_theme.R
PHENO_ALLOWLIST <- c("282.5","3006923","3007070","3009744","3013721","3022192",
  "3024929","3027114","3028288","3035995","BMI","CV_401","DE_668.1","EM_202.2",
  "EM_236.1","GI_522.1","GI_522.11","GI_522.12","MB_287.1","MB_290.1","MS_705",
  "NS_326.1","RE_475","height")
PHENO_LABELS <- c(BMI="BMI", height="Height", CV_401="Hypertension",
  EM_236.1="Obesity", RE_475="Asthma", NS_326.1="Multiple Sclerosis",
  "282.5"="Sickle Cell Anemia", "3007070"="HDL Cholesterol", "3009744"="MCHC",
  "3022192"="Triglycerides", "3024929"="Platelet Count", "3027114"="Total Cholesterol",
  "3028288"="LDL Cholesterol", DE_668.1="Atopic Dermatitis", EM_202.2="Type 2 Diabetes",
  MS_705="Rheumatoid arthritis", "3006923"="Alanine Aminotransferase",
  "3013721"="Aspartate Aminotransferase", "3035995"="Alkaline Phosphatase",
  GI_522.1="Inflammatory Bowel Disease", GI_522.11="Crohn's Disease",
  GI_522.12="Ulcerative Colitis", MB_287.1="Schizophrenia", MB_290.1="PTSD")

canon <- function(p) sub("^pheno_", "", as.character(p))
lbl   <- function(p){ p <- canon(p); ifelse(p %in% names(PHENO_LABELS), PHENO_LABELS[p], p) }

load_scatter <- function(name = "table_P_cct_admixed_c_vs_META.tsv") {
  dt <- fread(file.path(SCATTER, name))
  dt[canon(phenotype) %in% PHENO_ALLOWLIST]
}
get_locus <- function(dt, gene, ph, pos) {
  sub <- dt[canon(phenotype) == ph & grepl(gene, Gene, fixed = TRUE)]
  if (nrow(sub) == 0) return(NULL)
  sub[which.min(abs(as.numeric(Pos) - pos))]
}
p2chisq <- function(p) qchisq(pmax(as.numeric(p), 1e-300), 1, lower.tail = FALSE)
mlog    <- function(p, cap = 24) pmin(-log10(pmax(as.numeric(p), 1e-300)), cap)

SAIGE_SS <- file.path(ROOT, "saige_sumstat")
ABA_SS   <- file.path(ROOT, "aba_sumstat")
saige_region <- function(sp, chrom, pos, window = 5e5) {
  f <- file.path(SAIGE_SS, paste0("merged_saigetractor_", sp, ".txt.gz"))
  d <- fread(cmd = sprintf("gzcat %s", shQuote(f)))
  d[, CHR := sub("chr", "", as.character(CHR))][, POS := as.numeric(POS)]
  d[CHR == as.character(chrom) & POS >= pos - window & POS <= pos + window]
}
aba_region <- function(ph, anc, chrom, pos, window = 5e5) {
  f <- file.path(ABA_SS, sprintf("%s_%s.tsv", ph, anc))
  if (!file.exists(f)) return(data.table())
  d <- fread(f)
  d[, CHR := sub("chr", "", as.character(CHR))][, POS := as.numeric(POS)]
  d[CHR == as.character(chrom) & POS >= pos - window & POS <= pos + window]
}
