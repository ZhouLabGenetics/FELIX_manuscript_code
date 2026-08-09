#!/usr/bin/env Rscript
# LD locus-zooms. FELIXassoc CCT + All-by-All meta
source("scripts_R/00_theme.R"); source("scripts_R/01_data.R")
suppressMessages(library(viridisLite))

loci <- list(
  list(g="il23r",  sp="pheno_GI_522.11", ap="GI_522.11", chr=1,  pos=67240275,  trait="IL23R - Crohn's disease"),
  list(g="ikzf3",  sp="pheno_RE_475",    ap="RE_475",    chr=17, pos=39765489,  trait="IKZF3 - asthma"),
  list(g="adrb2",  sp="BMI",             ap="BMI",       chr=5,  pos=148898672, trait="ADRB2/SH3TC2 - BMI"),
  list(g="supt3h", sp="height",          ap="height",    chr=6,  pos=44831920,  trait="SUPT3H - height"))

for (L in loci) {
  ld <- fread(file.path(ROOT, "ref_1000g", "regions", paste0(L$g, ".ld")))
  r2map <- setNames(ld$R2, ld$BP_B)
  sg <- saige_region(L$sp, L$chr, L$pos)
  if (nrow(sg) == 0) { message("skip ", L$g); next }
  sg[, y := mlog(P_cct_admixed_c)]; sg <- sg[is.finite(y)]
  sg[, r2 := r2map[as.character(POS)]]
  lead <- sg[which.max(y)]; sg[POS == lead$POS, r2 := 1]
  ab <- aba_region(L$ap, "META", L$chr, L$pos)
  if (nrow(ab)) ab[, y := mlog(Pvalue)]
  p <- ggplot() +
    { if (nrow(ab)) geom_point(data = ab, aes(POS/1e6, y), shape = 1, color = "#0077BB",
                               size = 2, stroke = 0.7, alpha = 0.8) } +
    geom_point(data = sg[!is.na(r2)], aes(POS/1e6, y, fill = r2), shape = 21, size = 3, color = "white", stroke = 0.2) +
    geom_point(data = lead, aes(POS/1e6, y), shape = 18, size = 8, color = "#EE3377") +
    geom_hline(yintercept = GW_LOG, linetype = "dashed", color = "grey40", linewidth = 0.8) +
    scale_fill_viridis_c(option = "D", limits = c(0, 1), name = expression(r^2)) +
    labs(x = paste0("Chromosome ", L$chr, " position (Mb)"),
         y = expression("Association  " * -log[10](p)), title = L$trait)
  save_fig(p, paste0("fig4_lz_", L$g, "_LD"), width = 9.5, height = 5.6)
}
