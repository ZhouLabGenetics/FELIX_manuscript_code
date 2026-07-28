#!/usr/bin/env Rscript
# Signed Miami for the heterogeneity loci (ports scripts/20_signed_miami.py):
# per-ancestry sign(beta) * -log10(p) across the locus. AFR up, EUR down. Ancestry palette.
# No on-plot text; large fonts; PDF+PNG.
source("scripts_R/00_theme.R"); source("scripts_R/01_data.R")
suppressMessages(library(patchwork))

signed_panel <- function(sp, chrom, pos, trait, win = 2.5e5) {
  sg <- saige_region(sp, chrom, pos, window = win)
  if (nrow(sg) == 0) return(NULL)
  d <- rbindlist(lapply(seq_len(5), function(i) {
    b <- as.numeric(sg[[paste0("BETA_c_anc", i)]]); p <- as.numeric(sg[[paste0("p.value_c_anc", i)]])
    data.table(POS = sg$POS, anc = ANC_MAP$name[i], y = sign(b) * mlog(p))
  }))
  d <- d[is.finite(y)]; d[, anc := factor(anc, levels = ANC_MAP$name)]
  ggplot(d, aes(POS/1e6, y, color = anc)) +
    geom_hline(yintercept = 0, color = "grey30", linewidth = 0.6) +
    geom_hline(yintercept = c(-GW_LOG, GW_LOG), linetype = "dashed", color = "grey60", linewidth = 0.6) +
    geom_point(size = 2.2, alpha = 0.8) + scale_color_anc(name = "Local ancestry") +
    labs(x = paste0("Chromosome ", chrom, " position (Mb)"),
         y = expression("signed  " * sign(beta) %.% -log[10](p)), title = trait)
}
p1 <- signed_panel("pheno_3035995", 19, 44919689, "APOC1 - alkaline phosphatase")
p2 <- signed_panel("pheno_3028288", 16, 72080103, "HPR/TXNL4B - LDL cholesterol")
save_fig(p1 / p2, "SFig_mechB_signed_miami", width = 11, height = 11)
