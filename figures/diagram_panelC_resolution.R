#!/usr/bin/env Rscript
# Panel C "what FELIXassoc delivers" (rebuild of scripts/12_panel_c_resolution.py):
# a) tighter effect estimates (SE scatter at All-by-All GW loci); b) per-ancestry effect
# fingerprint for representative loci. Data-driven; no on-plot text; new palette.
source("scripts_R/00_theme.R"); source("scripts_R/01_data.R")
suppressMessages(library(patchwork))
dt <- load_scatter()

## a) SE scatter (All by All meta SE vs FELIXassoc combined-ancestry SE)
d <- dt[tophit_source == "ABA"]
d[, `:=`(x = as.numeric(ABA_META_SE), y = as.numeric(SAIGE_SE_c_ancALL))]
d <- d[is.finite(x) & is.finite(y) & x > 0 & y > 0]
d[, tighter := ifelse(y < x, "FELIXassoc tighter", "wider")]
lim <- range(c(d$x, d$y))
pA <- ggplot(d, aes(x, y, color = tighter)) +
  geom_abline(slope = 1, linetype = "dashed", color = "grey45") +
  geom_point(size = 2.4, alpha = 0.7) +
  scale_x_log10() + scale_y_log10() + coord_equal() +
  scale_color_manual(values = c("FELIXassoc tighter" = "#AA3377", "wider" = "#CCBB44"), name = NULL) +
  labs(x = "All by All SE (meta-analysis)", y = "FELIXassoc SE (combined-ancestry)",
       title = "a  Tighter effect estimates") + theme(legend.position = "top")

## b) per-ancestry effect fingerprint for representative loci
loci <- data.table(
  gene = c("IL23R","ADRB2","PNPLA3","APOC1","HPR","ZPR1"),
  ph   = c("GI_522.11","BMI","3013721","3035995","3028288","3022192"),
  pos  = c(67240275, 148898672, 43928850, 44919689, 72080103, 116778201))
fp <- rbindlist(lapply(seq_len(nrow(loci)), function(i) {
  r <- get_locus(dt, loci$gene[i], loci$ph[i], loci$pos[i]); if (is.null(r)) return(NULL)
  data.table(locus = loci$gene[i], anc = ANC_MAP$name,
             beta = sapply(paste0("SAIGE_BETA_c_", ANC_MAP$suf), function(cc) as.numeric(r[[cc]])))
}))
fp <- fp[is.finite(beta)]; fp[, `:=`(anc = factor(anc, levels = ANC_MAP$name),
                                     locus = factor(locus, levels = loci$gene))]
pB <- ggplot(fp, aes(anc, beta, fill = anc)) +
  geom_hline(yintercept = 0, color = "grey40") +
  geom_col(width = 0.8) + facet_wrap(~ locus, nrow = 1, scales = "free_y") +
  scale_fill_anc(guide = "none") +
  labs(x = NULL, y = expression("per-ancestry effect  " * beta),
       title = "b  Per-ancestry effect fingerprint") +
  theme(strip.background = element_blank(), strip.text = element_text(size = 17, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 13))
save_fig(pA / pB, "panel_c_resolution", width = 14, height = 11)
