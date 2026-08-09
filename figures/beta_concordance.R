#!/usr/bin/env Rscript
# Beta concordance (All by All vs FELIXassoc) at All-by-All top hits
source("scripts_R/00_theme.R"); source("scripts_R/01_data.R")
dt <- load_scatter()
aba <- dt[tophit_source == "ABA"]

# Fig 12 — overall meta vs mega. ALL loci included; betas shown on the absolute-value
# scale (negative betas folded to |beta|) so effect magnitudes compare on one axis.
d12 <- data.table(x = abs(as.numeric(aba$ABA_META_BETA)),
                  y = abs(as.numeric(aba$SAIGE_BETA_c_ancALL)))
d12 <- d12[is.finite(x) & is.finite(y)]
lim <- max(d12$x, d12$y) * 1.05
r12 <- cor(d12$x, d12$y)
p12 <- ggplot(d12, aes(x, y)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.9, color = "grey40") +
  geom_point(size = 3, color = METHOD_COLORS[["FELIXassoc"]], alpha = 0.75,
             stroke = 0.3, shape = 21, fill = METHOD_COLORS[["FELIXassoc"]]) +
  annotate("text", x = 0.03 * lim, y = 0.95 * lim, hjust = 0, size = 6, fontface = "bold",
           label = sprintf("Pearson r = %.2f\nn = %d loci", r12, nrow(d12))) +
  coord_equal(xlim = c(0, lim), ylim = c(0, lim)) +
  labs(x = expression("All by All  " * group("|", beta, "|") * "  (meta-analysis)"),
       y = expression("FELIXassoc  " * group("|", beta, "|") * "  (mega-analysis)"),
       title = "Effect-size concordance")
save_fig(p12, "fig12_beta_concordance_overall", width = 7.5, height = 7.5)

# Fig 12b — signed, ALL 715 loci (same set as the chi-square scatter; no positive-beta filter)
d12s <- data.table(x = as.numeric(aba$ABA_META_BETA), y = as.numeric(aba$SAIGE_BETA_c_ancALL))
d12s <- d12s[is.finite(x) & is.finite(y)]
lim2 <- max(abs(c(d12s$x, d12s$y))) * 1.05
r12s <- cor(d12s$x, d12s$y)
p12s <- ggplot(d12s, aes(x, y)) +
  geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.9, colour = "grey40") +
  geom_point(size = 2.6, alpha = 0.6, shape = 21, colour = METHOD_COLORS[["FELIXassoc"]],
             fill = METHOD_COLORS[["FELIXassoc"]], stroke = 0.3) +
  annotate("text", x = -lim2 * 0.96, y = lim2 * 0.96, hjust = 0, vjust = 1, size = 6, fontface = "bold",
           label = sprintf("Pearson r = %.2f\nn = %d loci", r12s, nrow(d12s))) +
  coord_equal(xlim = c(-lim2, lim2), ylim = c(-lim2, lim2)) +
  labs(x = "All by All beta (meta-analysis)", y = "FELIXassoc beta (mega-analysis)",
       title = "Effect-size concordance (signed, all loci)")
save_fig(p12s, "fig12b_beta_concordance_signed", width = 7.5, height = 7.5)

# Fig 13 — per-ancestry facets
ANC_MAP2 <- data.table(name = c("AFR","EAS","EUR","AMR","SAS"),
                       aba  = c("AFR","EAS","EUR","AMR","SAS"),
                       suf  = paste0("anc", 1:5))
# ALL loci per ancestry, folded to |beta| so every locus (positive and negative) is shown.
d13 <- rbindlist(lapply(seq_len(5), function(i) {
  data.table(anc = ANC_MAP2$name[i],
             x = abs(as.numeric(aba[[paste0("ABA_", ANC_MAP2$aba[i], "_BETA")]])),
             y = abs(as.numeric(aba[[paste0("SAIGE_BETA_c_", ANC_MAP2$suf[i])]])))
}))
d13 <- d13[is.finite(x) & is.finite(y)]
d13[, anc := factor(anc, levels = ANC_MAP2$name)]
r13 <- d13[, .(r = cor(x, y), x = min(x), y = max(y)), by = anc]
p13 <- ggplot(d13, aes(x, y, color = anc)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.8, color = "grey40") +
  geom_point(size = 2.6, alpha = 0.75) +
  geom_text(data = r13, aes(x = x, y = y, label = sprintf("r = %.2f", r)),
            hjust = 0, vjust = 1, size = 5, fontface = "bold", show.legend = FALSE) +
  facet_wrap(~ anc, scales = "free", nrow = 2) +
  scale_color_anc(guide = "none") +
  labs(x = expression("All by All  " * group("|", beta, "|")),
       y = expression("FELIXassoc  " * group("|", beta, "|")),
       title = "Per-ancestry effect-size concordance") +
  theme(strip.text = element_text(size = 20, face = "bold"),
        strip.background = element_blank())
save_fig(p13, "fig13_beta_concordance_per_ancestry", width = 13.5, height = 9)
