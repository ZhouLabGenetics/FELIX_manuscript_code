#!/usr/bin/env Rscript
# Chi-square comparison at the All-by-All baseline top hits (the SAME loci as the fig12 beta-
# concordance scatter). For each locus we compute the effect chi-square (beta/SE)^2 under All by All
# and under FELIXassoc. Two outputs:
#   (1) a plain boxplot of the two distributions (log y; heavy chi-square tail is intrinsic), and
#   (2) a per-locus scatter (FELIXassoc vs All by All), which mirrors the beta-concordance scatter:
#       tightly concordant, most points above the y=x line (FELIXassoc stronger), a few outliers.
source("scripts_R/00_theme.R"); source("scripts_R/01_data.R")
dt <- load_scatter()
b <- dt[tophit_source == "ABA"]
b[, ea := (as.numeric(ABA_META_BETA) / as.numeric(ABA_META_SE))^2]
b[, es := (as.numeric(SAIGE_BETA_c_ancALL) / as.numeric(SAIGE_SE_c_ancALL))^2]
b <- b[is.finite(ea) & is.finite(es) & ea > 0 & es > 0]

# (1) plain boxplot
d <- rbind(data.table(method = "All by All", x2 = b$ea),
           data.table(method = "FELIXassoc", x2 = b$es))
d[, method := factor(method, levels = c("All by All", "FELIXassoc"))]
pbox <- ggplot(d, aes(method, x2, fill = method)) +
  geom_boxplot(width = 0.6, alpha = 0.55, linewidth = 0.7, outlier.size = 1.2, outlier.alpha = 0.4) +
  geom_hline(yintercept = GW_CHISQ, linetype = "dashed", linewidth = 0.8, color = "grey50") +
  scale_y_log10() + scale_fill_manual(values = METHOD_COLORS, guide = "none") +
  labs(x = NULL, y = expression(chi^2 == (beta/SE)^2), title = "Signal strength at baseline loci")
save_fig(pbox, "SFig_chisq_boxplot_aba_vs_saige", width = 6.5, height = 6.2)

# (1b) SAME boxplot but the chi-square is derived from the reported p-value
#      (qchisq(p, df = 1, lower.tail = FALSE)) instead of (beta/SE)^2. This reflects each
#      method's headline test significance: All by All meta p vs FELIX Cauchy-combined p,
#      which includes the heterogeneous component. Juxtapose against panel (1).
bp <- dt[tophit_source == "ABA"]
bp[, pa := p2chisq(ABA_META_Pvalue)]
bp[, ps := p2chisq(SAIGE_P_cct_admixed_c)]
bp <- bp[is.finite(pa) & is.finite(ps) & pa > 0 & ps > 0]
dp <- rbind(data.table(method = "All by All", x2 = bp$pa),
            data.table(method = "FELIXassoc", x2 = bp$ps))
dp[, method := factor(method, levels = c("All by All", "FELIXassoc"))]
pboxp <- ggplot(dp, aes(method, x2, fill = method)) +
  geom_boxplot(width = 0.6, alpha = 0.55, linewidth = 0.7, outlier.size = 1.2, outlier.alpha = 0.4) +
  geom_hline(yintercept = GW_CHISQ, linetype = "dashed", linewidth = 0.8, color = "grey50") +
  scale_y_log10() + scale_fill_manual(values = METHOD_COLORS, guide = "none") +
  labs(x = NULL, y = expression(chi^2 == "qchisq(P, df=1)"),
       title = "Signal strength from reported P")
save_fig(pboxp, "SFig_chisq_boxplot_fromP_aba_vs_saige", width = 6.5, height = 6.2)
message(sprintf("from-P boxplot medians ABA=%.0f FELIX=%.0f (n=%d)",
                median(bp$pa), median(bp$ps), nrow(bp)))

# (2) chi-square concordance scatter (mirror of the beta scatter)
r <- cor(log(b$ea), log(b$es)); lim <- range(c(b$ea, b$es))
psc <- ggplot(b, aes(ea, es)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.9, color = "grey40") +
  geom_point(size = 2.8, alpha = 0.6, color = METHOD_COLORS[["FELIXassoc"]]) +
  annotate("text", x = lim[1], y = lim[2], hjust = 0, vjust = 1, size = 6, fontface = "bold",
           label = sprintf("Pearson r = %.2f", r)) +
  scale_x_log10() + scale_y_log10() + coord_equal() +
  labs(x = expression("All by All  " * chi^2), y = expression("FELIXassoc  " * chi^2),
       title = "Chi-square concordance (baseline loci)")
save_fig(psc, "SFig_chisq_scatter_aba_vs_saige", width = 7, height = 7)
message(sprintf("boxplot medians ABA=%.0f FELIX=%.0f ; scatter log-r=%.3f n=%d",
                median(b$ea), median(b$es), r, nrow(b)))
