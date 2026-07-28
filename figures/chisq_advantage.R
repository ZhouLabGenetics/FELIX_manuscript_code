#!/usr/bin/env Rscript
# Chi-square advantage of FELIXassoc over All by All.
#   (a) discovery breadth: each method's own genome-wide loci, detection chi-square (boxplot).
#   (b) signal depth (paired): per baseline locus, FELIXassoc vs All by All effect chi-square
#       (beta/SE)^2 on a log-log scale with the y=x line. Points above the line = FELIXassoc has
#       the stronger signal. The few far-out points are real large-effect loci (HBB sickle cell,
#       HLA) — labeled — not error; the median ratio quantifies the typical advantage.
source("scripts_R/00_theme.R"); source("scripts_R/01_data.R")
suppressMessages(library(patchwork))
dt <- load_scatter()

# (a) breadth — each method's own genome-wide loci, detection chi-square
dt[, det_aba := p2chisq(ABA_META_Pvalue)]
dt[, det_st  := p2chisq(SAIGE_P_cct_admixed_c)]
breadth <- rbind(
  data.table(method = "All by All", x2 = dt[tophit_source == "ABA",   det_aba]),
  data.table(method = "FELIXassoc", x2 = dt[tophit_source == "SAIGE", det_st]))
breadth <- breadth[is.finite(x2)]
breadth[, method := factor(method, levels = c("All by All", "FELIXassoc"))]

pa <- ggplot(breadth, aes(method, x2, fill = method)) +
  geom_boxplot(varwidth = TRUE, width = 0.75, alpha = 0.55, linewidth = 0.7, outlier.shape = NA) +
  geom_hline(yintercept = GW_CHISQ, linetype = "dashed", linewidth = 0.8, color = "grey50") +
  scale_y_log10() + scale_fill_manual(values = METHOD_COLORS, guide = "none") +
  labs(x = NULL, y = expression(chi^2~(1~df)), title = "a  Discovery breadth")

# (b) depth (paired) — baseline genome-wide loci, effect chi-square (beta/SE)^2
db <- dt[tophit_source == "ABA"]
db[, eff_aba := (as.numeric(ABA_META_BETA) / as.numeric(ABA_META_SE))^2]
db[, eff_st  := (as.numeric(SAIGE_BETA_c_ancALL) / as.numeric(SAIGE_SE_c_ancALL))^2]
db <- db[is.finite(eff_aba) & is.finite(eff_st) & eff_aba > 0 & eff_st > 0]
medratio <- median(db$eff_st / db$eff_aba)
frac_up  <- mean(db$eff_st > db$eff_aba)
lab <- db[order(-eff_st)][1:3]
lim <- range(c(db$eff_aba, db$eff_st))

pb <- ggplot(db, aes(eff_aba, eff_st)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.9, color = "grey40") +
  geom_point(size = 2.4, alpha = 0.55, color = METHOD_COLORS[["FELIXassoc"]]) +
  geom_text(data = lab, aes(label = Gene), size = 3.4, hjust = 1.1, vjust = -0.4, check_overlap = TRUE) +
  annotate("text", x = lim[1], y = lim[2], hjust = 0, vjust = 1, size = 5, fontface = "bold",
           label = sprintf("median ratio = %.2f\n%.0f%% of loci above the line", medratio, 100 * frac_up)) +
  scale_x_log10() + scale_y_log10() + coord_equal() +
  labs(x = expression("All by All  " * chi^2 == (beta/SE)^2),
       y = expression("FELIXassoc  " * chi^2 == (beta/SE)^2),
       title = "b  Signal depth (per locus)")

save_fig(pa | pb, "SFig_chisq_advantage_aba_vs_saige", width = 12, height = 6)
message(sprintf("median chisq ratio (FELIXassoc/ABA) = %.3f ; frac above line = %.2f", medratio, frac_up))
