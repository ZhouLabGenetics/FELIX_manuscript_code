#!/usr/bin/env Rscript
# plot_eur_global_local.R — EUR P+T global-vs-local, in the EXACT style of
# prs/plot_partial_prs.R (AFR) so the two figures sit side-by-side identically:
# same colors, legend labels, dodge/width, theme, alphabetical trait order, ΔR² (proportion).
#   global = All-by-All EUR (matched_EUR) ; local = FELIX EUR-tract (felix_tract_EUR)
# Also writes a standalone high-res legend PNG usable for BOTH the AFR and EUR figures.
#   Rscript plot_eur_global_local.R [summary_5e-8.tsv] [OUT_PREFIX]
suppressMessages({library(ggplot2); library(grid)})
args <- commandArgs(trailingOnly = TRUE)
inf <- ifelse(length(args) >= 1, args[1], "summary_5e-8.tsv")
OUT <- ifelse(length(args) >= 2, args[2], "fig_EUR_global_vs_local")

traits <- c(
  "3006923"="Alanine aminotransferase","3007070"="HDL cholesterol","3009744"="MCHC",
  "3013721"="Aspartate aminotransferase","3022192"="Triglycerides","3024929"="Platelet count",
  "3027114"="Total cholesterol","3028288"="LDL cholesterol","3035995"="Alkaline phosphatase",
  "BMI"="BMI","height"="Height")

# EUR rows: matched_EUR -> global, felix_tract_EUR -> local
d <- read.table(inf, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
d <- d[d$ancestry == "EUR" & d$score %in% c("matched_EUR__eur10k","felix_tract_EUR__eur10k"), ]
d$score <- ifelse(d$score == "matched_EUR__eur10k", "global", "local")
d$trait <- ifelse(d$trait %in% names(traits), traits[d$trait], d$trait)
df <- data.frame(trait = d$trait, score = d$score, estimate = d$incR2,
                 CI_low = d$CI_low, CI_high = d$CI_high, stringsAsFactors = FALSE)

# IDENTICAL to plot_partial_prs.R: fixed order + labels + palette (global, local subset)
ord <- c("global","local","partial","afrDS1")
lab <- c(global="Global (All by All)", local="Local (FELIXassoc)",
         partial="Partial (b_afr*G_afr + b_eur*G_eur)", afrDS1="AFR-segment (b_afr*G_afr)")
pal <- c(global="#6699C4", local="#C77BA6", partial="#5AA5A0", afrDS1="#E1A140")
present <- c(ord[ord %in% df$score], setdiff(unique(df$score), ord))
df$score <- factor(df$score, levels = present)
labs <- ifelse(present %in% names(lab), lab[present], present); names(labs) <- present
cols <- ifelse(present %in% names(pal), pal[present], "#999999"); names(cols) <- present

g <- ggplot(df, aes(trait, estimate, fill = score)) +
  geom_col(position = position_dodge(0.8), width = 0.72) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high),
                position = position_dodge(0.8), width = 0.25, linewidth = 0.5) +
  scale_fill_manual(values = cols, labels = labs, name = NULL) +
  labs(x = NULL, y = expression(Delta*R^2~"(incremental, EUR validation)"),
       title = "EUR polygenic prediction: global vs local PRS",
       subtitle = "bars = point estimate; whiskers = bootstrap 95% CI") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        legend.position = "top", panel.grid.major.x = element_blank())

ggsave(paste0(OUT, ".pdf"), g, width = 12, height = 6)
ggsave(paste0(OUT, ".png"), g, width = 12, height = 6, dpi = 300)

# no-legend version (for using the shared standalone legend instead)
ggsave(paste0(OUT, "_nolegend.pdf"), g + theme(legend.position = "none"), width = 12, height = 6)
ggsave(paste0(OUT, "_nolegend.png"), g + theme(legend.position = "none"), width = 12, height = 6, dpi = 300)

# standalone shared legend (same colors/labels -> use for BOTH AFR and EUR)
gt <- ggplotGrob(g)
leg <- gt$grobs[[which(vapply(gt$grobs, function(x) x$name, "") == "guide-box")]]
png("legend_global_local.png", width = 2600, height = 320, res = 300, bg = "white")
grid.newpage(); grid.draw(leg); dev.off()
pdf("legend_global_local.pdf", width = 8.7, height = 1.07); grid.newpage(); grid.draw(leg); dev.off()
cat("wrote", paste0(OUT, ".pdf/.png"), "+ _nolegend + legend_global_local.png/.pdf\n")
