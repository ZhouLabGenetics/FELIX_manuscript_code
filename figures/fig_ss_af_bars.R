#!/usr/bin/env Rscript
# Ancestry-coloured grouped bars: All by All = shaded (alpha), FELIXassoc = solid,
# same ancestry hue. Value labels kept ON the plot (exact N / AF). "not tested" in grey
# where an ancestry has no value. Regenerates fig4_ss_* and builds the PNPLA3 AF barplot.
source("scripts_R/00_theme.R"); source("scripts_R/01_data.R")
dt <- load_scatter()
ssL <- fread(file.path(REPLIC, "sample_size_table_long.tsv"))
ANCX <- c("AFR","EAS","EUR","NatAm","SAS")

bars_anc <- function(d, codings, ylab, title, file, fmt, ymax = NULL) {
  # explicit x positions so All by All is ALWAYS the left slot and FELIXassoc the
  # right slot, even when one is absent (dodge would otherwise centre the lone bar).
  d[, ancestry := factor(ancestry, levels = ANCX)]
  d[, coding := factor(coding, levels = codings)]
  off <- 0.205; bw <- 0.38
  d[, xi := as.integer(ancestry)]
  d[, xpos := xi + ifelse(coding == codings[1], -off, off)]
  d[, tested := is.finite(value) & value > 0]
  p <- ggplot() +
    geom_col(data = d[tested == TRUE], aes(xpos, value, fill = ancestry, alpha = coding),
             width = bw, colour = "grey35", linewidth = 0.3) +
    geom_text(data = d[tested == TRUE], aes(xpos, value, label = fmt(value)),
              vjust = -0.4, size = 4.1) +
    scale_fill_anc(guide = "none") +
    scale_alpha_manual(values = setNames(c(0.45, 1.0), codings), name = NULL,
                       guide = guide_legend(override.aes = list(fill = "grey35"))) +
    scale_x_continuous(breaks = seq_along(ANCX), labels = ANCX) +
    labs(x = NULL, y = ylab, title = title) + theme(legend.position = "top")
  if (any(!d$tested))
    p <- p + geom_text(data = d[tested == FALSE], aes(xpos, 0, label = "not tested"),
                       colour = "grey60", angle = 90, hjust = -0.03, size = 3.6, fontface = "italic")
  p <- p + scale_y_continuous(expand = expansion(mult = c(0, 0.16)),
                              limits = if (!is.null(ymax)) c(0, ymax) else NULL,
                              labels = scales::label_number(scale_cut = scales::cut_short_scale()))
  cols <- ANC_COLORS[ANCX]
  p <- p + theme(axis.text.x = element_text(colour = cols, face = "bold", size = 15))
  save_fig(p, file, width = 8.4, height = 5.8)
}

## ---- sample-size bars (exact N kept) ----
comma <- function(v) scales::comma(round(v))
mk_ss <- function(pheno, title, file) {
  d <- ssL[phenotype == pheno & ancestry %in% ANCX,
           .(ancestry, `All by All` = as.numeric(ABA_N), `FELIXassoc` = as.numeric(TRACTOR_N))]
  m <- melt(d, id.vars = "ancestry", variable.name = "coding", value.name = "value")
  bars_anc(m, c("All by All","FELIXassoc"), "Effective sample size", title, file, comma)
}
mk_ss("GI_522.11", "IL23R - Crohn's disease", "fig4_ss_il23r")
mk_ss("BMI",       "ADRB2 - BMI",             "fig4_ss_adrb2")
mk_ss("RE_475",    "IKZF3 - asthma",          "fig4_ss_ikzf3")
mk_ss("height",    "SUPT3H - height",         "fig4_ss_supt3h")

## ---- PNPLA3 I148M (rs738409, chr22:43,928,850) allele frequency by ancestry ----
r <- get_locus(dt, "PNPLA3", "3013721", 43928850)
af <- rbindlist(lapply(1:5, function(i) {
  data.table(ancestry = ANC_MAP$name[i],
    `All by All (global cluster)` = as.numeric(r[[paste0(ANC_MAP$aba[i], "_AF_Allele2")]]),
    `FELIXassoc (local ancestry)` = as.numeric(r[[paste0("SAIGE_AF_Allele2_", ANC_MAP$suf[i])]]))
}))
maf <- melt(af, id.vars = "ancestry", variable.name = "coding", value.name = "value")
bars_anc(maf, c("All by All (global cluster)","FELIXassoc (local ancestry)"),
         "Alternate-allele frequency", "PNPLA3 I148M (rs738409, chr22:43,928,850)",
         "fig_PNPLA3_rs738409_AF", function(v) sprintf("%.2f", v), ymax = 1.0)
