#!/usr/bin/env Rscript
# Nature-standard per-locus "triptych" main figures backing the two discovery mechanisms.
#   Fig 3 (Mechanism A - inclusion + AF refinement): IL23R-Crohn's, ADRB2/SH3TC2-BMI, PNPLA3 I148M-AST
#   Fig 4 (Mechanism B - cross-ancestry heterogeneity): APOC1-ALP, HPR/TXNL4B-LDL
# Each locus = 3 panels: (i) regional ancestry-resolved association, (ii) per-ancestry
# effect forest (FELIXassoc vs All by All), (iii) mechanism "money" panel.
# Locked CB-safe palette; no free-floating on-plot text (labels only in title/axis/legend,
# except the exact AF numbers the user asked to keep on the PNPLA3 bars).
source("scripts_R/00_theme.R"); source("scripts_R/01_data.R")
suppressMessages(library(patchwork))
dt <- load_scatter()
AN  <- ANC_MAP$name; SUF <- ANC_MAP$suf; ABAP <- ANC_MAP$aba

## ---------- (i) regional, ancestry-resolved (stacked SAIGE tracks + All by All meta) ----------
panel_regional <- function(sp, ph, chrom, pos, gene, window = 2.5e5) {
  sg <- saige_region(sp, chrom, pos, window)
  long <- rbindlist(lapply(1:5, function(i)
    data.table(POS = sg$POS, anc = AN[i], y = mlog(sg[[paste0("p.value_c_anc", i)]]))))
  long <- long[is.finite(y)]
  keep <- long[, .N, by = anc][N >= 3]$anc
  keep <- intersect(AN, keep)
  long <- long[anc %in% keep]
  ab <- aba_region(ph, "META", chrom, pos, window)
  aba_dt <- if (nrow(ab)) data.table(POS = ab$POS, anc = "All by All (meta)", y = mlog(ab$Pvalue))[is.finite(y)] else NULL
  both <- rbind(long, aba_dt)
  lev  <- c(keep, if (!is.null(aba_dt)) "All by All (meta)")
  both[, anc := factor(anc, levels = lev)]
  cols <- c(ANC_COLORS[keep]); cols["All by All (meta)"] <- META_COLOR
  ggplot(both, aes(POS / 1e6, y, color = anc)) +
    geom_vline(xintercept = pos / 1e6, linetype = "dotted", color = "grey55") +
    geom_hline(yintercept = GW_LOG, linetype = "dashed", color = "grey45") +
    geom_point(size = 1.4, alpha = 0.8) +
    facet_grid(anc ~ ., switch = "y") +
    scale_color_manual(values = cols, guide = "none") +
    scale_y_continuous(n.breaks = 3, expand = expansion(mult = c(0.06, 0.12))) +
    labs(x = paste0("chr", chrom, " (Mb)"), y = expression(-log[10] * p), title = gene) +
    theme(strip.text.y.left = element_text(angle = 0, face = "bold", size = 11),
          strip.background = element_blank(), strip.placement = "outside",
          axis.text.y = element_text(size = 8), panel.spacing = unit(2, "pt"))
}

## ---------- (ii) per-ancestry effect forest: FELIXassoc vs All by All ----------
panel_forest <- function(r, gene) {
  d <- rbindlist(lapply(1:5, function(i) data.table(anc = AN[i],
    st = as.numeric(r[[paste0("SAIGE_BETA_c_", SUF[i])]]), sse = as.numeric(r[[paste0("SAIGE_SE_c_", SUF[i])]]),
    ab = as.numeric(r[[paste0(ABAP[i], "_BETA")]]),   abse = as.numeric(r[[paste0(ABAP[i], "_SE")]]))))
  d <- d[is.finite(st) | is.finite(ab)]
  ord <- d[order(st, na.last = TRUE)]$anc
  L <- rbind(d[is.finite(st), .(anc, method = "FELIXassoc", b = st, se = sse)],
             d[is.finite(ab), .(anc, method = "All by All",    b = ab, se = abse)])
  L[, anc := factor(anc, levels = ord)]
  cct <- as.numeric(r$SAIGE_P_cct_admixed_c); het <- as.numeric(r$SAIGE_P_het_admixed_c)
  ggplot(L, aes(b, anc, color = anc, shape = method)) +
    geom_vline(xintercept = 0, color = "grey60") +
    geom_errorbarh(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se), height = 0,
                   linewidth = 1.0, position = position_dodge(0.55)) +
    geom_point(size = PTBIG, fill = "white", stroke = 1.1, position = position_dodge(0.55)) +
    scale_shape_manual(values = c("FELIXassoc" = 16, "All by All" = 22), name = NULL) +
    scale_color_anc(guide = "none") +
    labs(x = expression(beta * "  (95% CI)"), y = NULL, title = "Per-ancestry effect",
         subtitle = sprintf("CCT p = %.0e   HET p = %.0e", cct, het)) +
    theme(legend.position = "top")
}

## ---------- generic grouped bars (All by All shaded/left, FELIXassoc solid/right) ----------
gbars <- function(d, ylab, title, subtitle = NULL, fmt = NULL, gwline = FALSE, ymax = NULL) {
  ancs <- levels(d$ancestry); meths <- levels(d$method); off <- 0.205; bw <- 0.38
  d[, xi := as.integer(ancestry)][, xpos := xi + ifelse(method == meths[1], -off, off)]
  d[, tested := is.finite(value)]
  p <- ggplot() +
    geom_col(data = d[tested == TRUE], aes(xpos, value, fill = ancestry, alpha = method),
             width = bw, colour = "grey35", linewidth = 0.3) +
    scale_fill_anc(guide = "none") +
    scale_alpha_manual(values = setNames(c(0.45, 1.0), meths), name = NULL,
                       guide = guide_legend(override.aes = list(fill = "grey40"))) +
    scale_x_continuous(breaks = seq_along(ancs), labels = ancs) +
    labs(x = NULL, y = ylab, title = title, subtitle = subtitle) + theme(legend.position = "top")
  if (!is.null(fmt)) p <- p + geom_text(data = d[tested == TRUE], aes(xpos, value, label = fmt(value)),
                                        vjust = -0.4, size = 3.5)
  if (any(!d$tested)) p <- p + geom_text(data = d[tested == FALSE], aes(xpos, 0, label = "not tested"),
                                         colour = "grey60", angle = 90, hjust = -0.03, size = 3.1, fontface = "italic")
  if (gwline) p <- p + geom_hline(yintercept = GW_LOG, linetype = "dashed", colour = "grey45")
  p <- p + scale_y_continuous(expand = expansion(mult = c(0, 0.16)),
                              limits = if (!is.null(ymax)) c(0, ymax) else NULL)
  p + theme(axis.text.x = element_text(colour = ANC_COLORS[ancs], face = "bold"))
}

## ---------- (iii-A) discovery-gain: per-ancestry -log10 p, All by All vs FELIXassoc ----------
panel_gain <- function(r, gene) {
  d <- rbindlist(lapply(1:5, function(i) data.table(ancestry = AN[i],
    `All by All`    = mlog(r[[paste0(ABAP[i], "_Pvalue")]]),
    `FELIXassoc` = mlog(r[[paste0("SAIGE_p.value_c_", SUF[i])]]))))
  m <- melt(d, id.vars = "ancestry", variable.name = "method", value.name = "value")
  m[, ancestry := factor(ancestry, levels = AN)][, method := factor(method, levels = c("All by All", "FELIXassoc"))]
  gbars(m, expression("association  " * -log[10] * p), "Discovery gain", gwline = TRUE)
}

## ---------- (iii-B) PNPLA3 allele-frequency refinement (global cluster vs local ancestry) ----------
panel_af <- function(r, gene, title) {
  d <- rbindlist(lapply(1:5, function(i) data.table(ancestry = AN[i],
    `All by All` = as.numeric(r[[paste0(ABAP[i], "_AF_Allele2")]]),
    `FELIXassoc` = as.numeric(r[[paste0("SAIGE_AF_Allele2_", SUF[i])]]))))
  m <- melt(d, id.vars = "ancestry", variable.name = "method", value.name = "value")
  m[, ancestry := factor(ancestry, levels = AN)][, method := factor(method, levels = c("All by All", "FELIXassoc"))]
  gbars(m, "alternate-allele frequency", title, fmt = function(v) sprintf("%.2f", v), ymax = 1.0)
}

## ---------- (iii-C) pooling-cancels: per-ancestry effects vs one naive pooled effect ----------
panel_pool <- function(r, gene, hou = NULL) {
  d <- rbindlist(lapply(1:5, function(i) data.table(anc = AN[i],
    b = as.numeric(r[[paste0("SAIGE_BETA_c_", SUF[i])]]), se = as.numeric(r[[paste0("SAIGE_SE_c_", SUF[i])]]))))
  d <- d[is.finite(b) & is.finite(se) & se > 0]
  w <- 1 / d$se^2; pb <- sum(d$b * w) / sum(w); pse <- sqrt(1 / sum(w))
  D <- rbind(d[, .(anc, b, se, grp = "per-ancestry")],
             data.table(anc = "Pooled (1 effect)", b = pb, se = pse, grp = "pooled"))
  D[, anc := factor(anc, levels = c("Pooled (1 effect)", rev(d[order(b)]$anc)))]
  het <- as.numeric(r$SAIGE_P_het_admixed_c)
  fillv <- c(ANC_COLORS); fillv["Pooled (1 effect)"] <- META_COLOR
  ggplot(D, aes(b, anc)) +
    geom_vline(xintercept = 0, color = "grey60") +
    geom_errorbarh(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se, color = grp), height = 0, linewidth = 1.1) +
    geom_point(aes(fill = anc), size = PTBIG, shape = 21, stroke = 1.1, color = "grey25") +
    scale_color_manual(values = c("per-ancestry" = "grey35", "pooled" = META_COLOR), guide = "none") +
    scale_fill_manual(values = fillv, guide = "none") +
    labs(x = expression(beta * "  (95% CI)"), y = NULL, title = "Pooling cancels",
         subtitle = sprintf("HET p = %.0e%s", het, if (!is.null(hou)) paste0("\n", hou) else ""))
}

## ============================ FIGURE 3 - Mechanism A ============================
il <- get_locus(dt, "IL23R",  "GI_522.11", 67240275)
ad <- get_locus(dt, "ADRB2",  "BMI",       148898672)
pn <- get_locus(dt, "PNPLA3", "3013721",   43928850)

r_il <- panel_regional("pheno_GI_522.11", "GI_522.11", 1,  67240275, "IL23R")
r_ad <- panel_regional("BMI",             "BMI",       5,  148898672, "ADRB2/SH3TC2")
r_pn <- panel_regional("pheno_3013721",   "3013721",   22, 43928850, "PNPLA3")

fig3 <- (r_il | panel_forest(il, "IL23R")        | panel_gain(il, "IL23R")) /
        (r_ad | panel_forest(ad, "ADRB2/SH3TC2") | panel_gain(ad, "ADRB2/SH3TC2")) /
        (r_pn | panel_forest(pn, "PNPLA3")        | panel_af(pn, "PNPLA3", "AF refinement (rs738409)")) +
  plot_annotation(tag_levels = "a", title = "Mechanism A - inclusion and allele-frequency refinement",
                  theme = theme(plot.title = element_text(size = 22, face = "bold")))
save_fig(fig3, "fig3_mechanismA_triptych", width = 18, height = 15)

## ============================ FIGURE 4 - Mechanism B ============================
ap <- get_locus(dt, "APOC1", "3035995", 44919689)
hp <- get_locus(dt, "HPR",   "3028288", 72080103)

r_ap <- panel_regional("pheno_3035995", "3035995", 19, 44919689, "APOC1/APOE")
r_hp <- panel_regional("pheno_3028288", "3028288", 16, 72080103, "HPR/TXNL4B")
hou_note <- "Hou S11: HET=2.9e-8"

fig4 <- (r_ap | panel_forest(ap, "APOC1/APOE") | panel_pool(ap, "APOC1/APOE", hou_note)) /
        (r_hp | panel_forest(hp, "HPR/TXNL4B") | panel_pool(hp, "HPR/TXNL4B")) +
  plot_annotation(tag_levels = "a", title = "Mechanism B - cross-ancestry effect heterogeneity",
                  theme = theme(plot.title = element_text(size = 22, face = "bold")))
save_fig(fig4, "fig4_mechanismB_triptych", width = 18, height = 10.5)

message("done: fig3_mechanismA_triptych, fig4_mechanismB_triptych")
