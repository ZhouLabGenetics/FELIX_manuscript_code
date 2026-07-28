#!/usr/bin/env Rscript
# Mechanism B (cross-ancestry heterogeneity): per-ancestry forest + 4-test ladder.
# Ports scripts/19_mechB_heterogeneity.py. No on-plot text; large fonts; locked palette.
source("scripts_R/00_theme.R"); source("scripts_R/01_data.R")
suppressMessages(library(patchwork))
dt <- load_scatter()

forest_panel <- function(r, title) {
  d <- rbindlist(lapply(seq_len(5), function(i) {
    suf <- ANC_MAP$suf[i]; nm <- ANC_MAP$name[i]; aba <- ANC_MAP$aba[i]
    data.table(anc = nm,
      `FELIXassoc` = as.numeric(r[[paste0("SAIGE_BETA_c_", suf)]]),
      se_ST = as.numeric(r[[paste0("SAIGE_SE_c_", suf)]]),
      `All by All` = as.numeric(r[[paste0(aba, "_BETA")]]),
      se_AB = as.numeric(r[[paste0(aba, "_SE")]]))
  }))
  d <- d[!is.na(`FELIXassoc`)]
  d[, anc := factor(anc, levels = d[order(`FELIXassoc`)]$anc)]
  L <- rbind(
    d[, .(anc, method = "FELIXassoc", b = `FELIXassoc`, se = se_ST)],
    d[, .(anc, method = "All by All",   b = `All by All`,   se = se_AB)])
  ggplot(L, aes(b, anc, color = anc, shape = method)) +
    geom_vline(xintercept = 0, color = "grey60", linewidth = 0.6) +
    geom_errorbarh(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se), height = 0,
                   linewidth = 1.1, position = position_dodge(width = 0.5)) +
    geom_point(size = PTBIG, fill = "white", stroke = 1.2,
               position = position_dodge(width = 0.5)) +
    scale_shape_manual(values = c("FELIXassoc" = 16, "All by All" = 22), name = NULL) +
    scale_color_anc(guide = "none") +
    labs(x = expression("Per-ancestry effect ("*beta*", 95% CI)"), y = NULL, title = title) +
    theme(legend.position = "top")
}

ladder_panel <- function(r) {
  T <- data.table(
    test = c("All by All\nmeta", "FELIXassoc\nhomogeneous",
             "FELIXassoc\nheterogeneous", "FELIXassoc\nCCT"),
    p = c(as.numeric(r$ABA_META_Pvalue), as.numeric(r$SAIGE_P_hom_admixed_c),
          as.numeric(r$SAIGE_P_het_admixed_c), as.numeric(r$SAIGE_P_cct_admixed_c)),
    grp = c("meta", "HOM", "HET", "CCT"))
  T[, test := factor(test, levels = rev(test))]
  T[, y := mlog(p)]
  ggplot(T, aes(y, test, fill = grp)) +
    geom_col(width = 0.68) +
    geom_vline(xintercept = GW_LOG, linetype = "dashed", linewidth = 0.9, color = "grey40") +
    scale_fill_manual(values = c(meta = META_COLOR, TEST_COLORS), guide = "none") +
    labs(x = expression(-log[10]~italic(p)), y = NULL) +
    coord_cartesian(xlim = c(0, NA))
}

examples <- list(
  B1 = list(gene = "APOC1", ph = "3035995", pos = 44919689,
            title = "APOC1 — alkaline phosphatase", file = "fig4_mechB1_apoc1"),
  B2 = list(gene = "HPR", ph = "3028288", pos = 72080103,
            title = "HPR/TXNL4B — LDL cholesterol", file = "fig4_mechB2_hpr"))

panels <- list()
for (e in examples) {
  r <- get_locus(dt, e$gene, e$ph, e$pos)
  if (is.null(r)) { message("missing ", e$gene); next }
  fp <- forest_panel(r, e$title); lp <- ladder_panel(r)
  save_fig(fp | lp, e$file, width = 13, height = 5.4)
  panels[[e$file]] <- (fp | lp)
}
save_fig(panels[[1]] / panels[[2]], "fig4_mechB_heterogeneity", width = 13, height = 10.6)
