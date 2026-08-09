#!/usr/bin/env Rscript
# Allele-frequency refinement: per-ancestry local vs cluster AF + per-ancestry effect. 
source("scripts_R/00_theme.R"); source("scripts_R/01_data.R")
suppressMessages(library(patchwork))
dt <- load_scatter()

af_panel <- function(r, title) {
  d <- rbindlist(lapply(seq_len(5), function(i) {
    suf <- ANC_MAP$suf[i]; nm <- ANC_MAP$name[i]; aba <- ANC_MAP$aba[i]
    data.table(anc = nm,
      local   = as.numeric(r[[paste0("SAIGE_AF_Allele2_", suf)]]),
      cluster = as.numeric(r[[paste0(aba, "_AF_Allele2")]]))
  }))
  d <- d[!is.na(local) | !is.na(cluster)]
  d[, anc := factor(anc, levels = ANC_MAP$name)]
  L <- melt(d, id.vars = "anc", variable.name = "coding", value.name = "AF")
  L[, coding := factor(coding, levels = c("cluster", "local"),
        labels = c("All by All (cluster)", "FELIXassoc (local)"))]
  ggplot(L, aes(anc, AF, fill = anc, alpha = coding)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7,
             color = "grey30", linewidth = 0.4) +
    scale_fill_anc(guide = "none") +
    scale_alpha_manual(values = c(0.45, 1.0), name = NULL) +
    labs(x = NULL, y = "Alternate-allele frequency", title = title) +
    theme(legend.position = "top") + ylim(0, 1)
}

effect_panel <- function(r) {
  d <- rbindlist(lapply(seq_len(5), function(i) {
    suf <- ANC_MAP$suf[i]; nm <- ANC_MAP$name[i]; aba <- ANC_MAP$aba[i]
    data.table(anc = nm,
      `FELIXassoc` = as.numeric(r[[paste0("SAIGE_BETA_c_", suf)]]),
      se_ST = as.numeric(r[[paste0("SAIGE_SE_c_", suf)]]),
      `All by All` = as.numeric(r[[paste0(aba, "_BETA")]]),
      se_AB = as.numeric(r[[paste0(aba, "_SE")]]))
  }))
  d <- d[!is.na(`FELIXassoc`)]
  d[, anc := factor(anc, levels = ANC_MAP$name)]
  L <- rbind(d[, .(anc, method = "FELIXassoc", b = `FELIXassoc`, se = se_ST)],
             d[, .(anc, method = "All by All",   b = `All by All`,   se = se_AB)])
  ggplot(L, aes(b, anc, color = anc, shape = method)) +
    geom_vline(xintercept = 0, color = "grey60", linewidth = 0.6) +
    geom_errorbarh(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se), height = 0,
                   linewidth = 1.0, position = position_dodge(width = 0.5)) +
    geom_point(size = PTBIG, fill = "white", stroke = 1.2,
               position = position_dodge(width = 0.5)) +
    scale_shape_manual(values = c("FELIXassoc" = 16, "All by All" = 22), name = NULL) +
    scale_color_anc(guide = "none") +
    labs(x = expression("Per-ancestry effect ("*beta*", 95% CI)"), y = NULL) +
    theme(legend.position = "top")
}

examples <- list(
  list(gene="PNPLA3", ph="3013721", pos=43928850, title="PNPLA3 - aspartate aminotransferase", file="fig_af_refinement_PNPLA3"),
  list(gene="ABO",    ph="3035995", pos=133266804, title="ABO - alkaline phosphatase",         file="fig_af_refinement_ABO"),
  list(gene="ZPR1",   ph="3022192", pos=116778201, title="ZPR1 - triglycerides",               file="fig_af_refinement_ZPR1"))

for (e in examples) {
  r <- get_locus(dt, e$gene, e$ph, e$pos)
  if (is.null(r)) { message("missing ", e$gene); next }
  save_fig(af_panel(r, e$title) | effect_panel(r), e$file, width = 13, height = 5.4)
}
