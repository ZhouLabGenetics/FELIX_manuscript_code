#!/usr/bin/env Rscript
# Unified discovery-mechanism example figures (FELIX). One consistent 3-panel template per example,
# each quantity a DISTINCT encoding:
#   (i)   per-ancestry effective sample size  -> grouped bars (All by All vs FELIXassoc)
#   (ii)  per-ancestry alternate-allele freq   -> Cleveland dot plot (missing value = simply absent)
#   (iii) per-ancestry effect size             -> forest (beta +/- 95% CI)

source("scripts_R/00_theme.R"); source("scripts_R/01_data.R")
suppressMessages(library(patchwork))
dt  <- load_scatter()
ssL <- fread(file.path(REPLIC, "sample_size_table_long.tsv"))

DISP <- c("AFR","EAS","EUR","AMR","SAS")                 # display order; NatAm/NAT shown as AMR
SUF  <- setNames(paste0("anc", 1:5), DISP)               # SAIGE ancestry suffix (anc4 = AMR)
ABAP <- setNames(c("ABA_AFR","ABA_EAS","ABA_EUR","ABA_AMR","ABA_SAS"), DISP)
SS_ANC <- setNames(c("AFR","EAS","EUR","NatAm","SAS"), DISP)  # sample-size table uses NatAm
COLS <- ANC_COLORS[DISP]
fmtk <- function(v) ifelse(!is.finite(v), "", ifelse(v >= 1000, sprintf("%.1fk", v/1000), sprintf("%.0f", v)))
BASE <- 22

## (i) effective sample size — grouped bars (ABA shaded, FELIXassoc solid)
panel_N <- function(ph) {
  d <- ssL[phenotype == ph & ancestry %in% SS_ANC,
           .(ancestry, `All by All` = as.numeric(ABA_N), `FELIXassoc` = as.numeric(TRACTOR_N))]
  d[, ancestry := names(SS_ANC)[match(ancestry, SS_ANC)]]
  m <- melt(d, id.vars = "ancestry", variable.name = "method", value.name = "N")
  m[, ancestry := factor(ancestry, levels = DISP)][, method := factor(method, levels = c("All by All","FELIXassoc"))]
  m <- m[is.finite(N) & N > 0]
  off <- 0.22; m[, x := as.integer(ancestry) + ifelse(method == "All by All", -off, off)]
  ggplot(m, aes(x, N, fill = ancestry, alpha = method)) +
    geom_col(width = 0.4, colour = "grey35", linewidth = 0.3) +
    geom_text(aes(label = fmtk(N)), vjust = -0.5, size = 4.4, alpha = 1) +
    scale_fill_manual(values = COLS, guide = "none") +
    scale_alpha_manual(values = c("All by All" = 0.42, "FELIXassoc" = 1.0), name = NULL,
                       guide = guide_legend(override.aes = list(fill = "grey40"))) +
    scale_x_continuous(breaks = seq_along(DISP), labels = DISP) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.20)),
                       labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
    labs(x = NULL, y = "Effective sample size", title = "Sample size") +
    theme(legend.position = "top", plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(colour = COLS, face = "bold"))
}

## (ii) allele frequency — Cleveland dot plot (hollow = All by All, solid = FELIXassoc)
panel_AF <- function(r) {
  w <- rbindlist(lapply(DISP, function(a) data.table(anc = a,
    aba = as.numeric(r[[paste0(ABAP[a], "_AF_Allele2")]]),
    fx  = as.numeric(r[[paste0("SAIGE_AF_Allele2_", SUF[a])]]))))
  w[, anc := factor(anc, levels = rev(DISP))]
  # draw FELIXassoc (solid) first, All by All (hollow) last so a coincident ring stays visible
  pts <- rbind(w[is.finite(fx),  .(anc, af = fx,  method = "FELIXassoc")],
               w[is.finite(aba), .(anc, af = aba, method = "All by All")])
  pts[, method := factor(method, levels = c("All by All", "FELIXassoc"))]
  ggplot() +
    geom_segment(data = w, aes(x = aba, xend = fx, y = anc, yend = anc, colour = anc),
                 linewidth = 0.8, na.rm = TRUE) +
    geom_point(data = pts, aes(af, anc, colour = anc, shape = method, size = method), stroke = 1.6) +
    geom_text(data = w[is.finite(fx)], aes(fx, anc, label = sprintf("%.2f", fx)),
              vjust = -1.2, size = 4.0, colour = "grey15") +
    scale_shape_manual(values = c("All by All" = 1, "FELIXassoc" = 19), name = NULL,
                       guide = guide_legend(override.aes = list(colour = "grey25", size = 4))) +
    scale_size_manual(values = c("All by All" = PTBIG + 2.5, "FELIXassoc" = PTBIG), guide = "none") +
    scale_color_manual(values = COLS, guide = "none") +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, .5, 1)) +
    labs(x = "Alt-allele frequency", y = NULL, title = "Allele frequency") +
    theme(legend.position = "top", plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.y = element_text(colour = rev(COLS), face = "bold"))
}

## (iii) effect size — forest (open square = All by All, filled circle = FELIXassoc)
panel_beta <- function(r) {
  d <- rbindlist(lapply(DISP, function(a) data.table(anc = a,
    st = as.numeric(r[[paste0("SAIGE_BETA_c_", SUF[a])]]), sse = as.numeric(r[[paste0("SAIGE_SE_c_", SUF[a])]]),
    ab = as.numeric(r[[paste0(ABAP[a], "_BETA")]]),        abse = as.numeric(r[[paste0(ABAP[a], "_SE")]]))))
  L <- rbind(d[is.finite(st), .(anc, method = "FELIXassoc", b = st, se = sse)],
             d[is.finite(ab), .(anc, method = "All by All",  b = ab, se = abse)])
  L[, anc := factor(anc, levels = rev(DISP))]
  ggplot(L, aes(b, anc, colour = anc, shape = method)) +
    geom_vline(xintercept = 0, colour = "grey60") +
    geom_errorbarh(aes(xmin = b - 1.96*se, xmax = b + 1.96*se), height = 0, linewidth = 1.1,
                   position = position_dodge(0.55)) +
    geom_point(size = PTBIG + 1, fill = "white", stroke = 1.3, position = position_dodge(0.55)) +
    scale_shape_manual(values = c("FELIXassoc" = 16, "All by All" = 22), name = NULL) +
    scale_color_manual(values = COLS, guide = "none") +
    labs(x = expression(beta * "  (95% CI)"), y = NULL, title = "Effect size") +
    theme(legend.position = "top", plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.y = element_text(colour = rev(COLS), face = "bold"))
}

## (iv) significance — horizontal bars of -log10 P for All by All (meta) and the three FELIX
## tests (HOM, HET, CCT). A dashed line marks genome-wide significance. This shows FELIX
## recovering signal at loci where All by All falls short (the discovery mechanism).
panel_P <- function(r) {
  META_LAB <- "Global-Ancestry\nMeta-analysis"                 # two lines so it fits the axis
  # Monochrome: bar fill runs light grey -> near-black with -log10 P, so the darker the bar the
  # more significant. A single dark-red dashed line marks genome-wide significance.
  GW_RED <- "#B2182B"
  d <- data.table(
    lab = c(META_LAB, "FELIX HOM", "FELIX HET", "FELIX CCT"),
    p   = as.numeric(c(r$ABA_META_Pvalue, r$SAIGE_P_hom_admixed_c,
                       r$SAIGE_P_het_admixed_c, r$SAIGE_P_cct_admixed_c)))
  d <- d[is.finite(p)]
  d[, mlp := pmin(-log10(pmax(p, 1e-300)), 100)]
  d[, lab := factor(lab, levels = c("FELIX CCT", "FELIX HET", "FELIX HOM", META_LAB))]
  xmax <- max(d$mlp, GW_LOG) * 1.28
  ggplot(d, aes(mlp, lab, fill = mlp)) +
    geom_col(width = 0.62, colour = "grey20", linewidth = 0.35) +
    geom_vline(xintercept = GW_LOG, linetype = "dashed", linewidth = 1.1, colour = GW_RED) +
    geom_text(aes(label = formatC(p, format = "e", digits = 1)),
              hjust = -0.1, size = 4.0, colour = "grey15") +
    scale_fill_gradient(low = "#DBDBDB", high = "#1A1A1A", guide = "none") +
    scale_x_continuous(limits = c(0, xmax), expand = expansion(mult = c(0, 0.02))) +
    labs(x = expression(-log[10] * " P"), y = NULL, title = "Significance") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.y = element_text(face = "bold", colour = "grey15"))
}

## examples: gene, phenotype id, pos, rsid, file  (positions/rsids match the manuscript)
EX <- list(
  list("IL23R","GI_522.11",67240275,"rs11209026","fig_ex_il23r"),
  list("ADRB2","BMI",148898672,"rs9325126","fig_ex_adrb2"),
  list("PNPLA3","3013721",43928850,"rs738409","fig_ex_pnpla3"),
  list("APOC1","3035995",44919689,"rs4420638","fig_ex_apoc1"),
  list("HPR","3028288",72080103,"rs217181","fig_ex_hpr"))
gene_lab <- c(IL23R="IL23R", ADRB2="ADRB2/SH3TC2", PNPLA3="PNPLA3 I148M",
              APOC1="APOC1/APOE", HPR="HPR/TXNL4B")

for (e in EX) {
  gene <- e[[1]]; ph <- e[[2]]; pos <- e[[3]]; rs <- e[[4]]; fn <- e[[5]]
  r <- get_locus(dt, gene, ph, pos); if (is.null(r)) { message("miss ", gene); next }
  ttl <- sprintf("%s  (%s, chr%s:%s)  -  %s", gene_lab[gene], rs, r$SAIGE_CHR, format(pos, big.mark=","), lbl(ph))
  fig <- (panel_N(ph) | panel_AF(r) | panel_beta(r) | panel_P(r)) +
    plot_annotation(title = ttl, theme = theme(plot.title = element_text(size = BASE, face = "bold", hjust = 0.5)))
  save_fig(fig, fn, width = 23, height = 5.6, dpi = 600)
}
message("done: fig_ex_{il23r,adrb2,pnpla3,apoc1,hpr}")
