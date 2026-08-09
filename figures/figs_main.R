#!/usr/bin/env Rscript

source("scripts_R/00_theme.R"); source("scripts_R/01_data.R")
suppressMessages(library(patchwork))
ml <- function(p) -log10(pmax(as.numeric(p), 1e-320))

## Fig 1 — sample size (ABA global stacked, SAIGE local stacked)
ss <- fread(file.path(REPLIC, "sample_size_table_long.tsv"))
ss[, `:=`(ABA_N = as.numeric(ABA_N), TRACTOR_N = as.numeric(TRACTOR_N))]
ord <- ss[, .(tot = sum(ABA_N, na.rm = TRUE)), by = label][order(tot)]$label
ss[, label := factor(label, levels = ord)]
trunc32 <- function(x) ifelse(nchar(x) > 32, paste0(substr(x, 1, 30), "..."), x)
mk_ss <- function(col, ancs, title, xlab) {
  d <- ss[ancestry %in% ancs & is.finite(get(col)) & get(col) > 0]
  d[, ancestry := factor(ancestry, levels = ancs)]
  ggplot(d, aes(get(col), label, fill = ancestry)) +
    geom_col(width = 0.78, color = "white", linewidth = 0.3) +
    scale_fill_manual(values = ANC_COLORS, name = "Ancestry") +
    scale_x_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()),
                       n.breaks = 4) +
    scale_y_discrete(labels = trunc32) +
    labs(x = xlab, y = NULL, title = title)
}
p1a <- mk_ss("ABA_N", c("AFR","EAS","EUR","NatAm","SAS","MID"),
             "a  All by All (global N)", "Sample size (per global stratum)")
p1b <- mk_ss("TRACTOR_N", c("AFR","EAS","EUR","NatAm","SAS"),
             "b  FELIXassoc (local effective N)", "Effective N (per local test)") +
       theme(axis.text.y = element_blank())
save_fig(p1a | p1b, "fig1_sample_size_panel", width = 16, height = 9)

## Fig 2 — discovery counts per phenotype
dt <- load_scatter()
cnt <- dt[, .(
  Shared = sum(tophit_source == "SAIGE" & locus_status == "shared"),
  `FELIXassoc only` = sum(tophit_source == "SAIGE" & locus_status == "SAIGE_only"),
  `All by All only` = sum(tophit_source == "ABA" & locus_status == "ABA_only")), by = phenotype]
cnt[, label := sapply(phenotype, lbl)][, total := Shared + `FELIXassoc only` + `All by All only`]
cnt <- cnt[total > 0]; cnt[, label := factor(label, levels = cnt[order(total)]$label)]
m2 <- melt(cnt, id.vars = "label", measure.vars = names(STATUS_COLORS),
           variable.name = "status", value.name = "n")
m2[, status := factor(status, levels = names(STATUS_COLORS))]
p2 <- ggplot(m2, aes(n, label, fill = status)) +
  geom_col(width = 0.78) + scale_fill_manual(values = STATUS_COLORS, name = NULL) +
  scale_y_discrete(labels = function(x) ifelse(nchar(x) > 32, paste0(substr(x,1,30),"..."), x)) +
  labs(x = "Number of independent loci", y = NULL, title = "Locus discovery per phenotype") +
  theme(legend.position = "top")
save_fig(p2, "fig2_discovery_counts", width = 10, height = 9)

## Fig 3 — shared-locus p concordance, CCT/HOM/HET
mk3 <- function(tbl, lab) {
  d <- load_scatter(tbl)[locus_status == "shared" & tophit_source == "SAIGE"]
  data.table(test = lab, x = ml(d$SAIGE_pvalue), y = ml(d$ABA_pvalue))
}
d3 <- rbind(mk3("table_P_cct_admixed_c_vs_META.tsv", "CCT (combined)"),
            mk3("table_P_hom_admixed_c_vs_META.tsv", "HOM (homogeneous)"),
            mk3("table_P_het_admixed_c_vs_META.tsv", "HET (heterogeneous)"))
d3 <- d3[is.finite(x) & is.finite(y)]
d3[, stronger := ifelse(x - y > 1, "FELIXassoc", ifelse(y - x > 1, "All by All", "similar"))]
d3[, test := factor(test, levels = c("CCT (combined)","HOM (homogeneous)","HET (heterogeneous)"))]
p3 <- ggplot(d3, aes(x, y, color = stronger)) +
  geom_abline(slope = 1, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = GW_LOG, color = "grey70", linewidth = 0.4) +
  geom_vline(xintercept = GW_LOG, color = "grey70", linewidth = 0.4) +
  geom_point(size = 2.4, alpha = 0.75) + facet_wrap(~ test, nrow = 1) +
  scale_color_manual(values = c("FELIXassoc" = "#AA3377", "All by All" = "#0077BB",
                                "similar" = "#CCBB44"), name = NULL) +
  coord_equal() +
  labs(x = expression("FELIXassoc  " * -log[10](p)),
       y = expression("All by All  " * -log[10](p)),
       title = "Shared-locus signal strength") +
  theme(strip.background = element_blank(), strip.text = element_text(size = 18, face = "bold"),
        legend.position = "top")
save_fig(p3, "fig3_shared_pvalue_scatter", width = 18, height = 7)

## Fig 4 — per-ancestry p concordance (5 facets)
amap <- data.table(tbl = paste0("table_p.value_c_anc", 1:5, "_vs_", c("AFR","EAS","EUR","AMR","SAS"), ".tsv"),
                   anc = c("AFR","EAS","EUR","NatAm","SAS"))
d4 <- rbindlist(lapply(seq_len(5), function(i) {
  d <- tryCatch(load_scatter(amap$tbl[i]), error = function(e) NULL); if (is.null(d)) return(NULL)
  d <- d[locus_status == "shared" & tophit_source == "SAIGE"]
  data.table(anc = amap$anc[i], x = ml(d$SAIGE_pvalue), y = ml(d$ABA_pvalue))
}))
d4 <- d4[is.finite(x) & is.finite(y)]; d4[, anc := factor(anc, levels = amap$anc)]
p4 <- ggplot(d4, aes(x, y, color = anc)) +
  geom_abline(slope = 1, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = GW_LOG, color = "grey70", linewidth = 0.4) +
  geom_vline(xintercept = GW_LOG, color = "grey70", linewidth = 0.4) +
  geom_point(size = 2.4, alpha = 0.75) + facet_wrap(~ anc, nrow = 1, scales = "free") +
  scale_color_anc(guide = "none") +
  labs(x = expression("FELIXassoc  " * -log[10](p)),
       y = expression("All by All  " * -log[10](p)), title = "Per-ancestry signal strength") +
  theme(strip.background = element_blank(), strip.text = element_text(size = 18, face = "bold"))
save_fig(p4, "fig4_per_ancestry_scatter", width = 22, height = 6)

## Fig 5 — LA-conditioning boost waterfall (top 25 SAIGE-only)
so <- dt[tophit_source == "SAIGE" & locus_status == "SAIGE_only"]
so[, boost := ml(SAIGE_P_cct_admixed_c) - ml(SAIGE_P_cct_admixed)]
so5 <- so[is.finite(boost)][order(-boost)][1:min(25, .N)]
so5[, lab := paste0(tstrsplit(Gene, ";")[[1]], " (", sapply(phenotype, lbl), ")")]
so5[, lab := factor(lab, levels = rev(lab))]
p5 <- ggplot(so5, aes(boost, lab)) + geom_col(fill = "#0077BB", width = 0.75) +
  labs(x = expression(Delta ~ -log[10](p) ~ "(local-ancestry conditioning boost)"),
       y = NULL, title = "Local-ancestry conditioning unmasks signals")
save_fig(p5, "fig5_la_conditioning_boost", width = 11, height = 9)

## Fig 6 — HET vs HOM differential (SAIGE-only)
so[, delta := ml(SAIGE_P_het_admixed_c) - ml(SAIGE_P_hom_admixed_c)]
p6 <- ggplot(so[is.finite(delta)], aes(delta)) +
  geom_histogram(bins = 40, fill = "#CCBB44", color = "white") +
  geom_vline(xintercept = c(-2, 0, 2), linetype = c("dashed","solid","dashed"),
             color = c("#0077BB","grey20","#EE7733")) +
  labs(x = expression(Delta ~ -log[10](p) ~ "(HET - HOM)"), y = "FELIXassoc-only loci",
       title = "Heterogeneous vs homogeneous combiner")
save_fig(p6, "fig6_het_vs_hom", width = 9, height = 6)

## Fig 7 — replication barcode (mirrored)
rep <- fread(file.path(REPLIC, "replication_saige_only_in_aba.tsv"))
rep[, `:=`(nl_s = ml(SAIGE_p), nl_a = ml(ABA_META_p_window))]
rep <- rep[is.finite(nl_s)][order(-nl_s)][1:min(60, .N)]
rep[, lab := make.unique(paste0(substr(Gene, 1, 14), " - ", sapply(phenotype, lbl)))]
rep[, lab := factor(lab, levels = rev(lab))]
m7 <- rbind(rep[, .(lab, method = "FELIXassoc", v = nl_s)],
            rep[, .(lab, method = "All by All",   v = -nl_a)])
p7 <- ggplot(m7, aes(v, lab, fill = method)) + geom_col(width = 0.72) +
  geom_vline(xintercept = c(-GW_LOG, GW_LOG), linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = METHOD_COLORS, name = NULL) +
  labs(x = expression(-log[10](p)), y = NULL,
       title = "Replication of FELIXassoc unique loci in All by All (+/-500 kb)") +
  theme(legend.position = "top", axis.text.y = element_text(size = 11))
save_fig(p7, "fig7_replication_barcode", width = 11, height = 11)

## Fig 8 — SAIGE-only reasons (pie; counts in legend, none on the pie)
getv <- function(r, cols) suppressWarnings(as.numeric(sapply(cols, function(cc) r[[cc]])))
classify <- function(r) {
  b <- getv(r, paste0("SAIGE_BETA_c_anc", 1:5)); se <- getv(r, paste0("SAIGE_SE_c_anc", 1:5))
  b <- b[is.finite(b) & is.finite(se) & se > 0]
  has_opp <- any(b > 0) && any(b < 0); rng <- if (length(b)) max(b) - min(b) else 0
  lab <- ml(r$SAIGE_P_cct_admixed_c) - ml(r$SAIGE_P_cct_admixed)
  het <- ml(r$SAIGE_P_het_admixed_c) - ml(r$SAIGE_P_hom_admixed_c)
  sig <- sum(getv(r, paste0("SAIGE_p.value_c_anc", 1:5)) < 1e-4, na.rm = TRUE)
  if (has_opp && rng > 0.03) return("Opposing directions")
  if (is.finite(lab) && lab >= 5) return("LA conditioning boost")
  if (is.finite(het) && het > 5) return("Effect heterogeneity")
  if (sig <= 1 || length(b) <= 1) return("Single-ancestry signal")
  "LA refined resolution"
}
so[, reason := sapply(seq_len(.N), function(i) classify(so[i]))]
rc <- so[, .N, by = reason][order(-N)]
rc[, reason := factor(reason, levels = reason)]
p8 <- ggplot(rc, aes(x = "", y = N, fill = reason)) +
  geom_col(width = 1, color = "white") + coord_polar("y") +
  scale_fill_manual(values = setNames(CAT_COLORS[seq_len(nrow(rc))], levels(rc$reason)),
                    labels = paste0(rc$reason, "  (", rc$N, ")"), name = NULL) +
  labs(title = "Why FELIXassoc finds these loci") +
  theme_void(base_size = 20) +
  theme(plot.title = element_text(face = "bold", size = 22), legend.text = element_text(size = 16))
save_fig(p8, "fig8_saige_only_reasons", width = 9, height = 6)
