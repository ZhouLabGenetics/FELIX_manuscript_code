#!/usr/bin/env Rscript
# Ports scripts/05_extra_figures.py: fig9 meta-anchored scatter, fig10 BMI haplotype bars,
# fig11 method overlap, supp hom/het overlap. On-plot count labels/stat boxes removed.
source("scripts_R/00_theme.R"); source("scripts_R/01_data.R")
suppressMessages(library(patchwork))
ml <- function(p) -log10(pmax(as.numeric(p), 1e-320))
STRONG <- c("FELIXassoc" = "#AA3377", "All by All" = "#0077BB", "similar" = "#CCBB44")

## Fig 9 — All-by-All META tophits vs SAIGE HOM/HET/CCT
dt <- load_scatter()[tophit_source == "ABA"]
mk9 <- function(col, lab) {
  d <- data.table(test = lab, x = ml(dt$ABA_pvalue), y = ml(dt[[col]]))
  d[is.finite(x) & is.finite(y)]
}
d9 <- rbind(mk9("SAIGE_P_hom_admixed_c", "HOM vs All by All"),
            mk9("SAIGE_P_het_admixed_c", "HET vs All by All"),
            mk9("SAIGE_P_cct_admixed_c", "CCT vs All by All"))
d9[, stronger := ifelse(y - x > 1, "FELIXassoc", ifelse(x - y > 1, "All by All", "similar"))]
d9[, test := factor(test, levels = c("HOM vs All by All","HET vs All by All","CCT vs All by All"))]
p9 <- ggplot(d9, aes(x, y, color = stronger)) +
  geom_abline(slope = 1, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = GW_LOG, color = "grey75", linewidth = 0.4) +
  geom_vline(xintercept = GW_LOG, color = "grey75", linewidth = 0.4) +
  geom_point(size = 2.4, alpha = 0.75) + facet_wrap(~ test, nrow = 1) + coord_equal() +
  scale_color_manual(values = STRONG, name = NULL) +
  labs(x = expression("All by All  " * -log[10](p)), y = expression("FELIXassoc  " * -log[10](p)),
       title = "All by All top hits vs FELIXassoc combiners") +
  theme(strip.background = element_blank(), strip.text = element_text(size = 18, face = "bold"),
        legend.position = "top")
save_fig(p9, "fig9_meta_anchored_scatter", width = 18, height = 7)

## Fig 10 — BMI haplotype counts, local (SAIGE) vs global (ABA x2)
bmi_gz <- file.path(SAIGE_SS, "merged_felix_BMI.txt.gz")
hcols <- paste0("N_haplo_anc", 1:5)
bh <- fread(cmd = sprintf("gzcat %s", shQuote(bmi_gz)), select = hcols)
loc <- sapply(hcols, function(c) { v <- as.numeric(bh[[c]]); median(v[v > 0], na.rm = TRUE) })
csv <- fread(file.path(ROOT, "All_by_All_phenotypes_v7_Analyzed_phenotypes.csv"))
csv[, ancestry := toupper(ancestry)]
gN <- function(a) { r <- csv[phenoname == "BMI" & ancestry == a]; if (!nrow(r)) return(NA_real_)
  (sum(as.numeric(r$n_cases), na.rm = TRUE) + sum(as.numeric(r$n_controls), na.rm = TRUE)) * 2 }
fmtk <- function(v) ifelse(!is.finite(v), "", ifelse(v >= 1000, sprintf("%.1fk", v/1000), sprintf("%.0f", v)))
LAB_FX <- "FELIX - median N_haplo per tested variant (local ancestry)"
LAB_AB <- "All by All - sample size x 2 (global ancestry, haplotype equivalent)"
d10 <- data.table(anc = c("AFR","EAS","EUR","AMR","SAS","MID"),
                  local = c(loc, NA), global = sapply(c("AFR","EAS","EUR","AMR","SAS","MID"), gN))
d10[, anc := factor(anc, levels = c("AFR","EAS","EUR","AMR","SAS","MID"))]
m10 <- melt(d10, id.vars = "anc", variable.name = "coding", value.name = "N")
m10[, coding := factor(coding, levels = c("local","global"), labels = c(LAB_FX, LAB_AB))]
p10 <- ggplot(m10[!is.na(N)], aes(anc, N, fill = coding)) +
  geom_col(position = position_dodge2(preserve = "single", padding = 0.1), width = 0.7) +
  geom_text(aes(label = fmtk(N)), position = position_dodge2(width = 0.7, preserve = "single", padding = 0.1),
            vjust = -0.4, size = 5, fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values = setNames(c("#AA3377", "#0077BB"), c(LAB_FX, LAB_AB)), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)),
                     labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
  labs(x = NULL, y = "Haplotypes (BMI)", title = "Haplotype counts per ancestry") +
  theme(legend.position = "top", legend.direction = "vertical")
save_fig(p10, "fig10_bmi_haplotype_bars", width = 13, height = 7.5)

## Fig 11 — method overlap (stacked totals + near-miss breakdown)
dc <- load_scatter()
s <- dc[tophit_source == "SAIGE"]; a <- dc[tophit_source == "ABA"]
n_sh <- sum(s$locus_status == "shared"); n_so <- sum(s$locus_status == "SAIGE_only")
n_mo <- sum(a$locus_status == "ABA_only")
d11a <- rbind(
  data.table(method = "FELIXassoc (CCT)", status = "Shared", n = n_sh),
  data.table(method = "FELIXassoc (CCT)", status = "FELIXassoc only", n = n_so),
  data.table(method = "All by All meta", status = "Shared", n = n_sh),
  data.table(method = "All by All meta", status = "All by All only", n = n_mo))
d11a[, status := factor(status, levels = c("Shared","FELIXassoc only","All by All only"))]
pa <- ggplot(d11a, aes(method, n, fill = status)) + geom_col(width = 0.6) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 4.6,
            colour = "white", fontface = "bold") +
  scale_fill_manual(values = STATUS_COLORS, name = NULL) +
  labs(x = NULL, y = "Independent loci", title = "a  Locus discovery") +
  theme(legend.position = "top", legend.direction = "vertical")
brk <- function(sub, oc) { o <- as.numeric(sub[[oc]])
  c(`Genome-wide` = sum(o < 5e-8, na.rm = TRUE),
    `Near-miss`   = sum(o >= 5e-8 & o < 5e-6, na.rm = TRUE),
    `Weak`        = sum(o >= 5e-6 & o < 0.01, na.rm = TRUE),
    `Genuine miss`= sum(is.na(o) | o >= 0.01)) }
b1 <- brk(dc[tophit_source=="SAIGE" & locus_status=="SAIGE_only"], "ABA_pvalue")
b2 <- brk(dc[tophit_source=="ABA" & locus_status=="ABA_only"], "SAIGE_pvalue")
d11b <- rbind(data.table(set = "FELIXassoc unique", cat = names(b1), n = b1),
              data.table(set = "All by All unique", cat = names(b2), n = b2))
d11b[, cat := factor(cat, levels = c("Genome-wide","Near-miss","Weak","Genuine miss"))]
d11b[, pct := 100 * n / sum(n), by = set]
pb <- ggplot(d11b, aes(pct, set, fill = cat)) + geom_col(width = 0.55) +
  geom_text(data = d11b[n > 0], aes(label = n), position = position_stack(vjust = 0.5),
            size = 4.4, colour = "white", fontface = "bold") +
  scale_fill_manual(values = setNames(CAT_COLORS[1:4], levels(d11b$cat)), name = NULL,
                    guide = guide_legend(nrow = 2)) +
  labs(x = "Percent of unique loci", y = NULL, title = "b  Status in the other method") +
  theme(legend.position = "top")
save_fig((pa | pb) + patchwork::plot_layout(widths = c(1, 1.15)),
         "fig11_method_overlap", width = 18, height = 7)

## supp — CCT/HOM/HET overlap with META
tabs <- list(CCT = "table_P_cct_admixed_c_vs_META.tsv", HOM = "table_P_hom_admixed_c_vs_META.tsv",
             HET = "table_P_het_admixed_c_vs_META.tsv")
ds <- rbindlist(lapply(names(tabs), function(nm) {
  d <- load_scatter(tabs[[nm]]); s <- d[tophit_source=="SAIGE"]; a <- d[tophit_source=="ABA"]
  rbind(data.table(panel=nm, method=paste0("FELIXassoc (",nm,")"), status="Shared", n=sum(s$locus_status=="shared")),
        data.table(panel=nm, method=paste0("FELIXassoc (",nm,")"), status="FELIXassoc only", n=sum(s$locus_status=="SAIGE_only")),
        data.table(panel=nm, method="All by All meta", status="Shared", n=sum(s$locus_status=="shared")),
        data.table(panel=nm, method="All by All meta", status="All by All only", n=sum(a$locus_status=="ABA_only"))) }))
ds[, status := factor(status, levels = c("Shared","FELIXassoc only","All by All only"))]
ds[, panel := factor(panel, levels = c("CCT","HOM","HET"))]
ds[, mlab := ifelse(grepl("SAIGE", method), "FELIXassoc", "All by All")]
psupp <- ggplot(ds, aes(mlab, n, fill = status)) + geom_col(width = 0.6) +
  facet_wrap(~ panel, nrow = 1) + scale_fill_manual(values = STATUS_COLORS, name = NULL) +
  labs(x = NULL, y = "Independent loci", title = "Locus overlap with All by All meta") +
  theme(legend.position = "top", strip.background = element_blank(),
        strip.text = element_text(size = 18, face = "bold"))
save_fig(psupp, "supp_fig_hom_het_overlap", width = 15, height = 6)
