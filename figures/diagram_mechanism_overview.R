#!/usr/bin/env Rscript
# grid of per-ancestry effect forests for worked loci, showing shared vs heterogeneous architecture.

source("scripts_R/00_theme.R"); source("scripts_R/01_data.R")
dt <- load_scatter()
loci <- data.table(
  gene  = c("IL23R","ADRB2","PNPLA3","APOC1","HPR","ZPR1"),
  ph    = c("GI_522.11","BMI","3013721","3035995","3028288","3022192"),
  pos   = c(67240275, 148898672, 43928850, 44919689, 72080103, 116778201),
  label = c("IL23R (Crohn's)","ADRB2 (BMI)","PNPLA3 (AST)","APOC1 (ALP)","HPR (LDL)","ZPR1 (TG)"))
d <- rbindlist(lapply(seq_len(nrow(loci)), function(i) {
  r <- get_locus(dt, loci$gene[i], loci$ph[i], loci$pos[i]); if (is.null(r)) return(NULL)
  data.table(label = loci$label[i], anc = ANC_MAP$name,
             b  = sapply(paste0("SAIGE_BETA_c_", ANC_MAP$suf), function(c) as.numeric(r[[c]])),
             se = sapply(paste0("SAIGE_SE_c_",   ANC_MAP$suf), function(c) as.numeric(r[[c]])))
}))
d <- d[is.finite(b)]; d[, `:=`(anc = factor(anc, levels = rev(ANC_MAP$name)),
                               label = factor(label, levels = loci$label))]
p <- ggplot(d, aes(b, anc, color = anc)) +
  geom_vline(xintercept = 0, color = "grey55", linewidth = 0.6) +
  geom_errorbarh(aes(xmin = b - 1.96*se, xmax = b + 1.96*se), height = 0, linewidth = 1.1) +
  geom_point(size = PT) +
  facet_wrap(~ label, ncol = 3, scales = "free_x") +
  scale_color_anc(guide = "none") +
  labs(x = expression("per-ancestry effect  " * beta * " (95% CI)"), y = NULL,
       title = "Per-ancestry effect architecture across worked loci") +
  theme(strip.background = element_blank(), strip.text = element_text(size = 17, face = "bold"))
save_fig(p, "fig_mechanism_overview_6loci", width = 15, height = 9)
