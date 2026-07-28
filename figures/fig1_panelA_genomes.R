#!/usr/bin/env Rscript
# Figure 1 Panel A - realistic admixed-genome karyograms for the 3 narrative individuals
# (drop-in replacement for the "same three individuals - their genomes" sub-panel).
# Local-ancestry segments simulated from a recombination model: breakpoints ~ Poisson along
# the genetic map (g generations since admixture), segment ancestries drawn from each
# individual's global proportions. 22 chromosomes drawn to scale, diploid (2 haplotypes each).
# Colours = the exact palette sampled from the current figure (locked ANC palette).
source("scripts_R/00_theme.R")
set.seed(2024)

CHR_MB <- c(249,242,198,190,182,171,159,145,138,134,135,133,114,107,102,90,83,80,59,64,47,51)
names(CHR_MB) <- 1:22
CM_PER_MB <- 1.15
ANC5  <- c("EUR","AMR","AFR","EAS","SAS")                       # legend order matches current fig
cols5 <- c(EUR = ANC_COLORS[["EUR"]], AMR = ANC_COLORS[["AMR"]], AFR = ANC_COLORS[["AFR"]],
           EAS = ANC_COLORS[["EAS"]], SAS = ANC_COLORS[["SAS"]])

# --- simulate one haplotype's local-ancestry mosaic along a chromosome ---
sim_hap <- function(len_mb, g, props) {
  nbp  <- if (g <= 0) 0 else rpois(1, g * len_mb * CM_PER_MB / 100)
  edges <- c(0, sort(runif(nbp, 0, len_mb)), len_mb)
  anc  <- sample(names(props), length(edges) - 1, replace = TRUE, prob = props)
  data.table(x0 = edges[-length(edges)], x1 = edges[-1], anc = anc)
}
sim_indiv <- function(label, props, g) {
  rbindlist(lapply(1:22, function(c) rbindlist(lapply(1:2, function(h) {
    s <- sim_hap(CHR_MB[c], g, props); s[, `:=`(chr = c, hap = h, indiv = label)]; s
  }))))
}

# --- three plausible individuals (2-4 major ancestries; recent-admixture segment sizes) ---
d <- rbindlist(list(
  sim_indiv("Indiv 1  -  Global ancestry = EUR",              c(EUR = 1.0), g = 0),
  sim_indiv("Indiv 2  -  Global cluster = AMR (NatAm)",        c(EUR = 0.50, AMR = 0.42, AFR = 0.08), g = 6),
  sim_indiv("Indiv 3  -  Between global clusters (excluded)",  c(AFR = 0.40, EUR = 0.35, AMR = 0.20, EAS = 0.05), g = 7)
))
d[, indiv := factor(indiv, levels = unique(indiv))]
d[, anc  := factor(anc, levels = ANC5)]

# --- y geometry: each chromosome = two stacked haplotype bars ---
ytop <- function(c) -(c - 1)
d[, `:=`(ymin = ifelse(hap == 1, ytop(chr) - 0.02, ytop(chr) - 0.40) - ifelse(hap == 1, 0.34, 0.34),
         ymax = ifelse(hap == 1, ytop(chr) - 0.02, ytop(chr) - 0.40))]
# haplotype outline frames (full chromosome extent) for a crisp karyogram look
frames <- d[, .(x0 = 0, x1 = CHR_MB[chr], ymin = min(ymin), ymax = max(ymax)),
            by = .(indiv, chr, hap)]
frames[, `:=`(ymin = ifelse(hap == 1, ytop(chr) - 0.36, ytop(chr) - 0.74),
              ymax = ifelse(hap == 1, ytop(chr) - 0.02, ytop(chr) - 0.40))]

p <- ggplot() +
  geom_rect(data = frames, aes(xmin = x0, xmax = x1, ymin = ymin, ymax = ymax),
            fill = "grey96", colour = "grey70", linewidth = 0.25) +
  geom_rect(data = d, aes(xmin = x0, xmax = x1, ymin = ymin, ymax = ymax, fill = anc),
            colour = "white", linewidth = 0.08) +
  facet_wrap(~ indiv, nrow = 1) +
  scale_fill_manual(values = cols5, breaks = ANC5, name = "Local ancestry",
                    drop = FALSE, guide = guide_legend(nrow = 1, override.aes = list(colour = "grey70"))) +
  scale_y_continuous(breaks = ytop(1:22) - 0.38, labels = 1:22, expand = expansion(mult = c(0.02, 0.03))) +
  scale_x_continuous(name = "chromosome position (Mb)", expand = expansion(mult = c(0.01, 0.02))) +
  labs(y = "chromosome") +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.text.y = element_text(size = 11),
        strip.text = element_text(size = 15, face = "bold"),
        strip.background = element_blank(),
        panel.spacing = unit(14, "pt"))

save_fig(p, "fig1_panelA_admixed_genomes", width = 15, height = 8.2)
message("done: fig1_panelA_admixed_genomes")
