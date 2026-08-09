#!/usr/bin/env Rscript
# create a grid of human silhouettes (BioArt icon)
source("scripts_R/00_theme.R"); source("scripts_R/00_people.R")
suppressMessages(library(grid))
A <- ANC_COLORS
set.seed(3)
ANCS <- c("AFR","EUR","AMR","EAS","SAS")
hx <- function(a) A[[a]]

# composition -> length-10 hex vector (stripe widths ~ proportions); k=1 is solid single-ancestry
mkcomp <- function() {
  k <- sample(1:4, 1, prob = c(0.34, 0.36, 0.22, 0.08))
  ancs <- sample(ANCS, k)
  if (k == 1) return(hx(ancs))
  p <- as.vector(rmultinom(1, 10, prob = runif(k, 1, 3))); p[p == 0] <- 1
  unlist(lapply(seq_len(k), function(i) rep(hx(ancs[i]), p[i])))
}

NR <- 5; NC <- 9; N <- NR * NC
xs <- rep(1:NC, times = NR); ys <- rep(NR:1, each = NC)
g <- ggplot() + coord_fixed(xlim = c(0.4, NC + 0.6), ylim = c(0.2, NR + 1.4), expand = FALSE) +
  theme_void(base_size = 20)
for (i in seq_len(N)) g <- g + person_layer(xs[i], ys[i], 0.98, mkcomp())

# manual ancestry legend (a row of swatches near the top)
lx <- seq(NC/2 - 3.2, NC/2 + 2.0, length.out = 5) + 0.5
for (j in seq_along(ANCS)) g <- g +
  annotate("rect", xmin = lx[j] - 0.22, xmax = lx[j] - 0.02, ymin = NR + 0.85, ymax = NR + 1.12,
           fill = hx(ANCS[j]), colour = "grey35", linewidth = 0.5) +
  annotate("text", x = lx[j] + 0.04, y = NR + 0.985, label = ANCS[j], hjust = 0, size = 5.5, fontface = "bold")

g <- g + ggtitle("An admixed cohort: each genome a mosaic of local ancestries") +
  theme(plot.title = element_text(size = 22, face = "bold", hjust = 0.5, margin = margin(b = 6)))
save_fig(g, "fig_admixed_cohort", width = 12, height = 8.5)
message("done: fig_admixed_cohort")
