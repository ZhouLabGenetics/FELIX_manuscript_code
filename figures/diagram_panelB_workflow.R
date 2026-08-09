#!/usr/bin/env Rscript
# local-ancestry->dosage inset, HOM/HET/CCT tests, and outputs incl. the homogeneous joint p.
source("scripts_R/00_theme.R"); source("scripts_R/00_people.R")
suppressMessages(library(grid))
A <- ANC_COLORS

# canvas
g <- ggplot() + coord_fixed(xlim = c(0, 40), ylim = c(0, 12), expand = FALSE) +
  theme_void(base_size = 18)

txt <- function(x, y, lab, size = 5.2, face = "plain", col = "black", h = 0.5)
  annotate("text", x = x, y = y, label = lab, size = size, fontface = face, colour = col, hjust = h)
box <- function(x0, y0, x1, y1, col, fill = NA, lwd = 1.1)
  annotate("rect", xmin = x0, ymin = y0, xmax = x1, ymax = y1, colour = col, fill = fill, linewidth = lwd)
arrow_h <- function(x0, x1, y) annotate("segment", x = x0, xend = x1, y = y, yend = y,
  arrow = arrow(length = unit(0.28, "cm"), type = "closed"), linewidth = 1.4, colour = "#444")

# ---- Stage 1: inclusive cohort (real icons; some admixed = striped) ----
cohort <- list(
  list(A[["EUR"]]), list(A[["AFR"]]), list(c(A[["AFR"]],A[["EUR"]])), list(A[["NatAm"]]),
  list(c(A[["EUR"]],A[["NatAm"]],A[["AFR"]])), list(A[["EUR"]]), list(c(A[["EUR"]],A[["AFR"]])),
  list(c(A[["NatAm"]],A[["EUR"]])), list(A[["EAS"]]), list(A[["SAS"]]),
  list(c(A[["AFR"]],A[["NatAm"]])), list(A[["EUR"]]))
px <- rep(c(1.4, 3.0, 4.6), times = 4); py <- rep(c(9.5, 7.9, 6.3, 4.7), each = 3)
for (i in seq_along(cohort)) g <- g + person_layer(px[i], py[i], 1.5, cohort[[i]][[1]])
g <- g +
  txt(3.0, 11.4, "1  Inclusive cohort", 6, "bold") +
  txt(3.0, 3.3, "every participant,\nincluding admixed", 3.8, "italic", "#555")

# ---- Stage 2: local ancestry -> ancestry-specific dosages (haplotype deconvolution) ----
g <- g + txt(11.5, 11.4, "2  Local-ancestry dosages", 6, "bold")
# two phased haplotypes painted by local ancestry; a focal variant with alt/ref allele on each
seg <- function(x, y, cols) { n <- length(cols); w <- 3.0/n
  lapply(seq_len(n), function(k) annotate("rect", xmin = x + (k-1)*w, xmax = x + k*w,
    ymin = y - 0.26, ymax = y + 0.26, fill = cols[k], colour = "white", linewidth = 0.4)) }
xh <- 8.9; xsnp <- xh + 3.0 * (3.5/4)   # focal SNP in segment 4: maternal AFR, paternal EUR
for (L in seg(xh, 9.5, c(A[["EUR"]],A[["EUR"]],A[["AFR"]],A[["AFR"]]))) g <- g + L   # maternal
for (L in seg(xh, 8.6, c(A[["AFR"]],A[["AFR"]],A[["AFR"]],A[["EUR"]]))) g <- g + L   # paternal
g <- g +
  txt(xh, 10.15, "phased haplotypes", 3.5, "italic", "#555", 0) +
  txt(xh - 0.15, 9.5, "M", 3.3, "bold", "#333", 1) + txt(xh - 0.15, 8.6, "P", 3.3, "bold", "#333", 1) +
  annotate("point", x = xsnp, y = 9.5, size = 2.6, colour = "#111") +                # alt allele (filled)
  annotate("point", x = xsnp, y = 8.6, size = 2.6, shape = 21, fill = "white", colour = "#111") + # ref (open)
  annotate("segment", x = xsnp, xend = xsnp, y = 10.0, yend = 9.75, linewidth = 0.5, colour = "#333") +
  txt(xsnp, 10.3, "SNP", 3.0, "plain", "#333") +
  txt(xh + 3.25, 9.5, "alt", 2.9, "plain", "#111", 0) + txt(xh + 3.25, 8.6, "ref", 2.9, "plain", "#111", 0)
# deconvolution arrows -> per-ancestry dosage counts
g <- g + annotate("segment", x = 10.4, xend = 10.4, y = 8.1, yend = 7.4, linewidth = 1.1, colour = "#444",
                  arrow = arrow(length = unit(0.22, "cm"), type = "closed"))
dsbox <- function(y, a, ds, anc) { g2 <- list(
  box(9.0, y - 0.34, 13.4, y + 0.34, A[[a]], scales::alpha(A[[a]], 0.32)),
  txt(9.25, y, a, 3.5, "bold", h = 0),
  txt(13.15, y, sprintf("DS=%s  ANC=%s", ds, anc), 3.1, "plain", "#222", 1)); g2 }
for (L in dsbox(7.0, "AFR", "1", "1")) g <- g + L     # maternal AFR haplotype carries the alt allele
for (L in dsbox(6.1, "EUR", "0", "1")) g <- g + L     # paternal EUR haplotype carries the ref allele
g <- g +
  txt(11.2, 5.35, "g_AFR = x_M * 1(anc_M=AFR) + x_P * 1(anc_P=AFR)", 3.3, "plain", "#222", 0.5) +
  txt(11.2, 4.7, "per-ancestry dosage (DS) & local-ancestry count (ANC)", 3.4, "italic", "#555", 0.5)

# ---- Stage 3: two score tests + CCT ----
g <- g + txt(22.5, 11.4, "3  Two tests, combined", 6, "bold")
g <- g + box(19.6, 7.6, 23.0, 9.4, TEST_COLORS[["HOM"]], scales::alpha(TEST_COLORS[["HOM"]], 0.30)) +
  txt(21.3, 8.9, "HOM", 4.6, "bold", TEST_COLORS[["HOM"]]) +
  txt(21.3, 8.05, "1 df, summed dosage", 3.3, "plain", "#333")
g <- g + box(23.6, 7.6, 27.0, 9.4, TEST_COLORS[["HET"]], scales::alpha(TEST_COLORS[["HET"]], 0.30)) +
  txt(25.3, 8.9, "HET", 4.6, "bold", "#9a7b16") +
  txt(25.3, 8.05, "K df, per-ancestry", 3.3, "plain", "#333")
g <- g + box(21.0, 5.0, 25.6, 6.4, TEST_COLORS[["CCT"]], scales::alpha(TEST_COLORS[["CCT"]], 0.30), 1.6) +
  txt(23.3, 5.7, "Cauchy combination", 4.2, "bold", TEST_COLORS[["CCT"]])
g <- g +
  annotate("segment", x = 21.3, xend = 22.6, y = 7.6, yend = 6.4, linewidth = 1.1, colour = TEST_COLORS[["HOM"]],
           arrow = arrow(length = unit(0.22,"cm"), type = "closed")) +
  annotate("segment", x = 25.3, xend = 24.0, y = 7.6, yend = 6.4, linewidth = 1.1, colour = "#9a7b16",
           arrow = arrow(length = unit(0.22,"cm"), type = "closed"))

# ---- Stage 4: outputs ----
g <- g + txt(34.5, 11.4, "4  Outputs", 6, "bold")
# per-ancestry effect table
for (k in seq_len(5)) { a <- c("AFR","EAS","EUR","NatAm","SAS")[k]; yk <- 9.6 - (k-1)*0.62
  g <- g + annotate("point", x = 29.4, y = yk, size = 4, colour = A[[a]]) +
    txt(29.9, yk, a, 3.4, "bold", A[[a]], 0) +
    txt(33.0, yk, "beta, SE, p", 3.2, "plain", "#333", 1) }
g <- g + txt(29.4, 10.35, "per-ancestry effects", 3.8, "italic", "#555", 0)
# homogeneous joint p (from summed dosage) + CCT p
g <- g + box(29.2, 4.7, 36.6, 5.6, TEST_COLORS[["HOM"]], scales::alpha(TEST_COLORS[["HOM"]],0.22)) +
  txt(32.9, 5.15, "summed dosage -> joint p (HOM)", 3.4, "plain", "#333") +
  box(29.2, 3.6, 36.6, 4.5, TEST_COLORS[["CCT"]], scales::alpha(TEST_COLORS[["CCT"]],0.22)) +
  txt(32.9, 4.05, "combined p (CCT)", 3.6, "bold", TEST_COLORS[["CCT"]])
# rg / h2 matrix
mx <- 29.2; my <- 1.2; cell <- 0.62; labs3 <- c("AFR","EUR","NatAm")
for (i in 1:3) for (j in 1:3) { fillc <- if (i==j) A[[labs3[i]]] else "#d9d2e0"
  g <- g + annotate("rect", xmin = mx+(j-1)*cell, xmax = mx+j*cell, ymin = my+(3-i)*cell, ymax = my+(4-i)*cell,
                    fill = fillc, colour = "white", linewidth = 0.6) }
g <- g + txt(mx + 1.5*cell, my + 3*cell + 0.35, "rg / h² matrix", 3.8, "italic", "#555")

# ---- inter-stage arrows ----
g <- g + arrow_h(5.9, 8.9, 6.8) + arrow_h(13.8, 19.2, 6.8) + arrow_h(27.4, 28.8, 6.8)

save_fig(g, "panel_b_workflow", width = 20, height = 6.4)
