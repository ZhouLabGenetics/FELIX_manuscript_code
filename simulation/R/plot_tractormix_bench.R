## Three manuscript figures from aggregate_tractormix.R outputs:
##   1. Runtime + RSS bar charts (3x3 grid, DNF shown explicitly)
##   2. QQ plots of Tractor-Mix joint P, 3x3 grid (trait x scenario)
##   3. Tractor-Mix vs SAIGE causal-SNP pvalue scatter (3x3 grid)
##
## Usage:  Rscript plot_tractormix_bench.R <common|lowfreq>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: plot_tractormix_bench.R <common|lowfreq>")
mode <- args[1]

source(file.path(Sys.getenv("FELIX_SIM_BASE", unset = stop("Set FELIX_SIM_BASE to the simulation directory before running this script")), "R", "config.R"))
set_mode(mode)
library(data.table)
library(ggplot2)

bench_dir <- file.path(BASE, "Bench", mode)
plot_dir  <- file.path(bench_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

TRAIT_LBL <- c(quant = "Quantitative",
               bin10 = "Binary, 10% prev",
               bin01 = "Binary, 1% prev")
SCEN_LBL  <- c(shared = "Shared",
               afr    = "AFR-only",
               hetero = "Hetero (EUR=0.5B, NAT=B)")

apply_factors <- function(d) {
  d[, trait    := factor(trait,    levels = names(TRAIT_LBL), labels = TRAIT_LBL)]
  d[, scenario := factor(scenario, levels = names(SCEN_LBL),  labels = SCEN_LBL)]
  d
}

## ---- 1. Runtime + RSS bar chart ---------------------------------------
rt <- fread(file.path(bench_dir, "tractormix_runtime.tsv"))
rt <- apply_factors(rt)
rt[, wall_min := wall_s / 60]
rt[, rss_gb   := rss_kb / 1e6]

## Per-condition summary (mean +- SD); DNFs counted separately.
sum_ok <- rt[status == "OK",
             .(mean_wall_min = mean(wall_min, na.rm = TRUE),
               sd_wall_min   = sd(wall_min,   na.rm = TRUE),
               mean_rss_gb   = mean(rss_gb,   na.rm = TRUE),
               sd_rss_gb     = sd(rss_gb,     na.rm = TRUE),
               n_ok = .N),
             by = .(trait, scenario)]
n_total <- rt[, .N, by = .(trait, scenario)]
sum_ok <- merge(sum_ok, n_total, by = c("trait", "scenario"), all.y = TRUE)
sum_ok[, n_dnf := N - ifelse(is.na(n_ok), 0L, n_ok)]
sum_ok[is.na(mean_wall_min),
       label := sprintf("DNF (%d/%d)", n_dnf, N)]
sum_ok[!is.na(mean_wall_min),
       label := sprintf("%.1f ± %.1f min", mean_wall_min, sd_wall_min)]

p_rt <- ggplot(sum_ok, aes(x = scenario, y = mean_wall_min, fill = trait)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.8) +
  geom_errorbar(aes(ymin = pmax(0, mean_wall_min - sd_wall_min),
                    ymax = mean_wall_min + sd_wall_min),
                position = position_dodge(width = 0.85), width = 0.25) +
  geom_text(aes(label = label, y = mean_wall_min),
            position = position_dodge(width = 0.85),
            vjust = -0.5, size = 3) +
  facet_wrap(~ trait, scales = "free_y", ncol = 3) +
  labs(title = sprintf("Tractor-Mix wall time (%s, mean ± SD over seeds)", mode),
       x = NULL, y = "Wall time (min)", fill = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(file.path(plot_dir, "tractormix_runtime.pdf"), p_rt,
       width = 11, height = 5, units = "in")

p_rss <- ggplot(sum_ok, aes(x = scenario, y = mean_rss_gb, fill = trait)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.8) +
  geom_errorbar(aes(ymin = pmax(0, mean_rss_gb - sd_rss_gb),
                    ymax = mean_rss_gb + sd_rss_gb),
                position = position_dodge(width = 0.85), width = 0.25) +
  facet_wrap(~ trait, scales = "free_y", ncol = 3) +
  labs(title = sprintf("Tractor-Mix peak RSS (%s, mean ± SD over seeds)", mode),
       x = NULL, y = "Peak RSS (GB)", fill = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(file.path(plot_dir, "tractormix_rss.pdf"), p_rss,
       width = 11, height = 5, units = "in")

cat("Wrote tractormix_runtime.pdf, tractormix_rss.pdf\n")

## ---- 2. QQ plots ------------------------------------------------------
qq_path <- file.path(bench_dir, "tractormix_qq_pool.tsv")
if (file.exists(qq_path)) {
  qq <- fread(qq_path)
  qq <- qq[!is.na(P) & P > 0 & P <= 1]
  qq <- apply_factors(qq)

  ## Compute expected/observed -log10 per (trait, scen) panel.
  qq_plot <- qq[, {
    n <- .N
    obs <- sort(-log10(P))
    exp <- -log10((n:1) / (n + 1))
    .(exp = exp, obs = obs)
  }, by = .(trait, scenario)]

  lambda <- qq[, .(lambda = round(median(qchisq(P, df = 1, lower.tail = FALSE),
                                          na.rm = TRUE) / qchisq(0.5, df = 1), 3)),
               by = .(trait, scenario)]
  qq_plot <- merge(qq_plot, lambda, by = c("trait", "scenario"))

  p_qq <- ggplot(qq_plot, aes(x = exp, y = obs)) +
    geom_abline(slope = 1, intercept = 0, color = "grey50", linetype = "dashed") +
    geom_point(size = 0.6, alpha = 0.5, color = "#1f77b4") +
    geom_text(data = unique(qq_plot[, .(trait, scenario, lambda)]),
              aes(label = sprintf("λ = %.2f", lambda)),
              x = 0.5, y = Inf, hjust = 0, vjust = 1.5, size = 3.2,
              inherit.aes = FALSE) +
    facet_grid(trait ~ scenario, scales = "free") +
    labs(title = sprintf("Tractor-Mix QQ (joint P, non-causal variants pooled across seeds, %s)", mode),
         x = expression(Expected~-log[10](p)),
         y = expression(Observed~-log[10](p))) +
    theme_bw(base_size = 12)
  ggsave(file.path(plot_dir, "tractormix_qq.pdf"), p_qq,
         width = 12, height = 10, units = "in")
  cat("Wrote tractormix_qq.pdf\n")
} else {
  cat("Skipping QQ: no tractormix_qq_pool.tsv\n")
}

## ---- 3. Tractor-Mix vs SAIGE causal-SNP scatter ----------------------
cs_path <- file.path(bench_dir, "tractormix_causal.tsv")
if (file.exists(cs_path)) {
  cs <- fread(cs_path)
  cs <- apply_factors(cs)
  ## Pair the joint Tractor-Mix P against SAIGE CCT (the closest analog
  ## for a one-number admixed test). Also report HOM and HET if desired.
  cs[, neglogP_tmix  := -log10(tractormix_P)]
  cs[, neglogP_saige := -log10(saige_p_cct)]

  p_sc <- ggplot(cs, aes(x = neglogP_saige, y = neglogP_tmix, color = scenario)) +
    geom_abline(slope = 1, intercept = 0, color = "grey50", linetype = "dashed") +
    geom_point(size = 2) +
    facet_wrap(~ trait, ncol = 3, scales = "free") +
    labs(title = sprintf("Causal SNP: Tractor-Mix vs FELIX (CCT)  --  %s", mode),
         x = expression(FELIX~CCT~-log[10](p)),
         y = expression(Tractor-Mix~joint~-log[10](p)),
         color = "Scenario") +
    theme_bw(base_size = 12) +
    theme(legend.position = "top")
  ggsave(file.path(plot_dir, "tractormix_vs_saige_causal.pdf"), p_sc,
         width = 12, height = 5, units = "in")
  cat("Wrote tractormix_vs_saige_causal.pdf\n")
} else {
  cat("Skipping causal scatter: no tractormix_causal.tsv\n")
}

cat("All plots in ", plot_dir, "\n", sep = "")
