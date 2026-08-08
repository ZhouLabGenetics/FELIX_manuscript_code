## Big-picture power plots from aggregate_power.R outputs.
##
## Two PDFs per mode:
##   test_power_<thresh>.pdf      : HOM / HET / CCT,    9 panels (trait x scenario)
##   ancestry_power_<thresh>.pdf  : anc1 / anc2 / anc3, 9 panels (trait x scenario)
##
## Each panel has its own x-axis (free_x) because per-(trait, scenario) beta
## grids vary widely now -- e.g. quant common runs 0.1-0.5 while bin01 afr
## runs 0.8-4.0. Sharing an x-axis would either compress the small-beta panels
## or stretch the large-beta panels with empty space.
##
## Y-axis is shared at [0, 1] across all panels (it's a proportion). A faint
## dashed reference line at power = 0.8 marks the conventional well-powered
## threshold so it's easy to read off "smallest beta achieving 80% power".
##
## Usage:  Rscript plot_power.R <common|lowfreq> [threshold]
##         Rscript plot_power.R common              # threshold = 5e-8
##         Rscript plot_power.R common 1e-5         # custom threshold

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: plot_power.R <common|lowfreq> [threshold]")
mode   <- args[1]
THRESH <- if (length(args) >= 2) as.numeric(args[2]) else 5e-8

source(file.path(Sys.getenv("FELIX_SIM_BASE", unset = stop("Set FELIX_SIM_BASE to the simulation directory before running this script")), "R", "config.R"))
set_mode(mode)
library(data.table)
library(ggplot2)

POWER_DIR <- file.path(BASE, "Power", mode)
PLOT_DIR  <- file.path(POWER_DIR, "plots")
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

summ <- fread(file.path(POWER_DIR, "power_summary.tsv"))

dat <- summ[abs(threshold - THRESH) < 1e-12]
if (!nrow(dat)) {
  stop(sprintf("No rows for threshold=%g; available: %s",
               THRESH, paste(sort(unique(summ$threshold)), collapse = ", ")))
}

## Per-(mode, scenario, trait) x-axis caps. With facet_wrap(scales="free_x"),
## filtering rows above the cap is enough to clip the panel's x range -- the
## scale auto-fits the remaining data. Panels not listed here use the natural
## data max. If a panel's actual max beta is below the cap, nothing is dropped.
X_MAX_OVERRIDES <- list(
  common = list(
    shared = list(quant = 1.5, bin10 = 1.5, bin01 = 1.5),
    afr    = list(quant = 4.0, bin10 = 2.0, bin01 = 4.0),
    hetero = list(quant = 3.0, bin10 = 3.0, bin01 = 3.0)
  ),
  lowfreq = list(
    shared = list(quant = 0.6)
  )
)
mode_overrides <- X_MAX_OVERRIDES[[mode]]
if (!is.null(mode_overrides)) {
  for (sc in names(mode_overrides)) {
    for (tr in names(mode_overrides[[sc]])) {
      lim <- mode_overrides[[sc]][[tr]]
      n_before <- nrow(dat)
      dat <- dat[!(scenario == sc & trait == tr & beta > lim)]
      n_drop <- n_before - nrow(dat)
      if (n_drop > 0)
        cat(sprintf("Cap %s/%s/%s at beta<=%g: dropped %d rows\n",
                    mode, sc, tr, lim, n_drop))
    }
  }
}

P_TEST <- c(P_hom_admixed_c = "HOM",
            P_het_admixed_c = "HET",
            P_cct_admixed_c = "CCT")
P_ANC  <- c(p.value_c_anc1  = "AFR (anc1)",
            p.value_c_anc2  = "EUR (anc2)",
            p.value_c_anc3  = "NAT (anc3)")
SCEN_LBL  <- c(shared = "Shared (all betas = B)",
               afr    = "AFR-only (others = 0)",
               hetero = "Hetero (EUR = 0.5B, NAT = B)")
TRAIT_LBL <- c(quant = "Quantitative",
               bin10 = "Binary, 10% prev",
               bin01 = "Binary, 1% prev")

dat[, scenario_lbl := factor(scenario, levels = names(SCEN_LBL),
                                       labels = SCEN_LBL)]
dat[, trait_lbl    := factor(trait,    levels = names(TRAIT_LBL),
                                       labels = TRAIT_LBL)]

plot_grid <- function(d, color_map, title_text, fname) {
  d <- copy(d)
  d[, pcol_lbl := factor(pcol, levels = names(color_map), labels = color_map)]

  p <- ggplot(d, aes(x = beta, y = power,
                     color = pcol_lbl, group = pcol_lbl)) +
    geom_hline(yintercept = 0.8, linetype = "dashed",
               color = "grey60", linewidth = 0.3) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    facet_wrap(vars(trait_lbl, scenario_lbl),
               ncol = 3, scales = "free_x",
               labeller = label_wrap_gen(width = 28, multi_line = TRUE)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(title = title_text,
         x = expression(beta),
         y = sprintf("Power  (fraction of seeds with p < %g)", THRESH),
         color = NULL) +
    theme_bw(base_size = 12) +
    theme(legend.position  = "top",
          panel.grid.minor = element_blank(),
          strip.text       = element_text(size = 9, lineheight = 0.9))

  ggsave(file.path(PLOT_DIR, fname), p,
         width = 13, height = 9, units = "in")
}

thresh_tag <- gsub("[^0-9a-zA-Z]", "", format(THRESH, scientific = TRUE))

plot_grid(dat[pcol %in% names(P_TEST)], P_TEST,
          sprintf("Test power curve: HOM / HET / CCT  (%s, alpha = %g)",
                  mode, THRESH),
          sprintf("test_power_%s.pdf", thresh_tag))
plot_grid(dat[pcol %in% names(P_ANC)], P_ANC,
          sprintf("Per-ancestry power curve: AFR / EUR / NAT  (%s, alpha = %g)",
                  mode, THRESH),
          sprintf("ancestry_power_%s.pdf", thresh_tag))

cat(sprintf("Wrote test_power_%s.pdf and ancestry_power_%s.pdf to %s\n",
            thresh_tag, thresh_tag, PLOT_DIR))
