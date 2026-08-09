#!/usr/bin/env Rscript
# Fig_prs_r2 -- PRS incremental-R2 of all three training strategies (Stable 9)

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr); library(scales)
})

## ---- palette ----
C_FELIX <- "#2c7fb8"   # FELIX ancALL (local-ancestry-aware PRS)
C_META  <- "#d95f02"   # All by All META (trans-ancestry meta-analysis)
C_ANC   <- "#7570b3"   # All by All ancestrally-matched (single-ancestry GWAS)
GREY    <- "#888888"

## ---- theme (matplotlib-like: white bg, no top/right spine, ticks, readable type) ----
BASE <- 13
theme_ms <- function(base = BASE) {
  theme_classic(base_size = base) +
    theme(
      plot.title       = element_text(size = rel(1.0), hjust = 0, margin = margin(b = 4)),
      plot.subtitle    = element_text(size = rel(0.92), hjust = 0),
      axis.title       = element_text(size = rel(1.05)),
      axis.text        = element_text(size = rel(0.92)),
      strip.background = element_blank(),
      strip.text       = element_text(size = rel(1.0), face = "plain"),
      strip.placement  = "outside",
      legend.title     = element_text(size = rel(0.95)),
      legend.text      = element_text(size = rel(0.9)),
      legend.key       = element_blank(),
      legend.background = element_blank(),
      panel.spacing    = unit(1.1, "lines"),
      plot.title.position = "plot"
    )
}

## ---- data: Stable 9 (PRS incremental R2), embedded so this script is self-contained ----
ANC_LV <- c("EUR", "CSA", "AFR", "EAS")
method_lv <- c("FELIX ancALL", "All by All META", "All by All ancestrally-matched")

d <- read_tsv((tsv_data), show_col_types = FALSE,
              col_types = cols(trait = col_character(), ancestry = col_character(),
                                method = col_character(), N = col_double(), r2 = col_double(),
                                ci_lo = col_double(), ci_hi = col_double(), delta_r2 = col_double(),
                                delta_ci_lo = col_double(), delta_ci_hi = col_double(),
                                bootstrap_p = col_double())) %>%
  mutate(ancestry = factor(ancestry, levels = ANC_LV),
         method = factor(method, levels = method_lv))

n_lab <- d %>% group_by(ancestry) %>% summarize(medN = median(N), .groups = "drop") %>%
  mutate(lab = paste0(as.character(ancestry), " (N~", comma(round(medN / 100) * 100), ")"))
n_subtitle <- paste(n_lab$lab, collapse = ",  ")

# trait column order: FELIX ancALL's R2 in EUR (largest, most stable estimate), descending;
# same order used in every ancestry panel.
trait_order <- d %>% filter(method == "FELIX ancALL", ancestry == "EUR") %>%
  arrange(desc(r2)) %>% pull(trait)

# full trait x ancestry x method grid: every trait always reserves 3 dodge slots, even
# where a GWAS could not be estimated at that ancestry's sample size.
full_grid <- expand_grid(trait = trait_order, ancestry = ANC_LV, method = method_lv)
d_full <- full_grid %>%
  left_join(d %>% mutate(trait = as.character(trait), ancestry = as.character(ancestry),
                          method = as.character(method)) %>%
              select(trait, ancestry, method, r2, ci_lo, ci_hi),
            by = c("trait" = "trait", "ancestry" = "ancestry", "method" = "method")) %>%
  # left_join coerces the join keys to character, dropping factor level order -- restore it
  mutate(trait = factor(trait, levels = trait_order),
         ancestry = factor(ancestry, levels = ANC_LV),
         method = factor(method, levels = method_lv),
         estimable = !is.na(r2),
         r2_plot = ifelse(estimable, r2, 0),
         predictor = factor(ifelse(estimable, as.character(method), NA), levels = method_lv))

## group = method (not predictor) so the dodge slot always matches the method's intended
## position, even for the invisible (fill = NA) missing rows -- otherwise ggplot groups by
## fill alone and the empty slot can collapse into the wrong position.
## geom_errorbar and geom_text both dodge over the SAME full d_full (not a filtered
## subset), so all three layers agree on dodge width/offset -- filtering any one layer to a
## subset before dodging desyncs it from the others (misaligned CI whiskers / mispositioned
## asterisk); na.rm suppresses drawing individual rows only after dodge positions are set.
DODGE <- position_dodge(width = 0.75)
pA <- ggplot(d_full, aes(trait, r2_plot, fill = predictor, group = method)) +
  geom_col(position = DODGE, width = 0.68, colour = NA) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), position = DODGE,
                width = 0.2, linewidth = 0.4, colour = "black", na.rm = TRUE) +
  geom_text(aes(y = ifelse(estimable, NA_real_, 0.4),
                label = ifelse(estimable, NA_character_, "*")),
            position = DODGE, size = 6, colour = "grey40", vjust = 0.3, na.rm = TRUE) +
  facet_wrap(~ancestry, ncol = 1) +
  scale_fill_manual(values = setNames(c(C_FELIX, C_META, C_ANC), method_lv),
                     na.value = NA, na.translate = FALSE) +
  labs(x = NULL, y = "incremental R\u00b2 (%)",
       title = "Incremental R\u00b2 of all three PRS training strategies, across 11 traits and 4 ancestries",
       subtitle = n_subtitle) +
  theme_ms() + theme(legend.position = "top", legend.justification = "center",
                      axis.text.x = element_text(angle = 30, hjust = 1))

## ---- save (PDF vector + PNG raster; ragg gives better glyph rendering if installed) ----
ggsave("Fig_prs_r2.pdf", pA, width = 9, height = 12, units = "in")
if (requireNamespace("ragg", quietly = TRUE)) {
  ggsave("Fig_prs_r2.png", pA, width = 9, height = 12, units = "in", dpi = 600, device = ragg::agg_png)
} else {
  ggsave("Fig_prs_r2.png", pA, width = 9, height = 12, units = "in", dpi = 600)
}
message("Wrote Fig_prs_r2.pdf and Fig_prs_r2.png to ", getwd())
