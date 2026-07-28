# =====================================================================
# 00_theme.R  —  shared palette, theme, and save helper for ALL R figures
# Source at the top of every figure script:  source("scripts_R/00_theme.R")
# Run scripts from the repo root (or set STROOT env var to the repo path).
#
# Locked design rules:
#   * NO text on the plot beyond title / axis / legend (no in-plot annotations).
#   * Large production fonts; larger points/lines; R default (sans) typeface.
#   * Low-saturation, distinct palette; red-green colorblind friendly
#     (no red AND green together within one figure's colour set).
#   * Every figure written as BOTH vector PDF and 400-dpi PNG.
# =====================================================================
suppressMessages({
  library(ggplot2); library(data.table); library(scales)
  if (requireNamespace("ragg", quietly = TRUE)) library(ragg)
})

ROOT    <- Sys.getenv("STROOT", unset = getwd())
FIG_PDF <- file.path(ROOT, "figures_R", "pdf")
FIG_PNG <- file.path(ROOT, "figures_R", "png")
dir.create(FIG_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_PNG, recursive = TRUE, showWarnings = FALSE)

# ---- Colour scheme (LOCKED, colorblind-verified) ---------------------------
# Low-saturation (desaturated Paul Tol muted). Verified with Machado-2009
# deuteranopia/protanopia simulation: min pairwise separation ~40 (ancestries),
# ~75 (tests) — well distinguishable to red-green colourblind viewers.
# NO green anywhere in the ancestry set, so no green/red pair (the prior EUR/SAS clash).
ANC_COLORS <- c(
  AFR   = "#9B4E8E",   # muted purple
  EAS   = "#6699CC",   # medium blue (deepened to separate from teal SAS under deutan)
  EUR   = "#D8CB8B",   # sand   (no green anywhere)
  NatAm = "#DB8F75",   # soft terracotta
  SAS   = "#8ACCC9"    # teal   (user-chosen; no grey/black)
)
ANC_COLORS["AMR"] <- ANC_COLORS[["NatAm"]]           # alias
ANC_COLORS["MID"] <- "#B08968"                        # Middle-Eastern (ABA global only), muted brown
# min pairwise deutan/protan separation ~65 (colorblind-verified, Machado-2009).

# Non-ancestry palette — DISTINCT from every ancestry hue, CB-safe (deutan ~59),
# no grey/black/green. Used for tests and method contrasts so they never reuse an
# ancestry colour. 4 slots cover the ladder (meta + HOM + HET + CCT).
TEST_COLORS <- c(HOM = "#EE7733", HET = "#CCBB44", CCT = "#AA3377")
META_COLOR  <- "#0077BB"                              # All-by-All meta (ladder + methods)
# FELIX = suite name; FELIXassoc = association module (colour reused for FELIX / FELIXrg).
METHOD_COLORS <- c("All by All" = "#0077BB", "FELIXassoc" = "#AA3377",
                   "FELIX" = "#AA3377", "FELIXrg" = "#AA3377")

# Locus-status (shared / method-only) — non-ancestry, CB-safe.
STATUS_COLORS <- c(Shared = "#CCBB44", "FELIXassoc only" = "#AA3377", "All by All only" = "#0077BB")
# Extended categorical (for >4 non-ancestry categories, e.g. pie of reasons).
CAT_COLORS <- c("#0077BB","#EE7733","#CCBB44","#AA3377","#332288","#EE3377")

scale_color_anc  <- function(...) scale_color_manual(values = ANC_COLORS, ...)
scale_fill_anc   <- function(...) scale_fill_manual(values = ANC_COLORS, ...)
scale_color_test <- function(...) scale_color_manual(values = TEST_COLORS, ...)
scale_fill_test  <- function(...) scale_fill_manual(values = TEST_COLORS, ...)

# ---- Production theme (large fonts, R default sans, clean) -----------------
theme_st <- function(base_size = 20) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title      = element_text(size = base_size * 1.25, face = "bold", hjust = 0),
      axis.title      = element_text(size = base_size * 1.10),
      axis.text       = element_text(size = base_size * 0.85, colour = "black"),
      legend.title    = element_text(size = base_size * 1.00),
      legend.text     = element_text(size = base_size * 0.90),
      legend.key.size = unit(1.1, "lines"),
      axis.line       = element_line(linewidth = 0.7),
      axis.ticks      = element_line(linewidth = 0.7),
      plot.margin     = margin(10, 14, 10, 10)
    )
}
theme_set(theme_st())

PT    <- 3.4     # default point size
PTBIG <- 5.2     # emphasized point size
LWD   <- 1.1     # default line width

# ---- Save helper: PDF (vector) + PNG (400 dpi) -----------------------------
save_fig <- function(plot, name, width = 8, height = 6, dpi = 400) {
  # base pdf() device (cairo_pdf unavailable without X11/libXrender on this host);
  # keep titles/labels ASCII (hyphens, plotmath) so Helvetica encodes them.
  ggsave(file.path(FIG_PDF, paste0(name, ".pdf")), plot,
         width = width, height = height, device = "pdf", useDingbats = FALSE)
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave(file.path(FIG_PNG, paste0(name, ".png")), plot,
           width = width, height = height, dpi = dpi, device = ragg::agg_png)
  } else {
    ggsave(file.path(FIG_PNG, paste0(name, ".png")), plot,
           width = width, height = height, dpi = dpi)
  }
  message("[R] wrote ", name, ".pdf + .png")
}

# ---- shared data paths + constants -----------------------------------------
SCATTER <- file.path(ROOT, "scatter_output")
MTABLES <- file.path(ROOT, "manuscript_tables")
REPLIC  <- file.path(ROOT, "replication_output")
GW <- 5e-8; GW_LOG <- -log10(GW); GW_CHISQ <- qchisq(GW, 1, lower.tail = FALSE)
ANC_MAP <- data.table(  # anc suffix -> label, ABA prefix
  suf = paste0("anc", 1:5), name = c("AFR","EAS","EUR","NatAm","SAS"),
  aba = c("ABA_AFR","ABA_EAS","ABA_EUR","ABA_AMR","ABA_SAS"))
