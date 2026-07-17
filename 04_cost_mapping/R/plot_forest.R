# =============================================================================
# Forest plot of incidence risk ratios by GOVERNMENT FUNDING reduction
# -----------------------------------------------------------------------------
# Main-text Figure 1. A funding-axis version of the stage-03 coverage forest
# plot: scenarios are defined in program-funding terms (the quantity elicited
# from stakeholders), then mapped to ART/PrEP coverage reductions via the linear
# funding-to-coverage relationship (paper Eqs. 1-2) before being fed to the
# surrogate. Three scenarios are compared at each funding-reduction level:
#   - Vary ART (PrEP = 0)
#   - Vary PrEP (ART = 0)
#   - Vary both equally
#
# A funding cut of fraction d on a program reduces its coverage by
#   d * (P_baseline - gamma) / P_baseline
# (linear until coverage hits the privately-financed floor gamma), so the
# per-unit-funding coverage-reduction factors are passed in from cost_mapping.R
# where the mapping parameters live.
#
# Styling: randplot::theme_rand (RAND house style), PT Sans with a graceful
# fallback, no plot title (the caption lives in the paper), no shaded panel
# background. Output is vector (PDF, font-embedded) sized for the main text.
# =============================================================================

library(ggplot2)
library(data.table)
library(hetGP)     # so predict() dispatches on the stored hetGP models
library(randplot)

# ---- Font selection with graceful fallback ----------------------------------
# Register `preferred` with showtext (so it embeds as vector outlines in the PDF
# without relying on a working cairo/X11 stack) and return its family name. If
# the font or showtext is unavailable, fall back to a generic family so
# rendering never errors on machines without PT Sans.
resolve_font <- function(preferred = "PT Sans", fallback = "sans") {
  available <- tryCatch(systemfonts::system_fonts()$family,
                        error = function(e) character(0))
  if (!preferred %in% available) return(fallback)

  ok <- requireNamespace("showtext", quietly = TRUE) &&
        requireNamespace("sysfonts", quietly = TRUE)
  if (!ok) return(preferred)  # font exists; let the device resolve it directly

  reg <- tryCatch({
    faces <- systemfonts::match_fonts(rep(preferred, 2),
                                      weight = c("normal", "bold"))
    sysfonts::font_add(preferred, regular = faces$path[1], bold = faces$path[2])
    showtext::showtext_auto()
    # render text at the device DPI so glyph sizing matches raster output
    showtext::showtext_opts(dpi = 300)
    TRUE
  }, error = function(e) FALSE)

  if (reg) preferred else fallback
}

#' Forest plot of incidence risk ratios by government funding reduction
#'
#' The incidence risk ratios are computed by compute_irr_draws_crn() (see
#' R/irr_common_random_numbers.R), which supports common random numbers.
#'
#' @param gp_incidence_fit Composite incidence surrogate (from surrogate.Rdata).
#' @param art_cov_per_funding Coverage-reduction fraction per unit ART funding
#'   reduction, i.e. (P_ART_baseline - gamma_ART) / P_ART_baseline.
#' @param prep_cov_per_funding Same, for PrEP.
#' @param funding_levels Funding-reduction fractions to plot (y-axis levels).
#' @param n_samples_per_checkpoint Surrogate draws per checkpoint.
#' @param common_random_numbers If TRUE, the incidence risk ratio is computed
#'   with common random numbers (baseline and scenario share predictive draws
#'   within each checkpoint via compute_irr_draws_crn(); shared noise cancels so
#'   the (0,0) self-ratio is exactly 1). If FALSE, baseline and scenario draws
#'   are independent (classic behaviour, which inflates intervals for tiny
#'   effects and can push lower bounds below 1). Defaults to TRUE.
#' @param font Preferred font family; falls back to "sans" if unavailable.
#' @param save_path Path without extension. Saves "<save_path>.pdf" (vector,
#'   font-embedded) and "<save_path>.png" (raster preview) when provided.
#' @param width,height Figure size in inches.
#'
#' @return A list with the ggplot object ($plot) and the underlying forest data
#'   table ($data, one row per scenario x funding level with mean/lower/upper
#'   risk ratios and a crosses_one flag). Returned invisibly when saved.
plot_funding_forest <- function(gp_incidence_fit,
                                art_cov_per_funding,
                                prep_cov_per_funding,
                                funding_levels = c(0.10, 0.25, 0.40, 0.55, 0.70),
                                n_samples_per_checkpoint = 50,
                                common_random_numbers = TRUE,
                                font = "PT Sans",
                                save_path = NULL,
                                width = 7.0,
                                height = 4.6,
                                dpi = 300) {

  set.seed(0)  # surrogate prediction draws are stochastic

  # ---- Build the scenario grid (funding terms) ------------------------------
  scenarios <- rbind(
    data.table(scenario = "PrEP reduction", art_fund = 0,              prep_fund = funding_levels),
    data.table(scenario = "ART reduction",  art_fund = funding_levels, prep_fund = 0),
    data.table(scenario = "PrEP and ART",   art_fund = funding_levels, prep_fund = funding_levels)
  )
  scenarios[, reduction := pmax(art_fund, prep_fund)]  # funding reduction being varied

  # Map funding reductions -> coverage reductions (paper Eqs. 1-2).
  scenarios[, art_cov  := art_fund  * art_cov_per_funding]
  scenarios[, prep_cov := prep_fund * prep_cov_per_funding]

  # ---- Compute risk ratios --------------------------------------------------
  # All scenario points in one call: baseline (0,0) plus each (art, prep)
  # coverage reduction, native scale (reductions are negative). The IRR is
  # averaged over the trajectory horizon per draw (tick = NULL), matching the
  # forest-plot convention, and computed with or without common random numbers.
  newX_native <- as.matrix(scenarios[, .(-art_cov, -prep_cov)])
  irr_draws <- compute_irr_draws_crn(
    newX_native, gp_incidence_fit,
    n_samples_per_checkpoint = n_samples_per_checkpoint,
    common_random_numbers = common_random_numbers, tick = NULL)
  irr_summary <- summarise_irr_draws(irr_draws, probs = c(0.05, 0.95))

  forest_dt <- copy(scenarios)
  forest_dt[, `:=`(mean  = irr_summary$mean_irr,
                   lower = irr_summary$ci_lower,
                   upper = irr_summary$ci_upper)]

  forest_dt[, scenario := factor(scenario, levels = c("PrEP reduction",
                                                      "ART reduction",
                                                      "PrEP and ART"))]

  # Flag rows whose 90% credible interval crosses (or dips below) 1: these are
  # the whiskers the main text discusses. crosses_one is TRUE whenever the lower
  # bound is at or below the null risk ratio of 1.
  forest_dt[, crosses_one := lower <= 1]

  # ---- Plot -----------------------------------------------------------------
  pd <- position_dodge(width = 0.7)
  x_hi <- max(forest_dt$upper) * 1.04  # small headroom past the widest whisker

  chosen_font <- resolve_font(font)

  p <- ggplot(forest_dt,
              aes(x = mean, y = factor(reduction * 100), color = scenario)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.5) +
    geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y",
                  width = 0.4, linewidth = 0.8, position = pd) +
    geom_point(size = 3, position = pd) +
    # RAND categorical palette: blue = ART, green = PrEP, purple = both
    scale_color_manual(values = c("PrEP reduction" = "#45aF84",
                                  "ART reduction"  = "#597cbe",
                                  "PrEP and ART"   = "#af61a7"),
                       name = NULL) +
    scale_x_log10(breaks = c(1, 1.25, 1.5, 1.75, 2, 2.25),
                  limits = c(1 / 1.01, x_hi),
                  expand = expansion(mult = 0)) +
    labs(x = "Mean incidence risk ratio",
         y = "Government funding reduction (%)") +
    theme_rand(font = chosen_font) +
    theme(legend.position = "top",
          axis.title.x = element_text(face = "bold"),
          axis.title.y = element_text(face = "bold"),
          # theme_rand blanks the plot background (transparent); force white so
          # the raster PNG preview is not transparent.
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA))

  if (!is.null(save_path)) {
    # Vector PDF for the paper. showtext (set up in resolve_font) embeds the
    # font as outlines through the base pdf() device, avoiding the cairo/X11
    # dependency that breaks cairo_pdf on some machines.
    grDevices::pdf(paste0(save_path, ".pdf"), width = width, height = height)
    print(p)
    grDevices::dev.off()
    # PNG preview via ragg (resolves system fonts by family) when available.
    png_device <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"
    ggsave(paste0(save_path, ".png"), p, width = width, height = height,
           units = "in", dpi = dpi, device = png_device)
    return(invisible(list(plot = p, data = forest_dt)))
  }
  list(plot = p, data = forest_dt)
}
