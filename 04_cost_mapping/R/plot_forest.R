# =============================================================================
# Two-panel forest plot of incidence risk ratios (main-text Figure 1)
# -----------------------------------------------------------------------------
# Panel A: scenarios defined directly in INTERVENTION USE (coverage) reduction
#          terms. The labelled percentage IS the coverage reduction fed to the
#          surrogate (no cost mapping).
# Panel B: the same labelled percentages interpreted as GOVERNMENT FUNDING
#          reductions, mapped to coverage reductions via the linear
#          funding-to-coverage relationship (paper Eqs. 1-2) before the
#          surrogate query.
#
# Putting the two side by side on a shared incidence-risk-ratio axis makes the
# cost model's dampening explicit: a given percentage funding reduction (Panel
# B) produces a SMALLER coverage reduction -- and hence a smaller incidence
# increase -- than the same percentage use reduction (Panel A), because a
# funding cut of fraction d lowers coverage by only
#   d * (P_baseline - gamma) / P_baseline
# (linear until coverage hits the privately-financed floor gamma). The per-unit
# coverage-reduction factors are passed in from run_cost_mapping.R. PrEP has a
# larger private (gamma) share, so its funding cuts barely move coverage.
#
# In each panel three scenarios are compared at each reduction level:
#   - Vary ART (PrEP = 0)
#   - Vary PrEP (ART = 0)
#   - Vary both equally
#
# Styling: randplot::theme_rand (RAND house style), PT Sans with a graceful
# fallback, no plot title (the caption lives in the paper), no shaded panel
# background. Output is vector (PDF, font-embedded) sized for the main text.
# =============================================================================

library(ggplot2)
library(data.table)
library(hetGP)     # so predict() dispatches on the stored hetGP models
library(randplot)
library(patchwork)

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

# ---- Forest data for one panel ----------------------------------------------
# Build the scenario grid (PrEP-only, ART-only, both) at each reduction level,
# map the labelled reduction to a coverage reduction via the supplied per-unit
# factors (1 for the direct-use panel; the cost-model factors for the funding
# panel), query the surrogate, and summarise the incidence risk ratio with
# common random numbers by default. The IRR is averaged over the trajectory
# horizon per draw (tick = NULL), the forest-plot convention.
.forest_panel_data <- function(gp_incidence_fit,
                               art_cov_per_unit,
                               prep_cov_per_unit,
                               reduction_levels,
                               n_samples_per_checkpoint,
                               common_random_numbers,
                               ci_probs,
                               panel_label) {
  scenarios <- rbind(
    data.table(scenario = "PrEP reduction", art_red = 0,                prep_red = reduction_levels),
    data.table(scenario = "ART reduction",  art_red = reduction_levels, prep_red = 0),
    data.table(scenario = "PrEP and ART",   art_red = reduction_levels, prep_red = reduction_levels)
  )
  scenarios[, reduction := pmax(art_red, prep_red)]  # level being varied

  # Map labelled reduction -> coverage reduction (identity for the use panel).
  scenarios[, art_cov  := art_red  * art_cov_per_unit]
  scenarios[, prep_cov := prep_red * prep_cov_per_unit]

  # Baseline (0,0) plus each (art, prep) coverage reduction, native scale
  # (reductions are negative to the surrogate).
  newX_native <- as.matrix(scenarios[, .(-art_cov, -prep_cov)])
  draws <- compute_scenario_draws_crn(
    newX_native, gp_incidence_fit,
    n_samples_per_checkpoint = n_samples_per_checkpoint,
    common_random_numbers = common_random_numbers, tick = NULL)
  irr_summary <- summarise_draws(draws$irr, probs = ci_probs)

  forest_dt <- copy(scenarios)
  forest_dt[, `:=`(panel = panel_label,
                   mean  = irr_summary$mean,
                   lower = irr_summary$ci_lower,
                   upper = irr_summary$ci_upper)]
  forest_dt[, scenario := factor(scenario, levels = c("PrEP reduction",
                                                      "ART reduction",
                                                      "PrEP and ART"))]
  # crosses_one is TRUE when the lower bound is at or below the null (RR = 1).
  forest_dt[, crosses_one := lower <= 1]
  forest_dt[]
}

# ---- One panel's ggplot ------------------------------------------------------
.forest_panel_plot <- function(forest_dt, title, ylab, x_limits, x_breaks,
                               show_legend, chosen_font) {
  pd <- position_dodge(width = 0.7)
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
    scale_x_log10(breaks = x_breaks, limits = x_limits,
                  expand = expansion(mult = 0)) +
    labs(title = title, x = "Mean incidence risk ratio", y = ylab) +
    theme_rand(font = chosen_font) +
    theme(legend.position = if (show_legend) "top" else "none",
          plot.title = element_text(face = "bold", hjust = 0),
          axis.title.x = element_text(face = "bold"),
          axis.title.y = element_text(face = "bold"),
          # theme_rand blanks the plot background (transparent); force white so
          # the raster PNG preview is not transparent.
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA))
  p
}

#' Two-panel forest plot of incidence risk ratios (Figure 1)
#'
#' Panel A shows scenarios defined as direct intervention-use (coverage)
#' reductions; Panel B shows the same labelled percentages as government funding
#' reductions, mapped to coverage through the cost model. Both share the
#' incidence-risk-ratio axis so the funding-to-coverage dampening is visible.
#' Risk ratios are computed by compute_irr_draws_crn() (see
#' R/irr_common_random_numbers.R) with common random numbers by default.
#'
#' @param gp_incidence_fit Composite incidence surrogate (from surrogate.Rdata).
#' @param art_cov_per_funding Coverage-reduction fraction per unit ART funding
#'   reduction, i.e. (P_ART_baseline - gamma_ART) / P_ART_baseline. Used only in
#'   Panel B; Panel A uses an identity mapping (1).
#' @param prep_cov_per_funding Same, for PrEP.
#' @param reduction_levels Reduction fractions to plot on the y-axis of both
#'   panels (defaults to 10/25/40%).
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
#' @return A list with the combined ggplot object ($plot) and the underlying
#'   forest data table ($data, one row per panel x scenario x reduction level
#'   with mean/lower/upper risk ratios and a crosses_one flag). Returned
#'   invisibly when saved.
plot_funding_forest <- function(gp_incidence_fit,
                                art_cov_per_funding,
                                prep_cov_per_funding,
                                reduction_levels = c(0.10, 0.25, 0.40),
                                n_samples_per_checkpoint = 50,
                                ci_probs = c(0.025, 0.975),
                                common_random_numbers = TRUE,
                                font = "PT Sans",
                                save_path = NULL,
                                width = 10.0,
                                height = 4.8,
                                dpi = 300) {

  set.seed(0)  # surrogate prediction draws are stochastic

  # Panel A: direct use reduction (identity mapping). Panel B: funding reduction
  # mapped to coverage via the cost-model factors.
  data_use <- .forest_panel_data(
    gp_incidence_fit, art_cov_per_unit = 1, prep_cov_per_unit = 1,
    reduction_levels = reduction_levels,
    n_samples_per_checkpoint = n_samples_per_checkpoint,
    common_random_numbers = common_random_numbers, ci_probs = ci_probs,
    panel_label = "A. Intervention use reduction")

  data_fund <- .forest_panel_data(
    gp_incidence_fit, art_cov_per_unit = art_cov_per_funding,
    prep_cov_per_unit = prep_cov_per_funding,
    reduction_levels = reduction_levels,
    n_samples_per_checkpoint = n_samples_per_checkpoint,
    common_random_numbers = common_random_numbers, ci_probs = ci_probs,
    panel_label = "B. Government funding reduction")

  forest_dt <- rbind(data_use, data_fund)

  # Shared x-axis across both panels so the dampening reads straight across.
  x_hi <- max(forest_dt$upper) * 1.04
  x_breaks <- c(1, 1.5, 2, 2.5, 3, 3.5, 4)
  x_breaks <- x_breaks[x_breaks <= x_hi]
  x_limits <- c(1 / 1.01, x_hi)

  chosen_font <- resolve_font(font)

  p_use <- .forest_panel_plot(
    data_use, title = "(A)",
    ylab = "Intervention use reduction (%)",
    x_limits = x_limits, x_breaks = x_breaks,
    show_legend = TRUE, chosen_font = chosen_font)
  p_fund <- .forest_panel_plot(
    data_fund, title = "(B)",
    ylab = "Government funding reduction (%)",
    x_limits = x_limits, x_breaks = x_breaks,
    show_legend = TRUE, chosen_font = chosen_font)

  p <- (p_use + p_fund) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

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
