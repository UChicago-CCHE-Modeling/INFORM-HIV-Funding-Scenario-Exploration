# =============================================================================
# Incidence-trajectory ribbon plot (main-text figure)
# -----------------------------------------------------------------------------
# Shows how the mean HIV incidence risk ratio (scenario / no-reduction baseline)
# evolves over the projection horizon (years 0-10 since the intervention) under
# three reduction scenarios that cut ART and PrEP simultaneously by 10%, 20% and
# 40% (default).
#
# The labelled percentage is a reduction fed to the surrogate after mapping to
# coverage via the supplied per-unit factors (identity by default, so the label
# IS the coverage reduction; passing the cost-model factors reinterprets the
# labels as government funding reductions, which is how the combined figure uses
# it). For each scenario the surrogate is queried over the whole trajectory and
# the ratio to the (0,0) baseline is summarised into a posterior mean plus 50%
# and 95% credible ribbons. Common random numbers (see
# R/irr_common_random_numbers.R) share the baseline and scenario predictive
# draws within each checkpoint, so the ribbons reflect genuine scenario effects
# rather than the surrogate's raw predictive noise -- the same variance-
# reduction trick used for the forest plot.
#
# We show only the mean and the ribbons (no individual sample lines). Each
# ribbon is coloured by the heatmap band its terminal-year risk ratio falls in
# (see .band_color_fn), so this panel and the risk-ratio heatmap share a colour
# language. Styling matches the forest plot: randplot::theme_rand, PT Sans with
# a graceful fallback, no plot title, white background, vector PDF output.
# =============================================================================

library(ggplot2)
library(data.table)
library(hetGP)     # so predict() dispatches on the stored hetGP models
library(randplot)
library(patchwork)

# resolve_font() is defined in R/plot_forest.R (sourced by run_cost_mapping.R).
# plot_contour_heatmap() is defined in R/plot_heatmap.R (also sourced first).

# ---- Heatmap band -> colour map ----------------------------------------------
# The heatmap fills its contour bands (breaks e.g. seq(1, 2.5, 0.25)) with the
# discrete viridis "plasma" palette -- exactly what scale_fill_viridis_d() /
# geom_contour_filled() assign. Return a function that maps any risk-ratio value
# to the colour of the band [breaks[i], breaks[i+1]) it falls in, so panel A's
# ribbons can be coloured to match the band their terminal risk ratio sits in on
# panel B. Pulling the colours from viridisLite (rather than hard-coding) keeps
# the two panels in lockstep if the breaks change. Values at/below the first
# break take the first band's colour; values above the last take the last.
.band_color_fn <- function(breaks) {
  n_bands <- length(breaks) - 1L
  band_cols <- viridisLite::viridis(n_bands, option = "plasma")
  function(x) {
    # left.open matches geom_contour_filled's (a, b] bands; all.inside clamps
    # values at/below the first break or above the last into the end bands.
    idx <- findInterval(x, breaks, left.open = TRUE, all.inside = TRUE)
    band_cols[idx]
  }
}

# ---- Ribbon data for one scenario set ----------------------------------------
# Build the two scenarios at the single reduction level, map the labelled
# reduction to a coverage reduction with the supplied per-unit factors (identity
# for the direct-use framing), query the surrogate trajectory with common random
# numbers, and summarise the (scenario - baseline) incidence difference at each
# trajectory tick into a mean plus 50% and 95% credible bounds.
.ribbon_data <- function(gp_incidence_fit,
                         art_cov_per_unit,
                         prep_cov_per_unit,
                         reduction_levels,
                         years,
                         n_samples_per_checkpoint,
                         common_random_numbers,
                         inner_probs,
                         outer_probs) {
  # One scenario per reduction level, each cutting ART and PrEP simultaneously
  # by that level (e.g. "10% reduction" = a 10% cut to both). Labels ordered
  # ascending so the legend and the palette line up.
  labs <- paste0(round(reduction_levels * 100), "% reduction")
  scenarios <- data.table(scenario  = labs,
                          art_red    = reduction_levels,
                          prep_red   = reduction_levels)
  # Map labelled reduction -> coverage reduction (identity for the direct-use
  # framing; the cost-model factors for the funding framing).
  scenarios[, art_cov  := art_red  * art_cov_per_unit]
  scenarios[, prep_cov := prep_red * prep_cov_per_unit]

  # Reductions are negative to the surrogate.
  newX_native <- as.matrix(scenarios[, .(-art_cov, -prep_cov)])
  draws <- compute_scenario_trajectory_draws_crn(
    newX_native, gp_incidence_fit,
    n_samples_per_checkpoint = n_samples_per_checkpoint,
    common_random_numbers = common_random_numbers)

  n_ticks <- nrow(draws$irr[[1]])
  stopifnot(length(years) == n_ticks)

  # Summarise the incidence risk ratio (scenario / baseline) at each tick.
  rows <- lapply(seq_len(nrow(scenarios)), function(k) {
    d <- draws$irr[[k]]  # n_ticks x n_draws
    data.table(
      scenario = scenarios$scenario[k],
      year     = years,
      mean     = rowMeans(d),
      lo95     = apply(d, 1, quantile, probs = outer_probs[1]),
      hi95     = apply(d, 1, quantile, probs = outer_probs[2]),
      lo50     = apply(d, 1, quantile, probs = inner_probs[1]),
      hi50     = apply(d, 1, quantile, probs = inner_probs[2]))
  })
  ribbon_dt <- rbindlist(rows)
  ribbon_dt[, scenario := factor(scenario, levels = labs)]
  ribbon_dt[]
}

#' Incidence-trajectory ribbon plot
#'
#' Plots the mean HIV incidence risk ratio (scenario / no-reduction baseline)
#' over the projection horizon for reduction scenarios that cut ART and PrEP
#' simultaneously by each supplied level. Each scenario is drawn as a coloured
#' mean line with a darker 50% credible ribbon and a lighter 95% credible
#' ribbon; individual posterior sample trajectories are not shown. Scenario
#' draws use common random numbers (see compute_scenario_trajectory_draws_crn()).
#'
#' @param gp_incidence_fit Composite incidence surrogate (from surrogate.Rdata).
#' @param art_cov_per_unit Coverage-reduction fraction per unit labelled ART
#'   reduction. Use 1 for a direct intervention-use reduction (the labelled % IS
#'   the coverage reduction), or the cost-model factor
#'   (P_ART_baseline - gamma_ART) / P_ART_baseline for a funding reduction.
#' @param prep_cov_per_unit Same, for PrEP.
#' @param reduction_levels Reduction fractions, one scenario each, cutting ART
#'   and PrEP simultaneously by that fraction (default c(0.10, 0.20, 0.40)).
#' @param years Numeric vector of the trajectory ticks' labels; length must equal
#'   the surrogate's tick count (defaults to 0:10, years since the intervention).
#' @param n_samples_per_checkpoint Surrogate draws per checkpoint.
#' @param inner_probs,outer_probs Lower/upper quantile probabilities for the
#'   inner (50%) and outer (95%) credible ribbons.
#' @param common_random_numbers If TRUE (default), baseline and scenario draws
#'   share predictive randomness within each checkpoint (see the CRN module).
#' @param band_color_fn Optional function mapping a numeric risk ratio to a fill
#'   colour, used to colour each ribbon by the heatmap band its terminal-year
#'   mean risk ratio falls in (so panels A and B share a colour language). When
#'   NULL, the RAND categorical palette is used instead.
#' @param font Preferred font family; falls back to "sans" if unavailable.
#' @param save_path Path without extension; saves "<save_path>.pdf" (vector,
#'   font-embedded) and "<save_path>.png" (raster preview) when provided.
#' @param width,height Figure size in inches.
#'
#' @return A list with the ggplot object ($plot) and the summarised ribbon data
#'   ($data, one row per scenario x year with mean and 50%/95% bounds). Returned
#'   invisibly when saved.
plot_incidence_ribbon <- function(gp_incidence_fit,
                                  art_cov_per_unit = 1,
                                  prep_cov_per_unit = 1,
                                  reduction_levels = c(0.10, 0.20, 0.40),
                                  years = 0:10,
                                  n_samples_per_checkpoint = 50,
                                  inner_probs = c(0.25, 0.75),
                                  outer_probs = c(0.025, 0.975),
                                  common_random_numbers = TRUE,
                                  band_color_fn = NULL,
                                  font = "PT Sans",
                                  save_path = NULL,
                                  width = 8.0,
                                  height = 5.0,
                                  dpi = 300) {

  set.seed(0)  # surrogate prediction draws are stochastic

  ribbon_dt <- .ribbon_data(
    gp_incidence_fit,
    art_cov_per_unit = art_cov_per_unit,
    prep_cov_per_unit = prep_cov_per_unit,
    reduction_levels = reduction_levels,
    years = years,
    n_samples_per_checkpoint = n_samples_per_checkpoint,
    common_random_numbers = common_random_numbers,
    inner_probs = inner_probs,
    outer_probs = outer_probs)

  chosen_font <- resolve_font(font)

  # Colour each scenario by the heatmap band its terminal-year (last tick) mean
  # risk ratio falls into, so a ribbon's colour matches the contour band its
  # endpoint sits in on panel B. band_color_fn is derived from the same breaks +
  # plasma palette as the heatmap (see .band_color_fn()); if absent, fall back
  # to a RAND categorical palette so the function still works standalone.
  levs <- levels(ribbon_dt$scenario)
  if (!is.null(band_color_fn)) {
    terminal <- ribbon_dt[year == max(year),
                          .(mean = mean[1]), by = scenario]
    setkey(terminal, scenario)
    pal <- setNames(band_color_fn(terminal[levs, mean]), levs)
  } else {
    rand_cat <- c("#45aF84", "#597cbe", "#af61a7")
    pal <- setNames(rand_cat[seq_along(levs)], levs)
  }

  p <- ggplot(ribbon_dt, aes(x = year, group = scenario)) +
    # Null reference: no change from baseline (risk ratio = 1).
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey40",
               linewidth = 0.5) +
    # Lighter 95% ribbon, then darker 50% ribbon on top.
    geom_ribbon(aes(ymin = lo95, ymax = hi95, fill = scenario), alpha = 0.15) +
    geom_ribbon(aes(ymin = lo50, ymax = hi50, fill = scenario), alpha = 0.30) +
    geom_line(aes(y = mean, color = scenario), linewidth = 1.0) +
    scale_color_manual(values = pal, name = NULL) +
    scale_fill_manual(values = pal, name = NULL) +
    scale_x_continuous(breaks = scales::breaks_width(2)) +
    scale_y_continuous() +
    labs(x = "Year",
         y = "Mean incidence risk ratio") +
    theme_rand(font = chosen_font) +
    theme(legend.position = "bottom",
          # Small margin so the long rotated y-axis title is not clipped.
          plot.margin = margin(t = 8, r = 8, b = 4, l = 8),
          # theme_rand blanks the plot background; force white so the raster
          # PNG preview is not transparent.
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA))

  if (!is.null(save_path)) {
    # Vector PDF for the paper; showtext (set up in resolve_font) embeds the
    # font as outlines via the base pdf() device, avoiding the cairo/X11
    # dependency that breaks cairo_pdf on some machines.
    grDevices::pdf(paste0(save_path, ".pdf"), width = width, height = height)
    print(p)
    grDevices::dev.off()
    png_device <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"
    ggsave(paste0(save_path, ".png"), p, width = width, height = height,
           units = "in", dpi = dpi, device = png_device)
    return(invisible(list(plot = p, data = ribbon_dt)))
  }
  list(plot = p, data = ribbon_dt)
}

#' Combined two-panel figure: ribbon trajectory (A) + risk-ratio heatmap (B)
#'
#' Assembles a single main-text figure with the incidence-trajectory ribbon plot
#' as panel A (left) and the government-funding mean-incidence-risk-ratio contour
#' heatmap as panel B (right). Panel B reuses plot_contour_heatmap() (the same
#' code behind heatmap_risk_ratio_main_year10.png) but with both funding-
#' reduction axes restricted to 0-50% and fixed 0.25-wide colour bands from a
#' risk ratio of 1 up to 3.
#'
#' @param gp_incidence_fit Composite incidence surrogate (from surrogate.Rdata).
#' @param grid_scenarios The deduplicated funding->coverage->incidence surface
#'   (model_proportions_mean_incidence) that plot_contour_heatmap() consumes.
#' @param art_cov_per_funding Coverage-reduction fraction per unit ART funding
#'   reduction, (P_ART_baseline - gamma_ART) / P_ART_baseline. Panel A's ribbon
#'   scenarios use this so both panels share the government-funding framing.
#' @param prep_cov_per_funding Same, for PrEP.
#' @param reduction_levels Funding-reduction fractions for the ribbon scenarios
#'   (default c(0.10, 0.20, 0.40)); each cuts ART and PrEP simultaneously by
#'   that fraction, mapped to coverage via the per-funding factors above.
#' @param heatmap_max Upper bound (percent) for both heatmap axes (default 50).
#' @param rr_breaks Contour/fill breaks for the risk-ratio bands (default
#'   seq(1, 2.25, 0.25)). The 0-50% window's IRR tops out just above 2, so this
#'   gives exactly five fully-populated bands. Every band must contain data:
#'   panel B renders with scale_fill_viridis_d(drop = TRUE), so an empty top band
#'   would be dropped and the remaining bands re-interpolated to a DIFFERENT set
#'   of plasma colours than .band_color_fn() (which colours the ribbons) assigns,
#'   breaking the shared colour language. Keep the top break just above the
#'   window's max IRR.
#' @param n_samples_per_checkpoint Surrogate draws per checkpoint.
#' @param common_random_numbers Passed to the ribbon panel (default TRUE).
#' @param font Preferred font family; falls back to "sans" if unavailable.
#' @param save_path Path without extension; saves "<save_path>.pdf" and
#'   "<save_path>.png" when provided.
#' @param width,height Figure size in inches.
#'
#' @return A list with the combined ggplot ($plot) and the ribbon data
#'   ($ribbon_data). Returned invisibly when saved.
plot_ribbon_heatmap_panel <- function(gp_incidence_fit,
                                      grid_scenarios,
                                      art_cov_per_funding,
                                      prep_cov_per_funding,
                                      reduction_levels = c(0.10, 0.20, 0.40),
                                      heatmap_max = 50,
                                      rr_breaks = seq(1, 2.25, 0.25),
                                      n_samples_per_checkpoint = 50,
                                      common_random_numbers = TRUE,
                                      font = "PT Sans",
                                      save_path = NULL,
                                      width = 14.0,
                                      height = 5.5,
                                      dpi = 300) {

  chosen_font <- resolve_font(font)

  # Single explicit theme applied IDENTICALLY to both panels so their fonts,
  # sizes and weights match exactly. theme_rand() is an incomplete theme, so on
  # its own panel B keeps plot_contour_heatmap()'s larger (size 14, bold) axis
  # text -- this override pins both panels to the same non-bold sizes and trims
  # theme_rand's large bottom legend.box.margin (b = 24) that padded the figure.
  base_size <- 12
  shared_theme <- theme(
    text          = element_text(family = chosen_font, size = base_size),
    axis.title    = element_text(family = chosen_font, size = base_size, face = "plain"),
    axis.text     = element_text(family = chosen_font, size = base_size, face = "plain"),
    legend.text   = element_text(family = chosen_font, size = base_size - 1),
    legend.title  = element_text(family = chosen_font, size = base_size, face = "plain"),
    legend.position   = "bottom",
    legend.box.margin = margin(t = 4, r = 0, b = 0, l = 0),
    plot.margin       = margin(t = 8, r = 8, b = 4, l = 8),
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA))

  # ---- Panel A: incidence-trajectory ribbon (build without saving) ----------
  # Government-funding framing (to pair with panel B): the labelled reduction is
  # a funding reduction mapped to coverage via the cost model, so pass the
  # per-unit coverage factors rather than the identity default. Colour each
  # ribbon by the heatmap band its terminal risk ratio falls in (same rr_breaks
  # + plasma palette as panel B) so the two panels share a colour language.
  ribbon <- plot_incidence_ribbon(
    gp_incidence_fit = gp_incidence_fit,
    art_cov_per_unit = art_cov_per_funding,
    prep_cov_per_unit = prep_cov_per_funding,
    reduction_levels = reduction_levels,
    n_samples_per_checkpoint = n_samples_per_checkpoint,
    common_random_numbers = common_random_numbers,
    band_color_fn = .band_color_fn(rr_breaks),
    font = font,
    save_path = NULL)
  p_a <- ribbon$plot + shared_theme

  # ---- Panel B: risk-ratio contour heatmap, 0-50% axes, 0.5-wide bands ------
  # Restrict to the plotted window (small epsilon so a floating-point boundary
  # row at exactly heatmap_max is not dropped) and let the data hull set the
  # axis extent -- passing hard scale limits would clip that boundary row and
  # leave a white strip. Fixed breaks give the requested 0.5-wide bands.
  eps <- 1e-6
  grid_window <- grid_scenarios[
    grid_scenarios$pct_delta_ART_fund  <= heatmap_max + eps &
    grid_scenarios$pct_delta_PrEP_fund <= heatmap_max + eps, ]

  heatmap <- plot_contour_heatmap(
    funding_reduction_scenarios = grid_window,
    mean_var = "incidence_risk_ratio",
    pct_change_var = "relative_incidence_risk_ratio_change_pct",
    var_name = "Mean incidence risk ratio at year 10",
    baseline_value = 1.0,
    breaks = rr_breaks)
  # No coord_fixed: let panel B stretch to fill its column so both panels
  # occupy the same horizontal space. Apply theme_rand first (RAND base), then
  # the shared override so panel B matches panel A's font/size/weight exactly.
  # Sentence-case axis titles to match panel A. Lay the discrete risk-ratio
  # bands out in a single horizontal row along the bottom (one-line title above)
  # so the legend consumes width rather than vertical plot height; centre the
  # legend so the long title is not clipped at the figure edge.
  p_b <- heatmap$p_sim_heatmap +
    labs(x = "ART govt. funding reduction (%)",
         y = "PrEP govt. funding reduction (%)") +
    # Re-apply the fill scale with drop = TRUE so only the bands actually
    # present in the 0-50% window get legend keys. plot_contour_heatmap() uses
    # drop = FALSE (to show every band elsewhere), which here leaves a phantom
    # trailing key with no swatch that overflows the row.
    scale_fill_viridis_d(option = "plasma", drop = TRUE,
                         name = "Mean incidence risk ratio at year 10") +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE,
                               title.position = "top", title.hjust = 0.5)) +
    theme_rand(font = chosen_font) +
    theme(panel.grid = element_blank(),
          legend.justification = "center") +
    shared_theme +
    # Tighten the swatches, text and inter-key spacing so all six risk-ratio
    # bands fit in a single centred row under panel B (its column is only half
    # the figure width, so the default key size overflows the right edge).
    theme(legend.text     = element_text(family = chosen_font, size = base_size - 2),
          legend.key.width = unit(0.45, "cm"),
          legend.key.height = unit(0.45, "cm"),
          legend.key.spacing.x = unit(0.15, "cm"))

  # ---- Assemble ------------------------------------------------------------
  p <- (p_a | p_b) +
    plot_layout(widths = c(1, 1)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold"))

  if (!is.null(save_path)) {
    grDevices::pdf(paste0(save_path, ".pdf"), width = width, height = height)
    print(p)
    grDevices::dev.off()
    png_device <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"
    ggsave(paste0(save_path, ".png"), p, width = width, height = height,
           units = "in", dpi = dpi, device = png_device)
    return(invisible(list(plot = p, ribbon_data = ribbon$data)))
  }
  list(plot = p, ribbon_data = ribbon$data)
}
