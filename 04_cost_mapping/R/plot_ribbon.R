# =============================================================================
# Incidence-trajectory ribbon plot (main-text figure)
# -----------------------------------------------------------------------------
# Shows how the HIV incidence increase from baseline evolves over the projection
# horizon (years 0-10 since the intervention) under two intervention-use
# reduction scenarios, each at a single reduction level (default 40%):
#   - ART use reduction only (PrEP held at baseline)
#   - PrEP use reduction only (ART held at baseline)
#
# The labelled percentage is the coverage (intervention-use) reduction fed
# directly to the surrogate (identity mapping; the per-unit factors default to
# 1). Passing the cost-model factors instead would reinterpret the labels as
# government funding reductions. For each scenario the surrogate is queried over
# the whole trajectory and the difference from the (0,0) no-reduction baseline
# is summarised into a posterior mean plus 50% and 95% credible ribbons. Common
# random numbers (see R/irr_common_random_numbers.R) share the baseline and
# scenario predictive draws within each checkpoint, so the ribbons reflect
# genuine scenario effects rather than the surrogate's raw predictive noise --
# the same variance-reduction trick used for the forest plot.
#
# The mechanism source is pred_incidence_diff.png (surrogate trajectory draws),
# but here we show only the mean and the ribbons (no individual sample lines),
# and we colour the scenarios with the forest-plot palette. Styling matches the
# forest plot: randplot::theme_rand, PT Sans with a graceful fallback, no plot
# title, white background, vector PDF output.
# =============================================================================

library(ggplot2)
library(data.table)
library(hetGP)     # so predict() dispatches on the stored hetGP models
library(randplot)
library(patchwork)

# resolve_font() is defined in R/plot_forest.R (sourced by run_cost_mapping.R).
# plot_contour_heatmap() is defined in R/plot_heatmap.R (also sourced first).

# ---- Ribbon data for one scenario set ----------------------------------------
# Build the two scenarios at the single reduction level, map the labelled
# reduction to a coverage reduction with the supplied per-unit factors (identity
# for the direct-use framing), query the surrogate trajectory with common random
# numbers, and summarise the (scenario - baseline) incidence difference at each
# trajectory tick into a mean plus 50% and 95% credible bounds.
.ribbon_data <- function(gp_incidence_fit,
                         art_cov_per_unit,
                         prep_cov_per_unit,
                         reduction_level,
                         years,
                         n_samples_per_checkpoint,
                         common_random_numbers,
                         inner_probs,
                         outer_probs) {
  # Scenario labels include the reduction level, e.g. "40% PrEP reduction".
  pct_lab <- paste0(round(reduction_level * 100), "%")
  prep_lab <- paste(pct_lab, "PrEP reduction")
  art_lab  <- paste(pct_lab, "ART reduction")
  scenarios <- rbind(
    data.table(scenario = prep_lab, art_red = 0,               prep_red = reduction_level),
    data.table(scenario = art_lab,  art_red = reduction_level, prep_red = 0)
  )
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
  ribbon_dt[, scenario := factor(scenario, levels = c(prep_lab, art_lab))]
  ribbon_dt[]
}

#' Incidence-trajectory ribbon plot
#'
#' Plots the increase in mean HIV incidence from the no-reduction baseline over
#' the projection horizon for two intervention-use reduction scenarios (ART-only
#' and PrEP-only), each at a single reduction level. Each scenario is drawn as a
#' coloured mean line with a darker 50% credible ribbon and a lighter 95%
#' credible ribbon; individual posterior sample trajectories are not shown.
#' Scenario draws use common random numbers (see
#' compute_scenario_trajectory_draws_crn()).
#'
#' @param gp_incidence_fit Composite incidence surrogate (from surrogate.Rdata).
#' @param art_cov_per_unit Coverage-reduction fraction per unit labelled ART
#'   reduction. Use 1 for a direct intervention-use reduction (the labelled % IS
#'   the coverage reduction), or the cost-model factor
#'   (P_ART_baseline - gamma_ART) / P_ART_baseline for a funding reduction.
#' @param prep_cov_per_unit Same, for PrEP.
#' @param reduction_level Reduction fraction applied in every scenario
#'   (default 0.40 = a 40% reduction).
#' @param years Numeric vector of the trajectory ticks' labels; length must equal
#'   the surrogate's tick count (defaults to 0:10, years since the intervention).
#' @param n_samples_per_checkpoint Surrogate draws per checkpoint.
#' @param inner_probs,outer_probs Lower/upper quantile probabilities for the
#'   inner (50%) and outer (95%) credible ribbons.
#' @param common_random_numbers If TRUE (default), baseline and scenario draws
#'   share predictive randomness within each checkpoint (see the CRN module).
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
                                  reduction_level = 0.40,
                                  years = 0:10,
                                  n_samples_per_checkpoint = 50,
                                  inner_probs = c(0.25, 0.75),
                                  outer_probs = c(0.025, 0.975),
                                  common_random_numbers = TRUE,
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
    reduction_level = reduction_level,
    years = years,
    n_samples_per_checkpoint = n_samples_per_checkpoint,
    common_random_numbers = common_random_numbers,
    inner_probs = inner_probs,
    outer_probs = outer_probs)

  chosen_font <- resolve_font(font)

  # RAND categorical palette keyed to the factor levels (PrEP first, then ART):
  # green = PrEP, blue = ART (matches the forest plot).
  pal <- setNames(c("#45aF84", "#597cbe"), levels(ribbon_dt$scenario))

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
    # Log y-axis (as in the forest plot) so ratios near 1 read clearly.
    scale_y_log10(breaks = c(1, 1.5, 2, 2.5, 3, 4, 5)) +
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
#' @param reduction_level Funding-reduction fraction for the ribbon scenarios
#'   (default 0.40), mapped to coverage via the per-funding factors above.
#' @param heatmap_max Upper bound (percent) for both heatmap axes (default 50).
#' @param rr_breaks Contour/fill breaks for the risk-ratio bands (default
#'   seq(1, 2.5, 0.25), which covers the risk-ratio range within the 0-50%
#'   window; higher bands would be empty and overflow the horizontal legend).
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
                                      reduction_level = 0.40,
                                      heatmap_max = 50,
                                      rr_breaks = seq(1, 2.5, 0.25),
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
  # per-unit coverage factors rather than the identity default.
  ribbon <- plot_incidence_ribbon(
    gp_incidence_fit = gp_incidence_fit,
    art_cov_per_unit = art_cov_per_funding,
    prep_cov_per_unit = prep_cov_per_funding,
    reduction_level = reduction_level,
    n_samples_per_checkpoint = n_samples_per_checkpoint,
    common_random_numbers = common_random_numbers,
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
