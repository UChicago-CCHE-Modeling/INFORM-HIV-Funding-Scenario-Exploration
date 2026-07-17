# =============================================================================
# Surrogate outcomes: grid prediction and posterior predictive intervals
# -----------------------------------------------------------------------------
# Queries the stage-03 composite incidence surrogate over the ART x PrEP
# coverage-reduction grid and derives, at the reporting horizon, the mean
# incidence surface and the posterior predictive intervals for incidence,
# relative change, and the incidence risk ratio (IRR). The computations match
# the original code; the constants they used are now passed in as arguments.
# =============================================================================

#' Predict the incidence surface over the coverage-reduction grid.
#'
#' @param p Parameter list from cost_mapping_params().
#' @param gp_incidence_fit Composite incidence surrogate.
#' @param predict_composite_fn predict_hiv_composite_surrogate helper.
#' @return List with newX, newX_grid, predict_incidence (full array),
#'   mean_incidence, sd_incidence, and the horizon-tick model_incidence frame.
predict_incidence_grid <- function(p, gp_incidence_fit, predict_composite_fn) {
  ngrid <- p$ngrid
  newX  <- seq(-0.75, 0, length.out = ngrid)
  newX_grid <- as.matrix(expand.grid(newX, newX))

  predict_incidence <- predict_composite_fn(
    newX_grid, gp_incidence_fit,
    n_samples_per_checkpoint = p$n_samples_per_checkpoint)

  mean_incidence <- apply(predict_incidence, c(1, 2), mean)
  sd_incidence   <- apply(predict_incidence, c(1, 2), sd)

  tick <- p$horizon_tick
  model_incidence <- data.frame(cbind(
    newX_grid,
    year10_mean_incidence = mean_incidence[tick, ],
    year10_sd_incidence   = sd_incidence[tick, ]))
  colnames(model_incidence) <- c("ART_pct_reduction", "PrEP_pct_reduction",
                                 "year10_mean_incidence", "year10_sd_incidence")

  list(newX = newX, newX_grid = newX_grid,
       predict_incidence = predict_incidence,
       mean_incidence = mean_incidence, sd_incidence = sd_incidence,
       model_incidence = model_incidence)
}

#' Posterior predictive intervals for incidence, relative change, and IRR.
#'
#' @param grid Output of predict_incidence_grid().
#' @param p Parameter list from cost_mapping_params().
#' @return List with incidence_ci_results, relative_change_ci_results,
#'   irr_ci_results (each a data frame indexed by original grid scenario).
compute_ppi <- function(grid, p) {
  predict_incidence <- grid$predict_incidence
  newX_grid <- grid$newX_grid
  ngrid <- p$ngrid
  tick  <- p$horizon_tick
  probs <- p$ci_probs

  # --- Incidence PPI ---
  year10_data <- predict_incidence[tick, , ]
  ci_lower <- apply(year10_data, 1, quantile, probs = probs[1])
  ci_upper <- apply(year10_data, 1, quantile, probs = probs[2])
  ci_mean  <- apply(year10_data, 1, mean)

  incidence_ci_results <- data.frame(
    grid_scenario = 1:ngrid * ngrid,
    mean = ci_mean,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    ci_width = ci_upper - ci_lower)

  # --- Relative change PPI (sample-matched vs. baseline) ---
  baseline_idx <- which(newX_grid[, 1] == 0 & newX_grid[, 2] == 0)
  baseline_samples <- predict_incidence[tick, baseline_idx, ]
  n_scenarios <- nrow(newX_grid)
  n_samples   <- dim(predict_incidence)[3]

  relative_changes_matrix <- matrix(NA, nrow = n_scenarios, ncol = n_samples)
  for (i in 1:n_scenarios) {
    scenario_samples <- predict_incidence[tick, i, ]
    relative_changes_matrix[i, ] <-
      ((scenario_samples - baseline_samples) / baseline_samples) * 100
  }

  relative_change_ci_results <- data.frame(
    grid_scenario = 1:n_scenarios,
    mean_relative_change   = apply(relative_changes_matrix, 1, mean),
    median_relative_change = apply(relative_changes_matrix, 1, median),
    ci_lower = apply(relative_changes_matrix, 1, quantile, probs = probs[1]),
    ci_upper = apply(relative_changes_matrix, 1, quantile, probs = probs[2]),
    ci_width = apply(relative_changes_matrix, 1, quantile, probs = probs[2]) -
               apply(relative_changes_matrix, 1, quantile, probs = probs[1]))

  # --- IRR PPI (sample-matched vs. baseline) ---
  irr_matrix <- matrix(NA, nrow = n_scenarios, ncol = n_samples)
  for (i in 1:n_scenarios) {
    scenario_samples <- predict_incidence[tick, i, ]
    irr_matrix[i, ] <- scenario_samples / baseline_samples
  }

  irr_ci_results <- data.frame(
    grid_scenario = 1:n_scenarios,
    mean_irr   = apply(irr_matrix, 1, mean),
    median_irr = apply(irr_matrix, 1, median),
    ci_lower = apply(irr_matrix, 1, quantile, probs = probs[1]),
    ci_upper = apply(irr_matrix, 1, quantile, probs = probs[2]),
    ci_width = apply(irr_matrix, 1, quantile, probs = probs[2]) -
               apply(irr_matrix, 1, quantile, probs = probs[1]))

  list(incidence_ci_results = incidence_ci_results,
       relative_change_ci_results = relative_change_ci_results,
       irr_ci_results = irr_ci_results)
}
