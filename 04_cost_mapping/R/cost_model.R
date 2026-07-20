# =============================================================================
# Cost model: funding <-> coverage mapping and iso-budget geometry
# -----------------------------------------------------------------------------
# Pure functions translating between government funding reductions and
# intervention coverage proportions (paper Eqs. 1-2), plus the iso-budget line
# geometry used by the heatmaps. No parameters are hard-coded here; callers pass
# them in from cost_mapping_params().
# =============================================================================

#' Convert a coverage proportion to the implied dollar funding reduction.
#'
#' Solves the coverage model for new spending and returns the reduction from
#' baseline (floored at zero coverage response below the private substitution
#' level gamma).
proportion_to_dollar_reduction <- function(P,
                                            G_baseline_billions,
                                            C_per_person,
                                            N_population,
                                            alpha_or_beta,
                                            gamma) {
  G_new <- pmax(0, (P - gamma) * alpha_or_beta * C_per_person * N_population)
  G_baseline_billions - G_new
}

#' Map a grid of ART x PrEP funding-reduction fractions to coverage outcomes.
#'
#' Reproduces the cost-mapping block of the original code: builds the funding
#' grid, computes new coverage proportions (clamped to [gamma, 1]), the implied
#' dollar and percentage funding reductions, deduplicates coordinate collisions
#' at the gamma floor, and joins in surrogate incidence to derive IRR and
#' relative change vs. the no-cut baseline.
#'
#' The incidence surface must be keyed to the SAME cost-model funding->coverage
#' mapping the forest plot, Table 1, and the ribbon panel use: a funding-cut
#' fraction f maps to a coverage reduction f * (P_baseline - gamma) / P_baseline
#' (paper Eqs. 1-2). When `gp_incidence_fit`/`predict_composite_fn` are supplied
#' the surrogate is re-queried at those cost-model-mapped native reduction points
#' and year10_mean_incidence/year10_sd_incidence are overwritten before the IRR
#' block, so every heatmap/bar/CSV consumer is cost-model-consistent. If they are
#' NULL the incoming `model_incidence` values are used as-is (the legacy flat
#' 0.75 * funding coverage proxy, kept for back-compatibility).
#'
#' @param model_incidence Data frame with ART_pct_reduction, PrEP_pct_reduction,
#'   year10_mean_incidence, year10_sd_incidence (from surrogate grid).
#' @param p Parameter list from cost_mapping_params().
#' @param gp_incidence_fit Composite incidence surrogate (optional; enables the
#'   cost-model re-query of the incidence surface).
#' @param predict_composite_fn predict_hiv_composite_surrogate helper (optional).
#' @return The augmented, deduplicated funding-scenario data frame.
build_funding_scenarios <- function(model_incidence, p,
                                    gp_incidence_fit = NULL,
                                    predict_composite_fn = NULL) {
  ngrid  <- p$ngrid
  fund_X <- seq(1, 0, length.out = ngrid)
  funding_reduction_pct <- expand.grid(fund_X, fund_X)

  # Re-key the incidence surface to the cost-model funding->coverage mapping so
  # the heatmap agrees with the forest plot, Table 1, and the ribbon panel. A
  # funding cut of fraction f reduces coverage by f * (P_baseline - gamma) / P
  # (the surrogate's native input is a fractional reduction, negative = cut).
  if (!is.null(gp_incidence_fit) && !is.null(predict_composite_fn)) {
    art_cov_per_funding  <- (p$P_ART_baseline  - p$gamma_ART)  / p$P_ART_baseline
    prep_cov_per_funding <- (p$P_PrEP_baseline - p$gamma_PrEP) / p$P_PrEP_baseline

    newX_native <- cbind(
      -funding_reduction_pct[, 1] * art_cov_per_funding,
      -funding_reduction_pct[, 2] * prep_cov_per_funding)

    pred <- predict_composite_fn(
      as.matrix(newX_native), gp_incidence_fit,
      n_samples_per_checkpoint = p$n_samples_per_checkpoint)

    tick <- p$horizon_tick
    model_incidence$year10_mean_incidence <- apply(pred[tick, , ], 1, mean)
    model_incidence$year10_sd_incidence   <- apply(pred[tick, , ], 1, sd)
  }

  model_incidence$P_ART_new <- pmin(1, pmax(
    p$gamma_ART,
    (1 - funding_reduction_pct[, 1]) * p$G_ART_baseline /
      (p$alpha_ART * p$C_ART * p$N_HIV_pos) + p$gamma_ART))

  model_incidence$P_PrEP_new <- pmin(1, pmax(
    p$gamma_PrEP,
    (1 - funding_reduction_pct[, 2]) * p$G_PrEP_baseline /
      (p$beta_PrEP * p$C_PrEP * p$N_PrEP_eligible) + p$gamma_PrEP))

  model_incidence$delta_G_ART <- proportion_to_dollar_reduction(
    P = model_incidence$P_ART_new,
    G_baseline_billions = p$G_ART_baseline,
    C_per_person = p$C_ART,
    N_population = p$N_HIV_pos,
    alpha_or_beta = p$alpha_ART,
    gamma = p$gamma_ART)

  model_incidence$delta_G_PrEP <- proportion_to_dollar_reduction(
    P = model_incidence$P_PrEP_new,
    G_baseline_billions = p$G_PrEP_baseline,
    C_per_person = p$C_PrEP,
    N_population = p$N_PrEP_eligible,
    alpha_or_beta = p$beta_PrEP,
    gamma = p$gamma_PrEP)

  model_incidence$pct_delta_ART_fund  <- model_incidence$delta_G_ART  / p$G_ART_baseline  * 100
  model_incidence$pct_delta_PrEP_fund <- model_incidence$delta_G_PrEP / p$G_PrEP_baseline * 100

  model_incidence$G_ART_fund_new  <- p$G_ART_baseline  - model_incidence$delta_G_ART
  model_incidence$G_PrEP_fund_new <- p$G_PrEP_baseline - model_incidence$delta_G_PrEP

  # Deduplicate coordinates: when coverage hits the gamma floor, multiple grid
  # points collapse to the same (pct_delta_ART_fund, pct_delta_PrEP_fund).
  model_incidence <- model_incidence %>%
    group_by(pct_delta_ART_fund, pct_delta_PrEP_fund) %>%
    slice(1) %>%
    ungroup()

  # Baseline incidence at (0, 0) coverage reduction.
  baseline_idx <- which(
    model_incidence$`ART_pct_reduction` == 0 &
    model_incidence$`PrEP_pct_reduction` == 0)
  baseline_incidence <- model_incidence$year10_mean_incidence[baseline_idx]

  model_incidence$incidence_risk_ratio <-
    model_incidence$year10_mean_incidence / baseline_incidence
  baseline_incidence_risk_ratio <-
    model_incidence$incidence_risk_ratio[baseline_idx]

  model_incidence$relative_incidence_change_pct <-
    ((model_incidence$year10_mean_incidence - baseline_incidence) /
       baseline_incidence) * 100
  model_incidence$relative_incidence_risk_ratio_change_pct <-
    (model_incidence$incidence_risk_ratio - baseline_incidence_risk_ratio) /
      baseline_incidence_risk_ratio * 100

  attr(model_incidence, "baseline_incidence") <- baseline_incidence
  model_incidence
}

#' Iso-budget lines: loci of (ART%, PrEP%) cuts summing to a fixed dollar total.
create_isobudget_lines <- function(total_reduction,
                                    G_ART_baseline,
                                    G_PrEP_baseline) {
  lines_data <- lapply(total_reduction, function(total) {
    pct_ART_seq  <- seq(0, 100, by = 1)
    delta_G_ART  <- (pct_ART_seq / 100) * G_ART_baseline
    delta_G_PrEP <- total - delta_G_ART
    pct_PrEP_seq <- (delta_G_PrEP / G_PrEP_baseline) * 100

    feasible <- pct_ART_seq >= 0 & pct_ART_seq <= 100 &
                pct_PrEP_seq >= 0 & pct_PrEP_seq <= 100

    data.frame(
      total_reduction = total,
      iso_pct_ART  = pct_ART_seq[feasible],
      iso_pct_PrEP = pct_PrEP_seq[feasible],
      delta_G_ART  = delta_G_ART[feasible],
      delta_G_PrEP = delta_G_PrEP[feasible],
      feasible = TRUE
    )
  })
  result <- do.call(rbind, lines_data)
  rownames(result) <- NULL
  result
}
