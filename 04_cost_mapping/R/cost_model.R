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
#' @param model_incidence Data frame with ART_pct_reduction, PrEP_pct_reduction,
#'   year10_mean_incidence, year10_sd_incidence (from surrogate grid).
#' @param p Parameter list from cost_mapping_params().
#' @return The augmented, deduplicated funding-scenario data frame.
build_funding_scenarios <- function(model_incidence, p) {
  ngrid  <- p$ngrid
  fund_X <- seq(1, 0, length.out = ngrid)
  funding_reduction_pct <- expand.grid(fund_X, fund_X)

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
