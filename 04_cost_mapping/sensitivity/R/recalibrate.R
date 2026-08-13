# =============================================================================
# Recalibrate gamma_PrEP / beta_PrEP for a revised (r, P0) pair
# -----------------------------------------------------------------------------
# r  = payer-mix fraction: share of PrEP users assumed privately insured
#      (and therefore immune to a government funding cut).
# P0 = observed baseline PrEP coverage the model must reproduce.
# gamma_PrEP = r * P0 (paper's calibration Eq. 4-5, see
# INFORM-HIV-Analysis/Raff_tmp/gamma_prep_meeting_followup.tex sec. 1).
# beta_PrEP must be recalibrated in the same step so that
# P_PrEP_baseline == P0 is still reproduced (parameters.R line 72):
#   P0 = G_PrEP_baseline / (beta_PrEP * C_PrEP * N_PrEP_eligible) + gamma_PrEP
#   => beta_PrEP = G_PrEP_baseline / ((P0 - gamma_PrEP) * C_PrEP * N_PrEP_eligible)
# =============================================================================

#' @param r Payer-mix fraction (0-1).
#' @param P0 Target baseline PrEP coverage (0-1).
#' @param p_base Parameter list from cost_mapping_params() -- supplies the
#'   fixed quantities (G_PrEP_baseline, C_PrEP, N_PrEP_eligible) that do not
#'   change with r/P0.
#' @return List: r, P0, gamma_PrEP, beta_PrEP, prep_cov_per_funding,
#'   P_PrEP_baseline_check (should equal P0 by construction).
recalibrate_prep <- function(r, P0, p_base) {
  gamma_PrEP <- r * P0
  beta_PrEP <- p_base$G_PrEP_baseline /
    ((P0 - gamma_PrEP) * p_base$C_PrEP * p_base$N_PrEP_eligible)

  prep_cov_per_funding <- (P0 - gamma_PrEP) / P0

  # Sanity check: (P0 - gamma_PrEP)/P0 = (P0 - r*P0)/P0 = 1 - r algebraically,
  # independent of P0. If this ever fails, the derivation above is wrong.
  stopifnot(isTRUE(all.equal(prep_cov_per_funding, 1 - r)))

  P_PrEP_baseline_check <- p_base$G_PrEP_baseline /
    (beta_PrEP * p_base$C_PrEP * p_base$N_PrEP_eligible) + gamma_PrEP
  stopifnot(isTRUE(all.equal(P_PrEP_baseline_check, P0)))

  list(
    r = r,
    P0 = P0,
    gamma_PrEP = gamma_PrEP,
    beta_PrEP = beta_PrEP,
    prep_cov_per_funding = prep_cov_per_funding,
    P_PrEP_baseline_check = P_PrEP_baseline_check
  )
}
