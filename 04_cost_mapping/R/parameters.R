# =============================================================================
# Cost-mapping parameters
# -----------------------------------------------------------------------------
# Raw input values live in ../cost_params.yml (one documented entry per value,
# with provenance). This file reads that YAML and computes the derived
# quantities (population subgroups, MSM funding shares, baseline coverage). No
# raw parameter is hard-coded here; change a value in the YAML, not in R.
# =============================================================================

library(yaml)

# Default location of the parameter file (relative to script/ working dir).
COST_PARAMS_YML <- "../cost_params.yml"

#' Assemble all cost-mapping parameters from the YAML file.
#'
#' @param yml Path to cost_params.yml.
#' @return A named list of every parameter the pipeline needs, with derived
#'   quantities (baseline MSM funding, baseline coverage) already computed.
cost_mapping_params <- function(yml = COST_PARAMS_YML) {
  raw <- yaml::read_yaml(yml)

  # Left-fold sum: matches the original `a + b + c` arithmetic order exactly
  # (guards against any last-bit drift vs. sum()'s accumulator).
  add <- function(x) Reduce(`+`, x)

  # ---- Raw scalars --------------------------------------------------------
   total_MSM_2021              <- raw$population$total_MSM_2021
  total_MSM_with_hiv_2022     <- raw$population$total_MSM_with_hiv_2022
  total_num_hiv_infected_2022 <- raw$population$total_num_hiv_infected_2022
  total_num_prep_indication_MSM_2018 <- raw$population$total_num_prep_indication_MSM_2018
  total_num_prep_indication_pop_2018 <- raw$population$total_num_prep_indication_pop_2018

  # Government funding. Scalar totals are written in scientific notation in the
  # YAML and may parse as character; coerce. Component-based totals reproduce the
  # original `(a + b + ...) * scale` arithmetic.
  govt_fund_art  <- as.numeric(raw$government_funding$govt_fund_art)
  govt_fund_prep <- add(raw$government_funding$govt_fund_prep$components_billion) * 1e9
  govt_fund_discretion_art <-
    add(raw$government_funding$govt_fund_discretion_art$components_million) * 1e6
  govt_fund_discretion_prep <- as.numeric(raw$government_funding$govt_fund_discretion_prep)

  C_ART  <- raw$per_person_costs$C_ART
  C_PrEP <- raw$per_person_costs$C_PrEP

  alpha_ART  <- raw$coverage_model$alpha_ART
  beta_PrEP  <- raw$coverage_model$beta_PrEP
  gamma_ART  <- raw$coverage_model$gamma_ART
  gamma_PrEP <- raw$coverage_model$gamma_PrEP

  ngrid                    <- raw$grid$ngrid
  n_samples_per_checkpoint <- raw$grid$n_samples_per_checkpoint
  horizon_tick             <- raw$grid$horizon_tick
  ci_probs                 <- unlist(raw$grid$ci_probs)

  # ---- Derived: population subgroups --------------------------------------
  N_HIV_pos       <- total_MSM_with_hiv_2022
  N_PrEP_eligible <- total_num_prep_indication_MSM_2018

  # ---- Derived: MSM shares and MSM-attributed funding ---------------------
  prop_MSM_hiv_in_pop <- total_MSM_with_hiv_2022 / total_num_hiv_infected_2022
  prop_MSM_prep_indication_in_pop <-
    total_num_prep_indication_MSM_2018 / total_num_prep_indication_pop_2018

  G_ART_baseline  <- prop_MSM_hiv_in_pop *
    (govt_fund_art + govt_fund_discretion_art)
  G_PrEP_baseline <- prop_MSM_prep_indication_in_pop *
    (govt_fund_prep + govt_fund_discretion_prep)

  # ---- Derived: baseline coverage proportions -----------------------------
  P_ART_baseline  <- G_ART_baseline  / (alpha_ART * C_ART  * N_HIV_pos)       + gamma_ART
  P_PrEP_baseline <- G_PrEP_baseline / (beta_PrEP * C_PrEP * N_PrEP_eligible) + gamma_PrEP

  list(
    total_MSM_with_hiv_2022 = total_MSM_with_hiv_2022,
    total_num_hiv_infected_2022 = total_num_hiv_infected_2022,
    total_num_prep_indication_MSM_2018 = total_num_prep_indication_MSM_2018,
    total_num_prep_indication_pop_2018 = total_num_prep_indication_pop_2018,
    govt_fund_art = govt_fund_art,
    govt_fund_prep = govt_fund_prep,
    govt_fund_discretion_art = govt_fund_discretion_art,
    govt_fund_discretion_prep = govt_fund_discretion_prep,
    C_ART = C_ART,
    C_PrEP = C_PrEP,
    alpha_ART = alpha_ART,
    beta_PrEP = beta_PrEP,
    gamma_ART = gamma_ART,
    gamma_PrEP = gamma_PrEP,
    ngrid = ngrid,
    n_samples_per_checkpoint = n_samples_per_checkpoint,
    horizon_tick = horizon_tick,
    ci_probs = ci_probs,
    N_HIV_pos = N_HIV_pos,
    N_PrEP_eligible = N_PrEP_eligible,
    prop_MSM_hiv_in_pop = prop_MSM_hiv_in_pop,
    prop_MSM_prep_indication_in_pop = prop_MSM_prep_indication_in_pop,
    G_ART_baseline = G_ART_baseline,
    G_PrEP_baseline = G_PrEP_baseline,
    P_ART_baseline = P_ART_baseline,
    P_PrEP_baseline = P_PrEP_baseline
  )
}
