library(ggplot2)
library(data.table)
library(hetGP)
library(viridis)
library(dplyr)

# Load the fitted surrogate produced by stage 03. This restores gp_incidence_fit,
# baseline_incidence_prevalence, and the predict_hiv_* helper functions -- so this
# stage consumes the surrogate rather than re-fitting or redefining it.
load("../../03_intervention_scenario_surrogate/output/surrogate.Rdata")
# variables in surrogate.Rdata: gp_prevalence_fit, gp_incidence_fit,
# baseline_incidence_prevalence, predict_hiv_single_surrogate,
# predict_hiv_composite_surrogate, predict_hiv_diff_composite_surrogate
plots_folder_name <- "../output/"

source('plot01_heatmap.R')
source('plot02_funding_bar.R')
source('funding_scenario_table.R')
 
if (!dir.exists(plots_folder_name)) {
  dir.create(plots_folder_name, recursive = TRUE)
}

#------------------------------------------------
# Parameters
#------------------------------------------------

calculate_cost_2025 = function(inflation_rates, earlier_val) {
  FV = earlier_val * prod(1 + inflation_rates / 100)
  return(FV)
}

# Inflation ratio for 2018 to 2025 values
inflation_rate_19to25 = c(2.2, 2.3, 1.4, 6.0, 5.6, 3.9, 3.3)

#----------
# Unit costs

annual_cumulative_clinical_lab_cost_2018 = 120.40 + 140.36 + 16.92 + 33.13 +
  28.54 + 66.20 + 1052.70
first_year_physician_visit_total_cost_2018 = 671.01
quarterly_physician_visit_cost_2018 = 92.47 * 4
weighted_average_prep_meds_30days_2018 = 349.70

annual_hiv_healthcare_costs_2025 = 43 + 5903 + 9912 + 8833 + 31257 + 14214
annual_nonhiv_healthcare_costs_2025 = 5804
annual_art_costs_2025 = 22288

# Annual Prep costs: for HIV negative group
annual_cost_per_person_onPREP_2018 = annual_cumulative_clinical_lab_cost_2018 +
  first_year_physician_visit_total_cost_2018 + weighted_average_prep_meds_30days_2018 * 12
annual_cost_per_person_onPREP_2025 = calculate_cost_2025(inflation_rate_19to25, earlier_val = annual_cost_per_person_onPREP_2018)

## Annual ART + healthcare + lab costs: for HIV positive group
annual_cumulative_lab_costs_art_2018 = annual_cumulative_clinical_lab_cost_2018 - 120.40 + 120.40 /
  5
annual_cumulative_lab_costs_art_2025 = calculate_cost_2025(inflation_rate_19to25, earlier_val =
                                                             annual_cumulative_lab_costs_art_2018)
annual_cost_per_person_onART_2025 = annual_art_costs_2025 + annual_hiv_healthcare_costs_2025 +
  annual_nonhiv_healthcare_costs_2025 + annual_cumulative_lab_costs_art_2025

calculate_cost_2025(inflation_rate_19to25, annual_cumulative_lab_costs_art_2018)

#---
# Calculate proportions
#---

total_MSM_2021 = 4230000 #2017-2021 MSM pop. aged 13-64 estimate
total_MSM_with_hiv_2022 = 638900  #https://www.hiv.gov/hiv-basics/overview/data-and-trends/statistics
total_MSM_uninfected_2022 = total_MSM_2021 - total_MSM_with_hiv_2022
hiv_prevalence_MSM_2022 = total_MSM_with_hiv_2022 / total_MSM_2021

# Population estimate of all diagnosed HIV positive people
total_num_hiv_infected_2022 = 1218400 #https://www.hiv.gov/hiv-basics/overview/data-and-trends/statistics
prop_MSM_hiv_in_pop = total_MSM_with_hiv_2022 / total_num_hiv_infected_2022

total_num_prep_indication_MSM_2018 = 851240
total_num_prep_indication_pop_2018 = 1216210
prop_MSM_prep_indication_in_pop = total_num_prep_indication_MSM_2018 / total_num_prep_indication_pop_2018

govt_fund_art = 28.71e9 #2022 funding
govt_fund_prep = (2.68 + 1.08 + 3.34) * 10^9 #2022 funding
#2025 discretionary funding:
govt_fund_discretion_art = (2581.04 + 157.25 + 3294.00 + 119.3 + 60.00 + 505.00) * 10^6 #2022 funding
govt_fund_discretion_prep = 1013.71e6 #2022 funding

govt_fund_art_MSM = prop_MSM_hiv_in_pop * (govt_fund_art + govt_fund_discretion_art)
govt_fund_prep_MSM = prop_MSM_prep_indication_in_pop * (govt_fund_prep + govt_fund_discretion_prep)  

#---
# Setup
#---

# Baseline government spending (in billions USD, 2022)
G_ART_baseline <- govt_fund_art_MSM # billion USD for MSM
G_PrEP_baseline <- govt_fund_prep_MSM # billion USD for MSM

# Finalized Cost per person per year (2025 USD)
C_ART <- 50000
C_PrEP <- 21000

# Population sizes
N_HIV_pos <- 638900 # HIV+ MSM
N_MSM_total <- 4230000 # Total MSM
N_HIV_neg <- N_MSM_total - N_HIV_pos # HIV- MSM
N_PrEP_eligible <- 851240 # MSM with PrEP indications

# Model parameters (base case from LaTeX)
alpha_ART <- 1.26 # ART inefficiency factor
beta_PrEP <- 4.0 # PrEP inefficiency factor
gamma_ART <- .15 # Private substitution for ART
gamma_PrEP <- .28 # Private substitution for PrEP



# Calculate baseline coverage proportions
P_ART_baseline <- (G_ART_baseline ) / (alpha_ART * C_ART * N_HIV_pos) + gamma_ART
P_PrEP_baseline <- (G_PrEP_baseline ) / (beta_PrEP * C_PrEP * N_PrEP_eligible) + gamma_PrEP

cat(sprintf("Baseline ART coverage: %.2f%%\n", P_ART_baseline * 100))
cat(sprintf("Baseline PrEP coverage: %.2f%%\n", P_PrEP_baseline * 100))

#------------------------------------------------
# Surrogate simulated mean incidence from LHS
# (predict_hiv_single_surrogate and predict_hiv_composite_surrogate are loaded
#  from stage 03's surrogate.Rdata above.)
#------------------------------------------------

ngrid <- 101
newX <- seq(-0.75, 0, length.out = ngrid)
newX_grid <- as.matrix(expand.grid(newX, newX))

predict_incidence <- predict_hiv_composite_surrogate(newX_grid, gp_incidence_fit, n_samples_per_checkpoint = 50)

mean_incidence <- apply(predict_incidence, c(1, 2), mean)
sd_incidence <- apply(predict_incidence, c(1, 2), sd)

par(mfrow = c(1, 2))
image(newX, newX, matrix(mean_incidence[10, ], nrow = ngrid),
  main = "Predicted total incidence per 100 person year",
  xlab = "ART_pct_reduction", ylab = "PrEP_pct_reduction", col = terrain.colors(100)
)
contour(newX, newX, matrix(mean_incidence[10, ], nrow = ngrid), add = TRUE)

image(newX, newX, matrix(sd_incidence[10, ], nrow = ngrid),
  main = "sd(Total incidence per 100 person year)",
  xlab = "ART_pct_reduction", ylab = "PrEP_pct_reduction", col = terrain.colors(100)
)
contour(newX, newX, matrix(sd_incidence[10, ], nrow = ngrid), add = TRUE)


model_proportions_mean_incidence <- data.frame(cbind(newX_grid, year10_mean_incidence = mean_incidence[10, ], year10_sd_incidence = sd_incidence[10, ]))
colnames(model_proportions_mean_incidence) <- c("ART_pct_reduction", "PrEP_pct_reduction", "year10_mean_incidence", "year10_sd_incidence")
 
#------------------------------------------------
# Surrogate simulated mean prevalence
#------------------------------------------------

# # prevalence
# ngrid <- 101
# newX <- seq(-0.75, 0, length.out = ngrid)
# newX_grid <- as.matrix(expand.grid(newX, newX))

# predict_prevalence <- predict_hiv_composite_surrogate(newX_grid, gp_prevalenece_fit, n_samples_per_checkpoint = 50)

# mean_prevalence <- apply(predict_prevalence, c(1, 2), mean)
# sd_prevalence <- apply(predict_prevalence, c(1, 2), sd)

# par(mfrow = c(1, 2))
# image(newX, newX, matrix(mean_prevalence[10, ], nrow = ngrid),
#   main = "Predicted mean prevalence",
#   xlab = "ART_pct_reduction", ylab = "PrEP_pct_reduction", col = terrain.colors(100)
# )
# contour(newX, newX, matrix(mean_prevalence[10, ], nrow = ngrid), add = TRUE)
# image(newX, newX, matrix(sd_prevalence[10, ], nrow = ngrid),
#   main = "sd(mean prevalence)",
#   xlab = "ART_pct_reduction", ylab = "PrEP_pct_reduction", col = terrain.colors(100)
# )
# contour(newX, newX, matrix(sd_prevalence[10, ], nrow = ngrid), add = TRUE)

# model_proportions_mean_prevalence <- data.frame(cbind(newX_grid, year10_mean_prevalence = mean_prevalence[10, ], year10_sd_prevalence = sd_prevalence[10, ]))
# colnames(model_proportions_mean_prevalence) <- c("ART_pct_reduction", "PrEP_pct_reduction", "year10_mean_prevalence", "year10_sd_prevalence")


#------------------------------------------------
# Posterior Predictive Interval of Incidence
#------------------------------------------------

# Extract year 10 data: ngrid*ngrid grid scenarios × 5000 samples
year10_data <- predict_incidence[10, , ]  # dim: [ngrid*ngrid, 5000]

# Calculate 95% CI for each grid scenario (across the 5000 samples)
ci_lower <- apply(year10_data, 1, quantile, probs = 0.025)   
ci_upper <- apply(year10_data, 1, quantile, probs = 0.975)   
ci_mean  <- apply(year10_data, 1, mean)                       

# Combine into a data frame
incidence_ci_results <- data.frame(
  grid_scenario = 1:ngrid*ngrid,
  mean = ci_mean,
  ci_lower = ci_lower,
  ci_upper = ci_upper,
  ci_width = ci_upper - ci_lower
)

#------------------------------------------------
# Posterior Predictive Interval of Relative Change in Incidence
#------------------------------------------------

# Step 1: Identify and extract baseline scenario samples
# Baseline is at ART_pct_reduction = 0 and PrEP_pct_reduction = 0
baseline_idx <- which(newX_grid[, 1] == 0 & newX_grid[, 2] == 0)
baseline_samples <- predict_incidence[10, baseline_idx, ]  # dim: [5000]

# Step 2: Calculate sample-matched relative changes for all grid scenarios
n_scenarios <- nrow(newX_grid)
n_samples <- dim(predict_incidence)[3]

# Initialize matrix to store relative changes: scenarios × samples
relative_changes_matrix <- matrix(NA, nrow = n_scenarios, ncol = n_samples)

# Step 3: For each scenario, pair samples with baseline and calculate relative change
for (i in 1:n_scenarios) {
  scenario_samples <- predict_incidence[10, i, ]
  # Calculate relative change for each paired iteration
  relative_changes_matrix[i, ] <- ((scenario_samples - baseline_samples) / baseline_samples) * 100
}

# Step 4: Compute confidence intervals across iterations for each scenario
relative_change_ci_lower <- apply(relative_changes_matrix, 1, quantile, probs = 0.025)
relative_change_ci_upper <- apply(relative_changes_matrix, 1, quantile, probs = 0.975)
relative_change_mean <- apply(relative_changes_matrix, 1, mean)
relative_change_median <- apply(relative_changes_matrix, 1, median)

# Combine into a data frame
relative_change_ci_results <- data.frame(
  grid_scenario = 1:n_scenarios,
  mean_relative_change = relative_change_mean,
  median_relative_change = relative_change_median,
  ci_lower = relative_change_ci_lower,
  ci_upper = relative_change_ci_upper,
  ci_width = relative_change_ci_upper - relative_change_ci_lower
)

#------------------------------------------------
# Posterior Predictive Interval of IRR
#------------------------------------------------

# Reuse baseline and dimensions from relative change section above
# baseline_idx, baseline_samples, n_scenarios, n_samples are already defined

# Initialize matrix to store IRR: scenarios × samples
irr_matrix <- matrix(NA, nrow = n_scenarios, ncol = n_samples)

# For each scenario, pair samples with baseline and calculate IRR
for (i in 1:n_scenarios) {
  scenario_samples <- predict_incidence[10, i, ]
  # Calculate IRR for each paired iteration
  irr_matrix[i, ] <- scenario_samples / baseline_samples
}

# Compute confidence intervals across iterations for each scenario
irr_ci_lower <- apply(irr_matrix, 1, quantile, probs = 0.025)
irr_ci_upper <- apply(irr_matrix, 1, quantile, probs = 0.975)
irr_mean <- apply(irr_matrix, 1, mean)
irr_median <- apply(irr_matrix, 1, median)

# Combine into a data frame
irr_ci_results <- data.frame(
  grid_scenario = 1:n_scenarios,
  mean_irr = irr_mean,
  median_irr = irr_median,
  ci_lower = irr_ci_lower,
  ci_upper = irr_ci_upper,
  ci_width = irr_ci_upper - irr_ci_lower
)



#------------------------------------------------
# Cost Mapping
#------------------------------------------------

proportion_to_dollar_reduction <- function(P,
                                          G_baseline_billions,
                                           C_per_person,
                                           N_population,
                                           alpha_or_beta,
                                           gamma) {
  # # New coverage proportion
  # P_new <- P_baseline * (1 + prop_reduction) # prop_reduction is negative

  # Solve for new spending
  G_new <- pmax(0, (P - gamma) * alpha_or_beta * C_per_person * N_population)

  # Dollar reduction
  delta_G_billions <- G_baseline_billions - G_new

  return(delta_G_billions)
}

#---
## Calculate Cost-Informed Population ART, PrEP Proportions 
#---
 
# Calculate what coverage is after proportion reduction
ngrid <- 101
fund_X <- seq(1,0, length.out = ngrid)
funding_reduction_pct <- expand.grid(fund_X, fund_X)

model_proportions_mean_incidence$P_ART_new <- pmin(1, pmax(gamma_ART, (1-funding_reduction_pct[,1])*G_ART_baseline/ (alpha_ART*C_ART*N_HIV_pos) + gamma_ART))

model_proportions_mean_incidence$P_PrEP_new <- pmin(1, pmax(gamma_PrEP, (1-funding_reduction_pct[,2])*G_PrEP_baseline/ (beta_PrEP*C_PrEP*N_PrEP_eligible) + gamma_PrEP))


model_proportions_mean_incidence$delta_G_ART <- proportion_to_dollar_reduction(
  P = model_proportions_mean_incidence$P_ART_new,
  G_baseline_billions = G_ART_baseline,
  C_per_person = C_ART,
  N_population = N_HIV_pos,
  alpha_or_beta = alpha_ART,
  gamma = gamma_ART
)

model_proportions_mean_incidence$delta_G_PrEP <- proportion_to_dollar_reduction(
  P = model_proportions_mean_incidence$P_PrEP_new,
  G_baseline_billions = G_PrEP_baseline,
  C_per_person = C_PrEP,
  N_population = N_PrEP_eligible,
  alpha_or_beta = beta_PrEP,
  gamma = gamma_PrEP
)

# Calculate percentage reductions
model_proportions_mean_incidence$pct_delta_ART_fund <- model_proportions_mean_incidence$delta_G_ART / G_ART_baseline * 100
model_proportions_mean_incidence$pct_delta_PrEP_fund <- model_proportions_mean_incidence$delta_G_PrEP / G_PrEP_baseline * 100

# Calculate new funding levels for each scenario
model_proportions_mean_incidence$G_ART_fund_new <- G_ART_baseline - model_proportions_mean_incidence$delta_G_ART
model_proportions_mean_incidence$G_PrEP_fund_new <- G_PrEP_baseline - model_proportions_mean_incidence$delta_G_PrEP

# Deduplicate coordinates: when P_PrEP_new/P_ART_new hit gamma floor, multiple grid points
model_proportions_mean_incidence <- model_proportions_mean_incidence %>%
  group_by(pct_delta_ART_fund, pct_delta_PrEP_fund) %>%
  slice(1) %>%
  ungroup()
#--- 
## Calculate baselines and incidence risk ratio (IRR)
#---

# Find baseline incidence (0 ART_pct_reduction, 0 PrEP_pct_reduction)
baseline_idx <- which(
  model_proportions_mean_incidence$`ART_pct_reduction` == 0 &
  model_proportions_mean_incidence$`PrEP_pct_reduction` == 0
)
baseline_incidence <- model_proportions_mean_incidence$year10_mean_incidence[baseline_idx]

# Calculate incidence risk ratio (IRR)
model_proportions_mean_incidence$incidence_risk_ratio <-
  model_proportions_mean_incidence$year10_mean_incidence / baseline_incidence

# Baseline IRR
baseline_incidence_risk_ratio <- model_proportions_mean_incidence$incidence_risk_ratio[baseline_idx]

cat(sprintf("Baseline incidence (0%% ART_pct_reduction, 0%% PrEP_pct_reduction): %.4f per 100 PY\n", baseline_incidence))

 
#---
## Calculate relative change compared to baseline (no cuts)
#---

# Calculate relative change in incidence (%)
model_proportions_mean_incidence$relative_incidence_change_pct <-
  ((model_proportions_mean_incidence$year10_mean_incidence - baseline_incidence) /
   baseline_incidence) * 100

# Calculate relatibe change in incidence risk ratio (IRR)
model_proportions_mean_incidence$relative_incidence_risk_ratio_change_pct <-
  (model_proportions_mean_incidence$incidence_risk_ratio - baseline_incidence_risk_ratio) / baseline_incidence_risk_ratio * 100

#------------------------------------------------
# Generate Visualizations
#------------------------------------------------

#---
## Discrete Heatmaps - Absolute Incidence
#---

tile_heatmap_incidence <- plot_tile_matrix_heatmap(
  funding_reduction_scenarios = model_proportions_mean_incidence,
  mean_var = "year10_mean_incidence",
  pct_change_var = "relative_incidence_change_pct",
  var_name = "Mean Incidence"
)
print(tile_heatmap_incidence$p_sim_heatmap)
print(tile_heatmap_incidence$p_sim_relative)
 

#---
## Discrete Heatmaps - Incidence Risk Ratio
#---

tile_heatmap_risk_ratio <- plot_tile_matrix_heatmap(
  funding_reduction_scenarios = model_proportions_mean_incidence,
  mean_var = "incidence_risk_ratio",
  pct_change_var = "relative_incidence_risk_ratio_change_pct",
  var_name = "Mean Incidence\nRisk Ratio"
)
print(tile_heatmap_risk_ratio$p_sim_heatmap)
print(tile_heatmap_risk_ratio$p_sim_relative)
 

#---
## Contour Heatmaps - Absolute Incidence
#---

contour_heatmap_incidence <- plot_contour_heatmap(
  funding_reduction_scenarios = model_proportions_mean_incidence,
  mean_var = "year10_mean_incidence",
  pct_change_var = "relative_incidence_change_pct",
  var_name = "Mean Incidence",
  baseline_value = baseline_incidence
)
print(contour_heatmap_incidence$p_sim_heatmap)
print(contour_heatmap_incidence$p_sim_relative)
 

#---
## Contour Heatmaps - IRR
#---

contour_heatmap_incidence_risk_ratio <- plot_contour_heatmap(
  funding_reduction_scenarios = model_proportions_mean_incidence,
  mean_var = "incidence_risk_ratio",
  pct_change_var = "relative_incidence_risk_ratio_change_pct",
  var_name = "Mean Incidence\nRisk Ratio",
  baseline_value = 1.0
)
print(contour_heatmap_incidence_risk_ratio$p_sim_heatmap)
print(contour_heatmap_incidence_risk_ratio$p_sim_relative)
 


 
#---
## Generate heatmaps with iso budget lines
#---

# Generate percentage-based isobudget lines (10%, 20%, ..., 80%)
total_baseline = G_ART_baseline + G_PrEP_baseline
total_pct_cuts_valid <- seq(0.1, 1, by = 0.1) * total_baseline
isobudget_data <- create_isobudget_lines(
  total_pct_cuts_valid,
  G_ART_baseline,
  G_PrEP_baseline
)

# Add iso-budget lines to absolute incidence heatmaps
iso_heatmap_incidence_pct <- add_iso_budget_lines_to_plot(
  isobudget_data,
  heatmap_plots_list = contour_heatmap_incidence,
  total_baseline = total_baseline
)
print(iso_heatmap_incidence_pct$p_iso_main)
print(iso_heatmap_incidence_pct$p_iso_relative)

# Add iso-budget lines to risk ratio heatmaps
iso_heatmap_risk_ratio_pct <- add_iso_budget_lines_to_plot(
  isobudget_data,
  heatmap_plots_list = contour_heatmap_incidence_risk_ratio,
  total_baseline = total_baseline
)
print(iso_heatmap_risk_ratio_pct$p_iso_main)
print(iso_heatmap_risk_ratio_pct$p_iso_relative)


#---
## Save heatmaps
#---


# Save discrete tile matrix incidence heatmaps
tile_heatmap_plots <- c(
  "discrete_heatmap_incidence_main_year10",
  "discrete_heatmap_incidence_relative_year10"
)

for (p_num in seq_along(tile_heatmap_plots)) {
  ggsave(
    paste0(plots_folder_name, tile_heatmap_plots[p_num], ".png"),
    plot = tile_heatmap_incidence[[p_num]],
    width = 12, height = 10, dpi = 300
  )
}

# Save discrete tile matrix IRR heatmaps 
tile_heatmap_plots <- c(
  "discrete_heatmap_IRR_main_year10",
  "discrete_heatmap_IRR_relative_year10"
)

for (p_num in seq_along(tile_heatmap_plots)) {
  ggsave(
    paste0(plots_folder_name, tile_heatmap_plots[p_num], ".png"),
    plot = tile_heatmap_risk_ratio[[p_num]],
    width = 12, height = 10, dpi = 300
  )
}

# Save absolute incidence contour heatmaps
incidence_plots_pct <- c(
  "heatmap_incidence_main_year10",
  "heatmap_incidence_relative_year10"
)

for (p_num in seq_along(incidence_plots_pct)) {
  ggsave(
    paste0(plots_folder_name, incidence_plots_pct[p_num], ".png"),
    plot = iso_heatmap_incidence_pct[[p_num]],
    width = 10, height = 8, dpi = 300
  )
}

# Save risk ratio contour heatmaps
risk_ratio_plots_pct <- c(
  "heatmap_risk_ratio_main_year10",
  "heatmap_risk_ratio_relative_year10"
)

for (p_num in seq_along(risk_ratio_plots_pct)) {
  ggsave(
    paste0(plots_folder_name, risk_ratio_plots_pct[p_num], ".png"),
    plot = iso_heatmap_risk_ratio_pct[[p_num]],
    width = 10, height = 8, dpi = 300
  )
}


#---
## Bar Plots
#---

# EXAMPLE: Automatically generates and saves both incidence and relative_increase plots
plot_funding_bar(
  key_scenarios=model_proportions_mean_incidence,
  art_funding_reduction_pct = c(0, 10, 0, 25, 0, 40, 0, 40),
  prep_funding_reduction_pct = c(0, 0, 10, 0, 25, 0, 40, 40),
  title = "Selected Policy Scenarios",
  save_path = paste0(plots_folder_name, "bar_plot"),
  y_metric='all',
  incidence_ci_results = incidence_ci_results,
  relative_change_ci_results = relative_change_ci_results,
  irr_ci_results = irr_ci_results
)


#---
## Funding scenario table
#---

# EXAMPLE: //TODO
# define key policy scenarios
key_scenarios <- data.frame(
  Scenario = c(
    'Baseline (No cuts)',
    '10% ART cut only',
    "10% PrEP cut only",
    "25% ART cut only",
    "25% PrEP cut only",
    "40% ART cut only",
    "40% PrEP cut only",
    "40% cut for both"
  ),
  art_funding_reduction_pct = c(
    0, 
    10,
    0,
    25, 
    0,
    40,
    0,
    40
  ),
  prep_funding_reduction_pct = c(
    0,
    0, 
    10,
    0, 
    25,
    0, 
    40,
    40
  )
)
create_funding_scenario_table(
  key_scenarios,
  funding_data = model_proportions_mean_incidence,
  G_ART_baseline_val = G_ART_baseline,
  G_PrEP_baseline_val = G_PrEP_baseline,
  incidence_ci_results = incidence_ci_results,
  relative_change_ci_results = relative_change_ci_results,
  irr_ci_results = irr_ci_results,
  output_dir = plots_folder_name,
  metric_type = "both"
)
