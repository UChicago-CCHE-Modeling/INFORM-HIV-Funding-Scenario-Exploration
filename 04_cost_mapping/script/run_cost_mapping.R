# =============================================================================
# Cost mapping
# -----------------------------------------------------------------------------
# Runs the cost-mapping pipeline: load parameters and the stage-03 surrogate,
# predict the incidence surface, compute posterior predictive intervals, map
# funding reductions to coverage outcomes, and write all figures and the
# scenario table. Parameters are defined in R/parameters.R; the reusable
# functions live in R/*.R. This script calls them in order.
#
# Run with script/ as the working directory:
#   cd 04_cost_mapping/script && Rscript run_cost_mapping.R
# =============================================================================

library(ggplot2)
library(data.table)
library(hetGP)     # so predict() dispatches on the stored surrogate models
library(viridis)
library(dplyr)
library(tidyr)

# ---- Modules --------------------------------------------------------------
R_DIR <- "../R"
source(file.path(R_DIR, "parameters.R"))
source(file.path(R_DIR, "cost_model.R"))
source(file.path(R_DIR, "outcomes.R"))
source(file.path(R_DIR, "irr_common_random_numbers.R"))
source(file.path(R_DIR, "plot_heatmap.R"))
source(file.path(R_DIR, "plot_bar.R"))
source(file.path(R_DIR, "plot_forest.R"))
source(file.path(R_DIR, "table_scenarios.R"))

# Common random numbers for the incidence risk ratio: when TRUE, the baseline
# and scenario surrogate draws share randomness within each checkpoint so shared
# predictive noise cancels in the ratio (the (0,0) self-ratio is then exactly 1
# and small-effect intervals no longer spuriously cross below 1). Used by BOTH
# the forest figure and the scenario table so the two stay consistent. Flip to
# FALSE to reproduce the classic independent-draw behaviour for comparison.
USE_COMMON_RANDOM_NUMBERS <- TRUE

# ---- Inputs / outputs -----------------------------------------------------
# Stage-03 surrogate bundle restores gp_incidence_fit, baseline_incidence_prevalence,
# and the predict_hiv_* helpers.
load("../../03_intervention_scenario_surrogate/output/surrogate.Rdata")
plots_folder_name <- "../output/"
if (!dir.exists(plots_folder_name)) dir.create(plots_folder_name, recursive = TRUE)

# ---- Parameters -----------------------------------------------------------
p <- cost_mapping_params()
cat(sprintf("Baseline ART coverage: %.2f%%\n",  p$P_ART_baseline  * 100))
cat(sprintf("Baseline PrEP coverage: %.2f%%\n", p$P_PrEP_baseline * 100))

# ---- Surrogate incidence surface + posterior predictive intervals ---------
# The grid prediction draws from the surrogate's predictive distribution, so pin
# a seed here for reproducible heatmaps, PPIs, and the tracked results CSVs.
# (The forest plot and table CRN paths set their own seed downstream.)
set.seed(12345)
grid <- predict_incidence_grid(p, gp_incidence_fit, predict_hiv_composite_surrogate)
ppi  <- compute_ppi(grid, p)
incidence_ci_results       <- ppi$incidence_ci_results
relative_change_ci_results <- ppi$relative_change_ci_results
irr_ci_results             <- ppi$irr_ci_results

# ---- Funding -> coverage -> incidence outcomes ----------------------------
model_proportions_mean_incidence <- build_funding_scenarios(grid$model_incidence, p)
baseline_incidence <- attr(model_proportions_mean_incidence, "baseline_incidence")
cat(sprintf("Baseline incidence (0%% ART, 0%% PrEP reduction): %.4f per 100 PY\n",
            baseline_incidence))

# Summary CSV of the full grid-based cost-mapping surface (the deduplicated
# funding->coverage->incidence table that feeds the heatmaps and the scenario
# table), so these values are tracked in git across changes.
write.csv(as.data.frame(model_proportions_mean_incidence),
          paste0(plots_folder_name, "cost_mapping_grid_results.csv"),
          row.names = FALSE)

# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

# ---- Heatmaps: discrete tile matrices -------------------------------------
tile_heatmap_incidence <- plot_tile_matrix_heatmap(
  funding_reduction_scenarios = model_proportions_mean_incidence,
  mean_var = "year10_mean_incidence",
  pct_change_var = "relative_incidence_change_pct",
  var_name = "Mean Incidence")

tile_heatmap_risk_ratio <- plot_tile_matrix_heatmap(
  funding_reduction_scenarios = model_proportions_mean_incidence,
  mean_var = "incidence_risk_ratio",
  pct_change_var = "relative_incidence_risk_ratio_change_pct",
  var_name = "Mean Incidence\nRisk Ratio")

# ---- Heatmaps: contour, with iso-budget lines -----------------------------
contour_heatmap_incidence <- plot_contour_heatmap(
  funding_reduction_scenarios = model_proportions_mean_incidence,
  mean_var = "year10_mean_incidence",
  pct_change_var = "relative_incidence_change_pct",
  var_name = "Mean Incidence",
  baseline_value = baseline_incidence)

contour_heatmap_incidence_risk_ratio <- plot_contour_heatmap(
  funding_reduction_scenarios = model_proportions_mean_incidence,
  mean_var = "incidence_risk_ratio",
  pct_change_var = "relative_incidence_risk_ratio_change_pct",
  var_name = "Mean Incidence\nRisk Ratio",
  baseline_value = 1.0)

total_baseline <- p$G_ART_baseline + p$G_PrEP_baseline
total_pct_cuts_valid <- seq(0.1, 1, by = 0.1) * total_baseline
isobudget_data <- create_isobudget_lines(
  total_pct_cuts_valid, p$G_ART_baseline, p$G_PrEP_baseline)

iso_heatmap_incidence_pct <- add_iso_budget_lines_to_plot(
  isobudget_data, heatmap_plots_list = contour_heatmap_incidence,
  total_baseline = total_baseline,
  G_ART_baseline = p$G_ART_baseline, G_PrEP_baseline = p$G_PrEP_baseline)

iso_heatmap_risk_ratio_pct <- add_iso_budget_lines_to_plot(
  isobudget_data, heatmap_plots_list = contour_heatmap_incidence_risk_ratio,
  total_baseline = total_baseline,
  G_ART_baseline = p$G_ART_baseline, G_PrEP_baseline = p$G_PrEP_baseline)

# ---- Save heatmaps --------------------------------------------------------
tile_incidence_files <- c("discrete_heatmap_incidence_main_year10",
                          "discrete_heatmap_incidence_relative_year10")
for (i in seq_along(tile_incidence_files)) {
  ggsave(paste0(plots_folder_name, tile_incidence_files[i], ".png"),
         plot = tile_heatmap_incidence[[i]], width = 12, height = 10, dpi = 300)
}

tile_irr_files <- c("discrete_heatmap_IRR_main_year10",
                    "discrete_heatmap_IRR_relative_year10")
for (i in seq_along(tile_irr_files)) {
  ggsave(paste0(plots_folder_name, tile_irr_files[i], ".png"),
         plot = tile_heatmap_risk_ratio[[i]], width = 12, height = 10, dpi = 300)
}

incidence_contour_files <- c("heatmap_incidence_main_year10",
                             "heatmap_incidence_relative_year10")
for (i in seq_along(incidence_contour_files)) {
  ggsave(paste0(plots_folder_name, incidence_contour_files[i], ".png"),
         plot = iso_heatmap_incidence_pct[[i]], width = 10, height = 8, dpi = 300)
}

risk_ratio_contour_files <- c("heatmap_risk_ratio_main_year10",
                              "heatmap_risk_ratio_relative_year10")
for (i in seq_along(risk_ratio_contour_files)) {
  ggsave(paste0(plots_folder_name, risk_ratio_contour_files[i], ".png"),
         plot = iso_heatmap_risk_ratio_pct[[i]], width = 10, height = 8, dpi = 300)
}

# ---- Forest plot by funding reduction (main-text Figure 1) ----------------
# Coverage-reduction per unit funding reduction (invert paper Eqs. 1-2).
art_cov_per_funding  <- (p$P_ART_baseline  - p$gamma_ART)  / p$P_ART_baseline
prep_cov_per_funding <- (p$P_PrEP_baseline - p$gamma_PrEP) / p$P_PrEP_baseline

forest <- plot_funding_forest(
  gp_incidence_fit = gp_incidence_fit,
  art_cov_per_funding = art_cov_per_funding,
  prep_cov_per_funding = prep_cov_per_funding,
  funding_levels = c(0.10, 0.25, 0.40, 0.55, 0.70),
  common_random_numbers = USE_COMMON_RANDOM_NUMBERS,
  save_path = paste0(plots_folder_name, "forest_incidence_risk_ratio_funding"))

# Summary CSV of the forest results (mean incidence risk ratio + 90% CI per
# scenario x funding level), so the values behind Figure 1 are tracked in git.
forest_summary <- as.data.frame(forest$data)
forest_summary <- forest_summary[order(forest_summary$scenario,
                                        forest_summary$reduction), ]
write.csv(forest_summary,
          paste0(plots_folder_name, "forest_incidence_risk_ratio_funding.csv"),
          row.names = FALSE)

# ---- Bar plots of selected policy scenarios -------------------------------
plot_funding_bar(
  key_scenarios = model_proportions_mean_incidence,
  art_funding_reduction_pct  = c(0, 10, 0, 25, 0, 40, 0, 40),
  prep_funding_reduction_pct = c(0, 0, 10, 0, 25, 0, 40, 40),
  title = "Selected Policy Scenarios",
  save_path = paste0(plots_folder_name, "bar_plot"),
  y_metric = "all",
  incidence_ci_results = incidence_ci_results,
  relative_change_ci_results = relative_change_ci_results,
  irr_ci_results = irr_ci_results)

# ---------------------------------------------------------------------------
# Table
# ---------------------------------------------------------------------------
key_scenarios <- data.frame(
  Scenario = c("Baseline (No cuts)", "10% ART cut only", "10% PrEP cut only",
               "25% ART cut only", "25% PrEP cut only", "40% ART cut only",
               "40% PrEP cut only", "40% cut for both"),
  art_funding_reduction_pct  = c(0, 10, 0, 25, 0, 40, 0, 40),
  prep_funding_reduction_pct = c(0, 0, 10, 0, 25, 0, 40, 40))

create_funding_scenario_table(
  key_scenarios,
  funding_data = model_proportions_mean_incidence,
  G_ART_baseline_val = p$G_ART_baseline,
  G_PrEP_baseline_val = p$G_PrEP_baseline,
  incidence_ci_results = incidence_ci_results,
  relative_change_ci_results = relative_change_ci_results,
  irr_ci_results = irr_ci_results,
  output_dir = plots_folder_name,
  metric_type = "both",
  common_random_numbers = USE_COMMON_RANDOM_NUMBERS,
  gp_incidence_fit = gp_incidence_fit,
  n_samples_per_checkpoint = p$n_samples_per_checkpoint,
  horizon_tick = p$horizon_tick,
  ci_probs = p$ci_probs)
