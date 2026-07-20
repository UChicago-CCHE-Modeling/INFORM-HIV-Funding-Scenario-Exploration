# =============================================================================
# Cost mapping — archived extras
# -----------------------------------------------------------------------------
# Regenerates the cost-mapping outputs that are NOT the three main-text outputs
# produced by ../../script/run_cost_mapping.R:
#   - contour heatmaps (mean incidence, IRR, and their relative-change variants,
#     with iso-budget lines)
#   - discrete tile-matrix heatmaps
#   - bar plots of selected policy scenarios
#   - cost_mapping_grid_results.csv (full deduplicated grid surface)
#
# These were moved out of the main pipeline; this script is self-contained and
# reproduces them into ../output. It sources the still-live modules from
# ../../R and the archived plot_bar.R from ../R.
#
# Run with archive/script/ as the working directory:
#   cd 04_cost_mapping/archive/script && Rscript run_cost_mapping_extras.R
# =============================================================================

library(ggplot2)
library(data.table)
library(hetGP)     # so predict() dispatches on the stored surrogate models
library(viridis)
library(dplyr)
library(tidyr)

# ---- Modules --------------------------------------------------------------
R_DIR         <- "../../R"   # live modules
ARCHIVE_R_DIR <- "../R"      # archived modules
source(file.path(R_DIR, "parameters.R"))
source(file.path(R_DIR, "cost_model.R"))
source(file.path(R_DIR, "outcomes.R"))
source(file.path(R_DIR, "plot_heatmap.R"))
source(file.path(ARCHIVE_R_DIR, "plot_bar.R"))

# ---- Inputs / outputs -----------------------------------------------------
load("../../../03_intervention_scenario_surrogate/output/surrogate.Rdata")
plots_folder_name <- "../output/"
if (!dir.exists(plots_folder_name)) dir.create(plots_folder_name, recursive = TRUE)

# ---- Parameters -----------------------------------------------------------
# cost_params.yml lives at the stage root; parameters.R's default path assumes
# the live script/ working dir, so pass the archive-relative path explicitly.
p <- cost_mapping_params("../../cost_params.yml")

# ---- Surrogate incidence surface + posterior predictive intervals ---------
set.seed(12345)
grid <- predict_incidence_grid(p, gp_incidence_fit, predict_hiv_composite_surrogate)
ppi  <- compute_ppi(grid, p)
incidence_ci_results       <- ppi$incidence_ci_results
relative_change_ci_results <- ppi$relative_change_ci_results
irr_ci_results             <- ppi$irr_ci_results

# ---- Funding -> coverage -> incidence outcomes ----------------------------
model_proportions_mean_incidence <- build_funding_scenarios(
  grid$model_incidence, p,
  gp_incidence_fit = gp_incidence_fit,
  predict_composite_fn = predict_hiv_composite_surrogate)
baseline_incidence <- attr(model_proportions_mean_incidence, "baseline_incidence")

# Summary CSV of the full grid-based cost-mapping surface.
write.csv(as.data.frame(model_proportions_mean_incidence),
          paste0(plots_folder_name, "cost_mapping_grid_results.csv"),
          row.names = FALSE)

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
