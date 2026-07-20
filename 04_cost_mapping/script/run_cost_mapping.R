# =============================================================================
# Cost mapping
# -----------------------------------------------------------------------------
# Runs the cost-mapping pipeline and writes the three main-text outputs:
#   - forest_incidence_risk_ratio_funding  (two-panel forest plot + CSV)
#   - incidence_ribbon_heatmap_panel        (ribbon + heatmap panel + ribbon CSV)
#   - funding_scenarios_table               (Table 1 .tex + CSV)
# Load parameters and the stage-03 surrogate, map funding reductions to coverage
# outcomes, and produce those figures and the table. Parameters are defined in
# R/parameters.R; the reusable functions live in R/*.R. This script calls them
# in order.
#
# The heatmaps, discrete tile matrices, bar plots, and the full grid CSV that
# this pipeline used to emit now live in ../archive (run
# archive/script/run_cost_mapping_extras.R to regenerate them).
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
source(file.path(R_DIR, "plot_heatmap.R"))  # plot_contour_heatmap() for the panel
source(file.path(R_DIR, "plot_forest.R"))
source(file.path(R_DIR, "plot_ribbon.R"))
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

# ---- Surrogate incidence surface ------------------------------------------
# The grid prediction draws from the surrogate's predictive distribution, so pin
# a seed here for reproducibility. (The forest plot, ribbon panel, and table CRN
# paths set their own seed downstream.) The grid feeds the panel-B heatmap.
set.seed(12345)
grid <- predict_incidence_grid(p, gp_incidence_fit, predict_hiv_composite_surrogate)

# ---- Funding -> coverage -> incidence outcomes ----------------------------
model_proportions_mean_incidence <- build_funding_scenarios(
  grid$model_incidence, p,
  gp_incidence_fit = gp_incidence_fit,
  predict_composite_fn = predict_hiv_composite_surrogate)
baseline_incidence <- attr(model_proportions_mean_incidence, "baseline_incidence")
cat(sprintf("Baseline incidence (0%% ART, 0%% PrEP reduction): %.4f per 100 PY\n",
            baseline_incidence))

# Coverage-reduction per unit funding reduction (invert paper Eqs. 1-2). Shared
# by the forest plot, the ribbon panel, and the table.
art_cov_per_funding  <- (p$P_ART_baseline  - p$gamma_ART)  / p$P_ART_baseline
prep_cov_per_funding <- (p$P_PrEP_baseline - p$gamma_PrEP) / p$P_PrEP_baseline

# ---------------------------------------------------------------------------
# Output 1 — Forest plot by funding reduction (main-text Figure 1)
# ---------------------------------------------------------------------------
forest <- plot_funding_forest(
  gp_incidence_fit = gp_incidence_fit,
  art_cov_per_funding = art_cov_per_funding,
  prep_cov_per_funding = prep_cov_per_funding,
  reduction_levels = c(0.10, 0.25, 0.40),
  n_samples_per_checkpoint = p$n_samples_per_checkpoint,
  ci_probs = p$ci_probs,
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

# ---------------------------------------------------------------------------
# Output 2 — Combined ribbon + heatmap panel (main-text figure)
# ---------------------------------------------------------------------------
# Two-panel figure. Panel A: mean incidence-risk-ratio trajectory over years
# 0-10 under three government-funding reduction scenarios that cut ART and PrEP
# simultaneously by 10%, 20% and 40%, each with 50% and 95% credible ribbons
# (same common random numbers as the forest plot). Each ribbon is coloured by
# the panel-B heatmap band its terminal risk ratio falls in, so the two panels
# share a colour language. Panel B: the government-funding mean-incidence-risk-
# ratio contour heatmap with both axes restricted to 0-50% and fixed 0.25-wide
# colour bands.
ribbon_heatmap <- plot_ribbon_heatmap_panel(
  gp_incidence_fit = gp_incidence_fit,
  grid_scenarios = model_proportions_mean_incidence,
  art_cov_per_funding = art_cov_per_funding,
  prep_cov_per_funding = prep_cov_per_funding,
  reduction_levels = c(0.10, 0.20, 0.40),
  heatmap_max = 50,
  rr_breaks = seq(1, 2.25, 0.25),
  n_samples_per_checkpoint = p$n_samples_per_checkpoint,
  common_random_numbers = USE_COMMON_RANDOM_NUMBERS,
  save_path = paste0(plots_folder_name, "incidence_ribbon_heatmap_panel"))

# Summary CSV of the ribbon panel results (mean + 50%/95% credible bounds per
# scenario x year), so the values behind panel A are tracked in git.
write.csv(as.data.frame(ribbon_heatmap$ribbon_data),
          paste0(plots_folder_name, "incidence_trajectory_ribbon.csv"),
          row.names = FALSE)

# ---------------------------------------------------------------------------
# Output 3 — Scenario table (main-text Table 1)
# -----------------------------------------------------------------------------
# Two stacked blocks, shared baseline:
#   A. Intervention use reduction (labelled % = coverage reduction)
#   B. Government funding reduction (labelled % mapped to coverage via cost model)
# Both blocks and both figure panels use the same 10/25/40% levels and common
# random numbers. Also writes a CSV of the same values for git tracking.
# ---------------------------------------------------------------------------
table_out <- create_funding_scenario_table(
  gp_incidence_fit = gp_incidence_fit,
  art_cov_per_funding = art_cov_per_funding,
  prep_cov_per_funding = prep_cov_per_funding,
  reduction_levels = c(0.10, 0.25, 0.40),
  n_samples_per_checkpoint = p$n_samples_per_checkpoint,
  horizon_tick = p$horizon_tick,
  ci_probs = p$ci_probs,
  common_random_numbers = USE_COMMON_RANDOM_NUMBERS,
  output_dir = plots_folder_name)

write.csv(as.data.frame(table_out),
          paste0(plots_folder_name, "funding_scenarios_table.csv"),
          row.names = FALSE)
