# =============================================================================
# Placeholder figures for the paper's main text, at the new baseline r
# -----------------------------------------------------------------------------
# The paper's real Figure 1 (forest) and Figure 2 (ribbon+heatmap) are built by
# 04_cost_mapping/R/plot_forest.R / plot_ribbon.R, which `library(randplot)` --
# the RAND house-style package, not installed anywhere on this machine. This
# script reuses the underlying, randplot-free DATA functions those files call
# (compute_scenario_draws_crn/compute_scenario_trajectory_draws_crn from
# irr_common_random_numbers.R; predict_incidence_grid/build_funding_scenarios
# from outcomes.R/cost_model.R) and draws plain-ggplot2 stand-ins instead, so
# the paper has *something* visually present and directionally correct at the
# new baseline r, clearly marked PLACEHOLDER pending randplot.
#
# The heatmap sub-panel uses a reduced ngrid=25 (vs. the real pipeline's 101)
# purely for speed -- this is a stand-in, not the final figure.
#
# Run with 04_cost_mapping/sensitivity/ as the working directory:
#   /usr/local/bin/Rscript run_paper_figures_placeholder.R
# =============================================================================

library(ggplot2)
library(hetGP)
library(dplyr)
library(viridis)

source("../R/parameters.R")
source("../R/cost_model.R")
source("../R/outcomes.R")
source("../R/irr_common_random_numbers.R")
source("R/recalibrate.R")

set.seed(42)

fig_dir <- "../../hiv-funding-paper-tmp/figures"
if (!dir.exists(fig_dir)) stop("Expected paper figures/ dir not found: ", fig_dir)

load("../../03_intervention_scenario_surrogate/output/surrogate.Rdata")
p_base <- cost_mapping_params(yml = "../cost_params.yml")

R_UPPER    <- 0.78
R_LOWER    <- 0.25
R_BASELINE <- sqrt(R_LOWER * R_UPPER)
P0 <- p_base$P_PrEP_baseline
baseline <- recalibrate_prep(R_BASELINE, P0, p_base)

art_cov_per_funding  <- (p_base$P_ART_baseline - p_base$gamma_ART) / p_base$P_ART_baseline
prep_cov_per_funding <- baseline$prep_cov_per_funding
reduction_levels <- c(0.10, 0.25, 0.40)
n_spc <- p_base$n_samples_per_checkpoint
ci_probs <- p_base$ci_probs

PLACEHOLDER_CAPTION_NOTE <- sprintf(
  "PLACEHOLDER pending randplot (house style) -- baseline r=%.1f%%", R_BASELINE * 100)

# =============================================================================
# Figure 1 stand-in: two-panel forest plot (A. use reduction, B. funding
# reduction), mirrors plot_forest.R's .forest_panel_data() construction.
# =============================================================================
build_forest_data <- function(art_cov_per_unit, prep_cov_per_unit, panel_label) {
  scenarios <- rbind(
    data.frame(scenario = "Reduce PrEP only",         art_red = 0,                prep_red = reduction_levels),
    data.frame(scenario = "Reduce ART only",          art_red = reduction_levels, prep_red = 0),
    data.frame(scenario = "Reduce both PrEP and ART", art_red = reduction_levels, prep_red = reduction_levels)
  )
  scenarios$reduction <- pmax(scenarios$art_red, scenarios$prep_red)
  scenarios$art_cov  <- scenarios$art_red  * art_cov_per_unit
  scenarios$prep_cov <- scenarios$prep_red * prep_cov_per_unit
  newX <- as.matrix(cbind(-scenarios$art_cov, -scenarios$prep_cov))
  draws <- compute_scenario_draws_crn(newX, gp_incidence_fit,
                                       n_samples_per_checkpoint = n_spc,
                                       common_random_numbers = TRUE, tick = NULL)
  irr <- summarise_draws(draws$irr, probs = ci_probs)
  scenarios$mean  <- irr$mean
  scenarios$lower <- irr$ci_lower
  scenarios$upper <- irr$ci_upper
  scenarios$panel <- panel_label
  scenarios
}

forest_A <- build_forest_data(1, 1, "A. Intervention use reduction")
forest_B <- build_forest_data(art_cov_per_funding, prep_cov_per_funding, "B. Government funding reduction")
forest_dt <- rbind(forest_A, forest_B)
forest_dt$scenario <- factor(forest_dt$scenario,
  levels = c("Reduce PrEP only", "Reduce ART only", "Reduce both PrEP and ART"))
forest_dt$panel <- factor(forest_dt$panel,
  levels = c("A. Intervention use reduction", "B. Government funding reduction"))

forest_plot <- ggplot(forest_dt, aes(x = mean, y = factor(reduction * 100), color = scenario)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.4, linewidth = 0.8,
                 position = position_dodge(width = 0.7)) +
  geom_point(size = 3, position = position_dodge(width = 0.7)) +
  scale_color_manual(values = c("Reduce PrEP only" = "#45aF84",
                                 "Reduce ART only" = "#597cbe",
                                 "Reduce both PrEP and ART" = "#af61a7"), name = NULL) +
  facet_wrap(~panel, ncol = 2) +
  labs(x = "Mean incidence rate ratio", y = "Reduction level (%)",
       title = "Figure 1 PLACEHOLDER (pending randplot house style)",
       subtitle = PLACEHOLDER_CAPTION_NOTE) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "forest_incidence_risk_ratio_funding_PLACEHOLDER.png"),
       forest_plot, width = 11, height = 5.5, dpi = 150)

# =============================================================================
# Figure 2 stand-in: ribbon (IRR trajectory, both-funding cuts) + heatmap
# (year-10 IRR over ART x PrEP funding-reduction grid, reduced ngrid=25 for
# speed -- this is a stand-in, not the final-resolution figure).
# =============================================================================
traj_levels <- c(0.10, 0.20, 0.40)
traj_list <- lapply(traj_levels, function(lvl) {
  pt <- matrix(c(-lvl * art_cov_per_funding, -lvl * prep_cov_per_funding), nrow = 1)
  tr <- compute_scenario_trajectory_draws_crn(pt, gp_incidence_fit,
                                               n_samples_per_checkpoint = n_spc,
                                               common_random_numbers = TRUE)
  s <- summarise_draws(tr$irr[[1]], probs = ci_probs)
  s$year <- seq_len(nrow(s)) - 1
  s$level <- sprintf("%.0f%%", lvl * 100)
  s
})
ribbon_data <- do.call(rbind, traj_list)
ribbon_data$level <- factor(ribbon_data$level, levels = sprintf("%.0f%%", traj_levels * 100))

ribbon_plot <- ggplot(ribbon_data, aes(x = year, y = mean, color = level, fill = level)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  labs(x = "Year", y = "Incidence rate ratio", color = "Funding cut", fill = "Funding cut",
       title = "A. IRR trajectory (both ART+PrEP funding cuts)") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")

# ---- Reduced-resolution heatmap (ngrid=25) ---------------------------------
p_heatmap <- p_base
p_heatmap$gamma_PrEP <- baseline$gamma_PrEP
p_heatmap$beta_PrEP  <- baseline$beta_PrEP
p_heatmap$ngrid <- 25L
p_heatmap$P_PrEP_baseline <- baseline$P_PrEP_baseline_check  # == P0, by construction

grid <- predict_incidence_grid(p_heatmap, gp_incidence_fit, predict_hiv_composite_surrogate)
scenarios_grid <- build_funding_scenarios(
  grid$model_incidence, p_heatmap,
  gp_incidence_fit = gp_incidence_fit,
  predict_composite_fn = predict_hiv_composite_surrogate)

heatmap_data <- scenarios_grid[scenarios_grid$pct_delta_ART_fund  >= 0 & scenarios_grid$pct_delta_ART_fund  <= 50 &
                                scenarios_grid$pct_delta_PrEP_fund >= 0 & scenarios_grid$pct_delta_PrEP_fund <= 50, ]

heatmap_plot <- ggplot(heatmap_data, aes(x = pct_delta_ART_fund, y = pct_delta_PrEP_fund, fill = incidence_risk_ratio)) +
  geom_tile() +
  scale_fill_viridis(name = "IRR", option = "C") +
  labs(x = "ART funding reduction (%)", y = "PrEP funding reduction (%)",
       title = "B. Year-10 IRR heatmap (ngrid=25 placeholder, real fig uses 101)") +
  theme_minimal(base_size = 11)

library(patchwork)
combined_full <- (ribbon_plot | heatmap_plot) +
  patchwork::plot_annotation(
    title = "Figure 2 PLACEHOLDER (pending randplot house style)",
    subtitle = PLACEHOLDER_CAPTION_NOTE)

ggsave(file.path(fig_dir, "incidence_ribbon_heatmap_panel_PLACEHOLDER.png"),
       combined_full, width = 11, height = 5.5, dpi = 150)

cat("\nWrote placeholder figures:\n")
cat(sprintf("  %s/forest_incidence_risk_ratio_funding_PLACEHOLDER.png\n", fig_dir))
cat(sprintf("  %s/incidence_ribbon_heatmap_panel_PLACEHOLDER.png\n", fig_dir))
