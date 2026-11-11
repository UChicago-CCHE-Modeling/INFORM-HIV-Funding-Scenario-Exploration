# =============================================================================
# Case comparison: gamma_PrEP payer-mix parameter (r) sensitivity
# -----------------------------------------------------------------------------
# Exploratory analysis following the 2026-08-10 meeting
# (INFORM-HIV-Analysis/Raff_tmp/gamma_prep_meeting_followup.tex): the national
# r=0.78 (share of PrEP users privately insured, hence immune to a government
# funding cut) is likely too high for INFORM-HIV's actual population. The
# meeting provisionally proposed r ~= 25%.
#
# Three cases, funding-reduction convention (labelled % = government funding
# cut, mapped to coverage reduction via the cost model):
#   Baseline: r=0.78, P0 = current model baseline (~36%) -- the paper's
#             currently committed values, unchanged.
#   Case 1:   r=0.25, P0 unchanged (~36%) -- only the payer-mix assumption
#             revised.
#   Case 2:   r=0.25, P0=0.24 -- payer-mix AND baseline-coverage-level both
#             revised.
#
# NOTE: (P0 - gamma_PrEP)/P0 = 1 - r algebraically (see R/recalibrate.R), so
# Case 1 and Case 2 are expected to be IDENTICAL on every incidence/IRR output
# below -- P0 only feeds beta_PrEP's recalibrated value, not the coverage-
# reduction-per-unit-funding factor that drives the surrogate query. That is
# not a bug; it is the direct consequence of the additive coverage model
# (memo sec. 1.2) and confirms the memo's own analytical note (sec. 1.4).
#
# This script does NOT touch cost_params.yml, the tracked output/ directory,
# or plot_forest.R/plot_ribbon.R (which require the unavailable `randplot`
# package). It only uses the pure-computation modules (parameters.R,
# irr_common_random_numbers.R) and produces its own plain ggplot2 output.
#
# Run with 04_cost_mapping/sensitivity/ as the working directory:
#   /usr/local/bin/Rscript run_case_comparison.R
# =============================================================================

library(ggplot2)
library(hetGP)  # so predict() dispatches on the stored surrogate GPs

source("../R/parameters.R")
source("../R/irr_common_random_numbers.R")
source("R/recalibrate.R")

set.seed(42)  # reproducibility for the surrogate's predictive draws

out_dir <- "output"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---- Load surrogate + base parameters --------------------------------------
load("../../03_intervention_scenario_surrogate/output/surrogate.Rdata")
p_base <- cost_mapping_params(yml = "../cost_params.yml")

cat(sprintf("Current committed baseline: r=0.78 (implied), gamma_PrEP=%.4f, beta_PrEP=%.4f, P_PrEP_baseline=%.4f\n",
            p_base$gamma_PrEP, p_base$beta_PrEP, p_base$P_PrEP_baseline))

reduction_levels <- c(0.10, 0.25, 0.40)
n_samples_per_checkpoint <- p_base$n_samples_per_checkpoint
ci_probs <- p_base$ci_probs

art_cov_per_funding <- (p_base$P_ART_baseline - p_base$gamma_ART) / p_base$P_ART_baseline

# ---- Build the funding-reduction scenario matrix (no query yet) -----------
# Mirrors 04_cost_mapping/R/plot_forest.R's .forest_panel_data() "funding
# reduction" panel, without the randplot-dependent plotting.
build_case_scenarios <- function(prep_cov_per_funding, case_label) {
  scenarios <- rbind(
    data.frame(scenario = "Reduce PrEP only",
               art_red = 0, prep_red = reduction_levels),
    data.frame(scenario = "Reduce ART only",
               art_red = reduction_levels, prep_red = 0),
    data.frame(scenario = "Reduce both PrEP and ART",
               art_red = reduction_levels, prep_red = reduction_levels)
  )
  scenarios$reduction <- pmax(scenarios$art_red, scenarios$prep_red)
  scenarios$art_cov  <- scenarios$art_red  * art_cov_per_funding
  scenarios$prep_cov <- scenarios$prep_red * prep_cov_per_funding
  scenarios$case <- case_label
  scenarios
}

# ---- Three cases ------------------------------------------------------------
prep_cov_baseline <- (p_base$P_PrEP_baseline - p_base$gamma_PrEP) / p_base$P_PrEP_baseline
case1 <- recalibrate_prep(r = 0.25, P0 = p_base$P_PrEP_baseline, p_base = p_base)
case2 <- recalibrate_prep(r = 0.25, P0 = 0.24, p_base = p_base)

cat(sprintf("Case 1 (r=0.25, P0=%.4f): gamma_PrEP=%.4f, beta_PrEP=%.4f, prep_cov_per_funding=%.4f\n",
            case1$P0, case1$gamma_PrEP, case1$beta_PrEP, case1$prep_cov_per_funding))
cat(sprintf("Case 2 (r=0.25, P0=%.4f): gamma_PrEP=%.4f, beta_PrEP=%.4f, prep_cov_per_funding=%.4f\n",
            case2$P0, case2$gamma_PrEP, case2$beta_PrEP, case2$prep_cov_per_funding))

data_baseline <- build_case_scenarios(prep_cov_baseline, sprintf("Baseline (r=0.78, P0=%.0f%%)", p_base$P_PrEP_baseline * 100))
data_case1    <- build_case_scenarios(case1$prep_cov_per_funding, sprintf("Case 1 (r=25%%, P0=%.0f%%)", case1$P0 * 100))
data_case2    <- build_case_scenarios(case2$prep_cov_per_funding, sprintf("Case 2 (r=25%%, P0=%.0f%%)", case2$P0 * 100))

# Query ALL cases in ONE compute_scenario_draws_crn() call so cases with
# identical (art_cov, prep_cov) query points (Case 1 and Case 2, since
# prep_cov_per_funding depends only on r) get literally the same surrogate
# draws (same shared z per checkpoint, same GP mean/sd at the same point) --
# not just approximately equal up to Monte Carlo noise from separate calls.
case_comparison <- rbind(data_baseline, data_case1, data_case2)
newX_native <- as.matrix(cbind(-case_comparison$art_cov, -case_comparison$prep_cov))
draws <- compute_scenario_draws_crn(
  newX_native, gp_incidence_fit,
  n_samples_per_checkpoint = n_samples_per_checkpoint,
  common_random_numbers = TRUE, tick = NULL)
irr_summary <- summarise_draws(draws$irr, probs = ci_probs)
case_comparison$irr_mean     <- irr_summary$mean
case_comparison$irr_ci_lower <- irr_summary$ci_lower
case_comparison$irr_ci_upper <- irr_summary$ci_upper
case_comparison$case <- factor(case_comparison$case, levels = unique(case_comparison$case))
case_comparison$scenario <- factor(case_comparison$scenario,
  levels = c("Reduce PrEP only", "Reduce ART only", "Reduce both PrEP and ART"))

write.csv(case_comparison, file.path(out_dir, "case_comparison.csv"), row.names = FALSE)

# ---- Verification: Case 1 vs Case 2 should be numerically identical --------
c1 <- case_comparison[case_comparison$case == levels(case_comparison$case)[2], "irr_mean"]
c2 <- case_comparison[case_comparison$case == levels(case_comparison$case)[3], "irr_mean"]
identical_check <- isTRUE(all.equal(c1, c2))
cat(sprintf("\nCase 1 vs Case 2 identical on irr_mean: %s (max abs diff = %.2e)\n",
            identical_check, max(abs(c1 - c2))))

# ---- Plot: forest-style comparison, plain ggplot2 (no randplot) ------------
forest_plot <- ggplot(case_comparison,
                       aes(x = irr_mean, y = factor(reduction * 100), color = case)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  geom_errorbarh(aes(xmin = irr_ci_lower, xmax = irr_ci_upper),
                 height = 0.3, position = position_dodge(width = 0.6)) +
  geom_point(size = 2.5, position = position_dodge(width = 0.6)) +
  facet_wrap(~scenario, ncol = 1) +
  labs(x = "Mean incidence rate ratio (vs. no-cut baseline)",
       y = "Government funding reduction (%)",
       color = NULL,
       title = "PrEP payer-mix sensitivity: r=0.78 (current) vs r=25% (proposed)",
       subtitle = "Case 1 and Case 2 overlap exactly — (P⁰-γ)/P⁰ = 1-r is algebraically independent of P⁰") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "case_comparison_forest.png"), forest_plot,
       width = 8, height = 9, dpi = 150)

# ---- Trajectory comparison: Baseline vs Case 1, "reduce both, 25% cut" ----
traj_point_baseline <- matrix(c(-0.25 * art_cov_per_funding, -0.25 * prep_cov_baseline), nrow = 1)
traj_point_case1    <- matrix(c(-0.25 * art_cov_per_funding, -0.25 * case1$prep_cov_per_funding), nrow = 1)

traj_baseline <- compute_scenario_trajectory_draws_crn(
  traj_point_baseline, gp_incidence_fit,
  n_samples_per_checkpoint = n_samples_per_checkpoint, common_random_numbers = TRUE)
traj_case1 <- compute_scenario_trajectory_draws_crn(
  traj_point_case1, gp_incidence_fit,
  n_samples_per_checkpoint = n_samples_per_checkpoint, common_random_numbers = TRUE)

summ_baseline <- summarise_draws(traj_baseline$irr[[1]], probs = ci_probs)
summ_baseline$year <- seq_len(nrow(summ_baseline)) - 1
summ_baseline$case <- "Baseline (r=0.78)"

summ_case1 <- summarise_draws(traj_case1$irr[[1]], probs = ci_probs)
summ_case1$year <- seq_len(nrow(summ_case1)) - 1
summ_case1$case <- "Case 1 (r=25%)"

trajectory_comparison <- rbind(summ_baseline, summ_case1)
write.csv(trajectory_comparison, file.path(out_dir, "trajectory_comparison.csv"), row.names = FALSE)

trajectory_plot <- ggplot(trajectory_comparison, aes(x = year, y = mean, color = case, fill = case)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  labs(x = "Year", y = "Incidence rate ratio (vs. no-cut baseline)",
       title = "25% funding cut (ART + PrEP): IRR trajectory",
       subtitle = "Lower r → larger PrEP funding pass-through → higher IRR",
       color = NULL, fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "trajectory_comparison.png"), trajectory_plot,
       width = 8, height = 5, dpi = 150)

cat("\nWrote:\n  output/case_comparison.csv\n  output/case_comparison_forest.png\n  output/trajectory_comparison.csv\n  output/trajectory_comparison.png\n")
