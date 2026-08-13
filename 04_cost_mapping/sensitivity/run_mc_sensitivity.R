# =============================================================================
# Monte Carlo sensitivity: r and P0 uncertainty, +/-50% uniform around the
# proposed mode (r=25%, P0=24%)
# -----------------------------------------------------------------------------
# For each draw, gamma_PrEP and beta_PrEP are recalibrated from (r, P0) via
# R/recalibrate.R, then the surrogate is queried at a single scenario point
# ("reduce ART + PrEP both by 25% government funding") using that draw's own
# prep_cov_per_funding. All N draws are queried in ONE compute_scenario_draws_crn()
# call (an N x 2 matrix), not a loop, since each checkpoint's surrogate call
# already batches over points.
#
# Expectation (see run_case_comparison.R header / recalibrate.R for the
# algebra): the resulting IRR should vary strongly with r and show ~no
# dependence on P0, since (P0-gamma_PrEP)/P0 = 1-r regardless of P0.
#
# Run with 04_cost_mapping/sensitivity/ as the working directory:
#   /usr/local/bin/Rscript run_mc_sensitivity.R
# =============================================================================

library(ggplot2)
library(hetGP)

source("../R/parameters.R")
source("../R/irr_common_random_numbers.R")
source("R/recalibrate.R")

out_dir <- "output"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

load("../../03_intervention_scenario_surrogate/output/surrogate.Rdata")
p_base <- cost_mapping_params(yml = "../cost_params.yml")

N <- 200
set.seed(20260813)

r_mode  <- 0.25
p0_mode <- 0.24
r_draws  <- runif(N, r_mode * 0.5, r_mode * 1.5)
p0_draws <- runif(N, p0_mode * 0.5, p0_mode * 1.5)

art_cov_per_funding <- (p_base$P_ART_baseline - p_base$gamma_ART) / p_base$P_ART_baseline
FUNDING_CUT <- 0.25  # headline scenario: 25% government funding cut, ART+PrEP both

recal <- lapply(seq_len(N), function(i) recalibrate_prep(r_draws[i], p0_draws[i], p_base))
gamma_PrEP_draws <- vapply(recal, function(x) x$gamma_PrEP, numeric(1))
beta_PrEP_draws  <- vapply(recal, function(x) x$beta_PrEP, numeric(1))
prep_cov_draws   <- vapply(recal, function(x) x$prep_cov_per_funding, numeric(1))

newX_native <- cbind(-FUNDING_CUT * art_cov_per_funding, -FUNDING_CUT * prep_cov_draws)

draws <- compute_scenario_draws_crn(
  newX_native, gp_incidence_fit,
  n_samples_per_checkpoint = p_base$n_samples_per_checkpoint,
  common_random_numbers = TRUE, tick = NULL)
irr_summary <- summarise_draws(draws$irr, probs = p_base$ci_probs)

mc_results <- data.frame(
  draw = seq_len(N),
  r = r_draws,
  P0 = p0_draws,
  gamma_PrEP = gamma_PrEP_draws,
  beta_PrEP = beta_PrEP_draws,
  prep_cov_per_funding = prep_cov_draws,
  irr_mean = irr_summary$mean,
  irr_ci_lower = irr_summary$ci_lower,
  irr_ci_upper = irr_summary$ci_upper
)
write.csv(mc_results, file.path(out_dir, "mc_sensitivity_draws.csv"), row.names = FALSE)

cor_r  <- cor(mc_results$irr_mean, mc_results$r)
cor_p0 <- cor(mc_results$irr_mean, mc_results$P0)

summary_lines <- c(
  sprintf("Monte Carlo sensitivity: N=%d draws, r ~ U(%.3f, %.3f), P0 ~ U(%.3f, %.3f)",
          N, r_mode * 0.5, r_mode * 1.5, p0_mode * 0.5, p0_mode * 1.5),
  sprintf("Scenario queried: reduce ART+PrEP both by %.0f%% government funding, year-10-averaged IRR", FUNDING_CUT * 100),
  sprintf("IRR range across draws: [%.4f, %.4f], mean=%.4f", min(mc_results$irr_mean), max(mc_results$irr_mean), mean(mc_results$irr_mean)),
  sprintf("cor(irr_mean, r)  = %.4f  (expected: strongly negative -- higher r = more privately-insured/cut-immune PrEP users = smaller funding pass-through)", cor_r),
  sprintf("cor(irr_mean, P0) = %.4f  (expected: ~0, P0 cancels out algebraically)", cor_p0)
)
writeLines(summary_lines, file.path(out_dir, "mc_sensitivity_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")

# ---- Diagnostic plots -------------------------------------------------------
plot_irr_vs_r <- ggplot(mc_results, aes(x = r, y = irr_mean)) +
  geom_point(alpha = 0.6, color = "#597cbe") +
  geom_smooth(method = "lm", se = TRUE, color = "grey30") +
  labs(x = "Payer-mix fraction r", y = "Mean incidence rate ratio (25% funding cut)",
       title = sprintf("IRR vs r  (cor = %.3f)", cor_r)) +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "mc_irr_vs_r.png"), plot_irr_vs_r, width = 6, height = 5, dpi = 150)

plot_irr_vs_p0 <- ggplot(mc_results, aes(x = P0, y = irr_mean)) +
  geom_point(alpha = 0.6, color = "#45aF84") +
  geom_smooth(method = "lm", se = TRUE, color = "grey30") +
  labs(x = "Baseline PrEP coverage P0", y = "Mean incidence rate ratio (25% funding cut)",
       title = sprintf("IRR vs P0  (cor = %.3f, expected ~0)", cor_p0)) +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "mc_irr_vs_p0.png"), plot_irr_vs_p0, width = 6, height = 5, dpi = 150)

plot_hist <- ggplot(mc_results, aes(x = irr_mean)) +
  geom_histogram(bins = 30, fill = "#af61a7", alpha = 0.7, color = "white") +
  labs(x = "Mean incidence rate ratio (25% funding cut)", y = "Draws",
       title = "Distribution of IRR under r/P0 parameter uncertainty (N=200)") +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "mc_irr_histogram.png"), plot_hist, width = 6, height = 5, dpi = 150)

cat("\nWrote:\n  output/mc_sensitivity_draws.csv\n  output/mc_sensitivity_summary.txt\n  output/mc_irr_vs_r.png\n  output/mc_irr_vs_p0.png\n  output/mc_irr_histogram.png\n")
