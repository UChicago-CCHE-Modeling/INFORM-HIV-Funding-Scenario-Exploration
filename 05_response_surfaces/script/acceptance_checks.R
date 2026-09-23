# =============================================================================
# acceptance_checks.R -- independent R side of the Codex acceptance criteria
# -----------------------------------------------------------------------------
# Required demonstrations (see RAF_temp/CODEX_INDEPENDENT_REVIEW.md):
#   (A) agreement with existing results under MATCHING parameter and time-summary
#       conventions -> reproduce the committed forest CSV exactly;
#   (B) reference values at ASYMMETRIC ART/PrEP points, and at the MAPPED
#       coverage coordinates for the funding panel, so the Python surface can be
#       checked against an independent query.
# Read-only w.r.t. the deliverables; writes /tmp only.
# =============================================================================

COST_DIR <- "/Users/r_vardavas/Documents/Projects_2026/HIV-UChicago/INFORM-HIV-Funding-Scenario-Exploration/04_cost_mapping"
setwd(file.path(COST_DIR, "script"))

# plot_forest.R needs randplot (not installed here); stub that one package only.
library <- function(package, ...) {
  nm <- as.character(substitute(package))
  if (identical(nm, "randplot")) return(invisible(NULL))
  base::library(nm, character.only = TRUE, ...)
}

source("../R/parameters.R")
source("../R/irr_common_random_numbers.R")
source("../R/plot_forest.R")
load("../../03_intervention_scenario_surrogate/output/surrogate.Rdata")

p <- cost_mapping_params()
k_art  <- (p$P_ART_baseline  - p$gamma_ART)  / p$P_ART_baseline
k_prep <- (p$P_PrEP_baseline - p$gamma_PrEP) / p$P_PrEP_baseline

# ---- (A) exact reproduction of the committed forest CSV ---------------------
set.seed(0)
d_use <- .forest_panel_data(gp_incidence_fit, art_cov_per_unit = 1,
                            prep_cov_per_unit = 1,
                            reduction_levels = c(0.10, 0.25, 0.40),
                            n_samples_per_checkpoint = p$n_samples_per_checkpoint,
                            common_random_numbers = TRUE, ci_probs = p$ci_probs,
                            panel_label = "A. Intervention use reduction")
d_fund <- .forest_panel_data(gp_incidence_fit, art_cov_per_unit = k_art,
                             prep_cov_per_unit = k_prep,
                             reduction_levels = c(0.10, 0.25, 0.40),
                             n_samples_per_checkpoint = p$n_samples_per_checkpoint,
                             common_random_numbers = TRUE, ci_probs = p$ci_probs,
                             panel_label = "B. Government funding reduction")

repro <- rbind(as.data.frame(d_use), as.data.frame(d_fund))
committed <- read.csv("../output/forest_incidence_risk_ratio_funding.csv")
key_r <- paste(repro$panel, repro$scenario, repro$reduction)
key_c <- paste(committed$panel, committed$scenario, committed$reduction)
m <- match(key_c, key_r)
stopifnot(!any(is.na(m)))

rep_maxdiff <- function(col) max(abs(committed[[col]] - repro[[col]][m]))
forest_checks <- c(mean = rep_maxdiff("mean"), lower = rep_maxdiff("lower"),
                   upper = rep_maxdiff("upper"))
both40 <- committed$panel == "B. Government funding reduction" &
  committed$scenario == "Reduce both PrEP and ART" & committed$reduction == 0.40
use40 <- committed$panel == "A. Intervention use reduction" &
  committed$scenario == "Reduce both PrEP and ART" & committed$reduction == 0.40

cat("=== (A) forest CSV reproduced under matching conventions (horizon-mean, CRN) ===\n")
cat(sprintf("  max|diff| vs committed CSV: mean %.3e  lower %.3e  upper %.3e\n",
            forest_checks["mean"], forest_checks["lower"], forest_checks["upper"]))
cat(sprintf("  committed A(40%%,40%%) = %.10f   B(40%%,40%%) = %.10f   B/A = %.6f\n",
            committed$mean[use40], committed$mean[both40],
            committed$mean[both40] / committed$mean[use40]))

# ---- (B) asymmetric points, and the SAME points mapped to coverage ----------
# native coordinates are (ART, PrEP) reductions, negative to the surrogate
probe <- rbind(c(-0.40, -0.10),   # ART 40%, PrEP 10%
               c(-0.10, -0.40),   # ART 10%, PrEP 40%
               c(-0.40, -0.25),
               c(-0.25, -0.40),
               c(-0.40, -0.40))
probe_map <- cbind(probe[, 1] * k_art, probe[, 2] * k_prep)

set.seed(0)
a <- summarise_draws(compute_scenario_draws_crn(
  probe, gp_incidence_fit, n_samples_per_checkpoint = p$n_samples_per_checkpoint,
  common_random_numbers = TRUE, tick = NULL)$irr, probs = p$ci_probs)
set.seed(0)
b <- summarise_draws(compute_scenario_draws_crn(
  probe_map, gp_incidence_fit, n_samples_per_checkpoint = p$n_samples_per_checkpoint,
  common_random_numbers = TRUE, tick = NULL)$irr, probs = p$ci_probs)

ref <- data.frame(art_red = -probe[, 1], prep_red = -probe[, 2],
                  art_cov = -probe_map[, 1], prep_cov = -probe_map[, 2],
                  A_ref = a$mean, B_ref = b$mean)
write.csv(ref, "/tmp/accept_r_values.csv", row.names = FALSE)

cat("\n=== (B) independent R reference (horizon-mean, CRN) ===\n")
print(ref, row.names = FALSE, digits = 8)
cat(sprintf("\n  A at (ART .40, PrEP .10) = %.6f  vs  A at (ART .10, PrEP .40) = %.6f  (asymmetric)\n",
            ref$A_ref[1], ref$A_ref[2]))
cat(sprintf("  B at mapped (ART %.4f, PrEP %.4f) = %.6f\n",
            ref$art_cov[1], ref$prep_cov[1], ref$B_ref[1]))
cat("acceptance_checks.R done.\n")
