# =============================================================================
# Surrogate-only IRR surface over the coverage-reduction domain (r-space)
# =============================================================================
# Contains NO survey information. Writes output/coverage_domain_grid.csv:
# 101 x 101 nodes over r_PrEP, r_VS in [0, 0.75] (the surrogate's support), for
#   when = "final": surrogate row 11 = post-intervention year 10 (final year)
#   when = "year5": surrogate row 6  = year 5
# mean, 95% interval over 100 checkpoints x 50 draws, constant z block
# (set.seed(0) per chunk) so IRR(0,0) is exactly 1. `tick` in the helper is a
# ONE-BASED ROW INDEX into the 11-row trajectory.
# Run: /usr/local/bin/Rscript 05_response_surfaces/script/export_coverage_grid.R
# =============================================================================

args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grepl("^--file=", args)][1])
SCRIPT_DIR <- dirname(normalizePath(script_file))
STAGE_DIR  <- dirname(SCRIPT_DIR)
REPO_DIR   <- dirname(STAGE_DIR)
setwd(file.path(REPO_DIR, "04_cost_mapping", "script"))

library(data.table)
library(hetGP)
load("../../03_intervention_scenario_surrogate/output/surrogate.Rdata")
source("../R/parameters.R")
source("../R/irr_common_random_numbers.R")
p <- cost_mapping_params()

ROWS <- c(final = 11L, year5 = 6L)
R_MAX <- 0.75
stopifnot(nrow(gp_incidence_fit[[1]]$Kmat) == 11L)

draws_at <- function(r_art, r_prep, row) {
  stopifnot(all(r_art >= 0), all(r_prep >= 0), all(r_art <= R_MAX + 1e-12), all(r_prep <= R_MAX + 1e-12))
  set.seed(0)
  compute_scenario_draws_crn(cbind(-r_art, -r_prep), gp_incidence_fit,
                             n_samples_per_checkpoint = p$n_samples_per_checkpoint,
                             common_random_numbers = TRUE, tick = row)$irr
}
summ <- function(m) cbind(mean = rowMeans(m),
  lower = apply(m, 1, quantile, probs = p$ci_probs[1], names = FALSE),
  upper = apply(m, 1, quantile, probs = p$ci_probs[2], names = FALSE))
chunked <- function(r_art, r_prep, row, chunk = 500) {
  n <- length(r_art); o <- matrix(NA_real_, n, 3)
  for (k in split(seq_len(n), ceiling(seq_len(n) / chunk))) o[k, ] <- summ(draws_at(r_art[k], r_prep[k], row))
  colnames(o) <- c("mean", "lower", "upper"); as.data.frame(o)
}

rn <- round(seq(0, R_MAX, length.out = 101), 10)
rg <- as.data.table(expand.grid(r_prep = rn, r_art = rn))
out <- rbindlist(lapply(names(ROWS), function(w) {
  cat(sprintf("coverage domain: %s ...\n", w))
  data.table(when = w, row_index = ROWS[[w]], rg, chunked(rg$r_art, rg$r_prep, ROWS[[w]]))
}))
stopifnot(all(is.finite(out$mean)), all(out$lower <= out$mean + 1e-9), all(out$mean <= out$upper + 1e-9),
          all(abs(out[r_prep == 0 & r_art == 0]$mean - 1) < 1e-10))
fwrite(out, file.path(STAGE_DIR, "output", "coverage_domain_grid.csv"))
cat("wrote output/coverage_domain_grid.csv\n")
