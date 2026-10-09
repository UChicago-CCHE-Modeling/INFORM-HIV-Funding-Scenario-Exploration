# =============================================================================
# Insulation-case response surfaces on the Delta (funding / Medicaid) scale
# =============================================================================
# Mapping (Anna's fig_main footnote):  r_j = Delta_j * (1 - gamma_j), where
# gamma_j is the insulated share of current users from the stakeholder survey,
# used directly (NOT the legacy absolute floor gamma_old / P).
#
# Cases (median +/- 25 percent, Raff's decision 2026-10-01):
#   lower_insulation   gamma_PrEP 0.30   gamma_ART 0.5625
#   median_insulation  gamma_PrEP 0.40   gamma_ART 0.75
#   higher_insulation  gamma_PrEP 0.50   gamma_ART 0.9375
#
# Time: surrogate row 11 = post-intervention year 10 (final year, sim year 32);
# row 6 = year 5 (supplement only). The helper's `tick` is a ONE-BASED ROW INDEX.
#
# Native query point = (-r_ART, -r_PrEP), order (ART, PrEP), r must be <= 0.75.
# Constant z block (set.seed(0) per chunk) => identical draws for every case,
# point and row, so comparisons are matched and IRR(0,0) is exactly 1.
#
# Outputs (../output/):
#   gamma_case_grid.csv       Delta-grid surfaces per case and row
#   gamma_case_points.csv     per-draw IRR at the named scenarios
#   coverage_domain_grid.csv  gamma-independent r-space grid (supplement)
# Run: /usr/local/bin/Rscript 05_response_surfaces/script/export_gamma_cases.R
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

CASES <- list(
  list(name = "lower_insulation",  gamma_prep = 0.30, gamma_art = 0.5625),
  list(name = "median_insulation", gamma_prep = 0.40, gamma_art = 0.75),
  list(name = "higher_insulation", gamma_prep = 0.50, gamma_art = 0.9375))
ROWS <- c(final = 11L, year5 = 6L)       # surrogate row index -> reporting year
YEARS <- c(final = 10L, year5 = 5L)       # years since intervention start
R_MAX <- 0.75                              # surrogate support
N_ROWS <- nrow(gp_incidence_fit[[1]]$Kmat)
stopifnot(N_ROWS == 11L)

NDRAWS_PER <- p$n_samples_per_checkpoint
N_CKPT <- length(gp_incidence_fit)

# named scenarios in Delta space (Delta_prev, Delta_Med); Delta_Med -> ART
SCEN <- data.table(
  scenario = c("P1", "P2", "P3", "elicited"),
  delta_prev = c(0.0, 0.5, 1.0, 0.5),
  delta_med  = c(0.35, 0.0, 0.35, 0.35))

draws_at <- function(r_art, r_prep, row) {
  stopifnot(all(r_art >= 0), all(r_prep >= 0),
            all(r_art <= R_MAX + 1e-12), all(r_prep <= R_MAX + 1e-12))
  native <- cbind(-r_art, -r_prep)
  set.seed(0)
  compute_scenario_draws_crn(native, gp_incidence_fit,
                             n_samples_per_checkpoint = NDRAWS_PER,
                             common_random_numbers = TRUE, tick = row)$irr
}

summ <- function(irr_mat) {
  cbind(mean = rowMeans(irr_mat),
        lower = apply(irr_mat, 1, quantile, probs = p$ci_probs[1], names = FALSE),
        upper = apply(irr_mat, 1, quantile, probs = p$ci_probs[2], names = FALSE))
}

chunked <- function(r_art, r_prep, row, chunk = 500) {
  n <- length(r_art); out <- matrix(NA_real_, n, 3)
  for (k in split(seq_len(n), ceiling(seq_len(n) / chunk)))
    out[k, ] <- summ(draws_at(r_art[k], r_prep[k], row))
  colnames(out) <- c("mean", "lower", "upper"); as.data.frame(out)
}

# ---- Delta grid, 101 x 101 on [0,1] (contains 0.12, 0.25, 0.35, 0.5, 0.8, 1.0)
nodes <- round(seq(0, 1, by = 0.01), 10)
stopifnot(all(c(0.12, 0.25, 0.35, 0.5, 0.8, 1) %in% nodes))
grid <- as.data.table(expand.grid(delta_prev = nodes, delta_med = nodes))

t0 <- Sys.time()
grid_out <- rbindlist(lapply(CASES, function(cs) {
  r_art  <- grid$delta_med  * (1 - cs$gamma_art)
  r_prep <- grid$delta_prev * (1 - cs$gamma_prep)
  rbindlist(lapply(names(ROWS), function(w) {
    cat(sprintf("grid: %s, %s ...\n", cs$name, w))
    s <- chunked(r_art, r_prep, ROWS[[w]])
    data.table(case = cs$name, gamma_prep = cs$gamma_prep, gamma_art = cs$gamma_art,
               when = w, row_index = ROWS[[w]], years_since_start = YEARS[[w]],
               grid, r_art = r_art, r_prep = r_prep, s)
  }))
}))
cat(sprintf("grid done in %.1f s\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

stopifnot(all(is.finite(grid_out$mean)), all(grid_out$mean > 0),
          all(grid_out$lower <= grid_out$mean + 1e-9),
          all(grid_out$mean <= grid_out$upper + 1e-9))
z <- grid_out[delta_prev == 0 & delta_med == 0]
stopifnot(all(abs(z$mean - 1) < 1e-10))
cat("max r_art / r_prep per case:\n")
print(grid_out[when == "final", .(max_r_art = max(r_art), max_r_prep = max(r_prep)), by = case])
write.csv(grid_out, file.path(STAGE_DIR, "output", "gamma_case_grid.csv"), row.names = FALSE)

# ---- per-draw IRR at named scenarios ----------------------------------------
pts_out <- rbindlist(lapply(CASES, function(cs) {
  r_art  <- SCEN$delta_med  * (1 - cs$gamma_art)
  r_prep <- SCEN$delta_prev * (1 - cs$gamma_prep)
  rbindlist(lapply(names(ROWS), function(w) {
    m <- draws_at(r_art, r_prep, ROWS[[w]])
    rbindlist(lapply(seq_len(nrow(SCEN)), function(i)
      data.table(case = cs$name, when = w, row_index = ROWS[[w]],
                 scenario = SCEN$scenario[i], delta_prev = SCEN$delta_prev[i],
                 delta_med = SCEN$delta_med[i], r_prep = r_prep[i], r_art = r_art[i],
                 draw = seq_len(ncol(m)), checkpoint = ceiling(seq_len(ncol(m)) / NDRAWS_PER),
                 irr = m[i, ])))
  }))
}))
stopifnot(all(is.finite(pts_out$irr)), all(pts_out$irr > 0))
write.csv(pts_out, file.path(STAGE_DIR, "output", "gamma_case_points.csv"), row.names = FALSE)

# ---- gamma-independent coverage-domain grid (r on [0, 0.75]) ----------------
rn <- round(seq(0, R_MAX, length.out = 101), 10)
rg <- as.data.table(expand.grid(r_prep = rn, r_art = rn))
cov_out <- rbindlist(lapply(names(ROWS), function(w) {
  cat(sprintf("coverage domain: %s ...\n", w))
  s <- chunked(rg$r_art, rg$r_prep, ROWS[[w]])
  data.table(when = w, row_index = ROWS[[w]], rg, s)
}))
stopifnot(all(is.finite(cov_out$mean)))
write.csv(cov_out, file.path(STAGE_DIR, "output", "coverage_domain_grid.csv"), row.names = FALSE)
cat("done\n")
