# =============================================================================
# Export the response-surface grid for plotting  (REWRITTEN 2026-09-23)
# =============================================================================
# WHY THIS FILE WAS REWRITTEN
#   The previous version (i) computed Panel B by the IDENTICAL code path as
#   Panel A -- the cost factors were only written into metadata columns, so the
#   "funding" panel had no insulation effect at all and came out bitwise equal
#   to the use-reduction panel -- and (ii) queried the surrogate with
#   independent (non-CRN) draws while setting the seed after the draws, so it
#   had neither common random numbers nor the Figure-1 sampling convention.
#
# CONVENTIONS (explicit; see RAF_temp/generated_response_surfaces/VALIDATION.md)
#
#   Panel A  "Intervention use reduction": Delta_j IS the coverage reduction.
#            native query point = -Delta_j                     (identity mapping)
#
#   Panel B  "Government funding reduction": Delta_j is a FUNDING cut, mapped to
#            a coverage reduction through the committed linear cost model
#            (paper Eqs. 1-2):
#               coverage reduction = Delta_j * (1 - gamma_Anna)
#               gamma_Anna = gamma_old / P_baseline
#            where gamma_old is the ABSOLUTE substitution floor stored in
#            cost_params.yml (gamma_ART = 0.15, gamma_PrEP = 0.28) and
#            gamma_Anna is the insulated fraction of CURRENT USERS
#            (= 0.245287 ART, 0.779032 PrEP at the committed parameters).
#            native query point = -(Delta_j * (1 - gamma_Anna))
#            => Panel B is a SECOND surrogate query at the mapped coordinates,
#               never a relabelling of Panel A's predictions.
#
#   Time summary: BOTH are exported.
#            "horizon_mean" = mean over ticks of the per-draw ratio (tick=NULL)
#                             -- the convention Figure 1 uses.
#            "tick10"       = the single horizon tick (10y) -- the convention the
#                             scenario table uses.
#            Figures must state which one they show.
#
#   Draws: common random numbers (CRN) with a CONSTANT standard-normal block
#          across every chunk, panel and case (set.seed before each chunk), so
#          every point shares the same z, the baseline is identical everywhere,
#          and IRR(0,0) is exactly 1. Stronger than the Figure-1 CRN trick and
#          it makes Panel B / Panel A ratios noise-cancelled.
#
# Run from 04_cost_mapping/script:
#   /usr/local/bin/Rscript export_response_surface_grid.R
# =============================================================================

library(data.table)
library(hetGP)

load("../../03_intervention_scenario_surrogate/output/surrogate.Rdata")

source("../R/parameters.R")
source("../R/irr_common_random_numbers.R")

p <- cost_mapping_params()

# ---- gamma conventions, printed so the CSV can never be misread -------------
k_art  <- (p$P_ART_baseline  - p$gamma_ART)  / p$P_ART_baseline
k_prep <- (p$P_PrEP_baseline - p$gamma_PrEP) / p$P_PrEP_baseline

cat("=== parameters ===\n")
cat(sprintf("  P_ART_baseline = %.10f   gamma_ART (absolute floor) = %.4f   gamma_Anna,ART = %.6f\n",
            p$P_ART_baseline, p$gamma_ART, p$gamma_ART / p$P_ART_baseline))
cat(sprintf("  P_PrEP_baseline= %.10f   gamma_PrEP(absolute floor) = %.4f   gamma_Anna,PrEP= %.6f\n",
            p$P_PrEP_baseline, p$gamma_PrEP, p$gamma_PrEP / p$P_PrEP_baseline))
cat(sprintf("  coverage-loss multiplier per unit funding cut: ART = %.10f  PrEP = %.10f\n",
            k_art, k_prep))
cat("  (0.7547 / 0.2210 are COVERAGE-LOSS MULTIPLIERS (1 - gamma_Anna), NOT insulated shares.)\n\n")

# ---- gamma cases: add entries here to vary the insulation assumption --------
# Only the committed case is enabled by default. To vary the insulation
# assumption add cases using the SAME denominator (the insulated fraction of
# current users -- NOT the absolute floor). Example, the legacy PrEP sensitivity
# values (recalibrate.R semantics: r = insulated share, pass-through = 1 - r):
#   list(name = "prep_insulation_0.25", k_prep = 1 - 0.25),
#   list(name = "prep_insulation_0.44", k_prep = 1 - sqrt(0.25 * 0.78)),
#   list(name = "prep_insulation_0.78", k_prep = 1 - 0.78))
# Any side left unset falls back to the committed multiplier, and the case name
# is written into the CSV so a figure can never silently mix denominators.
GAMMA_CASES <- list(list(name = "committed", k_art = k_art, k_prep = k_prep))

# ---- grid -------------------------------------------------------------------
# 101 evenly spaced nodes over the surrogate's domain [0, 0.75] (step 0.0075)
# PLUS the policy levels Anna asked to highlight (0.10 / 0.25 / 0.40). The even
# grid alone does NOT contain any policy level (nearest node to 0.40 is 0.3975),
# which would put "40%" markers 0.0025 away from the actual cut and make
# regression against Figure 1's 40% scenarios compare the wrong number.
ngrid <- p$ngrid
POLICY_NODES <- c(0.10, 0.25, 0.40)
nodes <- sort(unique(c(seq(0, 0.75, length.out = ngrid), POLICY_NODES)))
grid  <- expand.grid(art_red = nodes, prep_red = nodes)
stopifnot(all(POLICY_NODES %in% nodes))
cat(sprintf("grid: %d x %d = %d nodes per panel, Delta in [%.4f, %.4f]; policy nodes %s pinned exactly\n\n",
            length(nodes), length(nodes), nrow(grid), min(nodes), max(nodes),
            paste(POLICY_NODES, collapse = "/")))

# ---- chunked, CRN-consistent surface query ----------------------------------
.irr_chunked <- function(native, tick, chunk = 500) {
  n   <- nrow(native)
  out <- matrix(NA_real_, n, 3)
  for (k in split(seq_len(n), ceiling(seq_len(n) / chunk))) {
    set.seed(0)   # constant z for every chunk -> constant z for the whole surface
    d <- compute_scenario_draws_crn(
      native[k, , drop = FALSE], gp_incidence_fit,
      n_samples_per_checkpoint = p$n_samples_per_checkpoint,
      common_random_numbers = TRUE, tick = tick)
    s <- summarise_draws(d$irr, probs = p$ci_probs)
    out[k, ] <- cbind(s$mean, s$ci_lower, s$ci_upper)
  }
  colnames(out) <- c("mean", "lower", "upper")
  as.data.frame(out)
}

surface_for <- function(case, tick, time_summary) {
  # Panel A: identity mapping (Delta_j = coverage reduction)
  a <- .irr_chunked(cbind(-grid$art_red, -grid$prep_red), tick)
  # Panel B: mapped to coverage through the cost model, then queried
  b <- .irr_chunked(cbind(-grid$art_red * case$k_art,
                          -grid$prep_red * case$k_prep), tick)
  rbind(
    data.frame(panel = "A. Intervention use reduction", case = case$name,
               time_summary = time_summary, art_red = grid$art_red,
               prep_red = grid$prep_red, art_cov = grid$art_red,
               prep_cov = grid$prep_red, a),
    data.frame(panel = "B. Government funding reduction", case = case$name,
               time_summary = time_summary, art_red = grid$art_red,
               prep_red = grid$prep_red, art_cov = grid$art_red * case$k_art,
               prep_cov = grid$prep_red * case$k_prep, b)
  )
}

t0 <- Sys.time()
out <- rbindlist(lapply(GAMMA_CASES, function(cs)
  rbind(surface_for(cs, NULL, "horizon_mean"),
        surface_for(cs, p$horizon_tick, "tick10"))))
cat(sprintf("surface computed in %.1f s\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))

# ---- self-checks that would have caught the previous defects ----------------
A <- out[panel == "A. Intervention use reduction" & time_summary == "horizon_mean"]
B <- out[panel == "B. Government funding reduction" & time_summary == "horizon_mean"]
stopifnot(!isTRUE(all.equal(A$mean, B$mean)))          # the panels must DIFFER
stopifnot(all(A$mean >= 1 - 1e-12), all(B$mean >= 1 - 1e-12))
stopifnot(all(A$lower <= A$mean + 1e-9), all(A$mean <= A$upper + 1e-9))

pick <- function(d, a, p_) d[which.min(abs(d$art_red - a) + abs(d$prep_red - p_)), ]
pa <- pick(A, 0.40, 0.10); pb <- pick(A, 0.10, 0.40)
cat("\n=== orientation check (asymmetric nodes; equal-cut nodes cannot see an axis swap) ===\n")
cat(sprintf("  A(ART=%.4f, PrEP=%.4f) = %.6f  vs  A(ART=%.4f, PrEP=%.4f) = %.6f  (must differ)\n",
            pa$art_red, pa$prep_red, pa$mean, pb$art_red, pb$prep_red, pb$mean))
stopifnot(abs(pa$mean - pb$mean) > 0.05)
zero <- A[which.min(A$art_red + A$prep_red), ]
cat(sprintf("  self-ratio at (0,0): A = %.10f (CRN => exactly 1)\n", zero$mean))
aa <- pick(A, 0.3975, 0.3975); bb <- pick(B, 0.3975, 0.3975)
cat(sprintf("  at Delta=(0.3975,0.3975): A = %.6f, B = %.6f, B/A = %.6f (insulation buffer visible)\n",
            aa$mean, bb$mean, bb$mean / aa$mean))

# ---- write ------------------------------------------------------------------
output_path <- "../output/response_surface_grid_full.csv"
write.csv(out, output_path, row.names = FALSE)
cat(sprintf("\nwrote %s  (%d rows)\n", output_path, nrow(out)))
cat(sprintf("  Panel A horizon_mean range: %.4f - %.4f\n", min(A$mean), max(A$mean)))
cat(sprintf("  Panel B horizon_mean range: %.4f - %.4f\n",
            min(B$mean), max(B$mean)))
