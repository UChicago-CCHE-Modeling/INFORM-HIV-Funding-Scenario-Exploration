# =============================================================================
# Point-value r sensitivity for the actual paper (not a Monte Carlo sweep)
# -----------------------------------------------------------------------------
# Follows the guidance from the follow-up meeting documented in
# tmp/Parameter Uncertainty and Sensitivity Analysis_otter_ai_transcript.txt:
# Pedro (statistician) explicitly rejected treating r/P0 as a uniform/beta
# random sweep merged into the surrogate's own well-characterized predictive
# intervals. Instead: ONE baseline r in the main text, lower/upper bound
# scenarios as point-value comparisons in the supplementary appendix. P0 gets
# NO sensitivity table -- only an algebraic explanation (it cancels out; see
# R/recalibrate.R and the existing exploratory report for the derivation).
#
# R_LOWER = 0.25 is kept as an EXPLICITLY UNVALIDATED PLACEHOLDER. A targeted
# literature search for a principled lower bound (population subgroup with
# lowest PrEP-payer private-insurance share) came back data-sparse; the one
# population-matched point found (Kelley et al., Atlanta young Black MSM PrEP
# cohort, PMC6781266: ~54% privately insured) does NOT support a bound this
# low. Kept anyway per team decision, pending AFC's actual number -- flagged
# here, in the generated macros, and in the paper prose itself.
#
# R_BASELINE = geometric mean of R_LOWER and R_UPPER (Pedro's rule-of-thumb
# for "upper bound / lower bound -> geometric mean" when a parameter is
# poorly characterized) REPLACES r=0.78 as the paper's main-text scenario;
# r=0.78 moves to the appendix as the "upper bound" comparison.
#
# Run with 04_cost_mapping/sensitivity/ as the working directory:
#   /usr/local/bin/Rscript run_paper_sensitivity.R
# =============================================================================

library(hetGP)      # so predict() dispatches on the stored surrogate GPs
library(kableExtra)  # create_funding_scenario_table()'s LaTeX rendering
library(dplyr)

source("../R/parameters.R")
source("../R/irr_common_random_numbers.R")
source("../R/table_scenarios.R")   # create_funding_scenario_table() -- no randplot dependency
source("R/recalibrate.R")

set.seed(42)

out_dir <- "output"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---- Load surrogate + base parameters --------------------------------------
load("../../03_intervention_scenario_surrogate/output/surrogate.Rdata")
p_base <- cost_mapping_params(yml = "../cost_params.yml")

R_UPPER    <- 0.78   # current committed value; becomes the appendix "upper bound"
R_LOWER    <- 0.25   # UNVALIDATED PLACEHOLDER -- see header note; pending AFC
R_BASELINE <- sqrt(R_LOWER * R_UPPER)  # geometric mean; new main-text scenario

P0 <- p_base$P_PrEP_baseline  # held fixed across all three r values (never varied)

lower    <- recalibrate_prep(R_LOWER,    P0, p_base)
baseline <- recalibrate_prep(R_BASELINE, P0, p_base)
upper    <- list(r = R_UPPER, P0 = P0, gamma_PrEP = p_base$gamma_PrEP,
                  beta_PrEP = p_base$beta_PrEP,
                  prep_cov_per_funding = (P0 - p_base$gamma_PrEP) / P0)
stopifnot(isTRUE(all.equal(upper$r, R_UPPER)))

cat(sprintf("Lower    r=%.4f  gamma_PrEP=%.4f  beta_PrEP=%.4f  prep_cov_per_funding=%.4f\n",
            lower$r, lower$gamma_PrEP, lower$beta_PrEP, lower$prep_cov_per_funding))
cat(sprintf("Baseline r=%.4f  gamma_PrEP=%.4f  beta_PrEP=%.4f  prep_cov_per_funding=%.4f\n",
            baseline$r, baseline$gamma_PrEP, baseline$beta_PrEP, baseline$prep_cov_per_funding))
cat(sprintf("Upper    r=%.4f  gamma_PrEP=%.4f  beta_PrEP=%.4f  prep_cov_per_funding=%.4f\n",
            upper$r, upper$gamma_PrEP, upper$beta_PrEP, upper$prep_cov_per_funding))

# Algebraic check: (P0-gamma)/P0 == 1-r, independent of P0. Exact for lower/
# baseline (freshly derived via recalibrate_prep()). NOT exact for "upper",
# because it reuses the ALREADY-COMMITTED gamma_PrEP=0.28 literal from
# cost_params.yml, which was itself rounded to 2 decimals rather than stored
# as the exact product 0.78*P0 (0.78*0.3594... = 0.2803..., not 0.28) -- a
# pre-existing rounding artifact in the committed value, not a bug here.
for (x in list(lower, baseline)) {
  stopifnot(isTRUE(all.equal(x$prep_cov_per_funding, 1 - x$r)))
}
upper_algebra_diff <- abs(upper$prep_cov_per_funding - (1 - upper$r))
cat(sprintf("\nAlgebraic check OK for lower/baseline. Upper (uses committed rounded gamma_PrEP=%.2f): (P0-gamma)/P0=%.4f vs 1-r=%.4f (diff=%.4f, expected from rounding).\n\n",
            p_base$gamma_PrEP, upper$prep_cov_per_funding, 1 - upper$r, upper_algebra_diff))

art_cov_per_funding <- (p_base$P_ART_baseline - p_base$gamma_ART) / p_base$P_ART_baseline
reduction_levels <- c(0.10, 0.25, 0.40)
n_samples_per_checkpoint <- p_base$n_samples_per_checkpoint
ci_probs <- p_base$ci_probs

# =============================================================================
# Part A: appendix comparison table (lower / baseline / upper), point values
# =============================================================================
build_scenarios <- function(prep_cov_per_funding, label) {
  s <- rbind(
    data.frame(scenario = "Reduce PrEP only",         art_red = 0,                prep_red = reduction_levels),
    data.frame(scenario = "Reduce ART only",          art_red = reduction_levels, prep_red = 0),
    data.frame(scenario = "Reduce both PrEP and ART", art_red = reduction_levels, prep_red = reduction_levels)
  )
  s$reduction <- pmax(s$art_red, s$prep_red)
  s$art_cov  <- s$art_red  * art_cov_per_funding
  s$prep_cov <- s$prep_red * prep_cov_per_funding
  s$r_scenario <- label
  s
}

data_lower    <- build_scenarios(lower$prep_cov_per_funding,    "Lower bound (r=25%)")
data_baseline <- build_scenarios(baseline$prep_cov_per_funding, sprintf("Baseline (r=%.1f%%)", R_BASELINE * 100))
data_upper    <- build_scenarios(upper$prep_cov_per_funding,    "Upper bound (r=78%)")

appendix_scenarios <- rbind(data_lower, data_baseline, data_upper)
newX_native <- as.matrix(cbind(-appendix_scenarios$art_cov, -appendix_scenarios$prep_cov))
draws <- compute_scenario_draws_crn(
  newX_native, gp_incidence_fit,
  n_samples_per_checkpoint = n_samples_per_checkpoint,
  common_random_numbers = TRUE, tick = NULL)
irr_summary <- summarise_draws(draws$irr, probs = ci_probs)
appendix_scenarios$irr_mean     <- irr_summary$mean
appendix_scenarios$irr_ci_lower <- irr_summary$ci_lower
appendix_scenarios$irr_ci_upper <- irr_summary$ci_upper
appendix_scenarios$r_scenario <- factor(appendix_scenarios$r_scenario, levels = unique(appendix_scenarios$r_scenario))
appendix_scenarios$scenario <- factor(appendix_scenarios$scenario,
  levels = c("Reduce PrEP only", "Reduce ART only", "Reduce both PrEP and ART"))

write.csv(appendix_scenarios, file.path(out_dir, "paper_sensitivity_table.csv"), row.names = FALSE)

# ---- A small, ready-to-\input LaTeX table for the appendix ----------------
# Styled as a simple booktabs table (matches tab:supp-param-estimates'
# plainness) rather than kableExtra's landscape/scale_down machinery used for
# the main Table 1 -- this is a compact 9-row comparison, not a wide table.
fmt_irr <- function(m, lo, hi) sprintf("%.2f [%.2f, %.2f]", m, lo, hi)
tbl <- appendix_scenarios[order(appendix_scenarios$scenario, appendix_scenarios$reduction), ]
tex_lines <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Sensitivity of the incidence rate ratio (mean [95\\% CrI]) to the PrEP payer-mix parameter $r$, at lower/baseline/upper point values.}",
  "\\label{tab:supp-r-sensitivity}",
  "\\resizebox{\\ifdim\\width>\\linewidth\\linewidth\\else\\width\\fi}{!}{%",
  "\\begin{tabular}{llccc}",
  "\\toprule",
  sprintf("Scenario & Funding cut & Lower ($r=25\\%%$) & Baseline ($r=%.1f\\%%$) & Upper ($r=78\\%%$) \\\\", R_BASELINE * 100),
  "\\midrule"
)
for (scen in levels(tbl$scenario)) {
  rows <- tbl[tbl$scenario == scen, ]
  for (lvl in reduction_levels) {
    lo_row  <- rows[rows$reduction == lvl & rows$r_scenario == levels(rows$r_scenario)[1], ]
    bas_row <- rows[rows$reduction == lvl & rows$r_scenario == levels(rows$r_scenario)[2], ]
    up_row  <- rows[rows$reduction == lvl & rows$r_scenario == levels(rows$r_scenario)[3], ]
    label <- if (lvl == reduction_levels[1]) paste0("\\multirow{3}{*}{", scen, "}") else ""
    tex_lines <- c(tex_lines, sprintf(
      "%s & %.0f\\%% & %s & %s & %s \\\\",
      label, lvl * 100,
      fmt_irr(lo_row$irr_mean, lo_row$irr_ci_lower, lo_row$irr_ci_upper),
      fmt_irr(bas_row$irr_mean, bas_row$irr_ci_lower, bas_row$irr_ci_upper),
      fmt_irr(up_row$irr_mean, up_row$irr_ci_lower, up_row$irr_ci_upper)))
  }
  tex_lines <- c(tex_lines, "\\addlinespace")
}
tex_lines <- c(tex_lines, "\\bottomrule", "\\end{tabular}}", "\\end{table}")
writeLines(tex_lines, file.path(out_dir, "paper_sensitivity_table.tex"))

# =============================================================================
# Part B: REAL main-text Table 1 replacement, at the new baseline r
# -----------------------------------------------------------------------------
# create_funding_scenario_table() has NO randplot dependency (kableExtra +
# dplyr only) -- this is not a placeholder, it is the actual function the
# paper's real pipeline uses, just called with the new baseline's
# prep_cov_per_funding.
# =============================================================================
table_dir <- file.path(out_dir, "paper_table1_baseline")
table_baseline <- create_funding_scenario_table(
  gp_incidence_fit = gp_incidence_fit,
  art_cov_per_funding = art_cov_per_funding,
  prep_cov_per_funding = baseline$prep_cov_per_funding,
  reduction_levels = reduction_levels,
  n_samples_per_checkpoint = n_samples_per_checkpoint,
  horizon_tick = p_base$horizon_tick,
  ci_probs = c(0.025, 0.975),
  common_random_numbers = TRUE,
  output_dir = table_dir)
write.csv(table_baseline, file.path(table_dir, "funding_scenarios_table.csv"), row.names = FALSE)

# Cross-check: regenerate the SAME function at r=0.78 (upper/original) and
# confirm it reproduces the numbers already committed in paper.tex's existing
# tab:funding_scenarios (e.g. "both" row: incidence 7.87, IRR 1.88).
table_dir_upper <- file.path(out_dir, "paper_table1_upper_crosscheck")
table_upper <- create_funding_scenario_table(
  gp_incidence_fit = gp_incidence_fit,
  art_cov_per_funding = art_cov_per_funding,
  prep_cov_per_funding = upper$prep_cov_per_funding,
  reduction_levels = reduction_levels,
  n_samples_per_checkpoint = n_samples_per_checkpoint,
  horizon_tick = p_base$horizon_tick,
  ci_probs = c(0.025, 0.975),
  common_random_numbers = TRUE,
  output_dir = table_dir_upper)
write.csv(table_upper, file.path(table_dir_upper, "funding_scenarios_table.csv"), row.names = FALSE)

cat("\n--- Cross-check: r=0.78 regeneration vs paper.tex's existing committed table ---\n")
cat("(compare the '40% cut for both' Block B row to paper.tex's '7.87 ... 1.88')\n")
print(table_upper[nrow(table_upper), ])

# =============================================================================
# Part C: macros file for the paper prose (r values, new baseline gamma/beta,
# and the specific "40% both, funding" sentence numbers that change)
# =============================================================================
# Extract straight from table_baseline's own "40% cut for both" row (Block B,
# last row) rather than re-querying the surrogate separately -- a fresh call
# would consume a different RNG state than create_funding_scenario_table()'s
# internal set.seed(0), giving a slightly different (though not materially
# different) Monte Carlo draw. Parsing keeps the macros bit-consistent with
# what the table itself displays -- single source of truth.
parse_ci <- function(s) {
  m <- regmatches(s, regexec("^(-?[0-9.]+) \\[(-?[0-9.]+), (-?[0-9.]+)\\]$", s))[[1]]
  list(mean = as.numeric(m[2]), ci_lower = as.numeric(m[3]), ci_upper = as.numeric(m[4]))
}
both_row <- table_baseline[nrow(table_baseline), ]  # "40% cut for both", Block B, last row
both_inc <- parse_ci(both_row[["Mean HIV Incidence per 100 p.y. [95\\% CrI]"]])
both_abs <- parse_ci(both_row[["Absolute Increase per 100 p.y. [95\\% CrI]"]])
both_irr <- parse_ci(both_row[["Incidence Rate Ratio [95\\% CrI]"]])

# PrEP-only funding rows at the new baseline (for the appendix table / any
# future prose), 10/25/40%.
prep_only_pct_use <- reduction_levels * baseline$prep_cov_per_funding * 100

# "Impacts of funding reductions... accumulate over time" paragraph: both
# ART+PrEP funding cut equally at 10/20/40%, year-1/year-2/year-10 IRR.
traj_levels <- c(0.10, 0.20, 0.40)
traj_irr <- sapply(traj_levels, function(lvl) {
  pt <- matrix(c(-lvl * art_cov_per_funding, -lvl * baseline$prep_cov_per_funding), nrow = 1)
  tr <- compute_scenario_trajectory_draws_crn(pt, gp_incidence_fit,
                                               n_samples_per_checkpoint = n_samples_per_checkpoint,
                                               common_random_numbers = TRUE)
  s <- summarise_draws(tr$irr[[1]], probs = ci_probs)
  c(year1 = s$mean[2], year2 = s$mean[3], year10 = s$mean[11])  # ticks are 1-indexed, year0=tick1
})
colnames(traj_irr) <- sprintf("%.0f", traj_levels * 100)

# Qualitative check (not a macro, printed for the record): is the "largest
# year-10 IRR when ART funding is cut more than PrEP" claim (paper.tex,
# results section) still true at the new baseline? A 50%-ART-only vs.
# 50%-PrEP-only funding cut, year-10 IRR:
qual_pt <- rbind(c(-0.50 * art_cov_per_funding, 0), c(0, -0.50 * baseline$prep_cov_per_funding))
qual_draws <- compute_scenario_draws_crn(qual_pt, gp_incidence_fit,
                                          n_samples_per_checkpoint = n_samples_per_checkpoint,
                                          common_random_numbers = TRUE, tick = p_base$horizon_tick)
cat(sprintf("\nQualitative check -- 50%% ART-only funding cut year-10 IRR: %.3f; 50%% PrEP-only: %.3f (ART should stay higher for the 'largest when ART cut more' claim to hold)\n",
            mean(qual_draws$irr[1, ]), mean(qual_draws$irr[2, ])))

macro <- function(name, value) sprintf("\\newcommand{\\%s}{%s}", name, value)
macros <- c(
  "% Generated by 04_cost_mapping/sensitivity/run_paper_sensitivity.R -- do not edit by hand.",
  "% R_LOWER=0.25 is an UNVALIDATED PLACEHOLDER pending AFC confirmation; see memo.",
  macro("rLower",    sprintf("%.0f\\%%", R_LOWER * 100)),
  macro("rBaseline", sprintf("%.1f\\%%", R_BASELINE * 100)),
  macro("rUpper",    sprintf("%.0f\\%%", R_UPPER * 100)),
  macro("rLowerDecimal",    sprintf("%.2f", R_LOWER)),
  macro("rBaselineDecimal", sprintf("%.2f", R_BASELINE)),
  macro("rUpperDecimal",    sprintf("%.2f", R_UPPER)),
  macro("PZeroPrEP", sprintf("%.1f\\%%", P0 * 100)),
  macro("gammaPrEPBaseline", sprintf("%.1f", baseline$gamma_PrEP * 100)),
  macro("gammaPrEPBaselineDecimal", sprintf("%.2f", baseline$gamma_PrEP)),
  macro("gammaPrEPLowerDecimal",    sprintf("%.2f", lower$gamma_PrEP)),
  macro("betaPrEPBaseline",  sprintf("%.2f", baseline$beta_PrEP)),
  macro("gammaPrEPUpper",    sprintf("%.1f", upper$gamma_PrEP * 100)),
  macro("gammaPrEPUpperDecimal", sprintf("%.2f", upper$gamma_PrEP)),
  macro("PrEPCovPerFundingBaseline", sprintf("%.1f\\%%", baseline$prep_cov_per_funding * 100)),
  macro("PrEPCovPerFundingUpper",    sprintf("%.1f\\%%", upper$prep_cov_per_funding * 100)),
  "% '40% cut for both ART and PrEP funding' sentence, at the new baseline r:",
  macro("bothFortyIncidenceBaseline",   sprintf("%.2f", both_inc$mean)),
  macro("bothFortyIncidenceBaselineLo", sprintf("%.2f", both_inc$ci_lower)),
  macro("bothFortyIncidenceBaselineHi", sprintf("%.2f", both_inc$ci_upper)),
  macro("bothFortyAbsBaseline",   sprintf("%.2f", both_abs$mean)),
  macro("bothFortyAbsBaselineLo", sprintf("%.2f", both_abs$ci_lower)),
  macro("bothFortyAbsBaselineHi", sprintf("%.2f", both_abs$ci_upper)),
  macro("bothFortyIRRBaseline",   sprintf("%.2f", both_irr$mean)),
  macro("bothFortyIRRBaselineLo", sprintf("%.2f", both_irr$ci_lower)),
  macro("bothFortyIRRBaselineHi", sprintf("%.2f", both_irr$ci_upper)),
  "% 'impacts of funding reductions accumulate over time' trajectory sentence.",
  "% NOTE: LaTeX command names cannot contain digits -- spell years/percents out.",
  macro("trajIRRYearOneAtForty",   sprintf("%.2f", traj_irr["year1", "40"])),
  macro("trajIRRYearTwoAtForty",   sprintf("%.2f", traj_irr["year2", "40"])),
  macro("trajIRRYearTenAtForty",   sprintf("%.2f", traj_irr["year10", "40"])),
  macro("trajIRRYearTenAtTen",     sprintf("%.2f", traj_irr["year10", "10"])),
  macro("trajIRRYearTenAtTwenty",  sprintf("%.2f", traj_irr["year10", "20"]))
)
writeLines(macros, file.path(out_dir, "paper_macros.tex"))

cat("\nWrote:\n  output/paper_sensitivity_table.csv\n  output/paper_sensitivity_table.tex\n")
cat(sprintf("  %s/funding_scenarios_table.tex  (real Table 1, at new baseline r=%.1f%%)\n", table_dir, R_BASELINE*100))
cat(sprintf("  %s/funding_scenarios_table.tex  (cross-check, r=78%% -- should match existing paper.tex numbers)\n", table_dir_upper))
cat("  output/paper_macros.tex\n")
