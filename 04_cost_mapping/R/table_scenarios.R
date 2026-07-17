# =============================================================================
# Funding / use scenario table (main-text Table 1)
# -----------------------------------------------------------------------------
# Two stacked blocks sharing one baseline row:
#   A. Intervention use reduction  -- the labelled percentage IS the coverage
#      reduction fed to the surrogate (no cost mapping).
#   B. Government funding reduction -- the labelled percentage is a funding
#      reduction, mapped to a (smaller) coverage reduction via the cost model
#      (paper Eqs. 1-2) before the surrogate query.
#
# The same percentage produces a larger incidence increase under direct use
# reduction than under funding reduction, because a funding cut of fraction d
# lowers coverage by only d * (P_baseline - gamma) / P_baseline. Stacking the
# blocks makes that dampening explicit and matches the two-panel Figure 1.
#
# Every quantity is computed from the stage-03 surrogate with common random
# numbers (compute_scenario_draws_crn), so the mean incidence, absolute increase
# (cases per 100 p.y.), and incidence risk ratio all reference the SAME baseline
# predictive draws and are mutually consistent.
# =============================================================================

library(kableExtra)
library(dplyr)

#' Create the two-block scenario table (Table 1).
#'
#' @param gp_incidence_fit Composite incidence surrogate (from surrogate.Rdata).
#' @param art_cov_per_funding Coverage-reduction fraction per unit ART funding
#'   reduction, (P_ART_baseline - gamma_ART) / P_ART_baseline. Used in block B.
#' @param prep_cov_per_funding Same, for PrEP.
#' @param reduction_levels ART-only / PrEP-only reduction levels per block
#'   (defaults to 10/25/40%). A "both at max level" row is added per block.
#' @param n_samples_per_checkpoint Surrogate draws per checkpoint.
#' @param horizon_tick Reporting horizon tick.
#' @param ci_probs Length-2 lower/upper quantile probabilities for the PPIs.
#' @param common_random_numbers If TRUE (default), baseline and scenario draws
#'   share randomness within each checkpoint (see irr_common_random_numbers.R).
#' @param output_dir Directory for funding_scenarios_table.tex.
#' @return The formatted output table (data.frame), invisibly.
create_funding_scenario_table <- function(gp_incidence_fit,
                                           art_cov_per_funding,
                                           prep_cov_per_funding,
                                           reduction_levels = c(0.10, 0.25, 0.40),
                                           n_samples_per_checkpoint = 50,
                                           horizon_tick = 10,
                                           ci_probs = c(0.05, 0.95),
                                           common_random_numbers = TRUE,
                                           output_dir = "output/") {

  max_level <- max(reduction_levels)

  # ---- Scenario definitions (per block, in labelled-reduction terms) --------
  # For each block: ART-only at each level, PrEP-only at each level, and a
  # "both" row at the maximum level. The labelled reduction is the coverage
  # reduction in block A and the funding reduction in block B.
  pct <- function(x) sprintf("%g%%", x * 100)
  build_block <- function() {
    rbind(
      data.frame(label = paste(pct(reduction_levels), "ART cut only"),
                 art = reduction_levels, prep = 0),
      data.frame(label = paste(pct(reduction_levels), "PrEP cut only"),
                 art = 0, prep = reduction_levels),
      data.frame(label = paste(pct(max_level), "cut for both"),
                 art = max_level, prep = max_level)
    )
  }
  block <- build_block()
  n_block <- nrow(block)

  # Coverage-reduction points fed to the surrogate. Block A: identity mapping.
  # Block B: funding reduction scaled to coverage by the cost-model factors.
  cov_use  <- cbind(-block$art,                       -block$prep)
  cov_fund <- cbind(-block$art * art_cov_per_funding, -block$prep * prep_cov_per_funding)
  newX_native <- rbind(cov_use, cov_fund)

  # ---- Surrogate draws with common random numbers ---------------------------
  set.seed(0)  # match the surrogate/figure RNG convention
  draws <- compute_scenario_draws_crn(
    newX_native, gp_incidence_fit,
    n_samples_per_checkpoint = n_samples_per_checkpoint,
    common_random_numbers = common_random_numbers, tick = horizon_tick)

  q_lo <- ci_probs[1]; q_hi <- ci_probs[2]
  fmt_inc <- function(v) sprintf("%.2f [%.2f, %.2f]",
                                 mean(v), quantile(v, q_lo), quantile(v, q_hi))
  fmt_irr <- function(v) sprintf("%.2f [%.2f, %.2f]",
                                 mean(v), quantile(v, q_lo), quantile(v, q_hi))

  # Per-draw absolute increase (cases per 100 p.y.) vs the shared baseline draws.
  abs_increase <- sweep(draws$incidence, 2, draws$baseline, "-")

  # ---- Assemble rows --------------------------------------------------------
  # Shared baseline row first, then block A rows, then block B rows.
  make_rows <- function(offset) {
    data.frame(
      Scenario = block$label,
      art_pct  = sprintf("%.1f", block$art * 100),
      prep_pct = sprintf("%.1f", block$prep * 100),
      incidence = vapply(seq_len(n_block),
                         function(k) fmt_inc(draws$incidence[offset + k, ]), character(1)),
      abs_inc  = vapply(seq_len(n_block),
                        function(k) fmt_inc(abs_increase[offset + k, ]), character(1)),
      irr      = vapply(seq_len(n_block),
                        function(k) fmt_irr(draws$irr[offset + k, ]), character(1)),
      check.names = FALSE, stringsAsFactors = FALSE)
  }

  baseline_row <- data.frame(
    Scenario = "Baseline (no cuts)", art_pct = "0.0", prep_pct = "0.0",
    incidence = fmt_inc(draws$baseline),
    abs_inc = sprintf("%.2f [%.2f, %.2f]", 0, 0, 0),
    irr = "1.00 [1.00, 1.00]",
    check.names = FALSE, stringsAsFactors = FALSE)

  rows_use  <- make_rows(0)
  rows_fund <- make_rows(n_block)
  all_rows  <- rbind(baseline_row, rows_use, rows_fund)

  output_table <- data.frame(
    `Scenario` = all_rows$Scenario,
    `ART Reduction (\\%)`  = all_rows$art_pct,
    `PrEP Reduction (\\%)` = all_rows$prep_pct,
    `Mean HIV Incidence per 100 p.y. [95\\% PPI]` = all_rows$incidence,
    `Absolute Increase per 100 p.y. [95\\% PPI]`  = all_rows$abs_inc,
    `Incidence Risk Ratio [95\\% PPI]` = all_rows$irr,
    check.names = FALSE, stringsAsFactors = FALSE)

  # ---- LaTeX table with two grouped blocks, landscape -----------------------
  # Escape the literal "%" in scenario labels for LaTeX (raw % starts a comment);
  # the returned/CSV table keeps the plain "%".
  latex_df <- output_table
  latex_df$Scenario <- gsub("%", "\\\\%", latex_df$Scenario)
  latex_table <- kable(latex_df,
                       format = "latex",
                       booktabs = TRUE,
                       escape = FALSE,
                       align = c("l", "r", "r", "l", "l", "l"),
                       caption = "HIV incidence outcomes under intervention-use and government-funding reductions") %>%
    kable_styling(latex_options = c("scale_down", "hold_position"),
                  full_width = FALSE) %>%
    pack_rows("A. Intervention use reduction", 2, 1 + n_block,
              bold = TRUE, latex_gap_space = "0.6em") %>%
    pack_rows("B. Government funding reduction", 2 + n_block, 1 + 2 * n_block,
              bold = TRUE, latex_gap_space = "0.6em") %>%
    landscape()

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  tex_file <- file.path(output_dir, "funding_scenarios_table.tex")
  writeLines(latex_table, tex_file)

  invisible(output_table)
}
