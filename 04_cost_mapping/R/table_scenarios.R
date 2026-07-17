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
# The two percentage columns report the resulting ART/PrEP intervention-USE
# reduction. In block A that equals the labelled reduction; in block B it is the
# smaller, cost-mapped value (e.g. a 10% ART funding cut -> ~7.5% ART use
# reduction), because a funding cut of fraction d lowers coverage by only
# d * (P_baseline - gamma) / P_baseline. The funding level itself lives in the
# row label and the block header. Stacking the blocks makes the dampening
# explicit and matches the two-panel Figure 1.
#
# Every quantity is computed from the stage-03 surrogate with common random
# numbers (compute_scenario_draws_crn), so the mean incidence, absolute increase
# (cases per 100 p.y.), and incidence risk ratio all reference the SAME baseline
# predictive draws and are mutually consistent.
# =============================================================================

library(kableExtra)
library(dplyr)

# Align a numeric vector on its decimal point for LaTeX, by prefixing phantom
# boxes: one \phantom{0} per missing integer digit (relative to the widest value
# in the vector) and \phantom{-} for non-negative values when any value is
# negative. Body-font digits and the minus sign are fixed-width, so the phantoms
# reserve exactly the space of the characters they stand in for and every cell
# ends up the same width -- decimals line up. Returns LaTeX strings; used only
# for the rendered table, never for the tracked CSV.
.align_decimal <- function(x, digits = 2) {
  neg      <- any(x < 0, na.rm = TRUE)
  body     <- formatC(abs(x), format = "f", digits = digits)  # unsigned, e.g. "10.12"
  int_part <- sub("\\..*$", "", body)
  max_int  <- max(nchar(int_part))
  vapply(seq_along(x), function(i) {
    sign_box <- if (isTRUE(x[i] < 0)) "-" else if (neg) "\\phantom{-}" else ""
    n_pad    <- max_int - nchar(int_part[i])
    pad_box  <- if (n_pad > 0) strrep("\\phantom{0}", n_pad) else ""
    paste0(sign_box, pad_box, body[i])
  }, character(1))
}

# "mean [lo, hi]" composite with each of the three sub-columns aligned across
# the whole column. `align` toggles the LaTeX phantom padding (TRUE for the
# rendered table, FALSE for the plain CSV strings).
.format_ci_column <- function(mean_v, lo_v, hi_v, digits = 2, align = FALSE) {
  fmt <- if (align) function(v) .align_decimal(v, digits)
         else       function(v) formatC(v, format = "f", digits = digits)
  paste0(fmt(mean_v), " [", fmt(lo_v), ", ", fmt(hi_v), "]")
}

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

  # Coverage-reduction points fed to the surrogate. Block A: identity mapping
  # (labelled reduction IS the use reduction). Block B: funding reduction scaled
  # to a use reduction by the cost-model factors.
  art_use_A  <- block$art;  prep_use_A <- block$prep
  art_use_B  <- block$art  * art_cov_per_funding
  prep_use_B <- block$prep * prep_cov_per_funding

  cov_use  <- cbind(-art_use_A, -prep_use_A)
  cov_fund <- cbind(-art_use_B, -prep_use_B)
  newX_native <- rbind(cov_use, cov_fund)

  # ---- Surrogate draws with common random numbers ---------------------------
  set.seed(0)  # match the surrogate/figure RNG convention
  draws <- compute_scenario_draws_crn(
    newX_native, gp_incidence_fit,
    n_samples_per_checkpoint = n_samples_per_checkpoint,
    common_random_numbers = common_random_numbers, tick = horizon_tick)

  q_lo <- ci_probs[1]; q_hi <- ci_probs[2]

  # Per-draw absolute increase (cases per 100 p.y.) vs the shared baseline draws.
  abs_increase <- sweep(draws$incidence, 2, draws$baseline, "-")

  # ---- Full per-row draw matrices (baseline row first, then A, then B) ------
  # Baseline is exact by construction: incidence = its own draws, absolute
  # increase = 0, IRR = 1.
  n_draws  <- length(draws$baseline)
  inc_full <- rbind(draws$baseline, draws$incidence)
  abs_full <- rbind(rep(0, n_draws), abs_increase)
  irr_full <- rbind(rep(1, n_draws), draws$irr)

  row_summary <- function(M) list(
    mean = rowMeans(M),
    lo   = apply(M, 1, quantile, probs = q_lo),
    hi   = apply(M, 1, quantile, probs = q_hi))
  s_inc <- row_summary(inc_full)
  s_abs <- row_summary(abs_full)
  s_irr <- row_summary(irr_full)

  # ---- Assemble the display columns -----------------------------------------
  scenario <- c("Baseline (no cuts)", block$label, block$label)
  art_use  <- c(0, art_use_A,  art_use_B)  * 100
  prep_use <- c(0, prep_use_A, prep_use_B) * 100

  build_table <- function(align_ci) {
    data.frame(
      `Scenario` = scenario,
      `ART Use Reduction (\\%)`  = sprintf("%.1f", art_use),
      `PrEP Use Reduction (\\%)` = sprintf("%.1f", prep_use),
      `Mean HIV Incidence per 100 p.y. [95\\% PPI]` =
        .format_ci_column(s_inc$mean, s_inc$lo, s_inc$hi, align = align_ci),
      `Absolute Increase per 100 p.y. [95\\% PPI]` =
        .format_ci_column(s_abs$mean, s_abs$lo, s_abs$hi, align = align_ci),
      `Incidence Risk Ratio [95\\% PPI]` =
        .format_ci_column(s_irr$mean, s_irr$lo, s_irr$hi, align = align_ci),
      check.names = FALSE, stringsAsFactors = FALSE)
  }

  # Plain, un-padded strings for the returned value and the tracked CSV.
  output_table <- build_table(align_ci = FALSE)

  # ---- LaTeX table with two grouped blocks, landscape -----------------------
  # Phantom-padded CI columns so decimals line up; escape the literal "%" in
  # scenario labels for LaTeX (raw % starts a comment). The returned/CSV table
  # keeps the plain, un-padded strings.
  latex_df <- build_table(align_ci = TRUE)
  latex_df$Scenario <- gsub("%", "\\\\%", latex_df$Scenario)
  latex_table <- kable(latex_df,
                       format = "latex",
                       booktabs = TRUE,
                       escape = FALSE,
                       align = c("l", "r", "r", "r", "r", "r"),
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
