library(kableExtra)
library(dplyr)

#' Create funding scenario table with incidence outcomes
#'
#' @param key_scenarios Data frame with scenario definitions
#' @param funding_data Full funding reduction scenario data
#' @param G_ART_baseline_val Baseline ART funding value
#' @param G_PrEP_baseline_val Baseline PrEP funding value
#' @param output_dir Directory to save output files
#' @param metric_type Which relative metric to display: "both", "relative_increase", or "risk_ratio"
#'
#' @return The formatted output table (data frame)
#'
#' @details
#' When metric_type = "both": includes both relative increase % and risk ratio columns
#' When metric_type = "relative_increase": only includes relative increase % column
#' When metric_type = "risk_ratio": only includes risk ratio column
#'
create_funding_scenario_table <- function(key_scenarios,
                                          funding_data,
                                          G_ART_baseline_val,
                                          G_PrEP_baseline_val,
                                          incidence_ci_results,
                                          relative_change_ci_results,
                                          irr_ci_results,
                                          output_dir = "results/orig_params/",
                                          metric_type = c("both", "relative_increase", "risk_ratio"),
                                          common_random_numbers = TRUE,
                                          gp_incidence_fit = NULL,
                                          n_samples_per_checkpoint = 50,
                                          horizon_tick = 10,
                                          ci_probs = c(0.05, 0.95)) {

  metric_type <- match.arg(metric_type)

  # When common_random_numbers = TRUE, the incidence risk ratio and relative
  # change are recomputed per scenario with common random numbers (baseline and
  # scenario surrogate draws share randomness within each checkpoint) via
  # compute_irr_draws_crn(), so they match the forest figure and no longer
  # inherit the independent-draw noise floor that pushes small-effect intervals
  # below the null. This requires the surrogate (gp_incidence_fit); if it is not
  # supplied, we fall back to the independent-draw grid CI results.
  use_crn <- isTRUE(common_random_numbers) && !is.null(gp_incidence_fit)

  # For each scenario, find the closest matching point in funding_data
  # and extract the mean incidence and sd
  tolerance <- 1  # 1% tolerance for matching

  key_scenarios$mean_incidence <- NA
  key_scenarios$sd_incidence <- NA

  for (i in 1:nrow(key_scenarios)) {
    # Find matching row in funding_data
    matching_idx <- which(
      abs(funding_data$pct_delta_ART_fund - key_scenarios$art_funding_reduction_pct[i]) < tolerance &
      abs(funding_data$pct_delta_PrEP_fund - key_scenarios$prep_funding_reduction_pct[i]) < tolerance
    )

    if (length(matching_idx) > 0) {
      # Take first match (or closest if multiple)
      if (length(matching_idx) > 1) {
        # Find closest match
        distances <- abs(funding_data$pct_delta_ART_fund[matching_idx] - key_scenarios$art_funding_reduction_pct[i]) +
                    abs(funding_data$pct_delta_PrEP_fund[matching_idx] - key_scenarios$prep_funding_reduction_pct[i])
        idx <- matching_idx[which.min(distances)]
      } else {
        idx <- matching_idx[1]
      }

      key_scenarios$mean_incidence[i] <- funding_data$year10_mean_incidence[idx]
      key_scenarios$sd_incidence[i] <- funding_data$year10_sd_incidence[idx]
    } else {
      warning(sprintf("No match found for scenario '%s' (ART: %.1f%%, PrEP: %.1f%%)",
                     key_scenarios$Scenario[i],
                     key_scenarios$art_funding_reduction_pct[i],
                     key_scenarios$prep_funding_reduction_pct[i]))
    }
  }

  # Match posterior predictive intervals from CI results
  # Note: CI data frames are indexed by original grid scenarios, not deduplicated funding_data rows
  key_scenarios$incidence_ci_lower <- NA
  key_scenarios$incidence_ci_upper <- NA
  key_scenarios$relative_change_ci_lower <- NA
  key_scenarios$relative_change_ci_upper <- NA
  key_scenarios$relative_change_mean <- NA
  key_scenarios$irr_ci_lower <- NA
  key_scenarios$irr_ci_upper <- NA
  key_scenarios$irr_mean <- NA

  # Native coverage-reduction coordinates per scenario (for the CRN recompute).
  key_scenarios$art_cov_reduction <- NA
  key_scenarios$prep_cov_reduction <- NA

  # Loop through scenarios to match with CI data
  for (i in 1:nrow(key_scenarios)) {
    # Find matching row in funding_data first
    matching_idx <- which(
      abs(funding_data$pct_delta_ART_fund - key_scenarios$art_funding_reduction_pct[i]) < tolerance &
      abs(funding_data$pct_delta_PrEP_fund - key_scenarios$prep_funding_reduction_pct[i]) < tolerance
    )

    if (length(matching_idx) > 0) {
      # Get the ART/PrEP reduction values from funding_data
      if (length(matching_idx) > 1) {
        distances <- abs(funding_data$pct_delta_ART_fund[matching_idx] - key_scenarios$art_funding_reduction_pct[i]) +
                    abs(funding_data$pct_delta_PrEP_fund[matching_idx] - key_scenarios$prep_funding_reduction_pct[i])
        idx <- matching_idx[which.min(distances)]
      } else {
        idx <- matching_idx[1]
      }

      art_reduction <- funding_data$ART_pct_reduction[idx]
      prep_reduction <- funding_data$PrEP_pct_reduction[idx]

      key_scenarios$art_cov_reduction[i] <- art_reduction
      key_scenarios$prep_cov_reduction[i] <- prep_reduction

      # Now find this in the CI results by matching ART/PrEP reduction values
      # The CI data frames' grid_scenario corresponds to rows in the original grid
      # We need to reconstruct which grid scenario this is

      # Create the grid to find the scenario index
      ngrid <- 101
      newX <- seq(-0.75, 0, length.out = ngrid)
      newX_grid <- as.matrix(expand.grid(newX, newX))

      # Find matching grid scenario
      grid_tolerance <- 0.001  # Small tolerance for floating point comparison
      grid_idx <- which(
        abs(newX_grid[, 1] - art_reduction) < grid_tolerance &
        abs(newX_grid[, 2] - prep_reduction) < grid_tolerance
      )

      if (length(grid_idx) > 0) {
        grid_scenario_num <- grid_idx[1]

        # Match incidence CI
        key_scenarios$incidence_ci_lower[i] <- incidence_ci_results$ci_lower[grid_scenario_num]
        key_scenarios$incidence_ci_upper[i] <- incidence_ci_results$ci_upper[grid_scenario_num]

        # Match relative change CI
        key_scenarios$relative_change_ci_lower[i] <- relative_change_ci_results$ci_lower[grid_scenario_num]
        key_scenarios$relative_change_ci_upper[i] <- relative_change_ci_results$ci_upper[grid_scenario_num]
        key_scenarios$relative_change_mean[i] <- relative_change_ci_results$mean_relative_change[grid_scenario_num]

        # Match IRR CI
        key_scenarios$irr_ci_lower[i] <- irr_ci_results$ci_lower[grid_scenario_num]
        key_scenarios$irr_ci_upper[i] <- irr_ci_results$ci_upper[grid_scenario_num]
        key_scenarios$irr_mean[i] <- irr_ci_results$mean_irr[grid_scenario_num]
      }
    }
  }

  # ---- Common random numbers override for IRR and relative change -----------
  # Recompute the risk ratio and relative-change columns per scenario with
  # common random numbers so they match the forest figure. Absolute incidence
  # and its PPI are left untouched (CRN affects only the ratio-to-baseline
  # quantities). Each scenario is queried at its own native coverage-reduction
  # point at the reporting horizon (single tick).
  if (use_crn) {
    valid <- which(!is.na(key_scenarios$art_cov_reduction))
    if (length(valid) > 0) {
      newX_native <- cbind(key_scenarios$art_cov_reduction[valid],
                           key_scenarios$prep_cov_reduction[valid])
      set.seed(0)  # match the surrogate/figure RNG convention
      irr_draws <- compute_irr_draws_crn(
        newX_native, gp_incidence_fit,
        n_samples_per_checkpoint = n_samples_per_checkpoint,
        common_random_numbers = TRUE, tick = horizon_tick)
      irr_sum <- summarise_irr_draws(irr_draws, probs = ci_probs)

      key_scenarios$irr_mean[valid]     <- irr_sum$mean_irr
      key_scenarios$irr_ci_lower[valid] <- irr_sum$ci_lower
      key_scenarios$irr_ci_upper[valid] <- irr_sum$ci_upper

      # Relative change (%) is the same CRN draws expressed as (IRR - 1) * 100.
      rel_draws <- (irr_draws - 1) * 100
      key_scenarios$relative_change_mean[valid]     <- apply(rel_draws, 1, mean)
      key_scenarios$relative_change_ci_lower[valid] <- apply(rel_draws, 1, quantile, probs = ci_probs[1])
      key_scenarios$relative_change_ci_upper[valid] <- apply(rel_draws, 1, quantile, probs = ci_probs[2])
    }
  }

  # Use proper posterior predictive intervals from CI results
  key_scenarios$lower_95 <- key_scenarios$incidence_ci_lower
  key_scenarios$upper_95 <- key_scenarios$incidence_ci_upper

  # Find baseline incidence (0% reductions) for relative comparison
  baseline_idx <- which(
    abs(key_scenarios$art_funding_reduction_pct) < 0.01 &
    abs(key_scenarios$prep_funding_reduction_pct) < 0.01
  )

  if (length(baseline_idx) > 0) {
    baseline_incidence <- key_scenarios$mean_incidence[baseline_idx[1]]
  } else {
    # If baseline not in scenarios, find it in the data
    baseline_idx_data <- which(
      abs(funding_data$pct_delta_ART_fund) < 0.01 &
      abs(funding_data$pct_delta_PrEP_fund) < 0.01
    )
    if (length(baseline_idx_data) > 0) {
      baseline_incidence <- funding_data$year10_mean_incidence[baseline_idx_data[1]]
      cat(sprintf("Baseline incidence (from data): %.4f per 100 PY\n", baseline_incidence))
    } else {
      warning("Could not find baseline incidence")
      baseline_incidence <- NA
    }
  }

  # Calculate metrics based on metric_type
  if (metric_type %in% c("both", "relative_increase")) {
    # Use CI mean from relative_change_ci_results
    key_scenarios$relative_increase_pct <- key_scenarios$relative_change_mean
  }

  if (metric_type %in% c("both", "risk_ratio")) {
    # Use CI mean from irr_ci_results
    key_scenarios$incidence_risk_ratio <- key_scenarios$irr_mean
  }

  # Format for table output - build columns based on metric_type
  output_table <- data.frame(
    `Govt. Funding Reduction Scenarios (%)` = key_scenarios$Scenario,
    `ART Govt. Reduction (\\%)` = sprintf("%.1f", key_scenarios$art_funding_reduction_pct),
    `PrEP Govt. Reduction (\\%)` = sprintf("%.1f", key_scenarios$prep_funding_reduction_pct),
    `Mean HIV Incidence per 100 p.y. [95\\% Posterior Predictive Interval]` = sprintf(
      "%.2f [%.2f, %.2f]",
      key_scenarios$mean_incidence,
      key_scenarios$lower_95,
      key_scenarios$upper_95
    ),
    check.names = FALSE
  )

  # Add metric columns based on metric_type with confidence intervals
  if (metric_type == "both") {
    output_table$`Incidence Risk Ratio [95\\% PPI]` <- sprintf(
      "%.2f [%.2f, %.2f]",
      key_scenarios$incidence_risk_ratio,
      key_scenarios$irr_ci_lower,
      key_scenarios$irr_ci_upper
    )

    output_table$`Relative Increase (\\%) [95\\% PPI]` <- sprintf(
      "%.1f [%.1f, %.1f]",
      key_scenarios$relative_increase_pct,
      key_scenarios$relative_change_ci_lower,
      key_scenarios$relative_change_ci_upper
    )
  } else if (metric_type == "risk_ratio") {
    output_table$`Incidence Risk Ratio [95\\% PPI]` <- sprintf(
      "%.2f [%.2f, %.2f]",
      key_scenarios$incidence_risk_ratio,
      key_scenarios$irr_ci_lower,
      key_scenarios$irr_ci_upper
    )
  } else if (metric_type == "relative_increase") {
    output_table$`Relative Increase (\\%) [95\\% PPI]` <- sprintf(
      "%.1f [%.1f, %.1f]",
      key_scenarios$relative_increase_pct,
      key_scenarios$relative_change_ci_lower,
      key_scenarios$relative_change_ci_upper
    )
  }

  # Create LaTeX table
  latex_table <- kable(output_table,
                       format = "latex",
                       booktabs = TRUE,
                       escape = FALSE,
                       caption = "Govt. Funding Reduction Scenarios and HIV Incidence Outcomes") %>%
    kable_styling(latex_options = c("scale_down", "hold_position"))

  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Write to files
  tex_file <- file.path(output_dir, "funding_scenarios_table.tex")

  writeLines(latex_table, tex_file)

  # Return the output table
  return(invisible(output_table))
}
