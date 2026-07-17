library(ggplot2)
library(dplyr)
library(tidyr)

#' Create bar plot comparing funding reduction scenarios
#'
#' @param key_scenarios Data frame with funding scenario results
#' @param title Plot title
#' @param y_metric Which metric to display: "all", "incidence", "relative_increase", or "risk_ratio"
#' @param art_funding_reduction_pct Numeric vector of ART reduction percentages to filter (NULL = all)
#' @param prep_funding_reduction_pct Numeric vector of PrEP reduction percentages to filter (NULL = all)
#' @param incidence_ci_results Data frame with columns: grid_scenario, ci_lower, ci_upper (optional)
#' @param relative_change_ci_results Data frame with columns: grid_scenario, ci_lower, ci_upper (optional)
#' @param irr_ci_results Data frame with columns: grid_scenario, ci_lower, ci_upper (optional)
#' @param save_path Path to save plot (without extension). When y_metric="all",
#'   saves three files: "<save_path>_incidence.png", "<save_path>_relative_increase.png", and "<save_path>_risk_ratio.png"
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @param dpi Resolution for saved plot
#'
#' @return When y_metric="all": list with three ggplot objects (incidence, relative_increase, and risk_ratio).
#'   Otherwise: single ggplot object
#'
#' @details
#' Filtering logic:
#' - If both art_funding_reduction_pct and prep_funding_reduction_pct have same length:
#'   Parallel matching (pair-wise combinations)
#' - If only one is specified: show all combinations with those values
#' - If both are different lengths: Cartesian product (all combinations)
#' - Baseline (0%, 0%) is automatically included unless showing relative_increase only
#'
plot_funding_bar <- function(key_scenarios,
                              title = "Govt. Funding Reduction Comparisons",
                              y_metric = c("all", "incidence", "relative_increase", "risk_ratio"),
                              art_funding_reduction_pct = NULL,
                              prep_funding_reduction_pct = NULL,
                              incidence_ci_results = NULL,
                              relative_change_ci_results = NULL,
                              irr_ci_results = NULL,
                              save_path = NULL,
                              text_size = 3.5,
                              width = 10,
                              height = 8,
                              dpi = 300) {

  y_metric <- match.arg(y_metric)

  # Ensure required columns exist
  required_cols <- c("pct_delta_ART_fund", "pct_delta_PrEP_fund", "year10_mean_incidence")
  missing_cols <- setdiff(required_cols, names(key_scenarios))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }

  # Calculate relative change if not present
  if (!"relative_increase_pct" %in% names(key_scenarios)) {
    baseline_idx <- which(
      abs(key_scenarios$pct_delta_ART_fund) < 0.01 &
      abs(key_scenarios$pct_delta_PrEP_fund) < 0.01
    )

    if (length(baseline_idx) == 0) {
      stop("Cannot find baseline scenario (0%, 0%) to calculate relative change")
    }

    baseline_incidence <- key_scenarios$year10_mean_incidence[baseline_idx[1]]
    key_scenarios$relative_increase_pct <-
      ((key_scenarios$year10_mean_incidence - baseline_incidence) /
         baseline_incidence) * 100
  }

  # Calculate incidence risk ratio if not present
  if (!"incidence_risk_ratio" %in% names(key_scenarios)) {
    baseline_idx <- which(
      abs(key_scenarios$pct_delta_ART_fund) < 0.01 &
      abs(key_scenarios$pct_delta_PrEP_fund) < 0.01
    )

    if (length(baseline_idx) == 0) {
      stop("Cannot find baseline scenario (0%, 0%) to calculate risk ratio")
    }

    baseline_incidence <- key_scenarios$year10_mean_incidence[baseline_idx[1]]
    key_scenarios$incidence_risk_ratio <-
      key_scenarios$year10_mean_incidence / baseline_incidence
  }

  # Filter scenarios based on funding reduction percentages
  filtered_data <- key_scenarios
  baseline_data <- filtered_data[
    abs(filtered_data$pct_delta_ART_fund) < 0.01 &
    abs(filtered_data$pct_delta_PrEP_fund) < 0.01,
  ]

  if (!is.null(art_funding_reduction_pct) || !is.null(prep_funding_reduction_pct)) {
    # Determine matching strategy
    parallel_match <- FALSE
    if (!is.null(art_funding_reduction_pct) &&
        !is.null(prep_funding_reduction_pct) &&
        length(art_funding_reduction_pct) == length(prep_funding_reduction_pct)) {
      parallel_match <- TRUE
    }

    if (parallel_match) {
      # Parallel matching: pair-wise combinations
      selected_rows <- vector("list", length(art_funding_reduction_pct))
      for (i in seq_along(art_funding_reduction_pct)) {
        matching_idx <- which(
          abs(filtered_data$pct_delta_ART_fund - art_funding_reduction_pct[i]) < 0.5 &
          abs(filtered_data$pct_delta_PrEP_fund - prep_funding_reduction_pct[i]) < 0.5
        )
        if (length(matching_idx) > 0) {
          selected_rows[[i]] <- filtered_data[matching_idx[1], ]
        }
      }
      filtered_data <- do.call(rbind, Filter(Negate(is.null), selected_rows))

    } else {
      # Cartesian product or single-dimension filtering
      if (!is.null(art_funding_reduction_pct)) {
        art_matches <- sapply(filtered_data$pct_delta_ART_fund, function(x) {
          any(abs(x - art_funding_reduction_pct) < 0.5)
        })
        filtered_data <- filtered_data[art_matches, ]
      }

      if (!is.null(prep_funding_reduction_pct)) {
        prep_matches <- sapply(filtered_data$pct_delta_PrEP_fund, function(x) {
          any(abs(x - prep_funding_reduction_pct) < 0.5)
        })
        filtered_data <- filtered_data[prep_matches, ]
      }
    }

    # Always include baseline for incidence and risk_ratio comparisons
    if (y_metric != "relative_increase" && y_metric != "risk_ratio" && nrow(baseline_data) > 0) {
      # Check if baseline is already in filtered data
      baseline_in_filtered <- any(
        abs(filtered_data$pct_delta_ART_fund) < 0.01 &
        abs(filtered_data$pct_delta_PrEP_fund) < 0.01
      )

      if (!baseline_in_filtered) {
        filtered_data <- rbind(baseline_data, filtered_data)
      }
    }
  }

  # Remove baseline if showing only relative change
  if (y_metric == "relative_increase") {
    filtered_data <- filtered_data[
      !(abs(filtered_data$pct_delta_ART_fund) < 0.01 &
        abs(filtered_data$pct_delta_PrEP_fund) < 0.01),
    ]
  }

  if (nrow(filtered_data) == 0) {
    stop("No scenarios match the specified filters")
  }

  # Calculate total funding reduction for color mapping
  filtered_data$total_reduction_pct <- abs(filtered_data$pct_delta_ART_fund) +
                                       abs(filtered_data$pct_delta_PrEP_fund)

  # Convert to proportion (0 to 1) for gradient scale
  filtered_data$total_reduction_prop <- filtered_data$total_reduction_pct / 100

  # Create scenario labels
  filtered_data$scenario_label <- sprintf(
    "ART: %.0f%%, PrEP: %.0f%%",
    abs(filtered_data$pct_delta_ART_fund),
    abs(filtered_data$pct_delta_PrEP_fund)
  )

  # Sort for consistent ordering (by total reduction, then ART, then PrEP)
  filtered_data <- filtered_data %>%
    arrange(total_reduction_pct, pct_delta_ART_fund, pct_delta_PrEP_fund)

  filtered_data$scenario_label <- factor(
    filtered_data$scenario_label,
    levels = unique(filtered_data$scenario_label)
  )

  # Match CI data if provided
  if (!is.null(incidence_ci_results) || !is.null(relative_change_ci_results) || !is.null(irr_ci_results)) {

    # Initialize CI columns
    filtered_data$incidence_ci_lower <- NA
    filtered_data$incidence_ci_upper <- NA
    filtered_data$relative_change_ci_lower <- NA
    filtered_data$relative_change_ci_upper <- NA
    filtered_data$irr_ci_lower <- NA
    filtered_data$irr_ci_upper <- NA

    # Grid parameters (must match cost_mapping.R)
    ngrid <- 101
    newX <- seq(-0.75, 0, length.out = ngrid)
    newX_grid <- as.matrix(expand.grid(newX, newX))
    tolerance <- 0.001

    # Match each filtered scenario to grid
    for (i in 1:nrow(filtered_data)) {
      # Find grid scenario by ART/PrEP reduction proportions
      # Use existing ART_pct_reduction and PrEP_pct_reduction columns (proportions, not percentages)
      # These should exist in key_scenarios if it's model_proportions_mean_incidence

      if ("ART_pct_reduction" %in% names(filtered_data) && "PrEP_pct_reduction" %in% names(filtered_data)) {
        art_reduction_prop <- filtered_data$ART_pct_reduction[i]
        prep_reduction_prop <- filtered_data$PrEP_pct_reduction[i]

        # Find matching grid scenario
        grid_idx <- which(
          abs(newX_grid[, 1] - art_reduction_prop) < tolerance &
          abs(newX_grid[, 2] - prep_reduction_prop) < tolerance
        )

        if (length(grid_idx) > 0) {
          grid_scenario_num <- grid_idx[1]

          # Match CIs from provided data frames
          if (!is.null(incidence_ci_results)) {
            filtered_data$incidence_ci_lower[i] <- incidence_ci_results$ci_lower[grid_scenario_num]
            filtered_data$incidence_ci_upper[i] <- incidence_ci_results$ci_upper[grid_scenario_num]
          }

          if (!is.null(relative_change_ci_results)) {
            filtered_data$relative_change_ci_lower[i] <- relative_change_ci_results$ci_lower[grid_scenario_num]
            filtered_data$relative_change_ci_upper[i] <- relative_change_ci_results$ci_upper[grid_scenario_num]
          }

          if (!is.null(irr_ci_results)) {
            filtered_data$irr_ci_lower[i] <- irr_ci_results$ci_lower[grid_scenario_num]
            filtered_data$irr_ci_upper[i] <- irr_ci_results$ci_upper[grid_scenario_num]
          }
        }
      }
    }
  }

  # Prepare data for plotting
  if (y_metric == "all") {
    # Create incidence plot
    p_incidence <- ggplot(filtered_data, aes(x = scenario_label, y = year10_mean_incidence, fill = year10_mean_incidence)) +
      geom_bar(stat = "identity", width = 0.7) +
      geom_errorbar(
        aes(ymin = incidence_ci_lower, ymax = incidence_ci_upper),
        width = 0.3,
        linewidth = 0.5,
        color = "black",
        alpha = 0.4
      ) +
      geom_text(
        aes(label = sprintf("%.2f", year10_mean_incidence)),
        vjust = -0.5,
        size = 3.5
      ) +
      scale_fill_gradient(
        low = "#ffc6d9",
        high = "#7b2d7d",
        name = "Mean Incidence (per 100 p.y.)",
        guide = guide_colorbar(
          barwidth = 1.5,
          barheight = 15,
          ticks.colour = "black",
          ticks.linewidth = 0.5,
          frame.colour = "black",
          frame.linewidth = 0.5
        )
      ) +
      labs(
        title = paste0(title, " - Mean Incidence (per 100 p.y.)"),
        x = "ART x PrEP Govt. Funding Reduction (%)",
        y = "HIV Mean Incidence (per 100 p.y.)"
      ) +
      theme_minimal(base_size = 13) +
      theme(
        axis.title = element_text(face = "bold", size=text_size*4),
        axis.text = element_text(size=text_size*3),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 15),
        legend.position = "right"
      )

    # Create relative change plot (excluding baseline)
    rel_data <- filtered_data[
      !(abs(filtered_data$pct_delta_ART_fund) < 0.01 &
        abs(filtered_data$pct_delta_PrEP_fund) < 0.01),
    ]

    p_relative <- ggplot(rel_data, aes(x = scenario_label, y = relative_increase_pct, fill = relative_increase_pct)) +
      geom_bar(stat = "identity", width = 0.7) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
      geom_errorbar(
        aes(ymin = relative_change_ci_lower, ymax = relative_change_ci_upper),
        width = 0.3,
        linewidth = 0.5,
        color = "black",
        alpha = 0.4
      ) +
      geom_text(
        aes(label = sprintf("%+.1f%%", relative_increase_pct)),
        vjust = ifelse(rel_data$relative_increase_pct >= 0, -0.5, 1.5),
        size = 3.5
      ) +
      scale_fill_gradient(
        low = "#ffc6d9",
        high = "#7b2d7d",
        name = "Relative Change in Mean Incidence (%)",
        guide = guide_colorbar(
          barwidth = 1.5,
          barheight = 15,
          ticks.colour = "black",
          ticks.linewidth = 0.5,
          frame.colour = "black",
          frame.linewidth = 0.5
        )
      ) +
      labs(
        title = paste0(title, " - Relative Change"),
        x = "ART x PrEP Govt. Funding Reduction (%)",
        y = "Relative Change in Mean Incidence (%)"
      ) +
      theme_minimal(base_size = 13) +
      theme(
        axis.title = element_text(face = "bold", size=text_size*4),
        axis.text = element_text(size=text_size*3),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 15),
        legend.position = "right"
      )

    # Create risk ratio plot (excluding baseline)
    p_risk_ratio <- ggplot(rel_data, aes(x = scenario_label, y = incidence_risk_ratio, fill = incidence_risk_ratio)) +
      geom_bar(stat = "identity", width = 0.7) +
      geom_hline(yintercept = 1.0, linetype = "dashed", color = "grey40", linewidth = 0.8) +
      geom_errorbar(
        aes(ymin = irr_ci_lower, ymax = irr_ci_upper),
        width = 0.3,
        linewidth = 0.5,
        color = "black",
        alpha = 0.4
      ) +
      geom_text(
        aes(label = sprintf("%.2f", incidence_risk_ratio)),
        vjust = -0.5,
        size = 3.5
      ) +
      scale_fill_gradient(
        low = "#ffc6d9",
        high = "#7b2d7d",
        name = "Mean Incidence Risk Ratio",
        guide = guide_colorbar(
          barwidth = 1.5,
          barheight = 15,
          ticks.colour = "black",
          ticks.linewidth = 0.5,
          frame.colour = "black",
          frame.linewidth = 0.5
        )
      ) +
      labs(
        title = paste0(title, " - Mean Incidence Risk Ratio"),
        x = "ART x PrEP Govt. Funding Reduction (%)",
        y = "Mean Incidence Risk Ratio"
      ) +
      theme_minimal(base_size = 13) +
      theme(
        axis.title = element_text(face = "bold", size=text_size*4),
        axis.text = element_text(size=text_size*3),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 15),
        legend.position = "right"
      )

    # Save all three plots if path is provided
    if (!is.null(save_path)) {
      ggsave(
        filename = paste0(save_path, "_incidence.png"),
        plot = p_incidence,
        width = width,
        height = height,
        dpi = dpi
      )

      ggsave(
        filename = paste0(save_path, "_relative_increase.png"),
        plot = p_relative,
        width = width,
        height = height,
        dpi = dpi
      )

      ggsave(
        filename = paste0(save_path, "_risk_ratio.png"),
        plot = p_risk_ratio,
        width = width,
        height = height,
        dpi = dpi
      )
    }

    return(list(
      p_incidence = p_incidence,
      p_relative = p_relative,
      p_risk_ratio = p_risk_ratio
    ))

  } else if (y_metric == "incidence") {
    p <- ggplot(filtered_data, aes(x = scenario_label, y = year10_mean_incidence, fill = year10_mean_incidence)) +
      geom_bar(stat = "identity", width = 0.7) +
      geom_errorbar(
        aes(ymin = incidence_ci_lower, ymax = incidence_ci_upper),
        width = 0.3,
        linewidth = 0.5,
        color = "black",
        alpha = 0.4
      ) +
      geom_text(
        aes(label = sprintf("%.2f", year10_mean_incidence)),
        vjust = -0.5,
        size = 3.5
      ) +
      scale_fill_gradient(
        low = "#ffc6d9",
        high = "#7b2d7d",
        name = "Mean Incidence (per 100 p.y.)",
        guide = guide_colorbar(
          barwidth = 1.5,
          barheight = 15,
          ticks.colour = "black",
          ticks.linewidth = 0.5,
          frame.colour = "black",
          frame.linewidth = 0.5
        )
      ) +
      labs(
        title = title,
        x = "ART x PrEP Govt. Funding Reduction (%)",
        y = "HIV Mean Incidence (per 100 p.y.)"
      ) +
      theme_minimal(base_size = 13) +
      theme(
        axis.title = element_text(face = "bold", size=text_size*4),
        axis.text = element_text(size=text_size*3),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 15),
        legend.position = "right"
      )

  } else if (y_metric == "relative_increase") {
    p <- ggplot(filtered_data, aes(x = scenario_label, y = relative_increase_pct, fill = relative_increase_pct)) +
      geom_bar(stat = "identity", width = 0.7) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
      geom_errorbar(
        aes(ymin = relative_change_ci_lower, ymax = relative_change_ci_upper),
        width = 0.3,
        linewidth = 0.5,
        color = "black",
        alpha = 0.4
      ) +
      geom_text(
        aes(label = sprintf("%+.1f%%", relative_increase_pct)),
        vjust = ifelse(filtered_data$relative_increase_pct >= 0, -0.5, 1.5),
        size = 3.5
      ) +
      scale_fill_gradient(
        low = "#ffc6d9",
        high = "#7b2d7d",
        name = "Relative Change in Mean Incidence (%)",
        guide = guide_colorbar(
          barwidth = 1.5,
          barheight = 15,
          ticks.colour = "black",
          ticks.linewidth = 0.5,
          frame.colour = "black",
          frame.linewidth = 0.5
        )
      ) +
      labs(
        title = title,
        x = "ART x PrEP Govt. Funding Reduction (%)",
        y = "Relative Change in Mean Incidence (%)"
      ) +
      theme_minimal(base_size = 13) +
      theme(
        axis.title = element_text(face = "bold", size=text_size*4),
        axis.text = element_text(size=text_size*3),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 15),
        legend.position = "right"
      )
  } else {  # risk_ratio
    p <- ggplot(filtered_data, aes(x = scenario_label, y = incidence_risk_ratio, fill = incidence_risk_ratio)) +
      geom_bar(stat = "identity", width = 0.7) +
      geom_hline(yintercept = 1.0, linetype = "dashed", color = "grey40", linewidth = 0.8) +
      geom_errorbar(
        aes(ymin = irr_ci_lower, ymax = irr_ci_upper),
        width = 0.3,
        linewidth = 0.5,
        color = "black",
        alpha = 0.4
      ) +
      geom_text(
        aes(label = sprintf("%.2f", incidence_risk_ratio)),
        vjust = -0.5,
        size = 3.5
      ) +
      scale_fill_gradient(
        low = "#ffc6d9",
        high = "#7b2d7d",
        name = "Mean Incidence Risk Ratio",
        guide = guide_colorbar(
          barwidth = 1.5,
          barheight = 15,
          ticks.colour = "black",
          ticks.linewidth = 0.5,
          frame.colour = "black",
          frame.linewidth = 0.5
        )
      ) +
      labs(
        title = title,
        x = "ART x PrEP Govt. Funding Reduction (%)",
        y = "Mean Incidence Risk Ratio"
      ) +
      theme_minimal(base_size = 13) +
      theme(
        axis.title = element_text(face = "bold", size=text_size*4),
        axis.text = element_text(size=text_size*3),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 15),
        legend.position = "right"
      )

    # Save plot if path is provided
    if (!is.null(save_path)) {
      ggsave(
        filename = paste0(save_path, ".png"),
        plot = p,
        width = width,
        height = height,
        dpi = dpi
      )
    }

    return(p)
  }
}
