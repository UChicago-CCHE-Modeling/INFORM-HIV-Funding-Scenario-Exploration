library(ggplot2)
library(viridis)
library(dplyr)
library(shadowtext)


#------------------------------------------------
# Contour Heatmap
#------------------------------------------------

plot_contour_heatmap <- function(funding_reduction_scenarios, mean_var, pct_change_var, var_name, n_bins=9, baseline_value=NULL, text_size=3.5) {
  # Determine if we're plotting a risk ratio (which doesn't need "per 100 PY" units)
  is_risk_ratio <- grepl("ratio|rr", mean_var, ignore.case = TRUE) || grepl("ratio|rr", var_name, ignore.case = TRUE)

  # Calculate breaks for contour intervals
  # If baseline_value provided, use it as the starting point; otherwise use data min
  data_range <- range(funding_reduction_scenarios[[mean_var]], na.rm = TRUE)
  breaks_start <- if (!is.null(baseline_value)) baseline_value else data_range[1]
  breaks <- seq(breaks_start, data_range[2], length.out = n_bins + 1)

  # Heatmap visualization
  p_sim_heatmap <- ggplot(
    funding_reduction_scenarios,
    aes(
      x = pct_delta_ART_fund, y = pct_delta_PrEP_fund,
      z = !!sym(mean_var)
    )
  ) +
  geom_contour_filled(breaks = breaks) +
  geom_contour(color = "white", alpha = 0.3, breaks = breaks) +
    scale_fill_viridis_d(
      option = "plasma",
      name = if (is_risk_ratio) var_name else paste0(var_name, "\n(per 100 PY)")
    ) +
    scale_x_continuous(
      labels = function(x) paste0(x, "%"),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      labels = function(x) paste0(x, "%"),
      expand = c(0, 0)
    ) +
    labs(
      x = "ART Govt. Funding Reduction (%)",
      y = "PrEP Govt. Funding Reduction (%)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.title = element_text(face = "bold", size=text_size*4),
      axis.text = element_text(size = text_size*4),
      legend.position = "right",
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
    )


  # Calculate breaks for relative change plot (always starts at 0%)
  pct_range <- range(funding_reduction_scenarios[[pct_change_var]], na.rm = TRUE)
  pct_breaks <- seq(0, pct_range[2], length.out = n_bins + 1)

  p_sim_relative <- ggplot(
    funding_reduction_scenarios,
    aes(
      x = pct_delta_ART_fund, y = pct_delta_PrEP_fund,
      z = !!sym(pct_change_var)
    )
  ) +
    geom_contour_filled(breaks = pct_breaks) +
    geom_contour(color = "white", alpha = 0.3, breaks = pct_breaks) +
    scale_fill_viridis_d(
      option = "plasma",
      name = paste0("% Relative Change\nin ", var_name)
    ) +
    scale_x_continuous(
      labels = function(x) paste0(x, "%"),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      labels = function(x) paste0(x, "%"),
      expand = c(0, 0)
    ) +
    labs(
      x = "ART Govt. Funding Reduction (%)",
      y = "PrEP Govt. Funding Reduction (%)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.title = element_text(face = "bold", size=text_size*4),
      axis.text = element_text(size = text_size*4),
      legend.position = "right",
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
    )

   
  return(list(
    p_sim_heatmap = p_sim_heatmap,
    p_sim_relative = p_sim_relative  
    ))
}
 

#------------------------------------------------
# Matrix Heatmap 
#------------------------------------------------
#' Create discrete tile matrix heatmap with labeled cells
#'
#' @param data Data frame with funding reduction scenarios and metrics
#' @param x_var Column name for x-axis (default: "pct_delta_PrEP_fund")
#' @param y_var Column name for y-axis (default: "pct_delta_ART_fund")
#' @param fill_var Column name for fill/color values (default: "relative_incidence_change_pct")
#' @param title Plot title
#' @param x_label X-axis label (default: "PrEP Funding Level (%)")
#' @param y_label Y-axis label (default: "ART Funding Level (%)")
#' @param legend_title Legend title
#' @param text_size Size of text labels in cells (default: 3.5)
#' @param decimal_places Number of decimal places for labels (default: 1)
#' @param save_path Optional path to save plot (with extension)
#' @param width Plot width in inches (default: 12)
#' @param height Plot height in inches (default: 10)
#' @param dpi Resolution for saved images (default: 300)
#'
#' @return A ggplot object
#'
#' @examples
#' # Basic usage
#' plot_tile_matrix_heatmap(data)
#'
#' # Custom metric
#' plot_tile_matrix_heatmap(data,
#'   fill_var = "incidence_risk_ratio",
#'   title = "Risk Ratio Matrix")
#'
plot_tile_matrix_heatmap <- function(funding_reduction_scenarios, mean_var, pct_change_var, var_name) {

  # Detect if this is a risk ratio variable
  is_risk_ratio <- grepl("ratio|rr", mean_var, ignore.case = TRUE) || grepl("ratio|rr", var_name, ignore.case = TRUE)

  # Helper function to create a tile matrix plot
  create_tile_plot <- function(data, fill_var, legend_title, decimal_places = 1, text_size = 3.5) {
    # Filter to discrete 10% increments only for matrix cell display
    # This reduces from 10,201 continuous points to 121 discrete cells (11×11 grid)
    data <- data %>%
      mutate(
        art_rounded = round(pct_delta_ART_fund / 10) * 10,
        prep_rounded = round(pct_delta_PrEP_fund / 10) * 10
      ) %>%
      filter(
        abs(pct_delta_ART_fund - art_rounded) < 0.5,
        abs(pct_delta_PrEP_fund - prep_rounded) < 0.5
      ) %>%
      select(-art_rounded, -prep_rounded)

    # Create working copy with renamed columns
    plot_data <- data %>%
      select(pct_delta_PrEP_fund, pct_delta_ART_fund, all_of(fill_var)) %>%
      rename(x = pct_delta_PrEP_fund,
             y = pct_delta_ART_fund,
             fill = !!sym(fill_var))

    # Round values for display
    plot_data$label <- round(plot_data$fill, decimal_places)

    # Create the plot
    p <- ggplot(plot_data, aes(x = x, y = y, fill = fill)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = label),
                color = "grey50", size = text_size * 2, fontface = "plain") +  # Medium gray text
      scale_fill_distiller(
        palette = "YlGnBu",
        direction = 1,
        name = legend_title,
        na.value = "grey90",
        guide = guide_colorbar(
          barwidth = 2,
          barheight = 15,
          ticks.colour = "black",
          ticks.linewidth = 0.5,
          frame.colour = "black",
          frame.linewidth = 0.5
        )
      ) +
      scale_color_identity() +
      scale_x_continuous(
        breaks = seq(0, 100, 10),
        expand = c(0, 0)
      ) +
      scale_y_continuous(
        breaks = seq(0, 100, 10),
        expand = c(0, 0)
      ) +
      coord_fixed(ratio = 1) +
      labs(
        x = "PrEP Govt. Funding Reduction (%)",
        y = "ART Govt. Funding Reduction (%)"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        panel.grid = element_blank(),
        axis.text = element_text(size = text_size*5),
        axis.title = element_text(size = text_size*6, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)),
        legend.position = "right",
        legend.title = element_text(size = text_size*6, margin = margin(b = 10)),
        legend.text = element_text(size = text_size*5),
        legend.spacing.y = unit(0.5, "cm"),
        legend.key.height = unit(1.2, "cm"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
      )

    return(p)
  }

  # Plot 1: Mean tile matrix
  mean_legend_title <- if (is_risk_ratio) {
    var_name
  } else {
    paste0(var_name, "\n(per 100 PY)")
  }
  p_sim_heatmap <- create_tile_plot(funding_reduction_scenarios, mean_var, mean_legend_title)

  # Plot 2: Relative change tile matrix
  pct_change_legend_title <- paste0("% Relative Change\nin ", var_name)
  p_sim_relative <- create_tile_plot(funding_reduction_scenarios, pct_change_var, pct_change_legend_title)
 
  # Return named list of plots
  return(list(
    p_sim_heatmap = p_sim_heatmap,
    p_sim_relative = p_sim_relative
  ))
}




#---
## Functions: Add Iso-budget lines for each level of funding cut
#---

#  Iso-Budget Lines
# For each total reduction, create a line: delta_G_ART + delta_G_PrEP = Total
# G_ART_baseline/G_PrEP_baseline ~ 3.27:1
create_isobudget_lines <- function(total_reduction,
                                   G_ART_baseline,
                                   G_PrEP_baseline) {
  # For each total reduction, create a line: delta_G_ART + delta_G_PrEP = Total
  lines_data <- lapply(total_reduction, function(total) {
    # Generate ART percentage sequence (fine-grained for smooth lines)
    pct_ART_seq <- seq(0, 100, by = 1)

    # Calculate corresponding dollar reductions
    delta_G_ART <- (pct_ART_seq / 100) * G_ART_baseline
    delta_G_PrEP <- total - delta_G_ART

    # Convert PrEP dollar reduction to percentage
    pct_PrEP_seq <- (delta_G_PrEP / G_PrEP_baseline) * 100

    # Filter for feasibility: both percentages must be in [0, 100]
    feasible <- pct_ART_seq >= 0 & pct_ART_seq <= 100 &
                pct_PrEP_seq >= 0 & pct_PrEP_seq <= 100

    data.frame(
      total_reduction = total,
      iso_pct_ART = pct_ART_seq[feasible],
      iso_pct_PrEP = pct_PrEP_seq[feasible],
      delta_G_ART = delta_G_ART[feasible],
      delta_G_PrEP = delta_G_PrEP[feasible],
      feasible = TRUE
    )
  })

  result <- do.call(rbind, lines_data)
  rownames(result) <- NULL
  return(result)
}

 

# Add dollar iso-budget lines to heatmaps
add_iso_budget_lines_to_plot <- function(isobudget_data,
                                         heatmap_plots_list,
                                         total_baseline,
                                         label_type = "both") {
  # label_type: "percentage" (e.g., "10% cut"),
  #             "absolute" (e.g., "$2.4B"),
  #             "both" (e.g., "10% cut ($2.4B)")

  # Filter to feasible points only
  isobudget_data_feasible <- isobudget_data %>% filter(feasible)

  # Calculate percentage of total baseline and create labels
  isobudget_data_feasible <- isobudget_data_feasible %>%
    mutate(
      pct_of_total = (total_reduction / total_baseline) * 100,
      label_text = case_when(
        label_type == "percentage" ~ paste0(round(pct_of_total, 0), "% cut"),
        label_type == "absolute" ~ paste0("$", round(total_reduction/1e9, 1), "B"),
        label_type == "both" ~ paste0(round(pct_of_total, 0), "% cut ($", round(total_reduction/1e9, 1), "B)"),
        TRUE ~ paste0(round(pct_of_total, 0), "% cut")  # default
      )
    )

  # Plot 1: Main heatmap (percentage space with iso-budget lines)
  p1_with_isobudget_fixed <- heatmap_plots_list$p_sim_heatmap +
    geom_path(
      data = isobudget_data_feasible,
      aes(x = iso_pct_ART, y = iso_pct_PrEP, group = total_reduction),
      color = "white",
      linetype = "dashed",
      linewidth = 0.5,
      alpha = 0.6,
      inherit.aes = FALSE
    ) +
    geom_shadowtext(
      data = isobudget_data_feasible %>%
        group_by(total_reduction) %>%
        slice(n() %/% 2),
      aes(
        x = iso_pct_ART, y = iso_pct_PrEP,
        label = label_text
      ),
      angle = -atan(G_ART_baseline / G_PrEP_baseline) * 180 / pi,
      hjust = -0.1,
      size = 4.5,
      color = 'white',
      bg.color = "black",
      inherit.aes = FALSE
    )

  # Plot 2: Relative plot (percentage space with iso-budget lines)
  p2_with_isobudget_fixed <- heatmap_plots_list$p_sim_relative +
    geom_path(
      data = isobudget_data_feasible,
      aes(x = iso_pct_ART, y = iso_pct_PrEP, group = total_reduction),
      color = "white",
      linetype = "dashed",
      linewidth = 0.5,
      alpha = 0.7,
      inherit.aes = FALSE
    ) +
    geom_shadowtext(
      data = isobudget_data_feasible %>%
        group_by(total_reduction) %>%
        slice(n() %/% 2),
      aes(
        x = iso_pct_ART, y = iso_pct_PrEP,
        label = label_text
      ),
      angle = -atan(G_ART_baseline / G_PrEP_baseline) * 180 / pi,
      hjust = -0.1,
      size = 4.5,
      color = 'white',
      bg.color = "black",
      inherit.aes = FALSE
    )

  return(list(
    p_iso_main = p1_with_isobudget_fixed,
    p_iso_relative = p2_with_isobudget_fixed
  ))
}

 