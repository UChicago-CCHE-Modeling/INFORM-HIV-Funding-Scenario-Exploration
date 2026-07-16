# =============================================================================
# Morris elementary-effects screening
# -----------------------------------------------------------------------------
# Computes Morris sensitivity indices (mu, mu*, sigma) for each model parameter
# against four quantities of interest (QOIs), and saves diagnostic plots.
#
# Inputs  (in ../data/):
#   morris_design.rds          - saved Morris design objects created with
#                                sensitivity::morris(). Contains:
#                                  * morris_design     - design on the natural
#                                                        parameter scale
#                                  * morris_design_01  - design on the normalized
#                                                        [0, 1] scale (used here)
#   morris_design_outputs.csv  - one row per design point (220 rows) with the
#                                four simulated QOIs appended as columns.
#   sensitivity_params.csv     - the 21 screened parameters (name + bounds).
#
# Outputs (written to ../output/):
#   morris_scatter_plot_*.png              - mu* vs sigma scatter, faceted by QOI
#   morris_mustar_heatmap_std_*.png        - mu* heatmap (parameter x QOI)
#   morris_sigma_heatmap_*.png             - sigma heatmap
#   morris_musigma_heatmap_*.png           - mu* x sigma heatmap
#   morris_ranked_heatmap_single_plot_std_*.png - parameter ranking by mu*
#
# Run this script with its own folder as the working directory, e.g.
#   cd 01_morris_screening/script && Rscript morris.R
# =============================================================================

library(data.table)
library(ggplot2)
library(sensitivity)

# ---- Paths ------------------------------------------------------------------
morris_data_dir <- "../data/"
out_dir <- "../output"

# Make sure the output directory exists before we write to it.
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load inputs ------------------------------------------------------------
# morris_design.rds restores morris_design and morris_design_01 into the session.
load(paste0(morris_data_dir, "morris_design.rds"))
simout <- fread(paste0(morris_data_dir, "morris_design_outputs.csv"))
# Parameter definitions (names + bounds) used for labelling the results.
params <- fread(paste0(morris_data_dir, "sensitivity_params.csv"))

# ---- Standardize the simulation output --------------------------------------
# The four QOIs are on very different scales, so we z-standardize each one. This
# is combined with the normalized [0, 1] design (morris_design_01) below so that
# the resulting Morris indices are comparable across QOIs and parameters.
std <- function(x) (x - mean(x)) / sd(x)
simout_01 <- copy(simout)
simout_01[, `:=`(avg_incidence_10yr_std  = std(avg_incidence_10yr),
                 avg_incidence_5yr_std   = std(avg_incidence_5yr),
                 avg_prevalence_10yr_std = std(avg_prevalence_10yr),
                 avg_prevalence_5yr_std  = std(avg_prevalence_5yr))]

# ---- Compute elementary effects ---------------------------------------------
# For a given Morris design object and a vector of responses y, tell() computes
# the elementary effects (morris_obj$ee), from which we summarize:
#   mu      - mean elementary effect (signed; indicates direction of influence)
#   mu.star - mean of absolute elementary effects (overall influence)
#   sigma   - sd of elementary effects (non-linear / interaction effects)
get_ee <- function(morris_obj, y){
  tell(morris_obj, y)
  mu      <- apply(morris_obj$ee, 2, mean)
  mu.star <- apply(morris_obj$ee, 2, function(x) mean(abs(x)))
  sigma   <- apply(morris_obj$ee, 2, sd)

  return(cbind(mu, mu.star, sigma))
}

# Compute indices for each standardized QOI and stack them into one table.
m_dt <- data.table()
qois <- c("avg_incidence_10yr", "avg_incidence_5yr", "avg_prevalence_10yr", "avg_prevalence_5yr")
qois <- paste0(qois, "_std")
for (qoi in qois){
  m_dt <- rbind(m_dt, cbind(get_ee(morris_design_01,
                                   unlist(simout_01[, ..qoi], use.names = FALSE)), qoi))
}

# cbind() with the character qoi column coerced everything to character, so cast
# the numeric indices back and add derived/labelling columns.
m_dt[, `:=`(mu       = as.numeric(mu),
            mu.star  = as.numeric(mu.star),
            sigma    = as.numeric(sigma),
            mu.sigma = as.numeric(mu.star) * as.numeric(sigma),
            param    = rep(params$parameter, length(qois)))]

morris_results <- m_dt

# ---- Plots ------------------------------------------------------------------
# mu* vs sigma scatter (log-log), one panel per QOI.
p_scatter <- ggplot(morris_results, aes(x = mu.star, y = sigma, fill = qoi, color = qoi)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~qoi) +
  theme_bw()

ggsave(file.path(out_dir, "morris_scatter_plot.png"),
       p_scatter, width = 10, height = 8, units = "in", dpi = 300)

# mu* heatmap (parameter x QOI). theme_bw() is applied before the axis-text
# overrides so the custom text sizes are not reset by the theme.
p_mustar <- ggplot(morris_results, aes(x = qoi, y = param, fill = mu.star)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "darkblue") +
  geom_text(aes(label = round(mu.star, 2)), color = "black", size = 3) +
  labs(title = "Morris mu* Index", x = "QOI", y = "Parameter") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 8))

ggsave(file.path(out_dir, "morris_mustar_heatmap_std.png"),
       p_mustar, width = 14, height = 12, units = "in", dpi = 300)

# sigma heatmap.
p_sigma <- ggplot(morris_results, aes(x = qoi, y = param, fill = sigma)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "darkblue") +
  geom_text(aes(label = round(sigma, 2)), color = "black", size = 3) +
  labs(title = "Morris sigma Index", x = "QOI", y = "Parameter") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 8))

ggsave(file.path(out_dir, "morris_sigma_heatmap.png"),
       p_sigma, width = 14, height = 12, units = "in", dpi = 300)

# mu* x sigma heatmap.
p_musigma <- ggplot(morris_results, aes(x = qoi, y = param, fill = mu.sigma)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "darkblue") +
  geom_text(aes(label = round(mu.sigma, 2)), color = "black", size = 3) +
  labs(title = "Morris mu* x sigma Index", x = "QOI", y = "Parameter") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 8))

ggsave(file.path(out_dir, "morris_musigma_heatmap.png"),
       p_musigma, width = 14, height = 12, units = "in", dpi = 300)

# ---- Parameter ranking ------------------------------------------------------
# Rank parameters by mu* within each QOI (rank 1 = most influential).
morris_ranked_by_qoi <- morris_results[, .(param, mu.star,
                                           rank = rank(-mu.star, ties.method = "first")),
                                       by = qoi]

# Single tile plot: rows are ranks, cells labelled with the parameter name.
p_single_ranked_plot <- ggplot(morris_ranked_by_qoi,
                               aes(x = qoi, y = rank, label = param, fill = mu.star)) +
  geom_tile(color = "black", linewidth = 0.5) +
  geom_text(aes(color = ifelse(mu.star > 2, "white", "black")), size = 3) +
  scale_y_reverse(breaks = 1:nrow(params), name = "Rank") +
  scale_fill_distiller(palette = "YlGnBu", direction = 1, name = "mu*") +
  scale_color_identity() +  # use the literal color names from the ifelse above
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Parameter Importance Ranking by QOI",
    x = "Quantity of Interest (QOI)"
  )

ggsave(file.path(out_dir, "morris_ranked_heatmap_single_plot_std.png"),
       p_single_ranked_plot, width = 12, height = 12, units = "in", dpi = 300)
