# =============================================================================
# Forest plot of intervention-scenario incidence risk ratios
# -----------------------------------------------------------------------------
# Loads the fitted surrogate bundle produced by surrogate.R and summarizes a set
# of (ART, PrEP) funding-reduction scenarios as a forest plot of the mean HIV
# incidence risk ratio (scenario vs. no-reduction baseline), with 90% credible
# intervals. Three facets:
#   - vary ART reduction, PrEP fixed at 0
#   - vary PrEP reduction, ART fixed at 0
#   - vary both by the same amount
#
# Input  (from ../output/): surrogate.Rdata
# Output (to   ../output/): forest_incidence_risk_ratio.png
#
# Run with script/ as the working directory:
#   cd 03_intervention_scenario_surrogate/script && Rscript forest_plot.R
# =============================================================================

library(data.table)
library(hetGP)     # needed so predict() dispatches on the stored hetGP models
library(ggplot2)

set.seed(0)  # surrogate prediction draws are stochastic

out_dir <- "../output"

# ---- Load the stored surrogate object ---------------------------------------
# Restores gp_prevalence_fit, gp_incidence_fit, baseline_incidence_prevalence,
# and the predict_* helper functions.
load(file.path(out_dir, "surrogate.Rdata"))

# ---- Summarize one scenario as a risk ratio ---------------------------------
# art_red / prep_red are funding-reduction fractions in [0, 0.75]; the surrogate
# takes their negatives as native inputs. Returns the mean incidence risk ratio
# over the horizon (years 20-30) with a 90% credible interval.
incidence_risk_ratio <- function(art_red, prep_red, n_samples_per_checkpoint = 50){
    samples <- predict_hiv_diff_composite_surrogate(
        c(-art_red, -prep_red), gp_incidence_fit, baseline_incidence_prevalence,
        qoi = "hiv_total_incidence_per_100py",
        n_samples_per_checkpoint = n_samples_per_checkpoint, diff_type = "times")
    # samples: n_ticks x 1 x n_draws -> mean risk ratio over the horizon per draw
    rr_by_draw <- apply(samples[, 1, ], 2, mean)
    data.table(mean  = mean(rr_by_draw),
               lower = quantile(rr_by_draw, 0.05),
               upper = quantile(rr_by_draw, 0.95))
}

# ---- Define the scenarios ---------------------------------------------------
red_levels <- c(0, 0.10, 0.20, 0.30, 0.40)  # funding-reduction fractions
scenarios <- rbind(
    data.table(scenario = "Vary ART (PrEP = 0)", art = red_levels, prep = 0),
    data.table(scenario = "Vary PrEP (ART = 0)", art = 0,          prep = red_levels),
    data.table(scenario = "Vary both equally",   art = red_levels, prep = red_levels)
)
scenarios[, reduction := pmax(art, prep)]  # the quantity being varied in each facet

# ---- Compute risk ratios for every scenario ---------------------------------
forest_dt <- scenarios[, incidence_risk_ratio(art, prep),
                       by = .(scenario, reduction, art, prep)]
forest_dt[, scenario := factor(scenario, levels = c("Vary PrEP (ART = 0)",
                                                    "Vary ART (PrEP = 0)",
                                                    "Vary both equally"))]

# ---- Forest plot ------------------------------------------------------------
p_forest <- ggplot(forest_dt,
                   aes(x = mean, y = factor(reduction * 100), color = reduction * 100)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
    geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = 0.25) +
    geom_point(size = 3) +
    facet_wrap(~scenario, ncol = 1) +
    # warm sequential: deeper funding cut reads as more "danger"
    scale_color_gradient(low = "#fdd0a2", high = "#a63603", name = "Reduction (%)") +
    labs(title = "Mean Incidence Risk Ratio by Intervention-Reduction Scenario",
         subtitle = "Scenario vs. no-reduction baseline; points are means, whiskers are 90% CIs",
         x = "Mean incidence risk ratio (scenario vs. baseline)",
         y = "Reduction") +
    theme_bw() +
    theme(legend.position = "right")

ggsave(file.path(out_dir, "forest_incidence_risk_ratio.png"), p_forest,
       width = 8, height = 9, units = "in", dpi = 300)

# ---- Grouped forest plot (single panel) -------------------------------------
# Same data, but colour encodes the scenario type (ART only / PrEP only / both)
# and the three scenarios are dodged together at each reduction level, so they
# can be compared directly at a fixed funding reduction.
pd <- position_dodge(width = 0.7)

# drop the 0% row (baseline anchor: all three scenarios are 1.00)
forest_dt_grouped <- forest_dt[reduction > 0]

# fixed x-axis window (RR units). We use finite band bounds tied to these
# limits: on a log axis xmin = -Inf becomes NaN and the bands silently vanish.
x_lo <- min(1, min(forest_dt_grouped$lower)) / 1.01
x_hi <- max(forest_dt_grouped$upper) * 1.18

# one shaded band per funding-reduction level (with thin white gaps between).
# Warm sequential fill encodes severity: pale at low reduction, warmer (more
# "danger") as the funding cut deepens.
band_pos <- seq_len(length(unique(forest_dt_grouped$reduction)))
band_fill <- colorRampPalette(c("#fff5ec", "#fdc58a"))(length(band_pos))
bands <- data.frame(ymin = band_pos - 0.42, ymax = band_pos + 0.42, fill = band_fill)

p_forest_grouped <- ggplot(forest_dt_grouped,
                           aes(x = mean, y = factor(reduction * 100), color = scenario)) +
    geom_rect(data = bands, inherit.aes = FALSE,
              aes(ymin = ymin, ymax = ymax, xmin = x_lo, xmax = x_hi, fill = fill)) +
    scale_fill_identity() +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.7) +
    geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y",
                  width = 0.4, linewidth = 1, position = pd) +
    geom_point(size = 4, position = pd) +
    # value labels just past the right whisker; dark ink so they read on any band.
    # group = scenario keeps the dodge (a fixed colour would otherwise drop the
    # grouping and stack all labels on one line).
    geom_text(aes(x = upper, label = sprintf("%.2f", mean), group = scenario),
              position = pd, hjust = -0.3, size = 2.9, fontface = "bold",
              color = "#0b0b0b", show.legend = FALSE) +
    # cool categorical hues (validated CVD-safe) so the data pops off the warm bands
    scale_color_manual(values = c("Vary PrEP (ART = 0)" = "#008300",
                                  "Vary ART (PrEP = 0)" = "#2a78d6",
                                  "Vary both equally"   = "#4a3aa7"),
                       name = "Scenario") +
    # natural-log scale for the risk ratio, labelled in RR units;
    # fixed limits so the shaded bands fill the panel edge to edge
    scale_x_continuous(transform = "log", breaks = c(1, 1.2, 1.5, 2),
                       limits = c(x_lo, x_hi), expand = expansion(mult = 0)) +
    labs(title = "Mean Incidence Risk Ratio by Funding Reduction",
         subtitle = "Scenario vs. no-reduction baseline; points are means, whiskers are 90% CIs",
         x = "Mean incidence risk ratio (scenario vs. baseline, natural-log scale)",
         y = "Reduction (%)") +
    theme_bw() +
    # legend as a horizontal row at the bottom; drop gridlines (the shaded
    # bands already separate the reduction levels) and the surrounding box,
    # keeping clean axis lines instead
    theme(legend.position = "bottom",
          panel.grid = element_blank(),
          panel.border = element_blank(),
          # keep only the bottom (x) axis line; drop the y line but keep its
          # ticks and labels — the shaded bands carry the vertical structure
          axis.line.x = element_line(color = "grey30", linewidth = 0.4),
          axis.line.y = element_blank())

ggsave(file.path(out_dir, "forest_incidence_risk_ratio_grouped.png"), p_forest_grouped,
       width = 8, height = 5, units = "in", dpi = 300)

print(forest_dt)
