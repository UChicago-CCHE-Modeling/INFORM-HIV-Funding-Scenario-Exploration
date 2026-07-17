# =============================================================================
# Surrogate scenario draws with optional common random numbers (CRN)
# -----------------------------------------------------------------------------
# One function, compute_scenario_draws_crn(), produces everything the figure and
# the table need for a set of scenarios: absolute incidence draws, incidence
# risk ratio (IRR) draws, and the shared baseline incidence draws -- all from the
# same predictive samples.
#
# The IRR is scenario_incidence / baseline_incidence. When the surrogate is
# queried at the scenario point and the baseline point with INDEPENDENT
# predictive draws, the ratio inherits the surrogate's full predictive noise in
# BOTH numerator and denominator. For scenarios whose true effect is small (e.g.
# modest PrEP funding cuts, which map to only a few percent coverage reduction),
# that noise floor is larger than the effect, so the 90% interval dips below 1
# even though the mean IRR is monotonically above 1. The self-ratio at (0,0)
# makes this obvious: with independent draws it spreads over roughly
# [0.98, 1.02] instead of collapsing to exactly 1.
#
# Common random numbers fix this: draw ONE set of standard-normal variates per
# (checkpoint, SVD component, sample) and reuse it for BOTH the baseline and the
# scenario prediction. The shared draws cancel in the ratio, so (0,0)/(0,0) = 1
# exactly and the interval reflects only genuine scenario-vs-baseline
# differences (the mean IRR is unchanged). This is the standard trick for
# reducing Monte Carlo variance when comparing two simulation configurations.
#
# `common_random_numbers = FALSE` reproduces the classic independent-draw
# behaviour; toggling the flag lets us compare the two on the same footing for
# both the figure and the table.
# =============================================================================

#' Predict incidence trajectories for a set of points, optionally sharing draws.
#'
#' Single-checkpoint predictor mirroring predict_hiv_single_surrogate, with one
#' addition: if a matrix of standard-normal variates `z` (ngp x n_samples) is
#' supplied, every point uses those SAME variates (scaled by each point's
#' predictive sd), which is what makes the ratio share randomness across points.
#' When `z` is NULL each point draws its own independent variates (classic).
#'
#' @param newX Points on the [0, 1] GP scale (n_points x 2).
#' @param gp_list One checkpoint surrogate (gps, Kmat, Y_mean, Y_sd).
#' @param n_samples Draws per point.
#' @param z Optional ngp x n_samples standard-normal matrix shared across points.
#' @return Array n_ticks x n_points x n_samples of incidence trajectories.
.predict_trajectories_shared_z <- function(newX, gp_list, n_samples, z = NULL) {
  gps    <- gp_list$gps
  Kmat   <- gp_list$Kmat
  Y_mean <- gp_list$Y_mean
  Y_sd   <- gp_list$Y_sd

  ngp      <- length(gps)
  n_points <- nrow(newX)
  n_ticks  <- nrow(Kmat)

  W_samples <- array(NA, dim = c(ngp, n_points, n_samples))
  for (ii in 1:ngp) {
    pred     <- predict(gps[[ii]], newX)
    pred_sd  <- sqrt(pred$sd2)
    for (j in 1:n_points) {
      if (is.null(z)) {
        W_samples[ii, j, ] <- rnorm(n_samples, mean = pred$mean[j], sd = pred_sd[j])
      } else {
        # Shared variates: same z for every point, scaled by this point's sd.
        W_samples[ii, j, ] <- pred$mean[j] + z[ii, ] * pred_sd[j]
      }
    }
  }

  Y_samples <- array(NA, dim = c(n_ticks, n_points, n_samples))
  for (s in 1:n_samples) {
    w_sample <- matrix(W_samples[, , s], nrow = ngp)  # ngp x n_points (keep dims)
    Y_samples[, , s] <- (Kmat %*% w_sample) * Y_sd +
      matrix(Y_mean, nrow = n_ticks, ncol = n_points)
  }
  Y_samples
}

#' Surrogate scenario draws: absolute incidence, IRR, and shared baseline.
#'
#' For each query point, returns the surrogate posterior draws of the absolute
#' incidence and the incidence risk ratio relative to the (0,0) no-reduction
#' baseline, pooled over all checkpoints, plus the baseline incidence draws
#' themselves. The baseline is the surrogate prediction at (0,0) (not the
#' observed baseline) so that, under CRN, its draws share randomness with the
#' scenario draws and cancel in the ratio. The baseline draws are returned so
#' callers can form the absolute increase (scenario - baseline) sample-matched.
#'
#' @param newX_native Query points in native coverage-reduction scale
#'   (n_points x 2, e.g. c(-0.05, 0) for a 5% ART reduction). Negative values
#'   are reductions, matching the surrogate's native inputs.
#' @param gp_fit_list List of per-checkpoint incidence surrogates (gp_incidence_fit).
#' @param n_samples_per_checkpoint Draws per checkpoint.
#' @param baseline_point Native coordinates of the baseline (default c(0, 0)).
#' @param common_random_numbers If TRUE (default), baseline and scenario share
#'   the same predictive draws within each checkpoint so shared noise cancels.
#'   If FALSE, draws are independent (classic behaviour).
#' @param tick Reporting horizon. If NULL (default), quantities are averaged over
#'   all trajectory ticks per draw (the forest-plot convention); if an integer,
#'   the single tick is used (the scenario-table convention).
#' @return List: `incidence` and `irr` are n_points x (n_samples_per_checkpoint *
#'   n_checkpoints) matrices; `baseline` is the length-n_draws vector of baseline
#'   incidence draws.
compute_scenario_draws_crn <- function(newX_native,
                                       gp_fit_list,
                                       n_samples_per_checkpoint,
                                       baseline_point = c(0, 0),
                                       common_random_numbers = TRUE,
                                       tick = NULL) {
  if (!is.matrix(newX_native)) newX_native <- matrix(newX_native, nrow = 1)
  n_points <- nrow(newX_native)
  n_ckpt   <- length(gp_fit_list)
  n_draws  <- n_samples_per_checkpoint * n_ckpt

  # Baseline first, scenarios after; map to the [0, 1] GP scale.
  allX_native <- rbind(baseline_point, newX_native)
  allX <- apply(allX_native, 2, function(x) (x + 0.75) / 0.75)
  if (!is.matrix(allX)) allX <- matrix(allX, nrow = 1)

  incidence <- matrix(NA, nrow = n_points, ncol = n_draws)
  irr       <- matrix(NA, nrow = n_points, ncol = n_draws)
  baseline  <- rep(NA, n_draws)

  reduce_tick <- function(m) if (is.null(tick)) colMeans(m) else m[tick, ]

  for (s in 1:n_ckpt) {
    gp_list <- gp_fit_list[[s]]
    ngp     <- length(gp_list$gps)

    # One shared standard-normal block per checkpoint enables CRN across points.
    z <- if (common_random_numbers) {
      matrix(rnorm(ngp * n_samples_per_checkpoint), nrow = ngp)
    } else {
      NULL
    }

    Y <- .predict_trajectories_shared_z(allX, gp_list, n_samples_per_checkpoint, z = z)
    # Y: n_ticks x (1 + n_points) x n_samples. Column 1 is the baseline.
    n_ticks <- dim(Y)[1]
    base <- matrix(Y[, 1, ], nrow = n_ticks)  # n_ticks x n_samples

    idx <- ((s - 1) * n_samples_per_checkpoint + 1):(s * n_samples_per_checkpoint)
    baseline[idx] <- reduce_tick(base)
    for (k in 1:n_points) {
      scen <- matrix(Y[, k + 1, ], nrow = n_ticks)  # n_ticks x n_samples
      incidence[k, idx] <- reduce_tick(scen)
      irr[k, idx]       <- reduce_tick(scen / base)
    }
  }

  list(incidence = incidence, irr = irr, baseline = baseline)
}

#' Summarise draws (per row) into mean and a credible interval.
#'
#' @param draws Numeric matrix (n_points x n_draws) or a single vector.
#' @param probs Length-2 lower/upper quantile probabilities.
#' @return data.frame with mean, ci_lower, ci_upper, ci_width per row.
summarise_draws <- function(draws, probs = c(0.05, 0.95)) {
  if (!is.matrix(draws)) draws <- matrix(draws, nrow = 1)
  data.frame(
    mean     = apply(draws, 1, mean),
    ci_lower = apply(draws, 1, quantile, probs = probs[1]),
    ci_upper = apply(draws, 1, quantile, probs = probs[2]),
    ci_width = apply(draws, 1, quantile, probs = probs[2]) -
               apply(draws, 1, quantile, probs = probs[1]))
}
