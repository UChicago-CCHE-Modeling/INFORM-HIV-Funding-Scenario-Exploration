# =============================================================================
# Intervention-scenario surrogate model
# -----------------------------------------------------------------------------
# Builds a Gaussian-process (GP) surrogate that predicts the HIV epidemic
# trajectory (prevalence and incidence over years 20-30) as a function of two
# intervention levers:
#   - art.yearly.adjustment   (yearly change in ART coverage)
#   - prep.yearly.adjustment  (yearly change in PrEP coverage)
#
# The training data are 20,000 model simulations run in 4 batches of 5,000 runs
# each. Each batch lives in its own folder under ../data/ (batch_1 .. batch_4,
# in ascending run order) and contains three files:
#   proportion_summary.csv - prevalence-style summary outcomes per (instance, ycat)
#   parameters.csv         - the full parameter table, one row per instance
#   incidence_data.csv     - incidence outcomes per (instance, ycat)
#
# Each simulation instance combines a "checkpoint" (a draw of six calibrated
# epidemic parameters, giving 100 checkpoints total, 25 per batch) with an
# intervention design point (an (art, prep) adjustment pair). For every
# checkpoint we fit an SVD + GP emulator over the intervention space, and the
# 100 per-checkpoint emulators are combined into a composite surrogate that
# marginalizes over calibration uncertainty.
#
# Output (written to ../output/):
#   surrogate.Rdata            - the fitted surrogates, baseline data, and the
#                                predict_* helper functions needed to use them.
#   pred_prevalence.png        - example prevalence prediction at one scenario
#   pred_incidence.png         - example incidence prediction at one scenario
#   pred_incidence_diff.png    - example incidence-vs-baseline difference
#
# Run this script with its own folder as the working directory, e.g.
#   cd 03_intervention_scenario_surrogate/script && Rscript surrogate.R
# =============================================================================

library(data.table)
library(hetGP)

# Seed for reproducibility (surrogate prediction draws from rnorm/sample).
set.seed(0)

# ---- Paths ------------------------------------------------------------------
out_dir <- "../output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Read and aggregate one batch of simulations ----------------------------
read_data <- function(data_dir){

    sim_dt       <- fread(paste0(data_dir, "proportion_summary.csv"))
    param_dt     <- fread(paste0(data_dir, "parameters.csv"))
    incidence_dt <- fread(paste0(data_dir, "incidence_data.csv"))

    # Inputs: the two intervention levers plus the six calibrated parameters
    # that define each checkpoint (kept so we can group by checkpoint below).
    input_param <- c("instance", "transmission.probability.insertive.infected",
                                    "prop.casual.sex.acts",
                                    "prop.steady.sex.acts",
                                    "prep.bl.use.prop",
                                    "prep.fraction.sd",
                                    "transmission.probability.receptive.infected",
                                    "art.yearly.adjustment",
                                    "prep.yearly.adjustment")
    param_sub <- param_dt[, ..input_param]

    # Keep only the post-intervention horizon (ycat > 21).
    sim_prevalence_dt <- sim_dt[ycat > 21, .(instance, mean_hiv_prevalence, ycat)]
    setkeyv(sim_prevalence_dt, c("instance", "ycat"))

    sim_incidence_dt <- incidence_dt[ycat > 21, .(instance, hiv_total_incidence_per_100py, ycat)]
    setkeyv(sim_incidence_dt, c("instance", "ycat"))

    dt_all <- merge(param_sub, sim_prevalence_dt, by = "instance")
    dt_all <- merge(dt_all, sim_incidence_dt[, .(instance, hiv_total_incidence_per_100py, ycat)], by = c("instance", "ycat"))
    setDT(dt_all)

    # param_id: unique id for each intervention design point (art, prep).
    dt_all[, param_id := .GRP, by = .(art.yearly.adjustment,
                                    prep.yearly.adjustment)]

    # post_id: unique id for each checkpoint (draw of the 6 calibrated params).
    dt_all[, post_id := .GRP, by = .(transmission.probability.insertive.infected,
                                    prop.casual.sex.acts,
                                    prop.steady.sex.acts,
                                    prep.bl.use.prop,
                                    prep.fraction.sd,
                                    transmission.probability.receptive.infected)]

    # Average replicate runs within each (checkpoint, design point, year).
    mean_hiv_dt <- dt_all[, .(mean_hiv_prevalence = mean(mean_hiv_prevalence),
                                hiv_total_incidence_per_100py = mean(hiv_total_incidence_per_100py),
                                art.yearly.adjustment = mean(art.yearly.adjustment),
                                prep.yearly.adjustment = mean(prep.yearly.adjustment)),
                                by = c("post_id", "param_id", "ycat")]

    return(mean_hiv_dt)
}

mean_hiv_dt_1 <- read_data("../data/batch_1/")
mean_hiv_dt_2 <- read_data("../data/batch_2/")
mean_hiv_dt_3 <- read_data("../data/batch_3/")
mean_hiv_dt_4 <- read_data("../data/batch_4/")


# ---- Fit an SVD + GP surrogate for a single checkpoint -----------------------
# Reduces the multivariate output (trajectory over years) to a few SVD basis
# components, then fits an independent heteroscedastic GP to each component's
# weights over the (art, prep) intervention space.
gp_surrogate <- function(mean_hiv_dt, qoi){

    X <- unique(mean_hiv_dt[, c("art.yearly.adjustment", "prep.yearly.adjustment")])
    Y_vec <- mean_hiv_dt[[qoi]]
    Y_mat <- matrix(Y_vec, nrow = 11)

    ## normalize: subtract per-row mean, divide by overall sd
    Y_mean <- apply(Y_mat, 1, mean)
    Y_sd <- sd(Y_vec)

    Y_mat_std <- (Y_mat - matrix(Y_mean, nrow = 11, ncol = 20)) / Y_sd

    # Map the intervention inputs to the [0, 1] scale used for GP fitting.
    X01 <- as.matrix(apply(X, 2, function(x) (x + 0.75)/(0.75)))

    ## SVD of the standardized output matrix
    udv <- svd(Y_mat_std)

    K <- udv$u %*% diag(udv$d) * sqrt(ncol(Y_mat_std))

    ## keep enough components to explain 99% of the spectrum (cap at 5)
    nK <- which(cumsum(udv$d / sum(udv$d)) > .99)[1] - 1
    nK <- min(nK, 5)

    K <- K[, 1:nK]

    W <- udv$v[, 1:nK] / sqrt(ncol(Y_mat_std))

    ## fit a GP to each column of the basis weights W
    gps <- list()
    for (ii in 1:nK){
        gps[[ii]] <- mleHetGP(X01, W[, ii], known = list(beta0 = 0), covtype = "Matern5_2")
    }

    return(list(gps = gps, Kmat = K, Y_mean = Y_mean, Y_sd = Y_sd))
}

# Fit one surrogate per checkpoint (post_id) within a batch.
build_surrogates <- function(dt, qoi) {
  results <- list()
  for (pid in unique(dt$post_id)) {
    results[[as.character(pid)]] <- gp_surrogate(dt[post_id == pid], qoi)
  }
  return(results)
}

# Fit surrogates for both QOIs, for each batch.
gp_prevalence_1 <- build_surrogates(mean_hiv_dt_1, qoi = "mean_hiv_prevalence")
gp_prevalence_2 <- build_surrogates(mean_hiv_dt_2, qoi = "mean_hiv_prevalence")
gp_prevalence_3 <- build_surrogates(mean_hiv_dt_3, qoi = "mean_hiv_prevalence")
gp_prevalence_4 <- build_surrogates(mean_hiv_dt_4, qoi = "mean_hiv_prevalence")

gp_incidence_1 <- build_surrogates(mean_hiv_dt_1, qoi = "hiv_total_incidence_per_100py")
gp_incidence_2 <- build_surrogates(mean_hiv_dt_2, qoi = "hiv_total_incidence_per_100py")
gp_incidence_3 <- build_surrogates(mean_hiv_dt_3, qoi = "hiv_total_incidence_per_100py")
gp_incidence_4 <- build_surrogates(mean_hiv_dt_4, qoi = "hiv_total_incidence_per_100py")

# Combine the 4 batches into one list of 100 per-checkpoint surrogates per QOI.
gp_prevalence_fit <- c(gp_prevalence_1, gp_prevalence_2, gp_prevalence_3, gp_prevalence_4)
gp_incidence_fit  <- c(gp_incidence_1, gp_incidence_2, gp_incidence_3, gp_incidence_4)

# ---- Prediction: single checkpoint ------------------------------------------
predict_hiv_single_surrogate <- function(newX, gp_list, n_samples = 50){
    # gp_list: surrogate model for one post_id (contains gps, Kmat, Y_mean, Y_sd)
    # newX: matrix of new points to predict (n_points x 2), on the [0, 1] scale
    # n_samples: number of samples from the predictive distribution

    gps <- gp_list$gps
    Kmat <- gp_list$Kmat
    Y_mean <- gp_list$Y_mean
    Y_sd <- gp_list$Y_sd

    ngp <- length(gps)
    n_points <- nrow(newX)
    n_ticks <- nrow(Kmat)

    # W_samples: predictive draws of the basis weights (ngp x n_points x n_samples)
    W_samples <- array(NA, dim = c(ngp, n_points, n_samples))

    # Sample from the predictive distribution for each GP component.
    for (ii in 1:ngp) {
        pred_obj <- predict(gps[[ii]], newX)
        w_pred_mean <- pred_obj$mean
        w_pred_var <- pred_obj$sd2

        for (j in 1:n_points) {
            W_samples[ii, j, ] <- rnorm(n_samples, mean = w_pred_mean[j], sd = sqrt(w_pred_var[j]))
        }
    }

    # Reconstruct trajectories from the basis: Y = K %*% W (de-standardized).
    Y_samples <- array(NA, dim = c(n_ticks, n_points, n_samples))

    for (s in 1:n_samples) {
        w_sample <- W_samples[, , s]  # ngp x n_points
        Y_samples[, , s] <- (Kmat %*% w_sample) * Y_sd + matrix(Y_mean, nrow = length(Y_mean), ncol = n_points)
    }

    return(Y_samples)  # n_ticks x n_points x n_samples
}


# ---- Prediction: composite over all checkpoints -----------------------------
predict_hiv_composite_surrogate <- function(newX_native, gp_fit_list, n_samples_per_checkpoint) {
    # gp_fit_list: list of 100 surrogate models (each with gps, Kmat, Y_mean, Y_sd)
    # newX_native: matrix of new points to predict in native scale (n_points x 2)
    # n_samples_per_checkpoint: draws from the predictive distribution per checkpoint

    if(!is.matrix(newX_native)) newX_native <- matrix(newX_native, nrow = 1)

    n_surrogates <- length(gp_fit_list)
    n_points <- nrow(newX_native)
    n_ticks <- nrow(gp_fit_list[[1]]$Kmat)

    Y_samples <- array(NA, dim = c(n_ticks, n_points, n_samples_per_checkpoint * n_surrogates))

    # transform newX_native to the [0, 1] scale
    newX <- apply(newX_native, 2, function(x) (x + 0.75)/(0.75))
    if(!is.matrix(newX)) newX <- matrix(newX, nrow = 1)

    # Stack draws from every checkpoint surrogate.
    for (s in 1:n_surrogates) {
        samples <- predict_hiv_single_surrogate(newX, gp_fit_list[[s]], n_samples = n_samples_per_checkpoint)
        idx <- ((s - 1) * n_samples_per_checkpoint + 1) : (s * n_samples_per_checkpoint)
        Y_samples[, , idx] <- samples
    }

  return(Y_samples)
}


# ---- Example: predict prevalence trajectory at one scenario -----------------
# Scenario: (art, prep) yearly reduction = (50%, 40%), i.e. native (-0.5, -0.4).
newX_native <- c(-0.5, -0.4)
samples <- predict_hiv_composite_surrogate(newX_native, gp_prevalence_fit, n_samples_per_checkpoint = 50)

pred_mean <- apply(samples, c(1, 2), mean)
pred_ci <- apply(samples, c(1, 2), quantile, probs = c(0.05, 0.95))

png(file.path(out_dir, "pred_prevalence.png"), width = 8, height = 6, units = "in", res = 300)
matplot(20:30, samples[, 1, sample(1000, 500)], type = "l", col = "cyan", xlab = "years", ylab = "mean HIV prevalence", main = "(art, prep) reduction = (50%, 40%)")
lines(20:30, pred_mean, col = "blue", lwd = 1.5)
matlines(20:30, t(pred_ci[,,1]), lty = 2, col = "red", type = "l")
legend("topleft", c("mean", "90% CI", "samples"), col = c("blue", "red", "cyan"), lty = c(1, 2, 1))
dev.off()


# ---- Example: predict incidence trajectory at one scenario ------------------
newX_native <- c(-0.5, -0.4)
samples <- predict_hiv_composite_surrogate(newX_native, gp_incidence_fit, n_samples_per_checkpoint = 50)

pred_mean <- apply(samples, c(1, 2), mean)
pred_ci <- apply(samples, c(1, 2), quantile, probs = c(0.05, 0.95))

png(file.path(out_dir, "pred_incidence.png"), width = 8, height = 6, units = "in", res = 300)
matplot(20:30, samples[, 1, sample(1000, 500)], type = "l", col = "cyan", xlab = "years", ylab = "mean HIV incidence", main = "(art, prep) reduction = (50%, 40%)")
lines(20:30, pred_mean, col = "blue", lwd = 1.5)
matlines(20:30, t(pred_ci[,,1]), lty = 2, col = "red", type = "l")
legend("topleft", c("mean", "90% CI", "samples"), col = c("blue", "red", "cyan"), lty = c(1, 2, 1))
dev.off()


# ---- Prediction relative to the no-intervention baseline --------------------
# Gather the baseline (art = 0, prep = 0) trajectories from the actual simulations.
baseline_incidence_prevalence <- rbind(mean_hiv_dt_1[art.yearly.adjustment == 0 & prep.yearly.adjustment == 0],
                                        mean_hiv_dt_2[art.yearly.adjustment == 0 & prep.yearly.adjustment == 0],
                                        mean_hiv_dt_3[art.yearly.adjustment == 0 & prep.yearly.adjustment == 0],
                                        mean_hiv_dt_4[art.yearly.adjustment == 0 & prep.yearly.adjustment == 0])

# post_id resets to 1 within each batch, so re-index to a global 1..100.
baseline_incidence_prevalence[1:(11*25), post_id := rep(1:25, each = 11)]
baseline_incidence_prevalence[((11*25*1)+1):(11*25*2), post_id := rep(26:50, each = 11)]
baseline_incidence_prevalence[((11*25*2)+1):(11*25*3), post_id := rep(51:75, each = 11)]
baseline_incidence_prevalence[((11*25*3)+1):(11*25*4), post_id := rep(76:100, each = 11)]


predict_hiv_diff_composite_surrogate <- function(newX_native, gp_fit_list, baseline_data, qoi, n_samples_per_checkpoint, diff_type = "percent") {
    # gp_fit_list: list of 100 surrogate models (each with gps, Kmat, Y_mean, Y_sd)
    # newX_native: matrix of new points to predict in native scale (n_points x 2)
    # baseline_data: baseline trajectories keyed by post_id
    # qoi: "mean_hiv_prevalence" or "hiv_total_incidence_per_100py"
    # n_samples_per_checkpoint: draws from the predictive distribution per checkpoint
    # diff_type: "percent" -> (x - baseline)/baseline; otherwise x/baseline

    if(!is.matrix(newX_native)) newX_native <- matrix(newX_native, nrow = 1)

    n_surrogates <- length(gp_fit_list)
    n_points <- nrow(newX_native)
    n_ticks <- nrow(gp_fit_list[[1]]$Kmat)

    Y_samples <- array(NA, dim = c(n_ticks, n_points, n_samples_per_checkpoint * n_surrogates))

    # transform newX_native to the [0, 1] scale
    newX <- apply(newX_native, 2, function(x) (x + 0.75)/(0.75))
    if(!is.matrix(newX)) newX <- matrix(newX, nrow = 1)

    # For each checkpoint, express the prediction relative to its baseline.
    for (s in 1:n_surrogates) {
        samples <- predict_hiv_single_surrogate(newX, gp_fit_list[[s]], n_samples = n_samples_per_checkpoint)
        baseline <- unlist(baseline_data[post_id == s, ..qoi], use.names = F)
        if(diff_type == "percent"){
            qoi_diff <- sweep(samples, MARGIN = 1, STATS = baseline, FUN = function(x, y) (x - y)/y)
        } else {
            qoi_diff <- sweep(samples, MARGIN = 1, STATS = baseline, FUN = function(x, y) x/y)
        }
        idx <- ((s - 1) * n_samples_per_checkpoint + 1) : (s * n_samples_per_checkpoint)
        Y_samples[, , idx] <- qoi_diff
    }

  return(Y_samples)
}


# ---- Example: incidence difference from baseline at one scenario ------------
newX_native <- c(-0.5, -0.4)
samples <- predict_hiv_diff_composite_surrogate(newX_native, gp_incidence_fit, baseline_incidence_prevalence,
                                        qoi = "hiv_total_incidence_per_100py", n_samples_per_checkpoint = 50)

pred_mean <- apply(samples, c(1, 2), mean)
pred_ci <- apply(samples, c(1, 2), quantile, probs = c(0.05, 0.95))

png(file.path(out_dir, "pred_incidence_diff.png"), width = 8, height = 6, units = "in", res = 300)
matplot(20:30, samples[, 1, sample(1000, 500)], type = "l", col = "cyan", xlab = "years", ylab = "difference in mean HIV incidence from baseline", main = "(art, prep) reduction = (50%, 40%)")
lines(20:30, pred_mean, col = "blue", lwd = 1.5)
matlines(20:30, t(pred_ci[,,1]), lty = 2, col = "red", type = "l")
legend("topleft", c("mean", "90% CI", "samples"), col = c("blue", "red", "cyan"), lty = c(1, 2, 1))
dev.off()


# ---- Save the surrogate bundle ----------------------------------------------
# Items needed to use the surrogate later:
#   - gp_prevalence_fit / gp_incidence_fit : 100 multivariate-output GP surrogates each
#   - baseline_incidence_prevalence        : baseline incidence/prevalence data table
#   - predict_hiv_single_surrogate         : single-checkpoint prediction helper
#   - predict_hiv_composite_surrogate      : composite prediction helper
#   - predict_hiv_diff_composite_surrogate : baseline-difference prediction helper
save(gp_prevalence_fit, gp_incidence_fit, baseline_incidence_prevalence,
     predict_hiv_single_surrogate, predict_hiv_composite_surrogate, predict_hiv_diff_composite_surrogate,
     file = file.path(out_dir, "surrogate.Rdata"))
