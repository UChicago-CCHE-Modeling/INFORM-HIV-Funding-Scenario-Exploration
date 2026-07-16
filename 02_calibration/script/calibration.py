"""
Model calibration via Simulation-Based Inference (SBI).

Uses a neural posterior estimator (SNPE) to calibrate six HIV-model parameters
to an observed HIV prevalence target. The training set is a Latin-hypercube
design of simulations; the estimator learns the mapping from (normalized)
parameters to (standardized) prevalence, and the posterior is then conditioned
on the observed target to recover plausible parameter values.

Inputs (in ../data/):
  proportion_summary.csv - simulated summary outcomes, one row per
                           (instance, year-category). We use the final
                           year-category prevalence per instance.
  parameters.csv         - full parameter table, one row per simulation
                           instance (wide format); the calibration parameters
                           are a subset of its columns.

Outputs (written to ../output/):
  posterior_samples_mean_10000.csv - 10,000 posterior draws (native scale)
  posterior_marginals_mean.png     - 1D marginal histograms
  posterior_pairplot_mean.png      - seaborn pairplot of the draws
  prior_vs_posterior_mean.png      - training (prior) vs posterior per parameter
  posterior_pairplot_kde_mean.png  - KDE pairplot (1D densities + 2D contours)

Run with script/ as the working directory:
    cd 02_calibration/script && python calibration.py

Requires Python with torch, sbi, pandas, numpy, matplotlib, seaborn, and scipy.
"""

import os

import numpy as np
import pandas as pd
import torch
from sbi import utils as utils_sbi
from sbi.inference import SNPE

import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import gaussian_kde

# Seed for reproducibility (SNPE training and posterior sampling are stochastic).
SEED = 0
np.random.seed(SEED)
torch.manual_seed(SEED)

# ---- Paths ------------------------------------------------------------------
data_dir = "../data"
out_dir = "../output"
os.makedirs(out_dir, exist_ok=True)

# ---- Load data --------------------------------------------------------------
df_sim = pd.read_csv(os.path.join(data_dir, "proportion_summary.csv"))
params = pd.read_csv(os.path.join(data_dir, "parameters.csv"))

# Parameters being calibrated (a subset of the columns in parameters.csv).
calib_params = [
    "transmission.probability.insertive.infected",
    "prop.casual.sex.acts",
    "prop.steady.sex.acts",
    "prep.bl.use.prop",
    "transmission.probability.receptive.infected",
    "prep.fraction.sd",
]
n_params = len(calib_params)

# ---- Build the training set (X, y) ------------------------------------------
# Summary statistic per instance: HIV prevalence in the final year-category.
prevalence_last = (
    df_sim.sort_values(["instance", "ycat"])
    .groupby("instance")
    .tail(1)
    .set_index("instance")
)
y_series = prevalence_last["mean_hiv_prevalence"]

# Align the parameter rows to the same instance order as the outcomes.
instances = y_series.index
X = params.set_index("instance").loc[instances, calib_params].values
y = y_series.values

# Normalize parameters to [0, 1] (per column).
X_min, X_max = X.min(axis=0), X.max(axis=0)
X_norm = (X - X_min) / (X_max - X_min)

# Standardize the observation.
y_mean, y_std = y.mean(), y.std()
y_norm = (y - y_mean) / y_std

# ---- Train the neural posterior estimator -----------------------------------
# Prior is uniform on the normalized [0, 1] hypercube.
prior = utils_sbi.BoxUniform(low=torch.zeros(n_params), high=torch.ones(n_params))

inference = SNPE(prior=prior)
inference.append_simulations(
    torch.tensor(X_norm, dtype=torch.float32),
    torch.tensor(y_norm.reshape(-1, 1), dtype=torch.float32),
)
inference.train()
posterior = inference.build_posterior()

# ---- Condition on the observed target ---------------------------------------
# Observed HIV prevalence (%) to calibrate to.
target = 29.490043
target_norm = (target - y_mean) / y_std

samples = posterior.sample((10000,), x=torch.tensor([[target_norm]]))
samples_norm = samples.numpy()  # posterior on the normalized [0, 1] scale

# Transform posterior samples back to the native parameter scale.
samples_native = samples_norm * (X_max - X_min) + X_min
samples_native_df = pd.DataFrame(samples_native, columns=calib_params)
samples_native_df.to_csv(
    os.path.join(out_dir, "posterior_samples_mean_10000.csv"), index=False
)

print(f"Posterior mean: {samples_native.mean(axis=0)}")
print(f"Posterior std: {samples_native.std(axis=0)}")

# ---- Plot: 1D marginal histograms -------------------------------------------
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
axes = axes.flatten()
for i, param in enumerate(calib_params):
    ax = axes[i]
    ax.hist(samples_native[:, i], bins=50, alpha=0.7, edgecolor="black")
    ax.axvline(
        samples_native[:, i].mean(),
        color="red",
        linestyle="--",
        linewidth=2,
        label="Mean",
    )
    ax.set_xlabel(param)
    ax.set_ylabel("Frequency")
    ax.legend()
plt.tight_layout()
plt.savefig(os.path.join(out_dir, "posterior_marginals_mean.png"), dpi=300)
plt.close(fig)

# ---- Plot: seaborn pairplot -------------------------------------------------
df_samples = pd.DataFrame(samples_native, columns=calib_params)
g = sns.pairplot(df_samples, diag_kind="hist", plot_kws={"alpha": 0.5})
g.figure.suptitle("Posterior Samples - Pairplot", y=1.001)
plt.savefig(
    os.path.join(out_dir, "posterior_pairplot_mean.png"), dpi=300, bbox_inches="tight"
)
plt.close(g.figure)

# ---- Plot: training (prior) vs posterior per parameter ----------------------
# Compare the normalized training distribution with the normalized posterior.
fig, axes = plt.subplots(1, n_params, figsize=(18, 4))
for i, param in enumerate(calib_params):
    ax = axes[i]
    ax.hist(
        X_norm[:, i],
        bins=30,
        alpha=0.5,
        label="Prior (training sims)",
        edgecolor="black",
    )
    ax.hist(
        samples_norm[:, i], bins=30, alpha=0.5, label="Posterior", edgecolor="black"
    )
    ax.set_xlabel(param)
    ax.set_ylabel("Frequency")
    ax.legend()
plt.tight_layout()
plt.savefig(os.path.join(out_dir, "prior_vs_posterior_mean.png"), dpi=300)
plt.close(fig)

# ---- Plot: KDE pairplot (1D densities on diagonal, 2D contours off) ---------
fig, axes = plt.subplots(n_params, n_params, figsize=(16, 16))
for i in range(n_params):
    for j in range(n_params):
        ax = axes[i, j]

        if i == j:
            # Diagonal: 1D density.
            kde = gaussian_kde(samples_native[:, i])
            x_range = np.linspace(
                samples_native[:, i].min(), samples_native[:, i].max(), 200
            )
            density = kde(x_range)
            ax.fill_between(x_range, density, alpha=0.7, color="steelblue")
            ax.plot(x_range, density, "k-", linewidth=1.5)
            ax.set_ylabel("Density")
        else:
            # Off-diagonal: 2D density contour.
            xy = np.vstack([samples_native[:, j], samples_native[:, i]])
            kde = gaussian_kde(xy)
            x = np.linspace(samples_native[:, j].min(), samples_native[:, j].max(), 100)
            y = np.linspace(samples_native[:, i].min(), samples_native[:, i].max(), 100)
            X_grid, Y_grid = np.meshgrid(x, y)
            positions = np.vstack([X_grid.ravel(), Y_grid.ravel()])
            Z = kde(positions).reshape(X_grid.shape)
            ax.contourf(X_grid, Y_grid, Z, levels=12, cmap="viridis", alpha=0.8)
            ax.contour(X_grid, Y_grid, Z, levels=6, colors="black", alpha=0.3, linewidths=0.5)

        # Axis labels only on the outer edges.
        ax.set_xlabel(calib_params[j], fontsize=9) if i == n_params - 1 else ax.set_xlabel("")
        ax.set_ylabel(calib_params[i], fontsize=9) if j == 0 else ax.set_ylabel("")
        ax.tick_params(labelsize=8)
plt.tight_layout()
plt.savefig(
    os.path.join(out_dir, "posterior_pairplot_kde_mean.png"),
    dpi=300,
    bbox_inches="tight",
)
plt.close(fig)
