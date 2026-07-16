# INFORM-HIV-Funding-Scenario-Exploration

Analysis code and supporting datasets accompanying the INFORM-HIV funding
scenario exploration paper. The material here reproduces the published analyses;
it was extracted from a larger internal project, keeping only the pieces needed
for the paper.

## Structure

The analysis is a numbered pipeline. Each stage is a self-contained top-level
folder with `data/` (inputs), `script/` (analysis code), and `output/`
(generated figures/tables).

| Stage | Description |
|-------|-------------|
| [`01_morris_screening/`](01_morris_screening/) | Morris elementary-effects screening to rank the influence of model parameters on epidemic outcomes. |
| [`02_calibration/`](02_calibration/) | Calibrates model parameters to an observed HIV prevalence target via simulation-based inference. |
| [`03_intervention_scenario_surrogate/`](03_intervention_scenario_surrogate/) | Gaussian-process surrogate that emulates HIV trajectories under ART/PrEP intervention scenarios. |

See each stage's `README.md` for details.

## Requirements

Requirements are listed per stage:

- **Stage 01** (R): `data.table`, `ggplot2`, `sensitivity`.
- **Stage 02** (Python): `torch`, `sbi`, `pandas`, `numpy`, `matplotlib`,
  `seaborn`, `scipy`.
- **Stage 03** (R): `data.table`, `hetGP`.

## Running

Run each script from within its own `script/` folder so that the relative
`../data/` and `../output/` paths resolve correctly:

```sh
# Stage 01 (R)
cd 01_morris_screening/script
Rscript morris.R

# Stage 02 (Python)
cd 02_calibration/script
python calibration.py

# Stage 03 (R)
cd 03_intervention_scenario_surrogate/script
Rscript surrogate.R
```

Each script reads only from its stage's `data/` folder and writes results to its
`output/` folder.
