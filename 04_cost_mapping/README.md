# 04 — Cost Mapping

Translates ART/PrEP government **funding reductions** into intervention coverage,
queries the stage-03 surrogate for the resulting HIV incidence, and produces the
cost-mapped figures and scenario table for the paper. This is a downstream
consumer of stage 03: it does not re-fit or re-define the surrogate, it loads the
saved bundle and calls its prediction helpers.

## Method

1. **Cost model.** Unit costs (2018 → 2025 inflation-adjusted), MSM population
   shares, and 2022/2025 government funding levels are combined into baseline ART
   and PrEP spending for the MSM population, and a per-person cost-to-coverage
   mapping (`proportion_to_dollar_reduction()`).
2. **Funding grid.** A 101 × 101 grid of ART × PrEP funding-reduction fractions is
   mapped to coverage proportions, then to the surrogate's native `(art, prep)`
   intervention inputs.
3. **Surrogate query.** The stage-03 composite incidence surrogate is evaluated
   over the grid, giving mean incidence, incidence risk ratio (IRR), and relative
   change vs. the no-cut baseline, each with a 95% posterior predictive interval.
4. **Outputs.** Contour and discrete-tile heatmaps (with iso-budget lines), bar
   plots of selected policy scenarios, and a LaTeX scenario table.

The reported outcome is the year-10 tick of the surrogate trajectory, consistent
with the original cost-mapping analysis.

### Incidence risk ratio and common random numbers

The incidence risk ratio (IRR) is `scenario_incidence / baseline_incidence`.
Drawing the numerator and denominator from **independent** surrogate predictive
samples gives the ratio a noise floor of roughly ±2%: even the self-ratio at the
no-cut baseline `(0, 0)` — which must equal 1 — spreads over about `[0.98, 1.02]`.
For scenarios whose true effect is small (modest PrEP funding cuts, which map to
only a few percent coverage reduction because most PrEP financing is private),
that noise exceeds the effect, so the 90% interval spuriously crosses below 1
even though the mean IRR is monotonically above 1.

`R/irr_common_random_numbers.R` (`compute_irr_draws_crn()`) fixes this with
**common random numbers**: within each checkpoint the baseline and scenario
predictions reuse the same standard-normal draws, so shared predictive noise
cancels in the ratio. The `(0, 0)` self-ratio then collapses to exactly 1 and
the intervals reflect only genuine scenario-vs-baseline differences (the mean
IRR is unchanged). A single flag, `USE_COMMON_RANDOM_NUMBERS` in
`run_cost_mapping.R`, toggles this on/off for **both** the forest figure and the
scenario table so the two stay consistent; set it to `FALSE` to reproduce the
classic independent-draw intervals.

## Contents

Parameters, reusable functions, and the run script are separated:

### `script/`
- **`run_cost_mapping.R`** — the script you run. Loads
  `../../03_intervention_scenario_surrogate/output/surrogate.Rdata`, then calls
  the functions in `R/` in order: assemble parameters, predict the incidence
  surface, compute posterior predictive intervals, map funding reductions to
  coverage outcomes, and write all figures and the table. Run it with `script/`
  as the working directory:

  ```sh
  cd 04_cost_mapping/script
  Rscript run_cost_mapping.R
  ```

  Requires `ggplot2`, `data.table`, `hetGP` (so `predict()` dispatches on the
  stored surrogate models), `viridis`, `dplyr`, `tidyr`, `shadowtext`, and
  `kableExtra`. Run stage 03 first so `surrogate.Rdata` exists.

### `cost_params.yml`
Every raw input parameter (unit costs, populations, government funding,
coverage-model parameters, grid settings), one documented entry per value with a
provenance comment tracing it to the paper where possible. This is the single
place to change a value. Derived quantities are computed in R, not stored here.

### `R/`
Reusable functions, sourced by the script:
- **`parameters.R`** — `cost_mapping_params()`, reads `cost_params.yml` and
  computes the derived quantities (population subgroups, MSM funding shares,
  baseline funding and coverage).
- **`cost_model.R`** — funding <-> coverage mapping (`build_funding_scenarios()`,
  `proportion_to_dollar_reduction()`) and iso-budget geometry
  (`create_isobudget_lines()`).
- **`outcomes.R`** — surrogate grid prediction (`predict_incidence_grid()`) and
  posterior predictive intervals (`compute_ppi()`).
- **`irr_common_random_numbers.R`** — `compute_irr_draws_crn()` and
  `summarise_irr_draws()`, the incidence risk ratio with optional common random
  numbers (see "Incidence risk ratio and common random numbers" above). Used by
  both the forest plot and the scenario table.
- **`plot_heatmap.R`** — contour and discrete tile-matrix heatmaps, plus
  iso-budget line helpers.
- **`plot_forest.R`** — `plot_funding_forest()`, the main-text forest plot of
  incidence risk ratio by funding reduction.
- **`plot_bar.R`** — `plot_funding_bar()`, bar plots comparing selected
  funding-reduction scenarios (incidence, relative increase, risk ratio).
- **`table_scenarios.R`** — `create_funding_scenario_table()`, the LaTeX
  scenario table.

### `output/`
Generated by `run_cost_mapping.R`. The binary figures (`*.png`, `*.pdf`) are
git-ignored (they regenerate from code and churn on every render); the CSVs and
the `.tex` table are the tracked, diffable record of results. Force-add a figure
(`git add -f`) only when its results actually change.

| File | Description |
|------|-------------|
| `forest_incidence_risk_ratio_funding.csv` | mean IRR + 90% CI per scenario x funding level, with a `crosses_one` flag (tracked) |
| `cost_mapping_grid_results.csv`            | full deduplicated funding->coverage->incidence surface (tracked) |
| `forest_incidence_risk_ratio_funding.pdf`/`.png` | main-text forest plot of incidence risk ratio by funding reduction |
| `heatmap_incidence_main_year10.png`          | contour heatmap of mean incidence, with iso-budget lines |
| `heatmap_incidence_relative_year10.png`      | contour heatmap of relative incidence change |
| `heatmap_risk_ratio_main_year10.png`         | contour heatmap of incidence risk ratio |
| `heatmap_risk_ratio_relative_year10.png`     | contour heatmap of relative risk-ratio change |
| `discrete_heatmap_incidence_main_year10.png` | discrete tile-matrix heatmap of mean incidence |
| `discrete_heatmap_incidence_relative_year10.png` | discrete tile-matrix heatmap of relative incidence change |
| `discrete_heatmap_IRR_main_year10.png`       | discrete tile-matrix heatmap of incidence risk ratio |
| `discrete_heatmap_IRR_relative_year10.png`   | discrete tile-matrix heatmap of relative risk-ratio change |
| `bar_plot_incidence.png`                     | bar plot of mean incidence for selected policy scenarios |
| `bar_plot_relative_increase.png`             | bar plot of relative incidence increase |
| `bar_plot_risk_ratio.png`                    | bar plot of incidence risk ratio |
| `funding_scenarios_table.tex`                | LaTeX table of scenarios and incidence outcomes |
