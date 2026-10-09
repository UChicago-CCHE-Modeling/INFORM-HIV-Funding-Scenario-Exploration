# Insulation-case response surfaces, r-space interpretation (October 2026)

Survey items 1-2 elicit r (reduction in use), items 3-4 elicit gamma; Delta = r/(1-gamma) is derived. The IRR depends on r only, so stakeholder points are fixed in r-space and only re-expressed in Delta. Final year = surrogate row 11 (post-intervention year 10). Stakeholder values are PROVISIONAL placeholders in `script/stakeholder_inputs.csv` (replace with respondent-level survey values, check r <= 1 - gamma per respondent). Plan: `UPDATE_PLAN.md`.

```sh
/usr/local/bin/Rscript 05_response_surfaces/script/export_gamma_cases.R   # ~6 min; r-space grid + Delta-grid cases
python3 05_response_surfaces/script/acceptance_checks_gamma.py
python3 05_response_surfaces/script/generate_gamma_case_figures.py
cd 05_response_surfaces && PATH=$PATH:/Library/TeX/texbin latexmk -pdf -outdir=build proportional_funding_reductions_response_surfaces.tex
```

Outputs: `fig_r_surface_{final,year5}.png`, `fig_delta_mapping_final.png`, `fig_funding_slices_final.png`, `fig_points_final.png`, `tables_*.tex`, plus the CSVs. Note `tick` in the surrogate helpers is a one-based row index; older "tick10" outputs are row 10, not the final year. The sections below describe the earlier legacy-parameter version, kept for regression only.

---


# Proportional funding reductions: response surface update

This folder contains the legacy-parameter response surfaces and a meeting report. It does not yet include stakeholder-gamma sensitivity cases or named policy bundles.

## Reproduce from the repository root

```sh
Rscript 05_response_surfaces/script/export_response_surface_grid.R
python3 05_response_surfaces/script/generate_response_surfaces_fullgrid.py
Rscript 05_response_surfaces/script/acceptance_checks.R
python3 05_response_surfaces/script/acceptance_checks.py
cd 05_response_surfaces
mkdir -p build output/pdf
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build proportional_funding_reductions_response_surfaces.tex
cp build/proportional_funding_reductions_response_surfaces.pdf output/pdf/
```

On this Mac, LaTeX is installed at `/Library/TeX/texbin`; add it to PATH if needed. Python needs numpy, pandas, matplotlib and Pillow. R uses the existing stage-04 dependencies including hetGP, data.table and yaml; acceptance checks also load the forest-data helper and its dependencies.

## Files and conventions

- `script/export_response_surface_grid.R`: queries the surrogate, resolving repository paths from the script location; owns this stage's CSV.
- `script/generate_response_surfaces_fullgrid.py`: reads this stage's CSV and writes all three PNGs into `output/`.
- `script/acceptance_checks.R` and `.py`: reference queries, forest regression, mapping/orientation and output checks. Temporary reference values go to ignored `build/`.
- `output/response_surface_grid_full.csv`: 43,264 rows; 104 x 104 nodes, two panels and two time summaries for one `committed` case. The 101 regular nodes per axis have exact 0.10/0.25/0.40 nodes added.
- Three PNGs: direct-use surface/slices, funding surface/slices, and comparison. Current figures show horizon-mean IRR. Tick-10 predictions are also in the CSV.
- `proportional_funding_reductions_response_surfaces.tex`: meeting update, definitions, measured results, limitations and feedback requests. PDF and LaTeX temporary files are ignored.
- `email_draft.md`: draft only; not sent.

Direct-use Delta is coverage loss. Funding Delta maps to coverage loss using `Delta * (1 - gamma_Anna)`, where `gamma_Anna = gamma_old / P_baseline`. Legacy insulated shares are 0.245287 ART and 0.779032 PrEP, giving multipliers 0.754713 and 0.220968. These are not the new stakeholder ranges.

Horizon-mean means the mean across ticks of each draw's incidence ratio, then the mean across draws. It is not a cumulative-incidence ratio. Current maxima are 2.897383 (use) and 2.196442 (funding). At 40% both the means are 2.042997 and 1.645269. Bands are 95% surrogate predictive intervals under the existing common-random-number convention; the coupling reduces Monte Carlo variation, not all uncertainty.

Adding cases: supply a name and coverage-loss multipliers (`k_art`, `k_prep`) in GAMMA_CASES; an omitted multiplier defaults to its legacy value. Select a case explicitly in the plotting script when exporting multiple cases. The acceptance comparisons target `committed`. TIME_SUMMARY can be `horizon_mean` or `tick10`; titles follow that selection. Save separate outputs if retaining multiple presentations.

Timing note: the saved surrogate contains 11 trajectory rows; horizon_mean averages all 11. The existing table convention selects row 10. Confirm its calendar-year mapping before labeling tick10 as year 10.
