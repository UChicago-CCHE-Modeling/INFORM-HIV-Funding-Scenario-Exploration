# Data pipeline

Private folder (outside Git): `config.json` + respondent CSV (read at run time only).

1. `script/survey_prepare.py`: validates the CSV (units, ranges, gamma < 100%), computes per-respondent r = q(1-gamma) and the alternative r = q, the componentwise medians of the respondent-level values, P1/P2/P3 (both readings), and hypothetical what-if insulation levels (median respondent gamma and relative +/-25%) with named funding policies S1-S3 mapped to r. Writes `out/respondent_r.csv` (anonymous labels R1..Rn), `out/points_to_query.csv`, `out/cases.json` in the private folder.
2. `script/export_survey_points.R`: queries the stage-03 surrogate directly at every point (rows 11 and 6; constant z block; points with a coordinate above 0.75 get NA, never clipped). Writes `out/points_draws.csv` (per-draw IRR with checkpoint ids).
3. `script/generate_report_assets.py`: figures, tables and `values.tex` macros, written to `out/`. Surfaces come from the tracked, survey-free `output/coverage_domain_grid.csv`.
4. `script/acceptance_checks_survey.py`: construction, median-of-products, support, surrogate agreement, privacy checks.
5. LaTeX with `TEXINPUTS` pointing at `out/`; PDF stays in the private folder.

Tracked inputs with no respondent information: `output/coverage_domain_grid.csv` (from `script/export_coverage_grid.R`).

## Summary-statistics mode (no respondent-level data)
`script/survey_prepare_summary.py` reads `pre_mtg_survey_summary_gtz_20261019.csv` (Anna's summary; only median, minimum, maximum used) from the private folder and writes `out_summary/{summary_values.json, points_to_query.csv, cases.json}`: the median of the respondent-level r_i (rows r_i_prep, r_i_art), marginal min/max bars, the r = q sensitivity, and the hypothetical named policies. `export_survey_points.R` (with `INFORM_OUT_SUBDIR=out_summary`) queries the surrogate; `generate_summary_assets.py` builds figures, tables and `values.tex`; `acceptance_checks_summary.py` verifies them; the `_summary.tex` report is compiled with `TEXINPUTS` pointing at `out_summary/`.
