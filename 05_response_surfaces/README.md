# Response surfaces from respondent-level effective reductions (first draft)

Implements Anna's October 2026 construction. For each respondent, q (stated reduction in PrEP use / ART adherence) is read as a **gross disruption before alternative access** and gamma is the same respondent's proportion retaining access by other means:

    r_PrEP = q_PrEP * (1 - gamma_PrEP),   r_VS = q_ART * (1 - gamma_ART)

Summaries are medians of the respondent-level r values (never products of marginal medians). P1 = median ART-only, P2 = median PrEP-only, P3 = both. The survey wording is ambiguous; the alternative reading r = q is run as a sensitivity. Under Delta = r/(1-gamma), the implied Delta corresponds to q only per respondent, and the survey did not ask about funding cuts, so q is a gross/implied disruption, not an observed funding reduction; funding-space figures are hypothetical what-ifs. ART adherence is a proxy for viral suppression. The earlier interim rescaling of gamma (r + gamma = 1) is withdrawn.

Outcome: IRR in the final year = surrogate row 11 (post-intervention year 10); row 6 = year 5 (supplement). `tick` in the surrogate helpers is a one-based row index.

## Confidentiality (read first)
Raw respondent data, files derived from them, and individualized PDFs must **never** enter Git. They live in a folder you choose **outside any repository** (`INFORM_PRIVATE_DIR`; this machine uses `HIV-UChicago/_private_survey/`, which is not under any repo). Code in this repo reads that folder at run time; every script refuses to read or write there if the path is inside a Git repository. The tracked `.tex` has no respondent numbers (it reads `values.tex`, tables and figures from the private folder through `TEXINPUTS`). `script/acceptance_checks_survey.py` verifies that only the survey-free grid CSV is tracked and that no respondent-derived value appears in tracked text.

## Two run modes
1. **Summary-statistics mode (current deliverable).** Uses the median, minimum and maximum across the five respondents from Anna's 19 Oct 2026 summary (`pre_mtg_survey_summary_gtz_20261019.csv` in the private folder), including her respondent-level r_i summaries (r_i_prep, r_i_art). The main point is the median of the respondent-level r_i (no product of medians); the lower/upper points combine marginal minima/maxima (the pairing is unknown), drawn as bars, never as a filled rectangle. Run `script/run_summary_analysis.sh`; report: `proportional_funding_reductions_response_surfaces_summary.tex`; PDF in `$INFORM_PRIVATE_DIR/out_summary/pdf/`.
2. **Respondent-level mode (when the paired data are run locally).** The construction described above, respondent by respondent. Run `script/run_survey_analysis.sh`; report: `proportional_funding_reductions_response_surfaces.tex`; outputs in `$INFORM_PRIVATE_DIR/out/`.

Both modes share the surrogate-only grid, the direct surrogate queries (`export_survey_points.R`, output folder chosen by `INFORM_OUT_SUBDIR`), the hypothetical funding-space view, and the r = q sensitivity.

## Run (respondent-level mode)
```sh
export INFORM_PRIVATE_DIR=/path/outside/git      # holds config.json + the respondent CSV (see respondent_tools/config_template.json)
05_response_surfaces/script/run_survey_analysis.sh
```
which runs, in order: `survey_prepare.py` (respondent r, points to query) -> `export_survey_points.R` (direct surrogate queries, final year and year 5) -> `generate_report_assets.py` (figures, tables, `values.tex`) -> `acceptance_checks_survey.py` -> LaTeX (PDF written to `$INFORM_PRIVATE_DIR/out/pdf/`). Set `"data_status": "synthetic"` in the config to stamp every figure and the title page with a synthetic-data banner (used for testing).

`script/export_coverage_grid.R` regenerates `output/coverage_domain_grid.csv` (surrogate only, tracked, no survey information; ~6 min). Funding-space (Delta) surfaces are bilinear interpolations of this grid (agreement with direct surrogate queries was checked to about 1e-4); all reported IRRs and intervals at named points are direct surrogate queries.

## Respondent validation sheets
`respondent_tools/` makes one anonymous PDF per respondent (their own q, gamma, constructed r, r-space plot, P1/P2/P3 outcomes, and the question whether their stated reduction already accounted for alternative access). See `respondent_tools/README.md`. Tests use synthetic data only: `python3 respondent_tools/test_respondent_tools.py`.

## Files
- `proportional_funding_reductions_response_surfaces.tex` (respondent-level) and `..._summary.tex` (summary statistics): the reports (numbers via macros from the private folder).
- `script/`: pipeline and checks. `respondent_tools/`: shared core, validation utility, synthetic-data generator, tests.
- `output/coverage_domain_grid.csv`: surrogate-only IRR grid over r in [0, 0.75]^2.
- `UPDATE_PLAN.md`: decisions and open items. `DATA_PIPELINE.md`: data flow. `email_draft.md`: draft note to the team (not sent).

**Paper layout (Anna, 20 Oct 2026).** Main text: Tables 1-3 and Figure 1 (`fig1_main.png`). Supplement: Table S1 and Figures S1-S3 (`figS1_funding_space`, `figS2_slices`, `figS3_year5`). The r = q reading is not in any table or figure (explained in the text only); scenarios are called PrEP component, ART component and both. `figD_*` (Figure 1 with the what-if rectangles, 3-panel funding space, box plots) are for discussion and are not for the paper. Medians vs means: see the footnote in the report; the summary file also has mean r_i.
