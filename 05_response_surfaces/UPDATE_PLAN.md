# Response-surface update plan (stakeholder gamma cases, tick-10, Anna's fig_main layout)

Status: PLAN ONLY, nothing below is implemented yet. Supersedes the "legacy parameter" figures in the current PDF.
Sources: meeting transcripts in `RAF_temp/` ("Modeling ART and PrEP Coverage Scenarios", "Model Parameter and Plot Updates"), Anna's mock-ups `RAF_temp/fig_main.png` and `RAF_temp/fig_coverage.png` (both carry "PLACEHOLDER SURFACE": layout reference only, not numbers), and Raff's decisions in this session.

## 1. Decisions (fixed)

| Item | Decision | Origin |
|---|---|---|
| Time summary | IRR at the **final year = surrogate row 11** (post-intervention year 10; sim year 32). No horizon-mean. The old `tick10` = row 10 (ninth intervention year) is kept only as a legacy regression reference. Export both reporting year and row index. Optional 5-year comparison = row 6. | Transcript 1 (Anna/Pedro); row mapping verified in code (`m[tick, ]`, 11 rows, `ycat > 21`) and by Codex review |
| Baseline coverage | ART (viral suppression) 0.61; PrEP 0.23 as a **provisional display baseline requested in the meeting** (replaces 0.36). Not yet reconciled with the ensemble: saved training data give mean `prep.bl.use.prop` 0.196 over the 100 checkpoints and about 18.7% realized pre-intervention use; one workflow setting is 0.2257. | Transcript 2; fig_coverage footer; Codex data check |
| Mapping | `r_j = Delta_j * (1 - gamma_j)`, gamma_j = insulated share of current users (survey), used directly. The legacy code `gamma` (absolute floor, `gamma_old/P`) is NOT used. | fig_main footnote |
| Gamma cases | Median gamma_PrEP=0.40, gamma_ART=0.75 (survey medians). Raff's executive choice: relative +/-25% of each (a sensitivity exercise, not survey bounds, IQR or CI). | Raff, this session |
| | Label as **lower insulation** 0.30 / 0.5625, **survey-median insulation** 0.40 / 0.75, **higher insulation** 0.50 / 0.9375 (avoid "worst/best"). Low gamma = larger coverage loss. Low-low and high-high pairings only; mixed combinations are not explored. | |
| Display | Three heatmaps (one per case), no separate Panel A / Panel B. | Raff |
| Medicaid = ART | Delta_Medicaid maps one-to-one to Delta_ART (Medicaid treated as all treatment). | Transcript 1 |

## 2. What fig_main actually shows (reference for the new heatmaps)

- Panel A heatmap domain is **Delta space, both axes 0 to 1**: x = Delta_prevention (bottom), y = Delta_Medicaid (left).
- Top axis = r_PrEP = Delta_prev*(1-gamma_PrEP) (0 to 0.6 at median); right axis = r_VS = Delta_Med*(1-gamma_ART) (0 to 0.25 at median). So the top/right tick values change with the gamma case, and the colour surface must be evaluated at those r values.
- Overlays: P1 (Delta_prev=0, Delta_Med=0.35) Medicaid change only; P2 (0.5, 0) prevention funding halved; P3 (1.0, 0.35) prevention eliminated plus Medicaid change; "elicited estimate" dot at (0.5, 0.35) = r of about (0.30, 0.0875) at median gamma; shaded rectangle = stakeholder range (Delta_prev 0.25-0.80, Delta_Med 0.12-0.80); contour lines labelled with IRR.
- Panel B: IRR vs Delta_prev, lines at fixed Delta_Med = 0, 0.35 (elicited), 1.0, 10% increase reference line, elicited-range band.
- fig_coverage (separate figure) is the gamma-independent coverage-domain view: bottom/left = proportional r_PrEP, r_VS (0 to 0.75); top/right = absolute percentage points lost (PrEP 0-17 = 0.75*23; VS 0-46 = 0.75*61); status-quo marker; Panel B = IRR vs absolute points lost with additivity line.

## 3. Points to confirm (possible mismatches with how this was discussed)

1. **The IRR surface is the same in r space for all three cases.** Gamma only changes (a) the r range reached on a Delta in [0,1] domain (lower: r_PrEP<=0.70, r_VS<=0.4375; median: 0.60, 0.25; upper: 0.50, 0.0625) and (b) where P-points and the stakeholder box sit in r. The three heatmaps are therefore different crops/stretches of one surface. Use one shared colour scale across the three. Do not assume the higher-insulation case is flat: a direct query at its Delta=(1,1) corner gives final-year IRR about 1.49 (ART effect compressed, but PrEP can still lose 50%). In Delta space the P-points and stakeholder box stay fixed; only their r coordinates and outcomes change.
2. **Absolute-percentage-point axes live in fig_coverage, not fig_main.** The four-sided axes Raff described (Delta bottom/left, r top/right) are fig_main. Absolute points need baselines 0.61/0.23 and are gamma-independent. Plan: Figure 1 = fig_main style x3; Figure 2 = fig_coverage style x1. Confirm that is wanted.
3. **P3 is defined** in fig_main's caption (prevention eliminated plus Medicaid change); "P1 about 0.3 and just under 0.1" is the elicited estimate in r coordinates (0.30, 0.0875), not P1. P1 is Delta_prev=0, Delta_Med=0.35.
4. **Gamma range vs survey.** Survey table: gamma_PrEP 0.14-0.40, gamma_ART 0.15-0.75 (respondent minimum to median, n=5). +/-25% gives 0.30-0.50 and 0.5625-0.9375: the lower case does not reach the respondent minimum, and the upper case is above every elicited median. Keep Raff's choice but state this in the document. The survey-minimum case (0.14/0.15) is only allowed as a supplement with unsupported regions masked: it reaches r = 0.86/0.85 at Delta=1, beyond the surrogate's r<=0.75 support (supported limits Delta_PrEP<=0.872, Delta_ART<=0.882). Never silently clip r to 0.75.
5. **Baseline PrEP effect on IRR (resolved).** The surrogate takes only the two proportional adjustments (ART, PrEP); baseline PrEP is embedded in the checkpoint fits. So 0.23 vs 0.36 changes only the percentage-point axes (and the r implied by a fixed absolute loss), not IRR at fixed r. It does not recalibrate the population. Reconcile the source of 0.23 before calling it calibrated.
6. **Time mapping (resolved).** Row 1 = sim year 22 (pre-intervention), row 6 = year 27, row 10 = year 31, row 11 = year 32 = tenth intervention year. Final year = row 11.
7. **Heatmap method.** The previous work avoided interpolation (audits in `RAF_temp/audit*.py`). Evaluating the surface on a Delta grid needs r = Delta*(1-gamma) off the old r nodes, so re-query the surrogate on the Delta grid per case instead of interpolating the old CSV.

## 4. Tasks

### A. Parameters and data
- [ ] A1. Put the three gamma cases and baselines 0.61/0.23 in one config (new YAML or GAMMA_CASES block); remove dependence on legacy `gamma_old/P`.
- [ ] A2. Reconcile the source of the 0.23 PrEP baseline with Anna (point 3.5); label it provisional until then.
- [ ] A3. Extend `script/export_response_surface_grid.R`: per case, Delta grid 0-1 (101 x 101, plus exact nodes 0.12, 0.25, 0.35, 0.5, 0.8, 1.0) -> native x = (-Delta_Med*(1-gamma_ART), -Delta_prev*(1-gamma_PrEP)), order (ART, PrEP) -> surrogate at **row 11** (also row 6 optional), CRN with identical checkpoint/draw IDs across cases, mean and 95% interval. Columns: case, reporting_year, row_index. Reject any query with r > 0.75. Keep a legacy row-10 export for regression only. Keep the r-space 0-0.75 grid at row 11 for Figure 2.
- [ ] A4. Export per-draw IRR (with checkpoint and draw IDs) at P1, P2, P3, elicited per case for the box plot.

### B. Figures
- [ ] B1. Figure 1 (x3 cases): fig_main Panel A with four-sided axes (Delta bottom/left, r top/right), shared colour scale and contour levels, P1/P2/P3, elicited dot, stakeholder box. Panel B slices at Delta_Med = 0, 0.35, 1.0 with bands. Layout: three columns or rows on one page.
- [ ] B2. Figure 2: coverage-domain heatmap (r bottom/left, absolute points top/right, display baselines 0.61/0.23, r<=0.75 so max loss 45.75 pts ART, 17.25 pts PrEP). Gamma-independent. Right panel only if precisely defined: common absolute loss a gives r_ART=a/(100*P_ART), r_PrEP=a/(100*P_PrEP), so the equal-points curve stops at 17.25 pts; additive reference R_add(a)=R_ART(a)+R_PrEP(a)-1, interaction from matched draws. Otherwise replace with the table (interaction plots were optional in the meeting).
- [ ] B3. Box plot (supplement): matched row-11 per-draw IRR at P1, P2, P3, elicited, three insulation cases per group; median/IQR boxes, whiskers 2.5-97.5%. Describes model uncertainty conditional on gamma, not the five survey responses.
- [ ] B4. Bar chart: percentage increase 100*(IRR-1), mean with 95% interval, same draws, zero baseline.
- [ ] B5. Scenario table (main document): IRR [95%] at P1, P2, P3, elicited x 3 cases, plus matched ART-only, PrEP-only and both-reduction rows to allow an interaction read.

### C. Checks (extend `script/acceptance_checks.*`)
- [ ] C1. Axis contract: top tick = Delta*(1-gamma) per case, corner values 0.7/0.4375, 0.6/0.25, 0.5/0.0625.
- [ ] C2. Independent surrogate query at P1/P2/P3/elicited matches plotted values; baseline ratio = 1 at Delta=0.
- [ ] C3. Same IRR at equal r across cases (surface invariance); shared colour limits.
- [ ] C4. Orientation test with asymmetric points (ART vs PrEP swap).
- [ ] C5. Time mapping: final year -> row 11, year 5 -> row 6; reject r outside support.
- [ ] C6. Changing display baselines leaves fixed-r predictions unchanged.
- [ ] C7. Box plots, bars, tables and surface agree at identical scenarios; chunk size/query order do not change paired results beyond tolerance.
- [ ] C8. No pooling of gamma cases; unique case/time/scenario keys, finite draws, positive denominators, labels correct, figures visually inspected.

### D. Document and delivery
- [ ] D1. Update the .tex: new definitions (gamma = insulated share), cases table, tick-10 definition, new figures, interpretation, caveats (points 3.1, 3.4, 3.5); drop the "legacy cost-model" sections and stale feedback questions.
- [ ] D2. README update (commands, cases, tick10). Compile PDF locally; do not commit PDF or build files.
- [ ] D3. Commit scripts, CSV decision (size check), .tex, PNG decision per Raff, to `sensitivity/gamma-prep-r-p0`. No push without approval.
- [ ] D4. Email to Anna, Jonathan, Pedro: what changed, ask to confirm 0.23 value and the gamma +/-25% choice.

## 5. Open questions for Raff / Anna
1. Source of the 0.23 PrEP baseline: ensemble data show about 0.196 input propensity / 18.7% realized use; a workflow has 0.2257. Which does Anna mean?
2. Box plot design as in B3 (supplement), with the table in the main text: agreed?
3. Add the optional 5-year (row 6) comparison Anna floated?
4. Tell Anna and team that the earlier "tick 10" figures were row 10 (ninth intervention year), not the final year.
