# Response Surface Analysis: Proportional Funding Reductions

This folder contains the complete analysis of HIV incidence response to proportional funding reductions for ART and PrEP interventions.

## Contents

### Scripts
- **`export_response_surface_grid.R`** — Generates the full 101×101 response surface grid by querying the stage-03 surrogate at every Δ pair. Outputs to `../04_cost_mapping/output/response_surface_grid_full.csv`
- **`generate_response_surfaces_fullgrid.py`** — Reads the grid CSV and creates three publication-quality figures with contours and slices
- **`acceptance_checks.R`** — Validates Panel A vs Panel B asymmetry, surrogate query correctness, and agreement with forest CSV benchmarks
- **`acceptance_checks.py`** — Python-side validation of grid structure, orientation, and file output

### Data
- **`output/response_surface_grid_full.csv`** — 20,402 rows (10,201 per panel × 2 panels); contains art_red, prep_red, art_cov, prep_cov, mean IRR, and 95% credible interval bounds
- **`output/fig_response_surface_intervention.png`** — Panel A surface + slices (direct Δ mapping)
- **`output/fig_response_surface_funding.png`** — Panel B surface + slices (cost-model mapped)
- **`output/fig_response_surfaces_comparison.png`** — Side-by-side panels with uncertainty bands

### Documentation
- **`proportional_funding_reductions_response_surfaces.tex`** — Main publication document with full methodology, math, embedded figures, and feedback requests

## Quick Start

### To regenerate the grid (requires R):
```bash
cd ../../04_cost_mapping/script
/usr/local/bin/Rscript export_response_surface_grid.R
```

### To regenerate the figures (requires Python + pandas, numpy, matplotlib):
```bash
python3 generate_response_surfaces_fullgrid.py
# Outputs saved to output/
```

### To compile the LaTeX document:
```bash
pdflatex proportional_funding_reductions_response_surfaces.tex
# Or use your preferred LaTeX editor
```

### To run validation checks:
```bash
# R checks
/usr/local/bin/Rscript acceptance_checks.R

# Python checks
python3 acceptance_checks.py
```

## Key Parameters

**Base case (from `cost_params.yml`):**
- γ_ART = 0.15 (private insurance substitution floor)
- γ_PrEP = 0.28
- P_ART_baseline = 0.6115 (61.15% baseline coverage)
- P_PrEP_baseline = 0.3594 (35.94% baseline coverage)
- Effective factors: 0.7547 × Δ_ART, 0.2210 × Δ_PrEP

**Grid specification:**
- Δ range: 0.0 to 0.75 (0 to 75% reduction)
- Grid points: 101 × 101 = 10,201 per panel
- Step size: 0.0075 (uniformly spaced)

## Panel Descriptions

**Panel A: Direct Intervention Use Reduction**
- Assumes funding cuts directly reduce coverage: r = Δ
- Shows maximum potential health impact
- IRR range: 1.0–3.56

**Panel B: Government Funding Reduction (Cost-Model Mapped)**
- Maps funding cuts to coverage via insulation parameters
- Shows realistic impact accounting for alternative funding
- IRR range: 1.0–3.56 (lower slopes in moderate range due to buffer)

## Interpretation

- Strong asymmetry: ART cuts have 2–3× larger impact than equivalent PrEP cuts
- Non-linear response: effects roughly additive in log-risk space
- Panel B impact is 19–21% lower than Panel A for same Δ values (buffering effect)
- Response surfaces enable continuous policy exploration, not discrete scenarios

## Next Steps (Awaiting Feedback)

1. Gamma sensitivity ranges (γ_ART, γ_PrEP)
2. Proportional reduction range for publication (0–0.4 vs 0–0.75)
3. Named policy scenarios to highlight
4. Time horizon (horizon-mean vs year-10)

---

**Generated:** 2026-09-23  
**Status:** Ready for stakeholder feedback
