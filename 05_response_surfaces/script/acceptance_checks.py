"""
acceptance_checks.py -- executable acceptance tests for the corrected surface.

These are the demonstrations the independent review (CODEX_INDEPENDENT_REVIEW.md)
requires before the figures can be used:

  D1  correct results at ASYMMETRIC ART/PrEP points (equal-cut diagonal points
      cannot detect an axis swap);
  D2  funding predictions queried at the MAPPED coverage coordinates;
  D3  agreement with existing results under matching parameter and time-summary
      conventions (the committed forest CSV / Figure 1);
  D4  correct figures and documentation in the ACTUAL delivery folder.

Run:  cd RAF_temp/generated_response_surfaces && python3 acceptance_checks.py
Needs: /tmp/accept_r_values.csv from acceptance_checks.R (run that first).
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

HERE = Path(__file__).resolve().parent
CSV = Path("/Users/r_vardavas/Documents/Projects_2026/HIV-UChicago/"
           "INFORM-HIV-Funding-Scenario-Exploration/04_cost_mapping/output/"
           "response_surface_grid_full.csv")
FOREST_CSV = Path("/Users/r_vardavas/Documents/Projects_2026/HIV-UChicago/"
                  "INFORM-HIV-Funding-Scenario-Exploration/04_cost_mapping/output/"
                  "forest_incidence_risk_ratio_funding.csv")
RAF_TEMP = HERE.parent
R_VALUES = Path("/tmp/accept_r_values.csv")
TIME_SUMMARY = "horizon_mean"
PANEL_A = "A. Intervention use reduction"
PANEL_B = "B. Government funding reduction"

results = []


def check(name, ok, detail):
    results.append((name, bool(ok), detail))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {detail}")


d = pd.read_csv(CSV)
sub = d[d.time_summary == TIME_SUMMARY]
art_nodes = np.sort(sub.art_red.unique())
prep_nodes = np.sort(sub.prep_red.unique())


def grid(panel, col):
    z = sub[sub.panel == panel].pivot_table(index="art_red", columns="prep_red",
                                            values=col, aggfunc="first")
    return z.reindex(index=art_nodes, columns=prep_nodes).values


A, B = grid(PANEL_A, "mean"), grid(PANEL_B, "mean")


def at(z, a, p_):
    return z[int(np.argmin(abs(art_nodes - a))), int(np.argmin(abs(prep_nodes - p_)))]


print("=" * 78)
print("D1  asymmetric ART/PrEP points (axis-swap detector)")
print("=" * 78)
r = pd.read_csv(R_VALUES)
a1, a2 = at(A, 0.40, 0.10), at(A, 0.10, 0.40)
check("pole asymmetry", abs(a1 - a2) > 0.05,
      f"A(ART .40, PrEP .10)={a1:.6f} vs A(ART .10, PrEP .40)={a2:.6f}")
check("not symmetric about the diagonal", not np.allclose(A, A.T, atol=1e-6),
      "max|A - A.T| = %.4f IRR" % np.max(np.abs(A - A.T)))
worst = max(abs(at(A, r.art_red[i], r.prep_red[i]) - r.A_ref[i]) for i in range(len(r)))
check("surface == independent R query at asymmetric points", worst < 5e-3,
      f"max|diff| = {worst:.2e} IRR over {len(r)} points")
print()
print("=" * 78)
print("D2  funding panel queried at the MAPPED coverage coordinates")
print("=" * 78)
wb = max(abs(at(B, r.art_red[i], r.prep_red[i]) - r.B_ref[i]) for i in range(len(r)))
check("Panel B == independent R query at mapped coverage", wb < 5e-3,
      f"max|diff| = {wb:.2e} IRR")
ratio = at(B, 0.40, 0.40) / at(A, 0.40, 0.40)
check("insulation buffer present (B < A)", ratio < 0.99,
      f"A(0.40,0.40)={at(A,0.40,0.40):.6f} B={at(B,0.40,0.40):.6f} B/A={ratio:.6f}")
check("panels are not clones",
      not np.array_equal(np.sort(A.ravel()), np.sort(B.ravel())),
      "Panel A and Panel B are different surfaces")
row_a = d[(d.panel == PANEL_B) & (np.isclose(d.art_red, 0.40, atol=1e-12))]
k_art = row_a.art_cov.iloc[0] / 0.40
row_p = d[(d.panel == PANEL_B) & (np.isclose(d.prep_red, 0.40, atol=1e-12))]
k_prep = row_p.prep_cov.iloc[0] / 0.40
check("mapping column = 1 - gamma_Anna multiplier",
      abs(k_art - 0.7547131505) < 1e-9 and abs(k_prep - 0.2209678043) < 1e-9,
      f"k_art={k_art:.10f} k_prep={k_prep:.10f}")
# structural: B(a, p) must equal the Panel-A surface at the mapped point
def bilinear(Z, a, p_):
    i = int(np.clip(np.searchsorted(art_nodes, a) - 1, 0, len(art_nodes) - 2))
    j = int(np.clip(np.searchsorted(prep_nodes, p_) - 1, 0, len(prep_nodes) - 2))
    ta = (a - art_nodes[i]) / (art_nodes[i + 1] - art_nodes[i])
    tp = (p_ - prep_nodes[j]) / (prep_nodes[j + 1] - prep_nodes[j])
    return ((1 - ta) * (1 - tp) * Z[i, j] + ta * (1 - tp) * Z[i + 1, j]
            + (1 - ta) * tp * Z[i, j + 1] + ta * tp * Z[i + 1, j + 1])


mapped = np.array([bilinear(A, a * k_art, p_ * k_prep)
                   for a in art_nodes for p_ in prep_nodes])
check("B == A at mapped coordinates (whole plane)",
      np.max(np.abs(mapped - B.ravel())) < 0.02,
      f"max|diff| = {np.max(np.abs(mapped - B.ravel())):.4f} IRR "
      "(bilinear interpolation of the A grid)")

print()
print("=" * 78)
print("D3  agreement with committed results under matching conventions")
print("=" * 78)
f = pd.read_csv(FOREST_CSV)
fB = f[(f.panel == PANEL_B) & (f.scenario == "Reduce both PrEP and ART")
       & (f.reduction == 0.40)]["mean"].iloc[0]
fA = f[(f.panel == PANEL_A) & (f.scenario == "Reduce both PrEP and ART")
       & (f.reduction == 0.40)]["mean"].iloc[0]
check("A(0.40,0.40) matches Figure 1 (horizon-mean)", abs(at(A, 0.40, 0.40) - fA) < 1e-6,
      f"surface {at(A,0.40,0.40):.10f} vs forest CSV {fA:.10f}")
# Panel B differs from Figure 1 only by the CRN realisation: the surface uses ONE
# constant z block for both panels (so the buffer is noise-free), whereas the
# forest plot draws B's block after A's. Measured difference at the nine shared
# nodes: <= 1.9e-4 IRR.
check("B(0.40,0.40) matches Figure 1 (horizon-mean)", abs(at(B, 0.40, 0.40) - fB) < 1e-3,
      f"surface {at(B,0.40,0.40):.10f} vs forest CSV {fB:.10f} "
      f"(diff {abs(at(B,0.40,0.40) - fB):.1e} = CRN realisation)")
check("policy nodes are exact grid nodes",
      all(any(np.isclose(art_nodes, v, atol=1e-12)) for v in (0.10, 0.25, 0.40)),
      "0.10 / 0.25 / 0.40 present as nodes")
print("  (bit-level reproduction of the whole forest CSV: see acceptance_checks.R)")

print()
print("=" * 78)
print("D4  figures and documentation in the actual delivery folder")
print("=" * 78)
names = ["fig_response_surface_intervention.png", "fig_response_surface_funding.png",
         "fig_response_surfaces_comparison.png"]
for n in names:
    p = HERE / n
    if p.exists():
        w, h = Image.open(p).size
        fresh = p.stat().st_mtime >= CSV.stat().st_mtime
        ok = p.stat().st_size > 20_000 and h < 3000 and w < 6000 and fresh
        detail = (f"{w}x{h} px, {p.stat().st_size // 1024} KB, "
                  f"newer than the CSV: {fresh}")
    else:
        ok, detail = False, "MISSING"
    check(f"delivery folder holds {n}", ok, detail)
stray = [p.name for p in RAF_TEMP.glob("fig_response_surface*.png")]
check("no stray surface PNG left in RAF_temp/", not stray,
      f"found {stray}" if stray else "clean")

docs = (HERE / "README.md").read_text() + (HERE / "DATA_PIPELINE.md").read_text()
check("no false gamma provenance claim in the docs",
      "gamma_ART = 0.75" not in docs and "0.75 / 0.40" not in docs,
      "old 'gamma_ART = 0.75, gamma_PrEP = 0.40' provenance gone")
check("gamma denominator stated explicitly",
      "gamma_Anna" in docs or "γ_Anna" in docs,
      "gamma_Anna = gamma_old / P_baseline documented")
check("time-summary convention documented", "horizon" in docs.lower(),
      "horizon-mean convention documented")

print()
print("=" * 78)
fails = [n for n, ok, _ in results if not ok]
print(f"{len(results) - len(fails)}/{len(results)} checks passed")
if fails:
    print("FAILED: " + "; ".join(fails))
    sys.exit(1)
print("ALL ACCEPTANCE CHECKS PASSED")

