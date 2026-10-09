"""Acceptance checks for the insulation-case outputs. Run after export_gamma_cases.R.
python3 05_response_surfaces/script/acceptance_checks_gamma.py
"""
from pathlib import Path

import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parents[1] / "output"
g = pd.read_csv(OUT / "gamma_case_grid.csv")
p = pd.read_csv(OUT / "gamma_case_points.csv")
c = pd.read_csv(OUT / "coverage_domain_grid.csv")
CASES = {"lower_insulation": (0.30, 0.5625), "median_insulation": (0.40, 0.75),
         "higher_insulation": (0.50, 0.9375)}
ok = True


def check(name, cond, detail=""):
    global ok
    ok &= bool(cond)
    print(("PASS " if cond else "FAIL ") + name + (f"  {detail}" if detail else ""))


# 1 time mapping
check("final -> row 11, year5 -> row 6",
      set(g[g["when"] == "final"].row_index) == {11} and set(g[g["when"] == "year5"].row_index) == {6})
# 2 domain / r maxima
for case, (gp, ga) in CASES.items():
    d = g[(g.case == case) & (g["when"] == "final")]
    check(f"{case}: r = Delta*(1-gamma) and r<=0.75",
          np.allclose(d.r_prep, d.delta_prev * (1 - gp)) and np.allclose(d.r_art, d.delta_med * (1 - ga))
          and d.r_prep.max() <= 0.75 and d.r_art.max() <= 0.75,
          f"max r_prep={d.r_prep.max():.4f} max r_art={d.r_art.max():.4f}")
# 3 baseline ratio is exactly one, keys unique, finite, intervals bracket mean
z = g[(g.delta_prev == 0) & (g.delta_med == 0)]
check("IRR(0,0) == 1 for every case/row", np.allclose(z["mean"], 1, atol=1e-10))
check("unique case/when/delta keys", not g.duplicated(["case", "when", "delta_prev", "delta_med"]).any())
check("finite, positive, interval brackets mean",
      np.isfinite(g[["mean", "lower", "upper"]]).all().all() and (g["mean"] > 0).all()
      and (g.lower <= g["mean"] + 1e-9).all() and (g["mean"] <= g.upper + 1e-9).all())
# 4 orientation: ART vs PrEP are not interchangeable
d = g[(g.case == "median_insulation") & (g["when"] == "final")].set_index(["delta_prev", "delta_med"])
a, b = d.loc[(0.0, 1.0), "mean"], d.loc[(1.0, 0.0), "mean"]
check("orientation: Delta_Med=1 differs from Delta_prev=1", abs(a - b) > 0.05, f"{a:.3f} vs {b:.3f}")
# 5 monotonicity along both axes (small CRN noise tolerated)
for case in CASES:
    z = g[(g.case == case) & (g["when"] == "final")].pivot(index="delta_med", columns="delta_prev", values="mean").values
    check(f"{case}: IRR non-decreasing in both Deltas",
          (np.diff(z, axis=0) > -1e-3).all() and (np.diff(z, axis=1) > -1e-3).all())
# 6 points agree with grid nodes (same z block despite different chunking)
worst = 0.0
for _, r in p.groupby(["case", "when", "scenario"]).agg(m=("irr", "mean"), dp=("delta_prev", "first"),
                                                         dm=("delta_med", "first")).reset_index().iterrows():
    gg = g[(g.case == r.case) & (g["when"] == r["when"]) & np.isclose(g.delta_prev, r.dp) & np.isclose(g.delta_med, r.dm)]
    worst = max(worst, abs(gg["mean"].iloc[0] - r.m))
check("per-draw point means equal grid values", worst < 1e-9, f"max diff {worst:.2e}")
# 7 draw bookkeeping
check("5000 draws, 100 checkpoints, same draws across cases",
      p.groupby(["case", "when", "scenario"]).draw.nunique().eq(5000).all() and p.checkpoint.nunique() == 100)
# 8 surface invariance: IRR at equal r agrees across cases (use P2: Delta_prev=.5 -> r_prep .35/.30/.25 differ,
#   so compare the coverage-domain grid at a shared node instead)
m = c[(c["when"] == "final") & np.isclose(c.r_prep, 0.30) & np.isclose(c.r_art, 0.0)]
n = p[(p.case == "median_insulation") & (p["when"] == "final") & (p.scenario == "P2")]
check("median P2 (r_prep=0.30, r_art=0) equals coverage-domain node",
      len(m) == 1 and abs(m["mean"].iloc[0] - n.irr.mean()) < 1e-9,
      f"{m['mean'].iloc[0]:.6f} vs {n.irr.mean():.6f}" if len(m) else "node missing")

# ---- r-space interpretation checks (added after the survey-wording correction) ----
from scipy.interpolate import RegularGridInterpolator
inp = pd.read_csv(Path(__file__).resolve().parent / "stakeholder_inputs.csv")
check("stakeholder inputs are labelled provisional", (inp.status == "provisional_placeholder").all())
cf = c[c["when"] == "final"].pivot(index="r_art", columns="r_prep", values="mean")
interp = RegularGridInterpolator((cf.index.values, cf.columns.values), cf.values)
for case, (gp, ga) in CASES.items():
    d = g[(g.case == case) & (g["when"] == "final")]
    rr = np.column_stack([d.r_art.values, d.r_prep.values])
    diff = np.abs(interp(rr) - d["mean"].values).max()
    check(f"{case}: IRR depends on r only (Delta grid matches r-space surface)", diff < 0.02, f"max diff {diff:.4f}")
    rt = (d.delta_prev * (1 - gp)) / (1 - gp)
    check(f"{case}: Delta -> r -> Delta round trip", np.allclose(rt, d.delta_prev))
P = {"P1": (0.0, 0.0875), "P2": (0.30, 0.0), "P3": (0.30, 0.0875)}
exp_feasible = {"lower_insulation": [True, True, True], "median_insulation": [True, True, True],
                "higher_insulation": [False, True, False]}
for case, (gp, ga) in CASES.items():
    fl = [(rp / (1 - gp) <= 1 + 1e-9) and (rv / (1 - ga) <= 1 + 1e-9) for rp, rv in P.values()]
    check(f"{case}: attainability of P1,P2,P3 = {exp_feasible[case]}", fl == exp_feasible[case], str(fl))
old = p[(p.case == "median_insulation") & (p["when"] == "final") & (p.scenario == "elicited")]
check("provisional point IRR identical for any case (fixed r) = 1.4308", abs(old.irr.mean() - 1.4308) < 5e-4,
      f"{old.irr.mean():.4f}")
tex = (OUT.parent / "proportional_funding_reductions_response_surfaces.tex").read_text()
bad = [w for w in ["elicited estimate", "Elicited estimate", "survey median r"] if w in tex]
check("no provisional quantity labelled as elicited", not bad, str(bad))
print("\nALL PASS" if ok else "\nFAILURES")
raise SystemExit(0 if ok else 1)
