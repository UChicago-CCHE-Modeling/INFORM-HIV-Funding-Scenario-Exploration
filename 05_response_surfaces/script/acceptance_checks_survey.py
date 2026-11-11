#!/usr/bin/env python3
"""Acceptance checks for the respondent-level analysis (run after steps 1-3).
  python3 script/acceptance_checks_survey.py --private-dir /path/outside/git
Reads only the private output folder and tracked repository files. Prints no respondent values.
"""
import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

STAGE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(STAGE / "respondent_tools"))
import respondent_core as core  # noqa: E402

ap = argparse.ArgumentParser(); ap.add_argument("--private-dir")
a = ap.parse_args()
PDIR = core.assert_outside_git(a.private_dir or os.environ.get("INFORM_PRIVATE_DIR", ""), "private directory")
OUT = PDIR / "out"
ok = True


def check(name, cond, detail=""):
    global ok
    ok &= bool(cond)
    print(("PASS " if cond else "FAIL ") + name + (f"  {detail}" if detail else ""))


meta = json.loads((OUT / "cases.json").read_text())
resp = pd.read_csv(OUT / "respondent_r.csv")
pts = pd.read_csv(OUT / "points_to_query.csv")
dr = pd.read_csv(OUT / "points_draws.csv")
grid = pd.read_csv(core.GRID_CSV)

# ---- construction ------------------------------------------------------------
check("respondent r = q(1-gamma), per respondent",
      np.allclose(resp.r_prep, resp.q_prep * (1 - resp.gamma_prep)) and np.allclose(resp.r_vs, resp.q_art * (1 - resp.gamma_art)))
check("alternative reading r = q", np.allclose(resp.r_prep_alt, resp.q_prep) and np.allclose(resp.r_vs_alt, resp.q_art))
check("effective r never exceeds stated q", (resp.r_prep <= resp.q_prep + 1e-12).all() and (resp.r_vs <= resp.q_art + 1e-12).all())
mp, mv = float(np.median(resp.r_prep)), float(np.median(resp.r_vs))
m = pts[pts.set == "main"].set_index("name")[["r_prep", "r_vs"]]
check("P1/P2/P3 use the median of respondent-level r (componentwise)",
      np.allclose(m.loc["P1"], [0, mv], atol=1e-12) and np.allclose(m.loc["P2"], [mp, 0], atol=1e-12)
      and np.allclose(m.loc["P3"], [mp, mv], atol=1e-12))
prod = float(np.median(resp.q_prep) * (1 - np.median(resp.gamma_prep)))
print(f"INFO median of products vs product of medians (PrEP) differ by {abs(mp - prod):.4f} (the former is used)")
check("one anonymous point per respondent in both readings",
      (pts.set == "respondent").sum() == meta["n"] and ((pts.set == "alt") & pts.name.str.startswith("R")).sum() == meta["n"])
check("no real respondent ids in outputs", resp.columns.tolist().count("id") == 0)

# ---- surrogate queries ---------------------------------------------------------
surf = {w: core.Surface(w) for w in ("final", "year5")}
check("rows: final = 11, year5 = 6", set(dr[dr["when"] == "final"].row_index) == {11} and set(dr[dr["when"] == "year5"].row_index) == {6})
check("5000 draws over 100 checkpoints per point", dr.groupby(["set", "name", "case", "when"]).draw.nunique().eq(5000).all()
      and dr.checkpoint.nunique() == 100)
inside = (dr.r_prep <= core.R_MAX + 1e-12) & (dr.r_vs <= core.R_MAX + 1e-12)
check("NA only where r exceeds the surrogate support (never clipped)", dr[~inside].irr.isna().all() and dr[inside].irr.notna().all())
worst = 0.0
for (s_, n_, c_, w_), g in dr[inside].groupby(["set", "name", "case", "when"]):
    worst = max(worst, abs(g.irr.mean() - float(surf[w_].at(g.r_prep.iloc[0], g.r_vs.iloc[0]))))
check("surrogate draws agree with the tracked surface at every point", worst < 0.03, f"max diff {worst:.4f}")
check("IRR(0,0) = 1 in the tracked grid", np.allclose(grid[(grid.r_prep == 0) & (grid.r_art == 0)]["mean"], 1, atol=1e-10))

# ---- funding-space what-if -------------------------------------------------------
cases = {k: tuple(v) for k, v in meta["cases"].items()}
for c, (gp, ga) in cases.items():
    for k, (dp, dm) in meta["policies"].items():
        row = pts[(pts.set == "policy") & (pts.case == c) & (pts.name == k)].iloc[0]
        check(f"{c} {k}: r = Delta(1-gamma)", abs(row.r_prep - dp * (1 - gp)) < 1e-12 and abs(row.r_vs - dm * (1 - ga)) < 1e-12)
    rt = (mp / (1 - gp)) * (1 - gp)
    check(f"{c}: Delta <-> r round trip", abs(rt - mp) < 1e-12)

# ---- outputs exist ----------------------------------------------------------------
need = ["fig_rspace_final.png", "fig_rspace_year5.png", "fig_rspace_alt.png", "fig_delta_whatif.png", "fig_slices.png",
        "fig_points.png", "fig_policy_by_case.png", "tab_summary.tex", "tab_points.tex", "tab_alt.tex", "tab_policy.tex",
        "tab_implied.tex", "values.tex"]
check("all report assets present", all((OUT / f).exists() for f in need))
if meta["data_status"] == "synthetic":
    check("synthetic data are bannered", "SYNTHETIC" in (OUT / "values.tex").read_text())

# ---- repository privacy ------------------------------------------------------------
_ls = subprocess.run(["git", "ls-files"], cwd=STAGE, capture_output=True, text=True).stdout.split()
tracked = [str((STAGE / t).resolve()) for t in _ls]
allowed_csv = {str(STAGE / "output" / "coverage_domain_grid.csv")}
check("only the survey-free grid CSV is tracked", {t for t in tracked if t.endswith(".csv")} <= allowed_csv,
      str(sorted(Path(t).name for t in tracked if t.endswith(".csv") and t not in allowed_csv)))
check("no PDFs or per-respondent files are tracked", not any(t.endswith(".pdf") or "respondent_" in Path(t).name and not t.endswith(".py")
                                                          and "respondent_tools" not in t for t in tracked))
check("private directory is outside Git", core.find_git_root(PDIR) is None)
tex = (STAGE / "proportional_funding_reductions_response_surfaces.tex").read_text()
check("tracked .tex defines no respondent values", "\\newcommand" not in tex and "\\input{values.tex}" in tex)
leak = []
needles = {f"{v:.4f}" for v in list(resp.r_prep) + list(resp.r_vs) + [mp, mv]} - {"0.0000"}
for t in tracked:
    if t.endswith((".tex", ".md", ".R", ".py", ".json", ".yml")):
        txt = Path(t).read_text(errors="ignore")
        leak += [(Path(t).name, n_) for n_ in needles if n_ in txt]
check("no respondent-derived value appears in tracked text files", not leak, str(leak))
print("\nALL PASS" if ok else "\nFAILURES")
sys.exit(0 if ok else 1)
