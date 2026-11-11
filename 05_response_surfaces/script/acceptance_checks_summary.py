#!/usr/bin/env python3
"""Acceptance checks for the summary-statistics mode. Prints no survey values.
  python3 script/acceptance_checks_summary.py --private-dir /path/outside/git
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
OUT = PDIR / "out_summary"
ok = True


def check(name, cond, detail=""):
    global ok
    ok &= bool(cond)
    print(("PASS " if cond else "FAIL ") + name + (f"  {detail}" if detail else ""))


V = json.loads((OUT / "summary_values.json").read_text()); meta = json.loads((OUT / "cases.json").read_text())
ST, RP, RV = V["stats"], V["r_prep"], V["r_vs"]
pts = pd.read_csv(OUT / "points_to_query.csv"); dr = pd.read_csv(OUT / "points_draws.csv")
raw = pd.read_csv(PDIR / "pre_mtg_survey_summary_gtz_20261019.csv", sep="\t", index_col="variable")
for k, var in {"q_prep": "prep_reduction_bmsm", "g_prep": "prep_other_means_bmsm", "q_art": "art_adh_reduction_bmsm",
               "g_art": "art_other_means_bmsm"}.items():
    check(f"{k}: only median/min/max used, in [0,1]",
          abs(ST[k]["med"] - raw.loc[var, "numeric.p50"] / 100) < 1e-12 and abs(ST[k]["min"] - raw.loc[var, "numeric.p0"] / 100) < 1e-12
          and abs(ST[k]["max"] - raw.loc[var, "numeric.p100"] / 100) < 1e-12)
check("main r = median / min / max of the respondent-level r_i supplied by Anna (no product of medians formed)",
      all(abs(R_["med"] - raw.loc[v, "numeric.p50"]) < 1e-12 and abs(R_["low"] - raw.loc[v, "numeric.p0"]) < 1e-12
          and abs(R_["high"] - raw.loc[v, "numeric.p100"]) < 1e-12 for R_, v in ((RP, "r_i_prep"), (RV, "r_i_art"))))
check("effective r <= stated q (medians and extremes)", RP["med"] <= ST["q_prep"]["med"] and RV["med"] <= ST["q_art"]["med"] and RP["high"] <= ST["q_prep"]["max"] and RV["high"] <= ST["q_art"]["max"])
m = pts[pts.set == "main"].set_index("name")[["r_prep", "r_vs"]]
check("P1/P2/P3 coordinates", np.allclose(m.loc["P1"], [0, RV["med"]]) and np.allclose(m.loc["P2"], [RP["med"], 0]) and np.allclose(m.loc["P3"], [RP["med"], RV["med"]]))
al = pts[pts.set == "alt"].set_index("name")[["r_prep", "r_vs"]]
check("alternative reading r = q is computed (query file only)", np.allclose(al.loc["P3"], [ST["q_prep"]["med"], ST["q_art"]["med"]]))
mn = pts[pts.set == "mean"].set_index("name")[["r_prep", "r_vs"]]
check("mean-r footnote points", np.allclose(mn.loc["P3"], [RP["mean"], RV["mean"]]) and abs(RP["mean"] - raw.loc["r_i_prep", "numeric.mean"]) < 1e-12)
for f_ in ("tab_alt.tex", "fig_rspace_alt.png"):
    check(f"no r=q table/figure in the output folder: {f_}", not (OUT / f_).exists())
import re
_tex = (STAGE / "proportional_funding_reductions_response_surfaces_summary.tex").read_text()
check("no P1/P2/P3 scenario labels in the report source", not re.search(r"\bP[123]\b", _tex))
check("included files all exist", all((OUT / f_).exists() for f_ in re.findall(r"includegraphics(?:\[[^\]]*\])?\{([^}]+)\}", _tex)))
surf = {w: core.Surface(w) for w in ("final", "year5")}
check("rows: final = 11, year5 = 6", set(dr[dr["when"] == "final"].row_index) == {11} and set(dr[dr["when"] == "year5"].row_index) == {6})
check("5000 draws over 100 checkpoints per point", dr.groupby(["set", "name", "case", "when"]).draw.nunique().eq(5000).all() and dr.checkpoint.nunique() == 100)
check("all points inside surrogate support (r <= 0.75), nothing clipped", (dr.r_prep <= 0.75 + 1e-12).all() and (dr.r_vs <= 0.75 + 1e-12).all() and dr.irr.notna().all())
worst = max(abs(g.irr.mean() - float(surf[w].at(g.r_prep.iloc[0], g.r_vs.iloc[0]))) for (s_, n_, c_, w), g in dr.groupby(["set", "name", "case", "when"]))
check("surrogate draws agree with the tracked surface", worst < 0.03, f"max diff {worst:.4f}")
for c, (g1, g2) in meta["cases"].items():
    for k, (dp, dm) in meta["policies"].items():
        row = pts[(pts.set == "policy") & (pts.case == c) & (pts.name == k)].iloc[0]
        check(f"{c} {k}: r = Delta(1-gamma)", abs(row.r_prep - dp * (1 - g1)) < 1e-12 and abs(row.r_vs - dm * (1 - g2)) < 1e-12)
check("implied Delta(P3) = r_med / (1 - gamma_med) (not equal to q_med in general: medians of r_i differ from products of medians)",
      all(np.isfinite(r / (1 - g)) for r, g in ((RP["med"], ST["g_prep"]["med"]), (RV["med"], ST["g_art"]["med"]))))
need = ["fig1_main.png", "figS1_funding_space.png", "figS2_slices.png", "figS3_year5.png", "figD_rect_final.png", "figD_delta_panels.png",
        "figD_box_points.png", "figD_box_policies.png", "tab_stats.tex", "tab_points.tex", "tab_policy.tex", "tab_implied.tex", "values.tex"]
check("all report assets present", all((OUT / f).exists() for f in need))
_ls = subprocess.run(["git", "ls-files"], cwd=STAGE, capture_output=True, text=True).stdout.split()
tracked = [str((STAGE / t).resolve()) for t in _ls]
check("only the survey-free grid CSV is tracked", {t for t in tracked if t.endswith(".csv")} <= {str(STAGE / "output" / "coverage_domain_grid.csv")})
check("no PDFs tracked", not any(t.endswith(".pdf") for t in tracked))
check("private directory is outside Git", core.find_git_root(PDIR) is None)
for tex in STAGE.glob("proportional_funding_reductions_response_surfaces*.tex"):
    t_ = tex.read_text()
    check(f"{tex.name}: no values defined in the tracked .tex", "\\newcommand{\\" not in t_ and "\\input{values.tex}" in t_)
needles = {f"{v:.4f}" for v in [RP["med"], RV["med"], RP["low"], RP["high"], RV["low"], RV["high"]]} - {"0.0000"}
leak = []
for t in tracked:
    if t.endswith((".tex", ".md", ".R", ".py", ".json", ".yml", ".sh")):
        txt = Path(t).read_text(errors="ignore"); leak += [(Path(t).name, n_) for n_ in needles if n_ in txt]
check("no survey-derived value appears in tracked text files", not leak, str(leak))
print("\nALL PASS" if ok else "\nFAILURES")
sys.exit(0 if ok else 1)
