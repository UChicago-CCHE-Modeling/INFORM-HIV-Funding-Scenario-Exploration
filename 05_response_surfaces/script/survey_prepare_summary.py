#!/usr/bin/env python3
"""Summary-statistics mode, step 1 (used when respondent-level data are not available to this analysis).

Reads <private-dir>/pre_mtg_survey_summary_gtz_20261019.csv (tab-separated; AFC/Anna's summary across the 5 respondents with
columns numeric.mean, numeric.sd, numeric.p0 (minimum), numeric.p25, numeric.p50 (median), numeric.p75, numeric.p100 (maximum)).
Rows prep_reduction_bmsm, prep_other_means_bmsm, art_adh_reduction_bmsm, art_other_means_bmsm are in percent; rows r_i_prep and
r_i_art are the respondent-level effective reductions r_i = q_i (1 - gamma_i) computed by Anna on the individual observations and
summarized (proportions). Only the median, minimum and maximum are used. Writes
<private-dir>/out_summary/{summary_values.json, points_to_query.csv, cases.json}. The private directory must be outside Git.

  python3 script/survey_prepare_summary.py --private-dir /path/outside/git

NOTE: r_med is the median of the respondent-level r_i (componentwise, computed upstream); no product of medians is formed.
Because only marginal summaries of r_i are available, the joint pairing is unknown: the lower/upper r values are marginal
minima/maxima (drawn as bars, never as a filled joint range). q and gamma summaries are used only for the r = q sensitivity and
for the hypothetical insulation what-ifs.
"""
import argparse
import json
import os
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "respondent_tools"))
import respondent_core as core  # noqa: E402

VARS = {"q_prep": "prep_reduction_bmsm", "g_prep": "prep_other_means_bmsm",
        "q_art": "art_adh_reduction_bmsm", "g_art": "art_other_means_bmsm"}
POLICIES = {"S1": (0.5, 0.0), "S2": (0.0, 0.536), "S3": (0.5, 0.536)}   # hypothetical (Delta_prevention, Delta_Medicaid)

ap = argparse.ArgumentParser(); ap.add_argument("--private-dir")
a = ap.parse_args()
pdir = core.assert_outside_git(a.private_dir or os.environ.get("INFORM_PRIVATE_DIR", ""), "private directory")
sv = pd.read_csv(pdir / "pre_mtg_survey_summary_gtz_20261019.csv", sep="\t", index_col="variable")
get = lambda v, sc=1.0: {"min": float(sv.loc[v, "numeric.p0"]) / sc, "med": float(sv.loc[v, "numeric.p50"]) / sc,
                         "max": float(sv.loc[v, "numeric.p100"]) / sc, "mean": float(sv.loc[v, "numeric.mean"]) / sc}
st = {k: get(v, 100) for k, v in VARS.items()}
for k, v in st.items():
    assert 0 <= v["min"] <= v["med"] <= v["max"] <= 1, k
r_p, r_v = get("r_i_prep"), get("r_i_art")
r_p, r_v = ({"med": d["med"], "low": d["min"], "high": d["max"], "mean": d["mean"]} for d in (r_p, r_v))
for r_ in (r_p, r_v):
    assert 0 <= r_["low"] <= r_["med"] <= r_["high"] <= 1
rows = []
for set_, (pp, vv, lo, hi) in {"main": (r_p["med"], r_v["med"], (r_p["low"], r_v["low"]), (r_p["high"], r_v["high"])),
                               "alt": (st["q_prep"]["med"], st["q_art"]["med"],
                                       (st["q_prep"]["min"], st["q_art"]["min"]), (st["q_prep"]["max"], st["q_art"]["max"]))}.items():
    rows += [(set_, "P1", "all", 0.0, vv), (set_, "P2", "all", pp, 0.0), (set_, "P3", "all", pp, vv),
             (set_, "low", "all", *lo), (set_, "high", "all", *hi)]
rows += [("mean", "P1", "all", 0.0, r_v["mean"]), ("mean", "P2", "all", r_p["mean"], 0.0), ("mean", "P3", "all", r_p["mean"], r_v["mean"])]
gp, ga = st["g_prep"]["med"], st["g_art"]["med"]
cases = {"lower_insulation": (0.75 * gp, 0.75 * ga), "median_insulation": (gp, ga),
         "higher_insulation": (min(1.25 * gp, 0.999), min(1.25 * ga, 0.999))}
for c, (g1, g2) in cases.items():
    for nm, (dp, dm) in POLICIES.items():
        rows.append(("policy", nm, c, dp * (1 - g1), dm * (1 - g2)))
out = pdir / "out_summary"; out.mkdir(exist_ok=True)
pd.DataFrame(rows, columns=["set", "name", "case", "r_prep", "r_vs"]).to_csv(out / "points_to_query.csv", index=False)
(out / "summary_values.json").write_text(json.dumps({"stats": st, "r_prep": r_p, "r_vs": r_v}, indent=2))
(out / "cases.json").write_text(json.dumps({"mode": "summary", "n": 5, "data_status": "real", "cases": cases,
                                            "policies": POLICIES}, indent=2))
print(f"summary mode: {len(rows)} query points -> {out}")
