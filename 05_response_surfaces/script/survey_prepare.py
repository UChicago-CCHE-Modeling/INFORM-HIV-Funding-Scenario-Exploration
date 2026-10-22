#!/usr/bin/env python3
"""Step 1: respondent file -> effective reductions and the points to query.

Reads <private-dir>/config.json (and the CSV it names). The private directory
must be OUTSIDE any Git repository; everything written goes to <private-dir>/out.
Nothing from the respondent file is written anywhere else.

  python3 script/survey_prepare.py --private-dir /path/outside/git
(or set INFORM_PRIVATE_DIR)
"""
import argparse
import json
import os
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "respondent_tools"))
import respondent_core as core  # noqa: E402

S2_DELTA = 0.536   # Medicaid loss quoted in Anna's mock-up (AFC published estimate): hypothetical policy S2
POLICIES = {"S1": (0.5, 0.0), "S2": (0.0, S2_DELTA), "S3": (0.5, S2_DELTA)}   # (Delta_prevention, Delta_Medicaid)


def private_dir(arg):
    p = arg or os.environ.get("INFORM_PRIVATE_DIR")
    if not p:
        raise SystemExit("give --private-dir or set INFORM_PRIVATE_DIR (a folder outside Git)")
    return core.assert_outside_git(p, "private directory")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--private-dir")
    a = ap.parse_args()
    pdir = private_dir(a.private_dir)
    cfg = core.load_config(pdir / "config.json")
    d = core.effective_r(core.load_respondents(cfg))
    out = pdir / "out"; out.mkdir(exist_ok=True)
    d.assign(anon=[f"R{i+1}" for i in range(len(d))]).drop(columns=["id"]).to_csv(out / "respondent_r.csv", index=False)
    rows = []
    for i, r in enumerate(d.itertuples(index=False)):
        rows.append(("respondent", f"R{i+1}", "all", r.r_prep, r.r_vs))
        rows.append(("alt", f"R{i+1}", "all", r.r_prep_alt, r.r_vs_alt))
    for nm, (a_, b_) in core.p_points(d).items():
        rows.append(("main", nm, "all", a_, b_))
    for nm, (a_, b_) in core.p_points(d, alt=True).items():
        rows.append(("alt", nm, "all", a_, b_))
    cases = core.insulation_cases(d)
    for c, (gp, ga) in cases.items():
        for nm, (dp, dm) in POLICIES.items():
            rows.append(("policy", nm, c, dp * (1 - gp), dm * (1 - ga)))
    pd.DataFrame(rows, columns=["set", "name", "case", "r_prep", "r_vs"]).to_csv(out / "points_to_query.csv", index=False)
    (out / "cases.json").write_text(json.dumps({"n": int(len(d)), "data_status": cfg["data_status"],
                                                "cases": cases, "policies": POLICIES}, indent=2))
    print(f"prepared {len(d)} respondents, {len(rows)} query points -> {out}")


if __name__ == "__main__":
    main()
