#!/usr/bin/env python3
"""Write SYNTHETIC respondent data and a matching config to a directory outside Git.
Usage: python3 make_synthetic_data.py --out-dir /path/outside/git [--n 5] [--seed 1]
The numbers are random and mean nothing; they exist only to test the pipeline.
"""
import argparse
import json
from pathlib import Path

import numpy as np

import respondent_core as core

ap = argparse.ArgumentParser()
ap.add_argument("--out-dir", required=True)
ap.add_argument("--n", type=int, default=5)
ap.add_argument("--seed", type=int, default=1)
a = ap.parse_args()
out = core.assert_outside_git(a.out_dir, "synthetic output directory")
out.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(a.seed)
cfg = json.loads((Path(__file__).parent / "config_template.json").read_text())
cols = cfg["columns"]
rows = ["syn_id," + ",".join(cols.values())]
for i in range(a.n):
    v = rng.integers(10, 91, size=4)
    rows.append(f"S{i+1:02d}," + ",".join(str(int(x)) for x in v))
(out / cfg["input_csv"]).write_text("\n".join(rows) + "\n")
(out / "config.json").write_text(json.dumps(cfg, indent=2) + "\n")
print(f"synthetic data and config written to {out}")
