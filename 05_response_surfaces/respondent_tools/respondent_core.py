"""Shared, data-free helpers for the respondent-level analysis.

No respondent values live in this repository. Everything here reads data at run
time from a path supplied by the user, and refuses to touch sensitive files that
sit inside a Git repository. There are no network calls.

Construction (first draft, per Anna's October 2026 note), respondent by respondent:
    r_PrEP = q_PrEP * (1 - gamma_PrEP)
    r_ART  = q_ART  * (1 - gamma_ART)
where q is the stated reduction read as a GROSS disruption before alternative
access, and gamma the same respondent's proportion retaining access by other
means. Alternative reading (sensitivity): the stated reduction already includes
alternative access, so r = q. ART adherence is used as a proxy for viral
suppression. Summaries are computed from the respondent-level r values; separate
marginal medians are never multiplied.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator

STAGE_DIR = Path(__file__).resolve().parents[1]
GRID_CSV = STAGE_DIR / "output" / "coverage_domain_grid.csv"   # surrogate-only, tracked
R_MAX = 0.75                                                    # surrogate support


# ---- privacy guard ----------------------------------------------------------
def find_git_root(path):
    p = Path(path).expanduser().resolve()
    for q in [p] + list(p.parents):
        if (q / ".git").exists():
            return q
    return None


def assert_outside_git(path, what):
    root = find_git_root(path)
    if root is not None:
        raise SystemExit(f"REFUSED: {what} ({Path(path).expanduser().resolve()}) is inside the Git repository {root}. "
                         "Keep respondent data, configs, intermediate files and PDFs outside Git.")
    return Path(path).expanduser().resolve()


# ---- config and respondent data ---------------------------------------------
REQUIRED = ["q_prep", "gamma_prep", "q_art", "gamma_art"]


def load_config(config_path):
    cfg_path = assert_outside_git(config_path, "config file")
    cfg = json.loads(cfg_path.read_text())
    for k in ("input_csv", "columns"):
        if k not in cfg:
            raise SystemExit(f"config missing '{k}'")
    missing = [k for k in REQUIRED if k not in cfg["columns"]]
    if missing:
        raise SystemExit(f"config.columns missing {missing}")
    cfg.setdefault("units", "percent")
    cfg.setdefault("id_column", None)
    cfg.setdefault("sep", ",")
    cfg.setdefault("data_status", "real")
    if cfg["units"] not in ("percent", "proportion"):
        raise SystemExit("units must be 'percent' or 'proportion'")
    if cfg["data_status"] not in ("real", "synthetic"):
        raise SystemExit("data_status must be 'real' or 'synthetic'")
    csv = Path(cfg["input_csv"]).expanduser()
    csv = csv if csv.is_absolute() else cfg_path.parent / csv
    cfg["_csv_path"] = assert_outside_git(csv, "respondent CSV")
    cfg["_dir"] = cfg_path.parent
    return cfg


def load_respondents(cfg):
    df = pd.read_csv(cfg["_csv_path"], sep=cfg["sep"])
    cols = cfg["columns"]
    absent = [v for v in cols.values() if v not in df.columns]
    if absent:
        raise SystemExit(f"CSV is missing configured columns: {absent}")
    scale = 100.0 if cfg["units"] == "percent" else 1.0
    out = pd.DataFrame({k: pd.to_numeric(df[cols[k]], errors="raise") / scale for k in REQUIRED})
    if cfg["id_column"]:
        out.insert(0, "id", df[cfg["id_column"]].astype(str).values)
    else:
        out.insert(0, "id", [str(i + 1) for i in range(len(df))])
    if out[REQUIRED].isna().any().any():
        raise SystemExit("missing values in respondent data; resolve before use")
    if ((out[REQUIRED] < 0) | (out[REQUIRED] > 1)).any().any():
        raise SystemExit("values outside [0, 100]% (check 'units' in the config)")
    if (out[["gamma_prep", "gamma_art"]] >= 1).any().any():
        raise SystemExit("gamma = 100% leaves no reduction to model; check the data")
    if out["id"].duplicated().any():
        raise SystemExit("duplicate respondent ids")
    return out.reset_index(drop=True)


def effective_r(df):
    """Respondent-level effective reductions (main) and the direct alternative r = q."""
    d = df.copy()
    d["r_prep"] = d.q_prep * (1 - d.gamma_prep)
    d["r_vs"] = d.q_art * (1 - d.gamma_art)
    d["r_prep_alt"] = d.q_prep
    d["r_vs_alt"] = d.q_art
    return d


def median_point(d, alt=False):
    """Componentwise median of the respondent-level r values (never a product of marginal medians)."""
    a, b = ("r_prep_alt", "r_vs_alt") if alt else ("r_prep", "r_vs")
    return float(np.median(d[a])), float(np.median(d[b]))


def p_points(d, alt=False):
    mp, mv = median_point(d, alt)
    return {"P1": (0.0, mv), "P2": (mp, 0.0), "P3": (mp, mv)}      # (r_PrEP, r_VS)


def insulation_cases(d):
    """What-if insulation levels for the funding-space view: median gamma and relative +/-25%."""
    gp, ga = float(np.median(d.gamma_prep)), float(np.median(d.gamma_art))
    cases = {"lower_insulation": (0.75 * gp, 0.75 * ga), "median_insulation": (gp, ga),
             "higher_insulation": (min(1.25 * gp, 0.999), min(1.25 * ga, 0.999))}
    return cases


# ---- surrogate surface (tracked grid, bilinear interpolation) ----------------
class Surface:
    def __init__(self, when="final"):
        g = pd.read_csv(GRID_CSV)
        g = g[g["when"] == when]
        axes = {}
        for col in ("mean", "lower", "upper"):
            z = g.pivot(index="r_art", columns="r_prep", values=col)
            self.r_vs_axis, self.r_prep_axis = z.index.values, z.columns.values
            axes[col] = RegularGridInterpolator((z.index.values, z.columns.values), z.values,
                                                bounds_error=False, fill_value=np.nan)
        self._f = axes
        self.mean_grid = g.pivot(index="r_art", columns="r_prep", values="mean").values

    def at(self, r_prep, r_vs, col="mean"):
        r_prep, r_vs = np.broadcast_arrays(np.asarray(r_prep, float), np.asarray(r_vs, float))
        pts = np.column_stack([r_vs.ravel(), r_prep.ravel()])
        return self._f[col](pts).reshape(r_prep.shape)

    def on_delta_grid(self, gamma_prep, gamma_art, n=101, col="mean"):
        x = np.linspace(0, 1, n)
        dp, dm = np.meshgrid(x, x)           # rows: Delta_Medicaid, cols: Delta_prevention
        return x, x, self.at(dp * (1 - gamma_prep), dm * (1 - gamma_art), col)
