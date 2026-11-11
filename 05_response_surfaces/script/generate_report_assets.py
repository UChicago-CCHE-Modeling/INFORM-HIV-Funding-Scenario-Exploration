#!/usr/bin/env python3
"""Step 3: figures, tables and LaTeX macros for the report.

Reads <private-dir>/out/{respondent_r.csv, cases.json, points_draws.csv} and the
tracked, survey-free grid output/coverage_domain_grid.csv. Writes ONLY into
<private-dir>/out (outside Git). The tracked .tex contains no survey numbers; it
reads values.tex, tables and figures from there (TEXINPUTS).

  python3 script/generate_report_assets.py --private-dir /path/outside/git
"""
import argparse
import json
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "respondent_tools"))
import respondent_core as core  # noqa: E402

ap = argparse.ArgumentParser(); ap.add_argument("--private-dir")
args = ap.parse_args()
PDIR = core.assert_outside_git(args.private_dir or os.environ.get("INFORM_PRIVATE_DIR", ""), "private directory")
OUT = PDIR / "out"
meta = json.loads((OUT / "cases.json").read_text())
STATUS, N = meta["data_status"], meta["n"]
CASES = {k: tuple(v) for k, v in meta["cases"].items()}
POL = {k: tuple(v) for k, v in meta["policies"].items()}
LABEL = {"lower_insulation": "Lower insulation", "median_insulation": "Median insulation",
         "higher_insulation": "Higher insulation"}
CASE_COLOR = {"lower_insulation": "#0072B2", "median_insulation": "#E69F00", "higher_insulation": "#7B2D8E"}
TEAL, PURPLE, GREY = "#1d6a85", "#6a3d8f", "#777777"
plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10, "savefig.dpi": 200})

resp = pd.read_csv(OUT / "respondent_r.csv")
draws = pd.read_csv(OUT / "points_draws.csv")
SURF = {"final": core.Surface("final"), "year5": core.Surface("year5")}


def watermark(fig):
    if STATUS == "synthetic":
        fig.text(0.5, 0.5, "SYNTHETIC TEST DATA", fontsize=44, color="red", alpha=0.18, ha="center", va="center",
                 rotation=25, zorder=100)


def summ(set_, when, case="all"):
    d = draws[(draws.set == set_) & (draws["when"] == when) & (draws.case == case)]
    g = d.groupby("name").irr
    s = pd.DataFrame({"mean": g.mean(), "lo": g.quantile(0.025), "hi": g.quantile(0.975)})
    s["pct"], s["pct_lo"], s["pct_hi"] = 100 * (s["mean"] - 1), 100 * (s.lo - 1), 100 * (s.hi - 1)
    return s


med_p, med_v = core.median_point(resp.assign(r_prep=resp.r_prep, r_vs=resp.r_vs))
med_pa, med_va = core.median_point(resp, alt=True)


def contour_levels(zmax):
    return np.round(np.arange(1.1, zmax, 0.2), 2)


# ---- Figure 1: r-space surface with all respondents ---------------------------
def rspace_fig(when, fname, alt=False):
    s = SURF[when]
    zmax = float(np.nanmax(s.mean_grid))
    ca, cb = ("r_prep_alt", "r_vs_alt") if alt else ("r_prep", "r_vs")
    pts = resp[[ca, cb]].values
    mp, mv = (med_pa, med_va) if alt else (med_p, med_v)
    lim = max(core.R_MAX, float(pts.max()) + 0.03)
    fig, ax = plt.subplots(figsize=(8.6, 7.3), layout="constrained")
    m = ax.pcolormesh(s.r_prep_axis, s.r_vs_axis, s.mean_grid, cmap="YlOrRd", vmin=1.0, shading="gouraud")
    cs = ax.contour(s.r_prep_axis, s.r_vs_axis, s.mean_grid, levels=contour_levels(zmax), colors="#444444", linewidths=0.6)
    ax.clabel(cs, fmt="%.1f", fontsize=7)
    if lim > core.R_MAX:
        ax.axvline(core.R_MAX, color="#555", ls=":", lw=1); ax.axhline(core.R_MAX, color="#555", ls=":", lw=1)
        ax.text(core.R_MAX + 0.005, 0.01, "beyond surrogate range", rotation=90, fontsize=7, color="#555")
    ax.scatter(pts[:, 0], pts[:, 1], s=70, facecolor="white", edgecolor=TEAL, linewidth=2, alpha=0.9, zorder=6,
               label=f"respondents (n = {N}, anonymous pairs)")
    ax.plot(mp, mv, "o", ms=13, mfc=TEAL, mec="white", mew=2, zorder=8, label="componentwise median (P3)")
    ax.plot(0, mv, "^", ms=11, mfc="white", mec="black", mew=2, clip_on=False, zorder=8)
    ax.plot(mp, 0, "^", ms=11, mfc="white", mec="black", mew=2, clip_on=False, zorder=8)
    box = dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.8)
    ax.text(0.012, mv + 0.012, "P1  ART only", fontsize=9, fontweight="bold", zorder=9, bbox=box)
    ax.text(mp + 0.012, 0.012, "P2  PrEP only", fontsize=9, fontweight="bold", zorder=9, bbox=box)
    ax.text(mp + 0.012, mv + 0.012, "P3  both", fontsize=9, fontweight="bold", zorder=9, bbox=box)
    ax.set_xlim(0, lim); ax.set_ylim(0, lim)
    ax.set_xlabel("Proportional reduction in PrEP use  ($r_{PrEP}$)")
    ax.set_ylabel("Proportional reduction in viral suppression  ($r_{VS}$; ART adherence as proxy)")
    ttl = "Direct interpretation $r=q$ (sensitivity)" if alt else "Effective reductions $r=q(1-\\gamma)$, respondent by respondent"
    ax.set_title(f"{ttl}\n{'final year' if when == 'final' else 'year 5'}", loc="left", fontweight="bold")
    ax.legend(loc="upper left", fontsize=8, framealpha=0.92)
    cb = fig.colorbar(m, ax=ax, shrink=0.85, pad=0.02); cb.set_label("Incidence rate ratio vs. status quo")
    watermark(fig); fig.savefig(OUT / fname); plt.close(fig)


# ---- Figure: hypothetical funding-loss space ----------------------------------
def delta_fig(fname):
    zmax = max(float(np.nanmax(SURF["final"].on_delta_grid(*g)[2])) for g in CASES.values())
    fig, axes = plt.subplots(1, 3, figsize=(17, 6.5), sharey=True, layout="constrained")
    for ax, (c, (gp, ga)) in zip(axes, CASES.items()):
        x, y, z = SURF["final"].on_delta_grid(gp, ga)
        ax.set_facecolor("#dcdcdc")
        dpmax, dmmax = core.R_MAX / (1 - gp), core.R_MAX / (1 - ga)
        if dpmax < 1:
            ax.text((dpmax + 1) / 2, 0.5, "beyond surrogate range (r > 0.75)", rotation=90, ha="center", va="center",
                    fontsize=8, color="#444", zorder=3)
        if dmmax < 1:
            ax.text(0.5, (dmmax + 1) / 2, "beyond surrogate range (r > 0.75)", ha="center", va="center",
                    fontsize=8, color="#444", zorder=3)
        mm = ax.pcolormesh(x, y, z, cmap="YlOrRd", vmin=1.0, vmax=zmax, shading="gouraud", zorder=1)
        cs = ax.contour(x, y, z, levels=contour_levels(zmax), colors="#444444", linewidths=0.6)
        ax.clabel(cs, fmt="%.1f", fontsize=7)
        for k, (dp, dm) in POL.items():
            ax.plot(dp, dm, "s", ms=9, mfc=PURPLE, mec="white", mew=1.5, clip_on=False, zorder=8)
            ax.text(dp + (0.02 if dp < 0.9 else -0.02), dm + 0.025, k, fontsize=11, fontweight="bold", color=PURPLE,
                    ha="left" if dp < 0.9 else "right", zorder=9,
                    bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.75))
        dp3, dm3 = med_p / (1 - gp), med_v / (1 - ga)
        if max(dp3, dm3) <= 1 + 1e-9:
            ax.plot(dp3, dm3, "o", ms=11, mfc="none", mec=TEAL, mew=2.4, zorder=8, clip_on=False)
            ax.annotate("P3 (median r)\nre-expressed", (dp3, dm3), xytext=(max(dp3 - 0.5, 0.04), max(dm3 - 0.25, 0.04)),
                        fontsize=8.5, color=TEAL, arrowprops=dict(arrowstyle="-", color=TEAL), zorder=9,
                        bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.85))
        else:
            ax.plot(min(dp3, 1), min(dm3, 1), "^", ms=11, mfc="#d62728", mec="black", clip_on=False, zorder=9)
            ax.text(0.97, 0.12, f"Red triangle: P3 (median r) would need\n$\\Delta_{{prev}}$ = {dp3:.2f}, $\\Delta_{{Med}}$ = {dm3:.2f}\n"
                    "(above 1: not attainable here)", ha="right", va="bottom", fontsize=8.5, zorder=9,
                    bbox=dict(fc="white", ec="#d62728", alpha=0.92))
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        ax.set_xlabel("Hypothetical prevention funding loss  ($\\Delta_{prevention}$)")
        ax.set_title(f"{LABEL[c]} (what-if)\n$\\gamma_{{PrEP}}$ = {gp:.2f}, $\\gamma_{{ART}}$ = {ga:.2f}", loc="left", fontweight="bold")
        if c != list(CASES)[0]:
            ax.tick_params(labelleft=False)
    axes[0].set_ylabel("Hypothetical Medicaid coverage loss  ($\\Delta_{Medicaid}$)")
    cb = fig.colorbar(mm, ax=axes, location="bottom", shrink=0.4, pad=0.03)
    cb.set_label("Incidence rate ratio vs. status quo (shared scale)")
    watermark(fig); fig.savefig(OUT / fname); plt.close(fig)


def slices_fig(fname):
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.2), sharey=True, layout="constrained")
    st = [(0.0, "$\\Delta_{Medicaid}$ = 0", "#0072B2", "-", 2.4), (0.536, "$\\Delta_{Medicaid}$ = 0.536 (S2 level)", "#D55E00", "-", 3.0),
          (1.0, "$\\Delta_{Medicaid}$ = 1.00", "#7B2D8E", "--", 2.4)]
    xs = np.linspace(0, 1, 101)
    for ax, (c, (gp, ga)) in zip(axes, CASES.items()):
        for dm, lab, col, ls, lw in st:
            rp_, rv_ = xs * (1 - gp), np.full_like(xs, dm * (1 - ga))
            mean, lo, hi = (SURF["final"].at(rp_, rv_, k) for k in ("mean", "lower", "upper"))
            ax.fill_between(xs, lo, hi, color=col, alpha=0.08, lw=0)
            ax.plot(xs, lo, color=col, lw=0.9, ls=":"); ax.plot(xs, hi, color=col, lw=0.9, ls=":")
            ax.plot(xs, mean, color=col, ls=ls, lw=lw, label=lab)
        ax.axhline(1.1, color="#555", lw=0.9, ls="--"); ax.text(0.99, 1.105, "10% increase", ha="right", fontsize=8)
        ax.set_xlim(0, 1); ax.grid(alpha=0.25)
        ax.set_xlabel("Hypothetical prevention funding loss  ($\\Delta_{prevention}$)")
        ax.set_title(f"{LABEL[c]} (what-if)\n$\\gamma_{{PrEP}}$ = {gp:.2f}, $\\gamma_{{ART}}$ = {ga:.2f}", loc="left", fontweight="bold")
    axes[0].set_ylabel("Incidence rate ratio (final year)"); axes[0].legend(loc="upper left", fontsize=8)
    watermark(fig); fig.savefig(OUT / fname); plt.close(fig)


def points_fig(fname):
    sm, sa = summ("main", "final"), summ("alt", "final")
    ks = ["P1", "P2", "P3"]; labs = {"P1": "P1\nART only", "P2": "P2\nPrEP only", "P3": "P3\nboth"}
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5.2), layout="constrained")
    w = 0.32
    for j, (set_, col, nm) in enumerate([("main", TEAL, "$r=q(1-\\gamma)$ (main)"), ("alt", GREY, "$r=q$ (sensitivity)")]):
        data = []; pos = []
        for i, k in enumerate(ks):
            v = draws[(draws.set == set_) & (draws.name == k) & (draws["when"] == "final")].irr.dropna().values
            if len(v):
                data.append(v); pos.append(i + (j - 0.5) * (w + 0.04))
        if data:
            bp = a1.boxplot(data, positions=pos, widths=w, whis=(2.5, 97.5), showfliers=False, patch_artist=True,
                            medianprops=dict(color="black", lw=1.5))
            for b in bp["boxes"]:
                b.set(facecolor=col, alpha=0.8, edgecolor="#333")
        s = sm if set_ == "main" else sa
        x = np.arange(3) + (j - 0.5) * (w + 0.04)
        a2.bar(x, s.loc[ks, "pct"], w, color=col, alpha=0.85, label=nm)
        a2.errorbar(x, s.loc[ks, "pct"], yerr=[s.loc[ks, "pct"] - s.loc[ks, "pct_lo"], s.loc[ks, "pct_hi"] - s.loc[ks, "pct"]],
                    fmt="none", ecolor="#222", capsize=3)
    for a in (a1, a2):
        a.set_xticks(range(3)); a.set_xticklabels([labs[k] for k in ks]); a.grid(axis="y", alpha=0.25)
    a1.axhline(1, color="#555", lw=0.8); a1.set_ylabel("Incidence rate ratio, final year")
    a1.set_title("Distribution of draws (box: IQR; whiskers: 2.5 to 97.5%)", loc="left", fontsize=10)
    a2.set_ylabel("Increase in incidence vs. status quo (%)"); a2.set_ylim(bottom=0); a2.legend(loc="upper left", fontsize=8)
    a2.set_title("Mean increase with 95% intervals", loc="left", fontsize=10)
    watermark(fig); fig.savefig(OUT / fname); plt.close(fig)


def policy_fig(fname):
    ks = ["S1", "S2", "S3"]; w = 0.26
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5.0), layout="constrained")
    for j, c in enumerate(CASES):
        d = draws[(draws.set == "policy") & (draws.case == c) & (draws["when"] == "final")]
        pos = np.arange(3) + (j - 1) * (w + 0.03)
        bp = a1.boxplot([d[d.name == k].irr.values for k in ks], positions=pos, widths=w, whis=(2.5, 97.5),
                        showfliers=False, patch_artist=True, medianprops=dict(color="black", lw=1.4))
        for b in bp["boxes"]:
            b.set(facecolor=CASE_COLOR[c], alpha=0.8, edgecolor="#333")
        s = summ("policy", "final", c).loc[ks]
        a2.bar(np.arange(3) + (j - 1) * w, s.pct, w, color=CASE_COLOR[c], alpha=0.85, label=LABEL[c])
        a2.errorbar(np.arange(3) + (j - 1) * w, s.pct, yerr=[s.pct - s.pct_lo, s.pct_hi - s.pct], fmt="none", ecolor="#222", capsize=3)
    nm = {"S1": "S1\nprevention halved", "S2": "S2\nMedicaid loss 0.536", "S3": "S3\nboth"}
    for a in (a1, a2):
        a.set_xticks(range(3)); a.set_xticklabels([nm[k] for k in ks], fontsize=8); a.grid(axis="y", alpha=0.25)
    a1.axhline(1, color="#555", lw=0.8); a1.set_ylabel("Incidence rate ratio, final year")
    a1.set_title("Distribution of draws", loc="left", fontsize=10)
    a2.set_ylabel("Increase in incidence vs. status quo (%)"); a2.set_ylim(bottom=0); a2.legend(loc="upper left", fontsize=8)
    a2.set_title("Mean increase with 95% intervals", loc="left", fontsize=10)
    watermark(fig); fig.savefig(OUT / fname); plt.close(fig)


# ---- tables and macros ---------------------------------------------------------
def ci(m, lo, hi):
    return "n/a (outside surrogate range)" if not np.isfinite(m) else f"{m:.2f} [{lo:.2f}, {hi:.2f}]"


def write_tables_and_macros():
    pc = lambda v: f"{100*v:.0f}"
    rng = lambda col: (float(resp[col].min()), float(np.median(resp[col])), float(resp[col].max()))
    with open(OUT / "tab_summary.tex", "w") as fh:
        fh.write("\\begin{table}[htbp]\\centering\\small\n\\begin{tabular}{llccc}\\toprule\nChannel & Quantity & Median & Minimum & Maximum\\\\\\midrule\n")
        for ch, rows in [("PrEP", [("Stated reduction $q$", "q_prep"), ("Other-access share $\\gamma$", "gamma_prep"),
                                   ("Effective reduction $r=q(1-\\gamma)$", "r_prep")]),
                         ("ART adherence", [("Stated reduction $q$", "q_art"), ("Other-access share $\\gamma$", "gamma_art"),
                                            ("Effective reduction $r=q(1-\\gamma)$", "r_vs")])]:
            for i, (lab, col) in enumerate(rows):
                lo, md, hi = rng(col)
                fh.write(f"{ch if i == 0 else ''} & {lab} & {pc(md)}\\% & {pc(lo)}\\% & {pc(hi)}\\%\\\\\n")
            if ch == "PrEP":
                fh.write("\\midrule\n")
        fh.write(f"\\bottomrule\\end{{tabular}}\\caption{{Respondent-level summary ($n={N}$). The effective reductions are computed for each "
                 "respondent first and then summarized; the median of $r$ is not the product of the medians of $q$ and $\\gamma$.}\\label{tab:summary}\\end{table}\n")

    def pts_table(fname, set_, label, cap, coords):
        sf, s5 = summ(set_, "final"), summ(set_, "year5")
        with open(OUT / fname, "w") as fh:
            fh.write("\\begin{table}[htbp]\\centering\\small\n\\begin{tabular}{lcccc}\\toprule\n"
                     "Scenario & $r_{PrEP}$ & $r_{VS}$ & Final-year IRR [95\\% interval] & Year-5 IRR [95\\% interval]\\\\\\midrule\n")
            for k, nm in [("P1", "P1 ART only (median)"), ("P2", "P2 PrEP only (median)"), ("P3", "P3 both (median)")]:
                a, b = coords[k]
                fh.write(f"{nm} & {a:.3f} & {b:.3f} & {ci(*(sf.loc[k, c] for c in ('mean', 'lo', 'hi')))} & "
                         f"{ci(*(s5.loc[k, c] for c in ('mean', 'lo', 'hi')))}\\\\\n")
            fh.write(f"\\bottomrule\\end{{tabular}}\\caption{{{cap}}}\\label{{{label}}}\\end{{table}}\n")

    pts_table("tab_points.tex", "main", "tab:points",
              "Incidence rate ratio at the median effective reductions (main interpretation, $r=q(1-\\gamma)$). These values do not "
              "depend on any insulation assumption.", core.p_points(resp))
    pts_table("tab_alt.tex", "alt", "tab:alt",
              "Sensitivity: the same scenarios under the direct interpretation $r=q$ (the stated reductions already include alternative access).",
              core.p_points(resp, alt=True))
    with open(OUT / "tab_policy.tex", "w") as fh:
        fh.write("\\begin{table}[htbp]\\centering\\small\n\\begin{tabular}{llcccc}\\toprule\n"
                 "Insulation (what-if) & Policy & $(\\Delta_{prev},\\Delta_{Med})$ & $(r_{PrEP},r_{VS})$ & Final-year IRR [95\\% interval]\\\\\\midrule\n")
        for c, (gp, ga) in CASES.items():
            s = summ("policy", "final", c)
            for i, k in enumerate(["S1", "S2", "S3"]):
                dp, dm = POL[k]
                lead = f"{LABEL[c]} ($\\gamma_{{PrEP}}={gp:.2f}$, $\\gamma_{{ART}}={ga:.2f}$)" if i == 0 else ""
                fh.write(f"{lead} & {k} & ({dp:g}, {dm:g}) & ({dp*(1-gp):.3f}, {dm*(1-ga):.3f}) & "
                         f"{ci(*(s.loc[k, x] for x in ('mean', 'lo', 'hi')))}\\\\\n")
            if c != list(CASES)[-1]:
                fh.write("\\midrule\n")
        fh.write("\\bottomrule\\end{tabular}\\caption{Hypothetical named policies defined on the funding axes (not survey quantities): S1 prevention funding "
                 "halved, S2 Medicaid coverage loss of $0.536$ (figure cited from AFC's published estimate), S3 both. Their $r$ coordinates, and hence "
                 "their IRR, depend on the assumed insulation.}\\label{tab:policy}\\end{table}\n")
    with open(OUT / "tab_implied.tex", "w") as fh:
        fh.write("\\begin{table}[htbp]\\centering\\small\\begin{tabular}{lccc}\\toprule\nMedian-based scenario & Lower & Median & Higher\\\\\n"
                 " & $(\\Delta_{prev},\\Delta_{Med})$ & $(\\Delta_{prev},\\Delta_{Med})$ & $(\\Delta_{prev},\\Delta_{Med})$\\\\\\midrule\n")
        for k, nm in [("P1", "P1 ART only"), ("P2", "P2 PrEP only"), ("P3", "P3 both")]:
            a, b = core.p_points(resp)[k]
            cells = []
            for c, (gp, ga) in CASES.items():
                dp, dm = a / (1 - gp), b / (1 - ga); txt = f"({dp:.2f}, {dm:.2f})"
                cells.append(txt if max(dp, dm) <= 1 + 1e-9 else f"\\textbf{{{txt}}}$^{{\\dagger}}$")
            fh.write(f"{nm} & " + " & ".join(cells) + "\\\\\n")
        fh.write("\\bottomrule\\end{tabular}\\caption{Hypothetical losses $\\Delta=r/(1-\\gamma)$ that would produce the median effective reductions "
                 "at each assumed insulation level. $^{\\dagger}$A component exceeds 1: not attainable at that insulation.}\\label{tab:implied}\\end{table}\n")

    sf, s5, sa_ = summ("main", "final"), summ("main", "year5"), summ("alt", "final")
    M = {"NResp": str(N), "RpMed": pc(med_p) + "\\%", "RvMed": pc(med_v) + "\\%",
         "RpMin": pc(resp.r_prep.min()) + "\\%", "RpMax": pc(resp.r_prep.max()) + "\\%",
         "RvMin": pc(resp.r_vs.min()) + "\\%", "RvMax": pc(resp.r_vs.max()) + "\\%"}
    for k, nm in [("P1", "One"), ("P2", "Two"), ("P3", "Three")]:
        M[f"Irr{nm}"] = f"{sf.loc[k,'mean']:.2f}"; M[f"Irr{nm}Lo"] = f"{sf.loc[k,'lo']:.2f}"; M[f"Irr{nm}Hi"] = f"{sf.loc[k,'hi']:.2f}"
        v = sa_.loc[k, "mean"]
        M[f"IrrAlt{nm}"] = f"{v:.2f}" if np.isfinite(v) else "n/a"
    M["IrrThreePct"] = f"{sf.loc['P3','pct']:.0f}"
    ex1, ex2, ex3 = sf.loc["P1", "mean"] - 1, sf.loc["P2", "mean"] - 1, sf.loc["P3", "mean"] - 1
    M["ExcSum"] = f"{ex1 + ex2:.2f}"; M["ExcBoth"] = f"{ex3:.2f}"
    M["AddPhrase"] = ("more than additive" if ex3 > 1.05 * (ex1 + ex2) else "less than additive" if ex3 < 0.95 * (ex1 + ex2)
                      else "approximately additive")
    att = []
    for c, (gp, ga) in CASES.items():
        dp, dm = med_p / (1 - gp), med_v / (1 - ga)
        M[f"Dp{c[0].upper()}"] = f"{dp:.2f}"; M[f"Dm{c[0].upper()}"] = f"{dm:.2f}"
        att.append(f"{LABEL[c].lower()}: $\\Delta=({dp:.2f},{dm:.2f})$, " + ("attainable" if max(dp, dm) <= 1 + 1e-9 else "not attainable"))
    M["AttainSentence"] = "; ".join(att)
    for k, nm in [("S1", "One"), ("S2", "Two"), ("S3", "Three")]:
        for c in CASES:
            M[f"Pol{nm}{c[0].upper()}"] = f"{summ('policy', 'final', c).loc[k, 'mean']:.2f}"
    M["DataBanner"] = ("\\par\\medskip{\\centering\\color{red}\\bfseries SYNTHETIC TEST DATA: NOT SURVEY RESULTS\\par}"
                       if STATUS == "synthetic" else "")
    with open(OUT / "values.tex", "w") as fh:
        for k, v in M.items():
            fh.write(f"\\newcommand{{\\{k}}}{{{v}\\xspace}}\n")
    return M


rspace_fig("final", "fig_rspace_final.png")
rspace_fig("year5", "fig_rspace_year5.png")
rspace_fig("final", "fig_rspace_alt.png", alt=True)
delta_fig("fig_delta_whatif.png"); slices_fig("fig_slices.png"); points_fig("fig_points.png"); policy_fig("fig_policy_by_case.png")
M = write_tables_and_macros()
print("assets written to", OUT, "| status:", STATUS, "| n =", N)
