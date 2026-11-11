#!/usr/bin/env python3
"""Summary-statistics mode, step 3: figures, tables and LaTeX macros.

Reads <private-dir>/out_summary/{summary_values.json, cases.json, points_draws.csv} and the tracked, survey-free grid
output/coverage_domain_grid.csv. Writes ONLY into <private-dir>/out_summary (outside Git). The tracked .tex has no
survey numbers; it reads values.tex, tables and figures from there (TEXINPUTS).

  python3 script/generate_summary_assets.py --private-dir /path/outside/git

Paper figures follow the team's mock-ups: fig1_main (Figure 1), figS1_funding_space, figS2_slices, figS3_year5.
Figures for discussion only (NOT for the paper): figD_rect_final (Figure 1 with the insulation rectangles), figD_delta_panels,
figD_box_points, figD_box_policies. The r = q reading is not shown in any table or figure.
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
OUT = PDIR / "out_summary"
V = json.loads((OUT / "summary_values.json").read_text())
meta = json.loads((OUT / "cases.json").read_text())
ST, RP, RV = V["stats"], V["r_prep"], V["r_vs"]
CASES = {k: tuple(v) for k, v in meta["cases"].items()}
POL = {k: tuple(v) for k, v in meta["policies"].items()}
LABEL = {"lower_insulation": "Lower other-access share", "median_insulation": "Median other-access share",
         "higher_insulation": "Higher other-access share"}
CASE_COLOR = {"lower_insulation": "#0072B2", "median_insulation": "#E69F00", "higher_insulation": "#7B2D8E"}
TEAL, PURPLE, GREY, NOTE = "#1d6a85", "#6a3d8f", "#777777", "#3a4a55"
plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "axes.labelsize": 11, "savefig.dpi": 200})
draws = pd.read_csv(OUT / "points_draws.csv")
SURF = {"final": core.Surface("final"), "year5": core.Surface("year5")}
MED = "median_insulation"


def summ(set_, when, case="all"):
    d = draws[(draws.set == set_) & (draws["when"] == when) & (draws.case == case)]
    g = d.groupby("name").irr
    s = pd.DataFrame({"mean": g.mean(), "lo": g.quantile(0.025), "hi": g.quantile(0.975)})
    s["pct"], s["pct_lo"], s["pct_hi"] = 100 * (s["mean"] - 1), 100 * (s.lo - 1), 100 * (s.hi - 1)
    return s


def levels(zmax, step=0.2, start=1.1):
    return np.round(np.arange(start, zmax, step), 2)


# ---------------------------------------------------------------- paper figures (follow the mock-ups)
def paper_rspace(when, fname, title, notes, step, start, two_line_label, clean=False, outdir=None):
    s = SURF[when]; zmax = float(np.nanmax(s.mean_grid)); mp, mv = RP["med"], RV["med"]
    irr = summ("main", when).loc["P3", "mean"]
    fig = plt.figure(figsize=(10, 6.2 if clean else 7.6))
    ax, cax = ((fig.add_axes([0.095, 0.12, 0.76, 0.84]), fig.add_axes([0.875, 0.12, 0.022, 0.84])) if clean else
               (fig.add_axes([0.095, 0.30, 0.76, 0.63]), fig.add_axes([0.875, 0.30, 0.022, 0.63])))
    m = ax.pcolormesh(s.r_prep_axis, s.r_vs_axis, s.mean_grid, cmap="YlOrRd", vmin=1.0, shading="gouraud")
    cs = ax.contour(s.r_prep_axis, s.r_vs_axis, s.mean_grid, levels=levels(zmax, step, start), colors="#444444", linewidths=0.6)
    ax.clabel(cs, fmt="%.1f", fontsize=8)
    if clean:   # Anna, 28 Oct 2026: one outline rectangle around the marginal min-max bounds, median dot, no triangles
        ax.add_patch(Rectangle((RP["low"], RV["low"]), RP["high"] - RP["low"], RV["high"] - RV["low"], fill=False, ec=TEAL, lw=2.4, zorder=7))
        ax.plot(mp, mv, "o", ms=11, mfc=TEAL, mec="white", mew=2, zorder=8)
    else:
        ax.errorbar(mp, mv, xerr=[[mp - RP["low"]], [RP["high"] - mp]], yerr=[[mv - RV["low"]], [RV["high"] - mv]], fmt="o", ms=11,
                    color=TEAL, mfc=TEAL, mec="white", mew=2, capsize=5, elinewidth=2.2, capthick=2.2, zorder=8)
        ax.plot(0, mv, "^", ms=10, mfc="none", mec=TEAL, mew=2, clip_on=False, zorder=8)
        ax.plot(mp, 0, "^", ms=10, mfc="none", mec=TEAL, mew=2, clip_on=False, zorder=8)
        ax.text(0.015, mv + 0.022, "ART component", color=TEAL, fontsize=9.5, zorder=9)
        ax.text(mp + 0.015, 0.02, "PrEP component", color=TEAL, fontsize=9.5, zorder=9)
    lab = f"median $r_{{PrEP}}$, $r_{{VS}}$\nIRR {irr:.2f}"
    ax.annotate(lab, (mp, mv), xytext=(mp + 0.13, 0.02), fontsize=11, color=TEAL, fontweight="bold", zorder=9,
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.35), arrowprops=dict(arrowstyle="-", color=TEAL))
    ax.set_xlim(0, core.R_MAX); ax.set_ylim(0, core.R_MAX)
    ax.set_xlabel("Proportional reduction in PrEP use  ($r_{PrEP}$)")
    ax.set_ylabel("Proportional reduction in viral suppression  ($r_{VS}$; ART adherence as proxy)")
    cb = fig.colorbar(m, cax=cax); cb.set_label(f"Incidence rate ratio vs. status quo, {'final year' if when == 'final' else 'year 5'}")
    if not clean:
        fig.text(0.095, 0.958, title, fontsize=14, fontweight="bold")
        fig.text(0.095, 0.205, "\n".join(notes), fontsize=9.5, color=NOTE, va="top", linespacing=1.55)
    fig.savefig((outdir or OUT) / fname); plt.close(fig)


def fig_s1(fname, clean=False, outdir=None):
    gp, ga = CASES[MED]; x, y, z = SURF["final"].on_delta_grid(gp, ga); zmax = float(np.nanmax(z))
    pol = summ("policy", "final", MED)
    fig = plt.figure(figsize=(10, 6.2 if clean else 8.0))
    ax, cax = ((fig.add_axes([0.105, 0.13, 0.75, 0.83]), fig.add_axes([0.875, 0.13, 0.022, 0.83])) if clean else
               (fig.add_axes([0.105, 0.37, 0.75, 0.56]), fig.add_axes([0.875, 0.37, 0.022, 0.56])))
    m = ax.pcolormesh(x, y, z, cmap="YlOrRd", vmin=1.0, shading="gouraud")
    cs = ax.contour(x, y, z, levels=levels(zmax, 0.2, 1.1), colors="#444444", linewidths=0.6); ax.clabel(cs, fmt="%.1f", fontsize=8)
    for k, (dp, dm) in POL.items():
        ax.plot(dp, dm, "s", ms=10, mfc=PURPLE, mec="white", mew=1.6, clip_on=False, zorder=8)
        ax.text(dp + 0.02, dm + 0.035, k, fontsize=12, fontweight="bold", color=PURPLE, zorder=9)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_xlabel("Hypothetical reduction in prevention funding reaching Chicago  ($\\Delta_{prevention}$)")
    ax.set_ylabel("Hypothetical loss of Medicaid coverage among\npeople in HIV care  ($\\Delta_{Medicaid}$)")
    fig.colorbar(m, cax=cax).set_label("Incidence rate ratio vs. status quo, final year")
    if clean:
        fig.savefig((outdir or OUT) / fname); plt.close(fig); return
    fig.text(0.105, 0.965, "Figure S1.  Hypothetical funding space, at the median other-access shares", fontsize=14, fontweight="bold")
    qp, qa = ST["q_prep"]["med"], ST["q_art"]["med"]
    txt = (f"The survey did not ask about funding reductions, so this view is a what-if. Coverage reductions are derived as $r=\\Delta(1-\\gamma)$ using the\n"
           f"median other-access shares ($\\gamma_{{PrEP}}$ = {gp:.2f}, $\\gamma_{{ART}}$ = {ga:.2f}); unlike Figure 1, this panel therefore does depend on that assumption.\n\n"
           f"S1 prevention funding halved \u2014 IRR {pol.loc['S1','mean']:.2f}.  Illustrative only: 0.50 was carried forward from an earlier mock-up and has no external source.\n"
           f"S2 Medicaid coverage loss of 0.536 \u2014 IRR {pol.loc['S2','mean']:.2f}.  From AFC's published estimate, 10,000 of 18,640 people with HIV in Illinois.\n"
           f"S3 both \u2014 IRR {pol.loc['S3','mean']:.2f}.\n\n"
           f"The survey estimate is deliberately not plotted here. Because $\\Delta_i=q_i$ exactly, the median implied funding loss is the median stated\n"
           f"reduction, {qp:.2f} for prevention and {qa:.2f} for the Medicaid channel \u2014 but that is a different respondent from the median of $r_i$ shown in\n"
           f"Figure 1. A single marker cannot represent the survey consistently in both coordinate systems.")
    fig.text(0.105, 0.30, txt, fontsize=9.2, color=NOTE, va="top", linespacing=1.5)
    fig.savefig(OUT / fname); plt.close(fig)


def fig_s2(fname, clean=False, outdir=None):
    st = [(0.0, "$\\Delta_{Medicaid}$ = 0", "#4a90b8", "-", 2.6), (0.536, "$\\Delta_{Medicaid}$ = 0.536  (S2 level)", "#d17a33", "-", 2.6),
          (1.0, "$\\Delta_{Medicaid}$ = 1.00", "#6a3d8f", "--", 2.6)]
    xs = np.linspace(0, 1, 101)
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.7), sharey=True)
    fig.subplots_adjust(left=0.065, right=0.985, top=0.86 if clean else 0.80, bottom=0.14 if clean else 0.34, wspace=0.07)
    for ax, (c, (gp, ga)) in zip(axes, CASES.items()):
        for dm, lab, col, ls, lw in st:
            rp_, rv_ = xs * (1 - gp), np.full_like(xs, dm * (1 - ga))
            mean, lo, hi = (SURF["final"].at(rp_, rv_, k) for k in ("mean", "lower", "upper"))
            ax.fill_between(xs, lo, hi, color=col, alpha=0.12, lw=0)
            ax.plot(xs, mean, color=col, ls=ls, lw=lw, label=lab)
        pol = summ("policy", "final", c)
        for k, (dp, dm) in POL.items():
            ax.plot(dp, pol.loc[k, "mean"], "s", ms=8, mfc=PURPLE, mec="white", mew=1.4, clip_on=False, zorder=8)
            ax.text(dp + 0.025, pol.loc[k, "mean"] + 0.03, k, fontsize=9.5, fontweight="bold", color=PURPLE, zorder=9)
        ax.axhline(1.1, color="#555", lw=0.9, ls=":"); ax.text(0.99, 1.115, "10% increase", fontsize=8, ha="right")
        ax.set_xlim(0, 1); ax.grid(alpha=0.25)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_xlabel("Hypothetical prevention funding loss  ($\\Delta_{prevention}$)")
        ax.set_title(f"{LABEL[c]}\n$\\gamma_{{PrEP}}$ = {gp:.2f},  $\\gamma_{{ART}}$ = {ga:.2f}", loc="left", fontweight="bold")
    axes[0].set_ylabel("Incidence rate ratio, final year"); axes[0].legend(loc="upper left", fontsize=8.5)
    if clean:
        fig.savefig((outdir or OUT) / fname); plt.close(fig); return
    fig.text(0.065, 0.955, "Figure S2.  Incidence rate ratio against hypothetical prevention funding loss, at fixed hypothetical Medicaid coverage loss",
             fontsize=14, fontweight="bold")
    fig.text(0.065, 0.19, "A what-if view: the survey did not ask about funding reductions. Coverage reductions are derived as $r=\\Delta(1-\\gamma)$, so these curves depend on the assumed other-access shares, unlike Figure 1.\n"
             "Squares mark the named policies S1 (prevention funding halved \u2014 illustrative, no external source), S2 (Medicaid coverage loss of 0.536, from AFC's published estimate) and S3 (both).\n"
             "Shaded bands are 95% model-uncertainty intervals.", fontsize=9.5, color=NOTE, va="top", linespacing=1.55)
    fig.savefig(OUT / fname); plt.close(fig)


# ---------------------------------------------------------------- figures for discussion only (not for the paper)
def rect_fig(fname):
    s = SURF["final"]; zmax = float(np.nanmax(s.mean_grid)); mp, mv = RP["med"], RV["med"]
    fig, ax = plt.subplots(figsize=(8.6, 7.3), layout="constrained")
    m = ax.pcolormesh(s.r_prep_axis, s.r_vs_axis, s.mean_grid, cmap="YlOrRd", vmin=1.0, shading="gouraud")
    cs = ax.contour(s.r_prep_axis, s.r_vs_axis, s.mean_grid, levels=levels(zmax, 0.3, 1.2), colors="#444444", linewidths=0.6)
    ax.clabel(cs, fmt="%.1f", fontsize=7)
    for c, (gp, ga) in CASES.items():
        x0, x1 = ST["q_prep"]["min"] * (1 - gp), ST["q_prep"]["max"] * (1 - gp)
        y0, y1 = ST["q_art"]["min"] * (1 - ga), ST["q_art"]["max"] * (1 - ga)
        ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, ec=CASE_COLOR[c], lw=2.2, zorder=6,
                               label=f"{LABEL[c]} ($\\gamma$ = {gp:.2f}, {ga:.2f}): range of stated reductions $q$ as $r=\\Delta(1-\\gamma)$"))
        ax.plot(ST["q_prep"]["med"] * (1 - gp), ST["q_art"]["med"] * (1 - ga), "D", ms=7, mfc=CASE_COLOR[c], mec="white", zorder=7)
    ax.errorbar(mp, mv, xerr=[[mp - RP["low"]], [RP["high"] - mp]], yerr=[[mv - RV["low"]], [RV["high"] - mv]], fmt="o", ms=12, color=TEAL,
                mfc=TEAL, mec="white", mew=2, capsize=7, elinewidth=2.4, capthick=2.4, zorder=8,
                label="median of respondent-level $r_i$; bars = marginal min\u2013max (not a joint range)")
    ax.plot(0, mv, "^", ms=11, mfc="white", mec="black", mew=2, clip_on=False, zorder=8)
    ax.plot(mp, 0, "^", ms=11, mfc="white", mec="black", mew=2, clip_on=False, zorder=8)
    box = dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.8)
    ax.text(0.012, mv + 0.012, "ART component", fontsize=9, fontweight="bold", zorder=9, bbox=box)
    ax.text(mp + 0.012, 0.012, "PrEP component", fontsize=9, fontweight="bold", zorder=9, bbox=box)
    ax.text(mp + 0.015, mv + 0.015, "both", fontsize=9, fontweight="bold", zorder=9, bbox=box)
    ax.set_xlim(0, core.R_MAX); ax.set_ylim(0, core.R_MAX)
    ax.set_xlabel("Proportional reduction in PrEP use  ($r_{PrEP}$)")
    ax.set_ylabel("Proportional reduction in viral suppression  ($r_{VS}$; ART adherence as proxy)")
    ax.set_title("Variant of Figure 1 with the funding what-if rectangles (final year)", loc="left", fontsize=11, fontweight="bold")
    ax.legend(loc="upper left", fontsize=7.5, framealpha=0.92)
    fig.colorbar(m, ax=ax, shrink=0.85, pad=0.02).set_label("Incidence rate ratio vs. status quo")
    fig.savefig(OUT / fname); plt.close(fig)


def delta_fig(fname):
    zmax = max(float(np.nanmax(SURF["final"].on_delta_grid(*g)[2])) for g in CASES.values())
    fig, axes = plt.subplots(1, 3, figsize=(17, 6.5), sharey=True, layout="constrained")
    for ax, (c, (gp, ga)) in zip(axes, CASES.items()):
        x, y, z = SURF["final"].on_delta_grid(gp, ga)
        ax.set_facecolor("#dcdcdc")
        mm = ax.pcolormesh(x, y, z, cmap="YlOrRd", vmin=1.0, vmax=zmax, shading="gouraud", zorder=1)
        cs = ax.contour(x, y, z, levels=levels(zmax), colors="#444444", linewidths=0.6); ax.clabel(cs, fmt="%.1f", fontsize=7)
        for k, (dp, dm) in POL.items():
            ax.plot(dp, dm, "s", ms=9, mfc=PURPLE, mec="white", mew=1.5, clip_on=False, zorder=8)
            ax.text(dp + (0.02 if dp < 0.9 else -0.02), dm + 0.025, k, fontsize=11, fontweight="bold", color=PURPLE,
                    ha="left" if dp < 0.9 else "right", zorder=9, bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.75))
        ax.add_patch(Rectangle((ST["q_prep"]["min"], ST["q_art"]["min"]), ST["q_prep"]["max"] - ST["q_prep"]["min"],
                               ST["q_art"]["max"] - ST["q_art"]["min"], fill=False, ec=CASE_COLOR[c], lw=2.2, zorder=6))
        ax.plot(ST["q_prep"]["med"], ST["q_art"]["med"], "D", ms=7, mfc=CASE_COLOR[c], mec="white", zorder=7)
        dp3, dm3 = RP["med"] / (1 - gp), RV["med"] / (1 - ga)
        if max(dp3, dm3) <= 1 + 1e-9:
            ax.plot(dp3, dm3, "o", ms=11, mfc="none", mec=TEAL, mew=2.4, zorder=8, clip_on=False)
            ax.annotate("both components,\nre-expressed", (dp3, dm3), xytext=(max(dp3 - 0.5, 0.04), min(dm3 + 0.3, 0.9)), fontsize=8.5, color=TEAL,
                        arrowprops=dict(arrowstyle="-", color=TEAL), zorder=9, bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.85))
        else:
            ax.plot(min(dp3, 1), min(dm3, 1), "^", ms=11, mfc="#d62728", mec="black", clip_on=False, zorder=9)
            ax.text(0.97, 0.12, f"Red triangle: both components would need\n$\\Delta_{{prev}}$ = {dp3:.2f}, $\\Delta_{{Med}}$ = {dm3:.2f}\n(above 1: not attainable here)",
                    ha="right", va="bottom", fontsize=8.5, zorder=9, bbox=dict(fc="white", ec="#d62728", alpha=0.92))
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        ax.set_xlabel("Hypothetical prevention funding loss  ($\\Delta_{prevention}$)")
        ax.set_title(f"{LABEL[c]} (what-if)\n$\\gamma_{{PrEP}}$ = {gp:.2f}, $\\gamma_{{ART}}$ = {ga:.2f}", loc="left", fontweight="bold")
        if c != list(CASES)[0]:
            ax.tick_params(labelleft=False)
    axes[0].set_ylabel("Hypothetical Medicaid coverage loss  ($\\Delta_{Medicaid}$)")
    axes[-1].tick_params(labelright=True, right=True)
    axes[-1].yaxis.set_label_position("right")
    axes[-1].set_ylabel("Hypothetical Medicaid coverage loss  ($\\Delta_{Medicaid}$)")
    cb = fig.colorbar(mm, ax=axes, location="bottom", shrink=0.4, pad=0.03); cb.set_label("Incidence rate ratio vs. status quo (shared scale)")
    fig.savefig(OUT / fname); plt.close(fig)


def box_points(fname):
    ks = ["P2", "P1", "P3"]; labs = {"P2": "PrEP\ncomponent", "P1": "ART\ncomponent", "P3": "both"}
    fig, ax = plt.subplots(figsize=(7, 5), layout="constrained")
    data = [draws[(draws.set == "main") & (draws.name == k) & (draws["when"] == "final")].irr.values for k in ks]
    bp = ax.boxplot(data, widths=0.5, whis=(2.5, 97.5), showfliers=False, patch_artist=True, medianprops=dict(color="black", lw=1.5))
    for b in bp["boxes"]:
        b.set(facecolor=TEAL, alpha=0.8, edgecolor="#333")
    ax.set_xticks(range(1, 4)); ax.set_xticklabels([labs[k] for k in ks]); ax.grid(axis="y", alpha=0.25)
    ax.axhline(1, color="#555", lw=0.8); ax.set_ylabel("Incidence rate ratio, final year")
    ax.set_title("Model-uncertainty draws (box: IQR; whiskers: 2.5 to 97.5%)", loc="left", fontsize=10)
    fig.savefig(OUT / fname); plt.close(fig)


def box_policies(fname):
    ks = ["S1", "S2", "S3"]; w = 0.26
    fig, ax = plt.subplots(figsize=(8, 5), layout="constrained")
    for j, c in enumerate(CASES):
        d = draws[(draws.set == "policy") & (draws.case == c) & (draws["when"] == "final")]
        pos = np.arange(3) + (j - 1) * (w + 0.03)
        bp = ax.boxplot([d[d.name == k].irr.values for k in ks], positions=pos, widths=w, whis=(2.5, 97.5), showfliers=False,
                        patch_artist=True, medianprops=dict(color="black", lw=1.4))
        for b in bp["boxes"]:
            b.set(facecolor=CASE_COLOR[c], alpha=0.8, edgecolor="#333")
        ax.plot([], [], "s", color=CASE_COLOR[c], label=LABEL[c])
    nm = {"S1": "S1\nprevention halved", "S2": "S2\nMedicaid loss 0.536", "S3": "S3\nboth"}
    ax.set_xticks(range(3)); ax.set_xticklabels([nm[k] for k in ks], fontsize=8); ax.grid(axis="y", alpha=0.25)
    ax.axhline(1, color="#555", lw=0.8); ax.set_ylabel("Incidence rate ratio, final year"); ax.legend(loc="upper left", fontsize=8)
    ax.set_title("Model-uncertainty draws at the named policies, by assumed other-access share", loc="left", fontsize=10)
    fig.savefig(OUT / fname); plt.close(fig)


def ci(m, lo, hi):
    return f"{m:.2f} [{lo:.2f}, {hi:.2f}]"


def write_tables_and_macros():
    pc = lambda v: f"{100*v:.0f}"
    with open(OUT / "tab_stats.tex", "w") as fh:
        fh.write("\\begin{table}[htbp]\\centering\\small\n\\begin{tabular}{llccc}\\toprule\nChannel & Quantity & Median & Minimum & Maximum\\\\\\midrule\n")
        for ch, rows in [("PrEP", [("Stated reduction $q$", "q_prep"), ("Other-access share $\\gamma$", "g_prep")]),
                         ("ART adherence", [("Stated reduction $q$", "q_art"), ("Other-access share $\\gamma$", "g_art")])]:
            for i, (lab, k) in enumerate(rows):
                fh.write(f"{ch if i == 0 else ''} & {lab} & {pc(ST[k]['med'])}\\% & {pc(ST[k]['min'])}\\% & {pc(ST[k]['max'])}\\%\\\\\n")
            rk = RP if ch == "PrEP" else RV
            fh.write(f" & Effective reduction $r_i=q_i(1-\\gamma_i)$ & {100*rk['med']:.1f}\\% & {100*rk['low']:.1f}\\% & {100*rk['high']:.1f}\\%\\\\\n")
            qk, gk = ("q_prep", "g_prep") if ch == "PrEP" else ("q_art", "g_art")
            fh.write(f" & \\emph{{For comparison, not used:}} $q_{{med}}(1-\\gamma_{{med}})$ & {100*ST[qk]['med']*(1-ST[gk]['med']):.1f}\\% & & \\\\\n")
            if ch == "PrEP":
                fh.write("\\midrule\n")
        fh.write("\\bottomrule\\end{tabular}\\caption{Summary across the five respondents (median, minimum, maximum). The effective reductions $r_i$ were computed by "
                 "respondent and then summarized; they are not products of the marginal medians (see the section on medians below).}\\label{tab:stats}\\end{table}\n")

    sf, s5 = summ("main", "final"), summ("main", "year5")
    main_pts = {"P2": (RP["med"], 0), "P1": (0, RV["med"]), "P3": (RP["med"], RV["med"]),
                "low": (RP["low"], RV["low"]), "high": (RP["high"], RV["high"])}
    names = {"P2": "PrEP component (median)", "P1": "ART component (median)", "P3": "Both (median)",
             "low": "Both, marginal minima", "high": "Both, marginal maxima"}
    with open(OUT / "tab_points.tex", "w") as fh:
        fh.write("\\begin{table}[htbp]\\centering\\small\n\\begin{tabular}{lcccc}\\toprule\n"
                 "Scenario & $r_{PrEP}$ & $r_{VS}$ & Final-year IRR [95\\% interval] & Year-5 IRR [95\\% interval]\\\\\\midrule\n")
        for k in ["P2", "P1", "P3", "low", "high"]:
            a, b = main_pts[k]
            fh.write(f"{names[k]} & {a:.3f} & {b:.3f} & {ci(*(sf.loc[k, c] for c in ('mean', 'lo', 'hi')))} & {ci(*(s5.loc[k, c] for c in ('mean', 'lo', 'hi')))}\\\\\n")
        fh.write("\\bottomrule\\end{tabular}\\caption{Incidence rate ratio at the respondent-level effective reductions ($r=q(1-\\gamma)$). The median rows use the median of the "
                 "respondent-level $r_i$ in each channel; the last two rows combine the marginal minima and the marginal maxima of the two channels, a pairing that no single respondent "
                 "need have reported. These values do not depend on any insulation assumption.}\\label{tab:points}\\end{table}\n")
    with open(OUT / "tab_policy.tex", "w") as fh:
        fh.write("\\begin{table}[htbp]\\centering\\small\n\\resizebox{\\textwidth}{!}{\\begin{tabular}{llcccc}\\toprule\n"
                 "Other-access share (what-if) & Policy & $(\\Delta_{prev},\\Delta_{Med})$ & $(r_{PrEP},r_{VS})$ & Final-year IRR [95\\% interval]\\\\\\midrule\n")
        for c, (gp, ga) in CASES.items():
            s = summ("policy", "final", c)
            for i, k in enumerate(["S1", "S2", "S3"]):
                dp, dm = POL[k]
                lead = f"{LABEL[c]} ($\\gamma_{{PrEP}}={gp:.2f}$, $\\gamma_{{ART}}={ga:.2f}$)" if i == 0 else ""
                fh.write(f"{lead} & {k} & ({dp:g}, {dm:g}) & ({dp*(1-gp):.3f}, {dm*(1-ga):.3f}) & {ci(*(s.loc[k, x] for x in ('mean', 'lo', 'hi')))}\\\\\n")
            if c != list(CASES)[-1]:
                fh.write("\\midrule\n")
        fh.write("\\bottomrule\\end{tabular}}\\caption{Hypothetical named policies defined on the funding axes (not survey quantities): S1 prevention funding halved "
                 "(illustrative; no external source), S2 Medicaid coverage loss of $0.536$ (AFC's published estimate), S3 both. Their $r$ coordinates, and hence their IRR, depend "
                 "on the assumed other-access share.}\\label{tab:policy}\\end{table}\n")
    with open(OUT / "tab_implied.tex", "w") as fh:
        fh.write("\\begin{table}[htbp]\\centering\\small\\begin{tabular}{lccc}\\toprule\nScenario (median $r_i$) & Lower & Median & Higher\\\\\n"
                 " & $(\\Delta_{prev},\\Delta_{Med})$ & $(\\Delta_{prev},\\Delta_{Med})$ & $(\\Delta_{prev},\\Delta_{Med})$\\\\\\midrule\n")
        for k, nm in [("P2", "PrEP component"), ("P1", "ART component"), ("P3", "Both")]:
            a, b = main_pts[k]; cells = []
            for c, (gp, ga) in CASES.items():
                dp, dm = a / (1 - gp), b / (1 - ga); txt = f"({dp:.2f}, {dm:.2f})"
                cells.append(txt if max(dp, dm) <= 1 + 1e-9 else f"\\textbf{{{txt}}}$^{{\\dagger}}$")
            fh.write(f"{nm} & " + " & ".join(cells) + "\\\\\n")
        fh.write("\\bottomrule\\end{tabular}\\caption{Hypothetical losses $\\Delta=r/(1-\\gamma)$ that would produce the median effective reductions at each assumed "
                 "other-access share. $\\Delta$ is not the stated reduction $q$ here because the medians of $r_i$ are not products of the marginal medians. "
                 "$^{\\dagger}$A component exceeds 1: not attainable.}\\label{tab:implied}\\end{table}\n")
    sm = summ("mean", "final")
    M = {"QpMed": pc(ST["q_prep"]["med"]) + "\\%", "QvMed": pc(ST["q_art"]["med"]) + "\\%",
         "GpMed": pc(ST["g_prep"]["med"]) + "\\%", "GaMed": pc(ST["g_art"]["med"]) + "\\%",
         "OneMinusGpMed": pc(1 - ST["g_prep"]["med"]) + "\\%", "OneMinusGaMed": pc(1 - ST["g_art"]["med"]) + "\\%",
         "RpMed": f"{100*RP['med']:.1f}\\%", "RvMed": f"{100*RV['med']:.1f}\\%",
         "RpLow": f"{100*RP['low']:.1f}\\%", "RpHigh": f"{100*RP['high']:.1f}\\%", "RvLow": f"{100*RV['low']:.1f}\\%", "RvHigh": f"{100*RV['high']:.1f}\\%",
         "ProdMedP": f"{100*ST['q_prep']['med']*(1-ST['g_prep']['med']):.1f}\\%", "ProdMedA": f"{100*ST['q_art']['med']*(1-ST['g_art']['med']):.1f}\\%",
         "MeanQp": f"{100*ST['q_prep']['mean']:.1f}\\%", "MeanGp": f"{100*ST['g_prep']['mean']:.1f}\\%",
         "MeanQa": f"{100*ST['q_art']['mean']:.1f}\\%", "MeanGa": f"{100*ST['g_art']['mean']:.1f}\\%",
         "MeanRp": f"{100*RP['mean']:.1f}\\%", "MeanRa": f"{100*RV['mean']:.1f}\\%",
         "ProdMeanP": f"{100*ST['q_prep']['mean']*(1-ST['g_prep']['mean']):.1f}\\%", "ProdMeanA": f"{100*ST['q_art']['mean']*(1-ST['g_art']['mean']):.1f}\\%",
         "IrrMeanBoth": f"{sm.loc['P3','mean']:.2f}", "IrrMeanPrep": f"{sm.loc['P2','mean']:.2f}", "IrrMeanArt": f"{sm.loc['P1','mean']:.2f}"}
    for k, nm in [("P2", "Prep"), ("P1", "Art"), ("P3", "Both"), ("low", "Low"), ("high", "High")]:
        M[f"Irr{nm}"] = f"{sf.loc[k,'mean']:.2f}"; M[f"Irr{nm}Lo"] = f"{sf.loc[k,'lo']:.2f}"; M[f"Irr{nm}Hi"] = f"{sf.loc[k,'hi']:.2f}"
    M["IrrBothPct"] = f"{sf.loc['P3','pct']:.0f}"
    ex1, ex2, ex3 = sf.loc["P1", "mean"] - 1, sf.loc["P2", "mean"] - 1, sf.loc["P3", "mean"] - 1
    M["ExcSum"] = f"{ex1 + ex2:.2f}"; M["ExcBoth"] = f"{ex3:.2f}"
    M["AddPhrase"] = ("more than additive" if ex3 > 1.05 * (ex1 + ex2) else "less than additive" if ex3 < 0.95 * (ex1 + ex2) else "approximately additive")
    att = []
    for c, (gp, ga) in CASES.items():
        dp, dm = RP["med"] / (1 - gp), RV["med"] / (1 - ga)
        att.append(f"{LABEL[c].lower()}: $\\Delta=({dp:.2f},{dm:.2f})$, " + ("attainable" if max(dp, dm) <= 1 + 1e-9 else "not attainable"))
    M["AttainSentence"] = "; ".join(att)
    for k, nm in [("S1", "One"), ("S2", "Two"), ("S3", "Three")]:
        for c in CASES:
            M[f"Pol{nm}{c[0].upper()}"] = f"{summ('policy', 'final', c).loc[k, 'mean']:.2f}"
    M["DataBanner"] = ""
    with open(OUT / "values.tex", "w") as fh:
        for k, v in M.items():
            fh.write(f"\\newcommand{{\\{k}}}{{{v}\\xspace}}\n")
    return M


# ---------------------------------------------------------------- build
n1 = ["This surface depends only on the coverage reductions; no insulation assumption enters it.",
      "Filled circle: median of the respondent-level effective reductions, computed per respondent as $r_i=q_i(1-\\gamma_i)$ and then summarised.",
      "Bars: marginal minimum to maximum of each channel separately. NOT a joint respondent range \u2014",
      "no respondent need have reported a corner.",
      "Triangles: the single-channel components of the same estimate."]
paper_rspace("final", "fig1_main.png", "Figure 1.  Incidence rate ratio over the coverage-reduction domain", n1, 0.3, 1.2, True)
paper_rspace("year5", "figS3_year5.png", "Figure S3.  Incidence rate ratio over the coverage-reduction domain, year 5",
             ["As Figure 1 but at year 5 rather than the final year, on its own colour scale. Markers and bars are identical in position:",
              "the coverage reductions do not change with the horizon, only the incidence rate ratio they produce.",
              "Comparing the two shows how the effect accumulates over the intervention period."], 0.2, 1.2, False)
fig_s1("figS1_funding_space.png"); fig_s2("figS2_slices.png")
OV = OUT / "overleaf_figures"; OV.mkdir(exist_ok=True)   # titles/footnotes removed, descriptive names (Anna, 28 Oct 2026)
paper_rspace("final", "IRR_over_coverage_reduction_domain_yr10.png", "", [], 0.3, 1.2, True, clean=True, outdir=OV)
paper_rspace("year5", "IRR_over_coverage_reduction_domain_yr5.png", "", [], 0.2, 1.2, False, clean=True, outdir=OV)
fig_s1("funding_space_at_median_other_access.png", clean=True, outdir=OV)
fig_s2("IRR_vs_prevention_funding_loss_by_Medicaid_loss.png", clean=True, outdir=OV)
rect_fig("figD_rect_final.png"); delta_fig("figD_delta_panels.png"); box_points("figD_box_points.png"); box_policies("figD_box_policies.png")
M = write_tables_and_macros()
print("summary-mode assets written to", OUT)
print({k: M[k] for k in ["RpMed", "RvMed", "IrrPrep", "IrrArt", "IrrBoth", "IrrLow", "IrrHigh", "IrrMeanBoth", "ProdMeanP", "MeanRp", "AddPhrase"]})
