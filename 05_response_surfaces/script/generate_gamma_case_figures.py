"""Figures and LaTeX tables, r-space interpretation.

Survey items 1-2 elicit r (reduction in use); items 3-4 elicit gamma. Delta is
derived: Delta = r / (1 - gamma). The IRR depends on r only, so every outcome is
computed once in r-space; gamma only moves points between r and Delta.

Reads output/{gamma_case_grid,gamma_case_points,coverage_domain_grid}.csv and
script/stakeholder_inputs.csv. Run after export_gamma_cases.R:
  python3 05_response_surfaces/script/generate_gamma_case_figures.py
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "output"
CASES = ["lower_insulation", "median_insulation", "higher_insulation"]
GAMMA = {"lower_insulation": (0.30, 0.5625), "median_insulation": (0.40, 0.75),
         "higher_insulation": (0.50, 0.9375)}              # (gamma_PrEP, gamma_ART)
LABEL = {"lower_insulation": "Lower insulation", "median_insulation": "Survey-median insulation",
         "higher_insulation": "Higher insulation"}
DOM_COLOR = {"lower_insulation": "#0072B2", "median_insulation": "#000000",
             "higher_insulation": "#7B2D8E"}
DOM_STYLE = {"lower_insulation": "-", "median_insulation": "--", "higher_insulation": ":"}
P_COLOR = {"P1": "#CC79A7", "P2": "#009E73", "P3": "#4D4D4D"}
P_NAME = {"P1": "P1  treatment (viral suppression) only", "P2": "P2  prevention (PrEP) only",
          "P3": "P3  both (provisional illustrative point)"}
REF_BASE = dict(prep=23.0, vs=61.0)      # reference baselines for approximate pp axes only

plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
                     "savefig.dpi": 200})

inp = pd.read_csv(HERE / "stakeholder_inputs.csv").set_index("name")
assert (inp.status == "provisional_placeholder").all(), "update labels before using survey data"
PT = (float(inp.loc["point", "r_prep"]), float(inp.loc["point", "r_vs"]))
LO = (float(inp.loc["range_low", "r_prep"]), float(inp.loc["range_low", "r_vs"]))
HI = (float(inp.loc["range_high", "r_prep"]), float(inp.loc["range_high", "r_vs"]))
P_R = {"P1": (0.0, PT[1]), "P2": (PT[0], 0.0), "P3": PT}   # (r_PrEP, r_VS)

grid = pd.read_csv(OUT / "gamma_case_grid.csv")
pts_raw = pd.read_csv(OUT / "gamma_case_points.csv")
cov = pd.read_csv(OUT / "coverage_domain_grid.csv")

# per-draw outcomes at the r-space points (median-case rows hold exactly these r values)
MAP_OLD = {"P1": "P1", "P2": "P2", "P3": "elicited"}
parts = []
for new, old in MAP_OLD.items():
    d = pts_raw[(pts_raw.case == "median_insulation") & (pts_raw.scenario == old)].copy()
    assert np.allclose(d.r_prep, P_R[new][0]) and np.allclose(d.r_art, P_R[new][1]), \
        f"re-run export_gamma_cases.R: stored r for {new} differs from stakeholder_inputs.csv"
    d["scenario"] = new
    parts.append(d)
pts = pd.concat(parts)


def pts_summary(when):
    d = pts[pts["when"] == when].copy()
    d["pct"] = 100 * (d.irr - 1)
    g = d.groupby("scenario")
    s = pd.DataFrame({"mean": g.irr.mean(), "lo": g.irr.quantile(0.025), "hi": g.irr.quantile(0.975),
                      "pct": g.pct.mean(), "pct_lo": g.pct.quantile(0.025), "pct_hi": g.pct.quantile(0.975)})
    return s.loc[["P1", "P2", "P3"]], d


def r_surface(when, fname, title):
    d = cov[cov["when"] == when]
    z = d.pivot(index="r_art", columns="r_prep", values="mean")
    assert z.shape == (101, 101)
    x, y, zz = z.columns.values, z.index.values, z.values
    fig, ax = plt.subplots(figsize=(9.2, 7.6), layout="constrained")
    m = ax.pcolormesh(x, y, zz, cmap="YlOrRd", vmin=1.0, shading="gouraud")
    cs = ax.contour(x, y, zz, levels=np.round(np.arange(1.1, zz.max(), 0.2), 2),
                    colors="#444444", linewidths=0.6)
    ax.clabel(cs, fmt="%.1f", fontsize=7)
    for c in CASES:
        gp, ga = GAMMA[c]
        ax.add_patch(Rectangle((0, 0), 1 - gp, 1 - ga, fill=False, ec=DOM_COLOR[c],
                               ls=DOM_STYLE[c], lw=2.2, zorder=4))
    ax.add_patch(Rectangle(LO, HI[0] - LO[0], HI[1] - LO[1], fc="#ffffff", alpha=0.25,
                           ec="#222222", lw=1.3, hatch="//", zorder=3))
    for k, (px, py) in P_R.items():
        ax.plot(px, py, "o" if k == "P3" else "s", ms=10, mfc="white",
                mec=P_COLOR[k] if k != "P3" else "black", mew=2.2, clip_on=False, zorder=7)
        ax.text(px + 0.014, py + 0.014, k, fontsize=11, fontweight="bold", color="#111111", zorder=8,
                bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.75))
    ax.set_xlim(0, 0.75); ax.set_ylim(0, 0.75)
    ax.set_xlabel("Proportional reduction in PrEP use  ($r_{PrEP}$)")
    ax.set_ylabel("Proportional reduction in viral suppression  ($r_{VS}$)")
    t = ax.secondary_xaxis("top", functions=(lambda r: r * REF_BASE["prep"], lambda a: a / REF_BASE["prep"]))
    t.set_xlabel("Approximate PrEP coverage lost (percentage points)")
    s = ax.secondary_yaxis("right", functions=(lambda r: r * REF_BASE["vs"], lambda a: a / REF_BASE["vs"]))
    s.set_ylabel("Approximate viral suppression lost (percentage points)")
    cb = fig.colorbar(m, ax=ax, pad=0.13, shrink=0.9)
    cb.set_label("Incidence rate ratio vs. status quo")
    handles = [Line2D([0], [0], color=DOM_COLOR[c], ls=DOM_STYLE[c], lw=2.2,
                      label=f"{LABEL[c]}: reachable with $\\Delta\\leq1$ "
                            f"($r_{{PrEP}}\\leq{1-GAMMA[c][0]:g}$, $r_{{VS}}\\leq{1-GAMMA[c][1]:g}$)")
               for c in CASES]
    handles.append(Rectangle((0, 0), 1, 1, fc="white", ec="#222222", hatch="//",
                             label="Provisional illustrative range (placeholder)"))
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.09), fontsize=8, frameon=False)
    ax.set_title(title, loc="left", fontweight="bold", pad=44)
    fig.savefig(OUT / fname); plt.close(fig)
    return float(zz.max())


def delta_map(fname):
    zmax = grid[grid["when"] == "final"]["mean"].max()
    levels = np.round(np.arange(1.1, zmax, 0.2), 2)
    fig, axes = plt.subplots(1, 3, figsize=(17, 6.6), sharey=True, layout="constrained")
    for ax, c in zip(axes, CASES):
        gp, ga = GAMMA[c]
        d = grid[(grid.case == c) & (grid["when"] == "final")]
        z = d.pivot(index="delta_med", columns="delta_prev", values="mean")
        assert z.shape == (101, 101)
        x, y = z.columns.values, z.index.values
        m = ax.pcolormesh(x, y, z.values, cmap="YlOrRd", vmin=1.0, vmax=zmax, shading="gouraud")
        cs = ax.contour(x, y, z.values, levels=levels, colors="#444444", linewidths=0.6)
        ax.clabel(cs, fmt="%.1f", fontsize=7)
        bx = (LO[0] / (1 - gp), HI[0] / (1 - gp)); by = (LO[1] / (1 - ga), HI[1] / (1 - ga))
        w, h = min(bx[1], 1) - bx[0], min(by[1], 1) - by[0]
        ax.add_patch(Rectangle((bx[0], by[0]), w, h, fc="white", alpha=0.25, ec="#222222",
                               lw=1.3, hatch="//", zorder=3))
        infeasible = []
        for k, (rp, rv) in P_R.items():
            dp, dm = rp / (1 - gp), rv / (1 - ga)
            if dp <= 1 + 1e-9 and dm <= 1 + 1e-9:
                ax.plot(dp, dm, "o" if k == "P3" else "s", ms=9, mfc="white",
                        mec="black" if k == "P3" else P_COLOR[k], mew=2.2, clip_on=False, zorder=7)
                ax.text(dp + 0.015, dm + 0.015, k, fontsize=11, fontweight="bold", zorder=8,
                        bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.75))
            else:
                infeasible.append((k, dp, dm))
                xm = min(dp, 1)
                ax.plot(xm, 1.0, "^", ms=11, mfc="#d62728", mec="black", mew=1.2, clip_on=False, zorder=8)
                ax.text(xm + (0.02 if xm < 0.9 else -0.02), 0.94, k, fontsize=11, fontweight="bold",
                        ha="left" if xm < 0.9 else "right", va="top", zorder=9,
                        bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.8))
        if infeasible:
            msg = "; ".join(f"{k}: $\\Delta_{{Medicaid}}$ = {dm:.2f}" for k, _, dm in infeasible)
            ax.text(0.5, 0.05, f"Red triangles: not attainable ($\\Delta>1$)\n{msg}", ha="center",
                    va="bottom", fontsize=8.5, bbox=dict(fc="white", ec="#d62728", alpha=0.92), zorder=9)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        ax.set_xlabel("Reduction in prevention funding  ($\\Delta_{prevention}$)")
        ax.set_title(f"{LABEL[c]}\n$\\gamma_{{PrEP}}$ = {gp:g}, $\\gamma_{{ART}}$ = {ga:g}\n"
                     f"reachable $r_{{PrEP}}$ 0 to {1-gp:g}; $r_{{VS}}$ 0 to {1-ga:g}",
                     loc="left", fontweight="bold", fontsize=10)
    axes[0].set_ylabel("Loss of Medicaid coverage among people in HIV care  ($\\Delta_{Medicaid}$)")
    for ax in axes[1:]:
        ax.tick_params(labelleft=False)
    cb = fig.colorbar(m, ax=axes, location="bottom", shrink=0.4, pad=0.03)
    cb.set_label("Incidence rate ratio vs. status quo (shared scale)")
    fig.savefig(OUT / fname); plt.close(fig)


def slice_figure(fname):
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.2), sharey=True, layout="constrained")
    styles = [(0.0, "$\\Delta_{Medicaid}$ = 0", "#0072B2", "-", 2.4),
              (0.35, "$\\Delta_{Medicaid}$ = 0.35 (illustrative)", "#D55E00", "-", 3.0),
              (1.0, "$\\Delta_{Medicaid}$ = 1.00 (complete loss)", "#7B2D8E", "--", 2.4)]
    for ax, c in zip(axes, CASES):
        gp, ga = GAMMA[c]
        d = grid[(grid.case == c) & (grid["when"] == "final")]
        piv = lambda col: d.pivot(index="delta_med", columns="delta_prev", values=col)
        z, lo, hi = piv("mean"), piv("lower").values, piv("upper").values
        x, y = z.columns.values, z.index.values; z = z.values
        for dm, lab, col, ls, lw in styles:
            i = int(np.argmin(np.abs(y - dm)))
            ax.fill_between(x, lo[i], hi[i], color=col, alpha=0.08, lw=0)
            ax.plot(x, lo[i], color=col, lw=0.9, ls=":", alpha=0.9)
            ax.plot(x, hi[i], color=col, lw=0.9, ls=":", alpha=0.9)
            ax.plot(x, z[i], color=col, ls=ls, lw=lw, label=lab)
        ax.axhline(1.1, color="#555555", lw=0.9, ls="--")
        ax.text(0.99, 1.105, "10% increase", ha="right", fontsize=8)
        ax.set_xlim(0, 1); ax.grid(alpha=0.25)
        ax.set_xlabel("Reduction in prevention funding  ($\\Delta_{prevention}$)")
        ax.set_title(f"{LABEL[c]}\n$\\gamma_{{PrEP}}$ = {gp:g}, $\\gamma_{{ART}}$ = {ga:g}", loc="left",
                     fontweight="bold")
    axes[0].set_ylabel("Incidence rate ratio (final year)")
    axes[0].legend(loc="upper left", fontsize=8)
    fig.savefig(OUT / fname); plt.close(fig)


def points_figure(fname):
    s, d = pts_summary("final")
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.5, 5.0), layout="constrained")
    ks = ["P1", "P2", "P3"]
    bp = a1.boxplot([d[d.scenario == k].irr.values for k in ks], whis=(2.5, 97.5), showfliers=False,
                    patch_artist=True, widths=0.5, medianprops=dict(color="black", lw=1.6))
    for b, k in zip(bp["boxes"], ks):
        b.set(facecolor=P_COLOR[k], alpha=0.8, edgecolor="#333333")
    a1.axhline(1, color="#555555", lw=0.8)
    a1.set_xticklabels([P_NAME[k].replace("  ", "\n", 1).replace(" (", "\n(") for k in ks], fontsize=8)
    a1.set_ylabel("Incidence rate ratio, final year"); a1.grid(axis="y", alpha=0.25)
    a1.set_title("Distribution of draws (box: IQR; whiskers: 2.5 to 97.5%)", loc="left", fontsize=10)
    a2.bar(ks, s.pct, color=[P_COLOR[k] for k in ks], alpha=0.85, width=0.55)
    a2.errorbar(ks, s.pct, yerr=[s.pct - s.pct_lo, s.pct_hi - s.pct], fmt="none", ecolor="#222", capsize=4)
    a2.set_ylabel("Increase in incidence vs. status quo, final year (%)"); a2.set_ylim(bottom=0)
    a2.grid(axis="y", alpha=0.25)
    a2.set_title("Mean increase, 100 x (IRR - 1), with 95% intervals", loc="left", fontsize=10)
    fig.savefig(OUT / fname); plt.close(fig)


def fmt(m, lo, hi):
    return f"{m:.2f} [{lo:.2f}, {hi:.2f}]"


def tables():
    s5, _ = pts_summary("year5"); s, _ = pts_summary("final")
    name = {"P1": "P1 treatment only", "P2": "P2 prevention only", "P3": "P3 both (provisional point)"}
    for tag, ss, cap, lab in [("final", s, "Final-year (post-intervention year 10) incidence rate ratio at the "
                               "provisional $r$-space points, mean [95\\% interval]. The same values apply under every "
                               "insulation case.", "tab:outcome"),
                              ("year5", s5, "Year-5 incidence rate ratio at the same points, mean [95\\% interval].",
                               "tab:year5")]:
        with open(OUT / f"tables_outcome_{tag}.tex", "w") as fh:
            fh.write("\\begin{table}[htbp]\\centering\\small\n\\begin{tabular}{lccc}\\toprule\n"
                     "Scenario & $r_{PrEP}$ & $r_{VS}$ & IRR [95\\% interval]\\\\\\midrule\n")
            for k in ["P1", "P2", "P3"]:
                r = ss.loc[k]
                fh.write(f"{name[k]} & {P_R[k][0]:.2f} & {P_R[k][1]:.4f} & {fmt(r['mean'], r.lo, r.hi)}\\\\\n")
            fh.write(f"\\bottomrule\\end{{tabular}}\\caption{{{cap}}}\\label{{{lab}}}\\end{{table}}\n")
    with open(OUT / "tables_mapping.tex", "w") as fh:
        fh.write("\\begin{table}[htbp]\\centering\\small\n\\begin{tabular}{llccc}\\toprule\n"
                 "Insulation case & Scenario & $\\Delta_{prevention}$ & $\\Delta_{Medicaid}$ & Attainable ($\\Delta\\leq1$)\\\\\\midrule\n")
        for c in CASES:
            gp, ga = GAMMA[c]
            for i, k in enumerate(["P1", "P2", "P3"]):
                dp, dm = P_R[k][0] / (1 - gp), P_R[k][1] / (1 - ga)
                ok = dp <= 1 + 1e-9 and dm <= 1 + 1e-9
                lead = f"{LABEL[c]} ($\\gamma_{{PrEP}}={gp:g}$, $\\gamma_{{ART}}={ga:g}$)" if i == 0 else ""
                flag = "yes" if ok else "\\textbf{no}"
                fh.write(f"{lead} & {k} & {dp:.3f} & {dm:.3f} & {flag}\\\\\n")
            if c != CASES[-1]:
                fh.write("\\midrule\n")
        fh.write("\\bottomrule\\end{tabular}\\caption{Funding or coverage losses $\\Delta=r/(1-\\gamma)$ implied by the "
                 "provisional $r$-space points. The IRR at each point is identical in all three cases "
                 "(Table~\\ref{tab:outcome}); only the implied $\\Delta$ differs.}\\label{tab:mapping}\\end{table}\n")


v1 = r_surface("final", "fig_r_surface_final.png",
               "Incidence rate ratio over the coverage-reduction domain, final year")
v5 = r_surface("year5", "fig_r_surface_year5.png",
               "Incidence rate ratio over the coverage-reduction domain, year 5")
delta_map("fig_delta_mapping_final.png")
slice_figure("fig_funding_slices_final.png")
points_figure("fig_points_final.png")
tables()
print("max IRR:", round(v1, 3), round(v5, 3))
print(pts_summary("final")[0].round(3).to_string())
