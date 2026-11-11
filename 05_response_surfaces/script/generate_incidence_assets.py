"""Absolute-incidence tables (main Table 3, supplementary Table S4) and the incidence ribbon figure.

Reads <private>/out_summary/trajectory_summary.csv (export_incidence_trajectories.R).
Writes, all outside Git:
  out_summary/tab3_incidence_yr10.tex, out_summary/tabS4_incidence_yr5_yr10.tex
  out_summary/overleaf_figures/incidence_trajectories_elicited_medians.{png,pdf}
  out_summary/overleaf_figures/incidence_trajectories_elicited_medians_with_irr.{png,pdf}  (variant)
  out_summary/for_anna/incidence_additions_for_paper.tex  (paste-ready tables, figure, caption, text)
Year 10 = surrogate row 11 and year 5 = row 6, as in the rest of stage 05.
"""
import os, pathlib, sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

pdir = os.environ.get("INFORM_PRIVATE_DIR")
if not pdir:
    sys.exit("set INFORM_PRIVATE_DIR to a folder outside Git")
P = pathlib.Path(pdir).resolve()
if any((q / ".git").exists() for q in [P, *P.parents]):
    sys.exit(f"REFUSED: {P} is inside a Git repository")
O = P / "out_summary"
FIG = O / "overleaf_figures"; FIG.mkdir(exist_ok=True)
ANNA = O / "for_anna"; ANNA.mkdir(exist_ok=True)

T = pd.read_csv(O / "trajectory_summary.csv").set_index(["point", "year"])

ROWS = [  # (label, point id) in Anna's table order; None marks a section heading
    ("Baseline (no reduction)", "baseline"),
    (None, "Reference reductions"),
    ("Viral suppression $-$10\\%", "ref_vs_10"), ("Viral suppression $-$25\\%", "ref_vs_25"),
    ("Viral suppression $-$40\\%", "ref_vs_40"),
    ("PrEP use $-$10\\%", "ref_prep_10"), ("PrEP use $-$25\\%", "ref_prep_25"), ("PrEP use $-$40\\%", "ref_prep_40"),
    ("Both $-$10\\%", "ref_both_10"), ("Both $-$25\\%", "ref_both_25"), ("Both $-$40\\%", "ref_both_40"),
    (None, "Stakeholder-elicited reductions"),
    ("PrEP component (median)", "elic_P2"), ("ART component (median)", "elic_P1"),
    ("\\textbf{Both (median)}", "elic_P3"),
    ("Both, marginal minima", "elic_low"), ("Both, marginal maxima", "elic_high"),
]


def ci(pt, yr, q):
    r = T.loc[(pt, yr)]
    return f"{r[q + '_mean']:.2f} [{r[q + '_lo']:.2f}, {r[q + '_hi']:.2f}]"


def irr(pt, yr):
    return "1.00 (ref.)" if pt == "baseline" else ci(pt, yr, "irr")


def rvals(pt):
    r = T.loc[(pt, 10)]
    return f"{r.r_prep:.3f}", f"{r.r_vs:.3f}"


def bold(cells, on):
    return [f"\\textbf{{{c}}}" if on else c for c in cells]


NOTE_CI = ("Rates and IRRs are summarised from the same draws: the rate intervals also carry "
           "uncertainty in the baseline incidence level across calibrations, so they overlap "
           "the baseline even where the IRR interval excludes 1.")

body3 = []
for lab, pt in ROWS:
    if lab is None:
        body3 += ["\\addlinespace", f"\\multicolumn{{5}}{{@{{}}l}}{{\\textit{{{pt}}}}} \\\\"]
        continue
    c = bold([*rvals(pt), ci(pt, 10, "inc"), irr(pt, 10)], pt == "elic_P3")
    body3.append(f"{lab} & " + " & ".join(c) + " \\\\")
    if pt == "baseline":
        body3.append("\\addlinespace")
body3 = [b for i, b in enumerate(body3) if not (b == "\\addlinespace" and i and body3[i - 1] == "\\addlinespace")]

tab3 = r"""\begin{table}[htbp]
\centering
\begin{threeparttable}
\caption{Projected HIV incidence at 10 years under reductions in viral
  suppression and PrEP use.}
\label{tab:results}
\small
\setlength{\tabcolsep}{4pt}
\begin{tabular}{@{}l r r l l@{}}
\toprule
 & & & Incidence rate & \\
Scenario & $r_{\mathrm{PrEP}}$ & $r_{VS}$
  & per 100 PY & IRR  \\
  &  & & [95\% PPI] & [95\% PPI] \\
\midrule
""" + "\n".join(body3) + r"""
\bottomrule
\end{tabular}
\begin{tablenotes}[flushleft]\footnotesize
\item PY, person-years; IRR, incidence rate ratio; PPI, posterior prediction interval. Reference reductions are round-number values spanning the surrogate training domain and are not tied to any elicited quantity. Elicited rows use the median of the respondent-level effective reductions $r_i = q_i(1-\gamma_i)$ in each channel; the final two rows combine the marginal minima and the marginal maxima of the two channels, a pairing that no single respondent need have reported, and describe the width of the marginal ranges rather than an interval for the population. All values depend only on the coverage reductions and are identical under any assumption about the proportion able to retain access by other means. Intervals are 95\% intervals over 100 calibration checkpoints and 50 surrogate draws with common random numbers; they do not include uncertainty in the stakeholder responses. """ + NOTE_CI + r"""
\end{tablenotes}
\end{threeparttable}
\end{table}
"""

bodyS = []
for lab, pt in ROWS:
    if lab is None:
        bodyS += ["\\addlinespace", f"\\multicolumn{{7}}{{@{{}}l}}{{\\textit{{{pt}}}}} \\\\"]
        continue
    rate = lambda y: f"{T.loc[(pt, y)].inc_mean:.2f}"
    c = bold([*rvals(pt), rate(5), irr(pt, 5), rate(10), irr(pt, 10)], pt == "elic_P3")
    bodyS.append(f"{lab} & " + " & ".join(c) + " \\\\")
    if pt == "baseline":
        bodyS.append("\\addlinespace")
bodyS = [b for i, b in enumerate(bodyS) if not (b == "\\addlinespace" and i and bodyS[i - 1] == "\\addlinespace")]

tabS4 = r"""\begin{table}[htbp]
\centering
\begin{threeparttable}
\caption{Projected HIV incidence at years 5 and 10 under reductions in viral
  suppression and PrEP use.}
\label{tab:supp-results-5yr-10yr}
\footnotesize
\setlength{\tabcolsep}{4pt}
\begin{tabular}{@{}l r r r l r l@{}}
\toprule
 & & & \multicolumn{2}{c}{Year 5} & \multicolumn{2}{c}{Year 10} \\
\cmidrule(lr){4-5}\cmidrule(l){6-7}
Scenario & $r_{\mathrm{PrEP}}$ & $r_{VS}$
  & \multicolumn{1}{c}{Rate} & IRR [95\% interval]
  & \multicolumn{1}{c}{Rate} & IRR [95\% interval] \\
\midrule
""" + "\n".join(bodyS) + r"""
\bottomrule
\end{tabular}
\begin{tablenotes}[flushleft]\footnotesize
\item Rate, incidence per 100 person-years; IRR, incidence rate ratio relative to baseline at the same horizon. This table reproduces Table~\ref{tab:results} with the addition of year-5 values. Reference reductions are round-number values spanning the surrogate training domain. Elicited rows use the median of the respondent-level effective reductions $r_i = q_i(1-\gamma_i)$ in each channel; the final two rows combine the marginal minima and maxima of the two channels, a pairing that no single respondent need have reported. All values depend only on the coverage reductions and are identical under any assumption about the proportion able to retain access by other means.
\end{tablenotes}
\end{threeparttable}
\end{table}
"""
(O / "tab3_incidence_yr10.tex").write_text(tab3)
(O / "tabS4_incidence_yr5_yr10.tex").write_text(tabS4)

# ---- Ribbon figure: absolute incidence, baseline + three elicited medians -------------------------
SERIES = [("baseline", "Baseline (no reduction)", "#555555", "--"),
          ("elic_P2", "PrEP component (median)", "#2a9d8f", "-"),
          ("elic_P1", "ART component (median)", "#e07a1f", "-"),
          ("elic_P3", "Both (median)", "#6a3d8f", "-")]
plt.rcParams.update({"font.size": 11, "axes.labelsize": 12, "savefig.dpi": 300})


def draw(ax, q, ylabel):
    for pt, lab, col, ls in SERIES:
        d = T.loc[pt].reset_index()
        if not (q == "irr" and pt == "baseline"):
            ax.fill_between(d.year, d[q + "_lo"], d[q + "_hi"], color=col, alpha=0.13, lw=0)
        ax.plot(d.year, d[q + "_mean"], color=col, ls=ls, lw=2.2, label=lab)
    ax.set_xlabel("Years since reduction"); ax.set_ylabel(ylabel)
    ax.set_xticks(range(0, 11, 2)); ax.set_xlim(0, 10)
    ax.spines[["top", "right"]].set_visible(False); ax.grid(axis="y", alpha=0.25)


fig, ax = plt.subplots(figsize=(7.2, 5.0))
draw(ax, "inc", "HIV incidence per 100 person-years")
ax.set_ylim(bottom=0)
ax.legend(loc="lower left", frameon=False, fontsize=10, ncol=2)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(FIG / f"incidence_trajectories_elicited_medians.{ext}", bbox_inches="tight")
plt.close(fig)

fig, axes = plt.subplots(1, 2, figsize=(12, 5.0))
draw(axes[0], "inc", "HIV incidence per 100 person-years"); axes[0].set_ylim(bottom=0)
draw(axes[1], "irr", "Incidence rate ratio vs baseline")
for a, t in zip(axes, "AB"):
    a.text(-0.12, 1.02, f"({t})", transform=a.transAxes, fontsize=13, fontweight="bold")
axes[0].legend(loc="lower left", frameon=False, fontsize=10, ncol=2)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(FIG / f"incidence_trajectories_elicited_medians_with_irr.{ext}", bbox_inches="tight")
plt.close(fig)

# ---- Paste-ready LaTeX for Anna (separate from the paper) -----------------------------------------
v = lambda pt, y, q="inc": T.loc[(pt, y)]
b10, p3 = v("baseline", 10), v("elic_P3", 10)
fmt = lambda r, q="inc": f"{r[q + '_mean']:.2f} (95\\% PPI {r[q + '_lo']:.2f}--{r[q + '_hi']:.2f})"
text = (
    "Under no reduction, projected incidence remained roughly stable over the ten-year horizon, at "
    f"{fmt(b10)} per 100 person-years in year 10 (Table~\\ref{{tab:results}}). At the elicited medians, "
    f"year-10 incidence rose to {fmt(v('elic_P2', 10))} for the PrEP component alone, "
    f"{fmt(v('elic_P1', 10))} for the ART component alone, and {fmt(p3)} for both together, "
    f"an absolute increase of {p3.inc_mean - b10.inc_mean:.2f} per 100 person-years over baseline. "
    "The trajectories separate within the first two years and continue to diverge thereafter "
    "(Figure~\\ref{fig:incidence-trajectories}), so the effect of the reductions accumulates over time rather than appearing "
    "as a one-off step."
)
fig_env = r"""\begin{figure}[htbp]
\centering
\includegraphics[width=0.8\linewidth]{incidence_trajectories_elicited_medians}
\caption{Projected HIV incidence per 100 person-years over ten years with no reduction (baseline) and at the
  elicited median effective reductions for the PrEP component, the ART component, and both together. Lines are
  means and shaded bands 95\% posterior prediction intervals over 100 calibration checkpoints and 50 surrogate
  draws. The bands include uncertainty in the baseline incidence level, which is shared by all scenarios; the
  incidence rate ratios in Table~\ref{tab:results} are computed draw by draw and are more precise than the overlap
  of the bands suggests.}
\label{fig:incidence-trajectories}
\end{figure}
"""
doc = r"""\documentclass[11pt]{article}
\usepackage[margin=1in]{geometry}
\usepackage{booktabs,threeparttable,graphicx,caption}
\graphicspath{{../overleaf_figures/}}
\begin{document}
\section*{Additions for the paper: absolute incidence}
\raggedright
Copy each block between the \texttt{\%\% ---} markers into the paper. Tables need \texttt{booktabs} and
\texttt{threeparttable}. Year 10 is the final surrogate year; year 5 is the midpoint.

%% --- Suggested results text ----------------------------------------------------------------------
\subsection*{Suggested results text}
""" + text + r"""

%% --- Main text Table 3 ---------------------------------------------------------------------------
""" + tab3 + r"""
%% --- Main text or supplementary figure: incidence trajectories ------------------------------------
""" + fig_env + r"""
%% --- Supplementary Table S4 ----------------------------------------------------------------------
""" + tabS4 + r"""
\end{document}
"""
(ANNA / "incidence_additions_for_paper.tex").write_text(doc)
print("wrote tab3_incidence_yr10.tex, tabS4_incidence_yr5_yr10.tex, two ribbon figures, for_anna/incidence_additions_for_paper.tex")
