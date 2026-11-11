#!/usr/bin/env python3
"""Render one anonymous validation PDF per respondent.

Usage (data and outputs must be OUTSIDE any Git repository):
  python3 respondent_validation.py --config /path/outside/git/config.json --out-dir /path/outside/git/pdfs

The config maps CSV columns to q_prep, gamma_prep, q_art, gamma_art (see
config_template.json, which uses synthetic column names). Each PDF shows only
that respondent's own answers, their constructed effective reductions
r = q(1-gamma), an r-space plot, the outcomes at their own P1/P2/P3, and a
question about whether the stated reduction already accounted for alternative
access. No network access is used (sockets are disabled below); the surrogate
surface is read from the tracked, survey-free grid.
"""
import argparse
import socket
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

import respondent_core as core


def _no_network(*a, **k):
    raise RuntimeError("network access is disabled in this utility")


socket.socket = _no_network
socket.create_connection = _no_network
socket.getaddrinfo = _no_network


def fmt_irr(surf, rp, rv):
    if rp > core.R_MAX + 1e-9 or rv > core.R_MAX + 1e-9:
        return "outside model range (r > 0.75)"
    m, lo, hi = (float(surf.at(rp, rv, c)) for c in ("mean", "lower", "upper"))
    return f"{m:.2f} [{lo:.2f}, {hi:.2f}]"


def render(row, surf, out_pdf, status):
    q_p, g_p, q_a, g_a = row.q_prep, row.gamma_prep, row.q_art, row.gamma_art
    r_p, r_v = row.r_prep, row.r_vs
    fig = plt.figure(figsize=(8.27, 11.69))
    fig.text(0.07, 0.965, f"Respondent {row.id}: anonymous validation sheet", fontsize=15, fontweight="bold")
    if status == "synthetic":
        fig.text(0.07, 0.94, "SYNTHETIC TEST DATA, NOT SURVEY RESULTS", fontsize=11, color="#b00020", fontweight="bold")
    p = lambda v: f"{100*v:.0f}%"
    fig.text(0.07, 0.905,
             "Your answers\n"
             f"  PrEP:  reduction in use (q) = {p(q_p)};  access by other means (gamma) = {p(g_p)}\n"
             f"         constructed reduction  r = q(1 - gamma) = {p(r_p)}\n"
             f"  ART adherence:  reduction (q) = {p(q_a)};  access by other means (gamma) = {p(g_a)}\n"
             f"         constructed reduction  r = q(1 - gamma) = {p(r_v)}",
             fontsize=9, va="top", family="DejaVu Sans")
    ax = fig.add_axes([0.10, 0.50, 0.72, 0.31])
    x, y = surf.r_prep_axis, surf.r_vs_axis
    m = ax.pcolormesh(x, y, surf.mean_grid, cmap="YlOrRd", vmin=1.0, shading="gouraud")
    ax.contour(x, y, surf.mean_grid, levels=np.round(np.arange(1.1, np.nanmax(surf.mean_grid), 0.2), 2),
               colors="#444444", linewidths=0.5)
    for name, (a, b) in {"P1": (0.0, r_v), "P2": (r_p, 0.0), "P3": (r_p, r_v)}.items():
        ax.plot(a, b, "s" if name != "P3" else "o", ms=7, mfc="white", mec="black", mew=1.6, clip_on=False, zorder=6)
        ax.annotate(name, (a, b), xytext=(5, 5), textcoords="offset points", fontsize=9, fontweight="bold")
    ax.plot(q_p, q_a, "D", ms=7, mfc="none", mec="#1d6a85", mew=1.8, zorder=6)
    ax.annotate("if r = q", (q_p, q_a), xytext=(5, -12), textcoords="offset points", fontsize=8, color="#1d6a85")
    ax.set_xlim(0, core.R_MAX); ax.set_ylim(0, core.R_MAX)
    ax.set_xlabel("Proportional reduction in PrEP use (r_PrEP)"); ax.set_ylabel("Proportional reduction in viral suppression (r_VS)")
    cax = fig.add_axes([0.85, 0.50, 0.02, 0.31]); fig.colorbar(m, cax=cax).set_label("Incidence rate ratio, post-intervention year 10", fontsize=8)
    rows = [("P1  ART only", 0.0, r_v), ("P2  PrEP only", r_p, 0.0), ("P3  both", r_p, r_v)]
    lines = ["Model outcomes at your constructed reductions (incidence rate ratio vs status quo, year 10; 95% interval):", ""]
    for nm, a, b in rows:
        lines.append(f"  {nm:<14} r_PrEP = {100*a:4.0f}%  r_VS = {100*b:4.0f}%   IRR = {fmt_irr(surf, a, b)}")
    lines += ["", "Sensitivity, if your stated reductions already included alternative access (r = q):",
              f"  P3 both        r_PrEP = {100*q_p:4.0f}%  r_VS = {100*q_a:4.0f}%   IRR = {fmt_irr(surf, q_p, q_a)}"]
    fig.text(0.07, 0.455, "\n".join(lines), fontsize=8.2, va="top", family="DejaVu Sans Mono")
    q_txt = (f"Question for you: When you said PrEP use would fall by {p(q_p)} and ART adherence by {p(q_a)}, did those "
             f"figures already take into account the people who would keep access through other means (you said {p(g_p)} "
             f"for PrEP and {p(g_a)} for ART)? If yes, we would use your stated reductions directly (r = q). If no, we "
             f"apply the other-access shares, r = q(1 - gamma), as shown above.")
    fig.text(0.07, 0.30, "\n".join(textwrap.wrap(q_txt, 82)), fontsize=9.5, va="top", fontweight="bold")
    note = ("Notes. Incidence rate ratio = incidence relative to the modeled status quo at post-intervention year 10, from the "
            "surrogate model; interval = 95% model uncertainty. P1, P2 and P3 are your own ART-only, PrEP-only and combined "
            "reductions. ART adherence is used as a stand-in for viral suppression. This sheet shows only your own answers.")
    fig.text(0.07, 0.18, "\n".join(textwrap.wrap(note, 120)), fontsize=7.5, va="top", color="#444444")
    with PdfPages(out_pdf) as pdf:
        pdf.savefig(fig)
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--anonymize", action="store_true", help="replace ids by R01.. and write id_key.csv to --out-dir")
    a = ap.parse_args(argv)
    cfg = core.load_config(a.config)
    out_dir = core.assert_outside_git(a.out_dir, "output directory")
    out_dir.mkdir(parents=True, exist_ok=True)
    d = core.effective_r(core.load_respondents(cfg))
    if a.anonymize:
        key = d[["id"]].copy(); d["id"] = [f"R{i+1:02d}" for i in range(len(d))]
        key.insert(0, "anonymous_id", d["id"]); key.to_csv(out_dir / "id_key.csv", index=False)
    surf = core.Surface("final")
    for row in d.itertuples(index=False):
        safe = "".join(ch for ch in str(row.id) if ch.isalnum() or ch in "-_")
        render(row, surf, out_dir / f"respondent_{safe}.pdf", cfg["data_status"])
    print(f"wrote {len(d)} PDFs to {out_dir}")


if __name__ == "__main__":
    main()
